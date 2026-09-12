#!/usr/bin/env python3

import re
import csv
import bisect
import argparse
from collections import defaultdict

# === Precompiled Regex for GTF/GFF attributes ===
ATTR_RE = re.compile(r'(\S+) "([^"]+)"')

# Regex to extract genomic coordinates from gene_id
# e.g., ESPRESSO_chr1_0_13-chr1:16748-17310_F:2_P:1_M
COORD_RE = re.compile(r'(chr[\dXYM]+):(\d+)-(\d+)')

def parse_attributes(attr_str):
    return dict(ATTR_RE.findall(attr_str))

def make_ucsc_link(gene_id, db='hg38'):
    """
    Extract genomic coordinates from gene_id and return a clickable UCSC hyperlink
    formula for Excel/Google Sheets.
    """
    match = COORD_RE.search(gene_id)
    if match:
        chrom = match.group(1)
        start = match.group(2)
        end = match.group(3)
        position = f"{chrom}:{start}-{end}"
        url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={db}&position={position}"
        return f'=HYPERLINK("{url}", "{position}")'
    return ''

class smORFAnnotator:
    def __init__(self, args):
        self.intersect_file = args.intersect_output
        self.non_intersect_file = args.non_intersect_output
        self.output_file = args.output_file
        self.ensembl_gtf = args.ensembl_gtf
        self.genome = getattr(args, 'genome', 'hg38')

    def _build_cds_start_lookup(self):
        """
        Parse the Ensembl GTF to build lookups:
        1. transcript_cds: transcript_id → (CDS start position, strand)
        2. gene_cds: gene_id → (most-upstream CDS start across ALL isoforms, strand)
        3. gene_cds_end: gene_id → (most-downstream CDS stop across ALL isoforms, strand)
        4. gene_all_cds_starts: gene_id → set of ALL CDS start positions across isoforms
        5. locus_cds_starts: (chrom, position) set of ALL CDS start positions genome-wide
        6. locus_cds_boundaries: chrom → list of (gene_id, most_upstream_start, most_downstream_end, strand)

        For + strand: CDS start = min start, CDS end = max end.
        For - strand: CDS start = max end, CDS end = min start.
        """
        transcript_cds = {}  # transcript_id → (position, strand)
        gene_cds = {}        # gene_id → (most_upstream_cds_start, strand)
        gene_cds_end = {}    # gene_id → (most_downstream_cds_stop, strand)
        gene_all_cds_starts = defaultdict(set)  # gene_id → {all CDS start positions}
        locus_cds_starts = set()  # (chrom, position) for ALL CDS starts genome-wide
        gene_chrom = {}      # gene_id → chrom

        with open(self.ensembl_gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue

                strand = parts[6]
                start = int(parts[3])
                end = int(parts[4])
                attrs = parse_attributes(parts[8])
                transcript_id = attrs.get('transcript_id')
                gene_id = attrs.get('gene_id')
                if not transcript_id:
                    continue

                # --- Transcript-level: CDS start for this transcript ---
                cds_start_pos = start if strand == '+' else end
                if transcript_id not in transcript_cds:
                    transcript_cds[transcript_id] = (cds_start_pos, strand)
                else:
                    existing_pos, _ = transcript_cds[transcript_id]
                    if strand == '+':
                        transcript_cds[transcript_id] = (min(existing_pos, start), strand)
                    else:
                        transcript_cds[transcript_id] = (max(existing_pos, end), strand)

                # --- Gene-level: most UPSTREAM CDS start across all isoforms ---
                # + strand: most upstream = smallest start
                # - strand: most upstream = largest end
                if gene_id:
                    if gene_id not in gene_cds:
                        gene_cds[gene_id] = (cds_start_pos, strand)
                    else:
                        existing_pos, _ = gene_cds[gene_id]
                        if strand == '+':
                            gene_cds[gene_id] = (min(existing_pos, start), strand)
                        else:
                            gene_cds[gene_id] = (max(existing_pos, end), strand)

                    # --- Collect ALL CDS start positions per gene ---
                    gene_all_cds_starts[gene_id].add(cds_start_pos)

                    # --- Locus-level: collect ALL CDS starts genome-wide ---
                    chrom = parts[0]
                    locus_cds_starts.add((chrom, cds_start_pos))
                    gene_chrom[gene_id] = chrom

                    # --- Gene-level: most DOWNSTREAM CDS stop across all isoforms ---
                    # + strand: most downstream stop = largest end
                    # - strand: most downstream stop = smallest start
                    cds_end_pos = end if strand == '+' else start
                    if gene_id not in gene_cds_end:
                        gene_cds_end[gene_id] = (cds_end_pos, strand)
                    else:
                        existing_end, _ = gene_cds_end[gene_id]
                        if strand == '+':
                            gene_cds_end[gene_id] = (max(existing_end, end), strand)
                        else:
                            gene_cds_end[gene_id] = (min(existing_end, start), strand)

        # --- Build locus-level CDS intervals for range checking ---
        # Normalize to (lo, hi) per chromosome, merge overlapping intervals,
        # and precompute lo-array for O(log n) bisect lookups.
        # (handles readthrough genes where matched gene has narrow CDS but
        # a component gene at the same locus has a wider CDS range)
        raw_intervals = defaultdict(list)  # chrom → [(lo, hi), ...]
        for gid in gene_cds:
            if gid in gene_cds_end and gid in gene_chrom:
                gchrom = gene_chrom[gid]
                g_upstream = gene_cds[gid][0]
                g_downstream = gene_cds_end[gid][0]
                lo, hi = min(g_upstream, g_downstream), max(g_upstream, g_downstream)
                raw_intervals[gchrom].append((lo, hi))

        locus_cds_intervals = {}  # chrom → [(lo, hi), ...] merged
        locus_cds_lo = {}         # chrom → [lo1, lo2, ...] for bisect
        for chrom, ivs in raw_intervals.items():
            ivs.sort()
            merged = [ivs[0]]
            for lo, hi in ivs[1:]:
                if lo <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
                else:
                    merged.append((lo, hi))
            locus_cds_intervals[chrom] = merged
            locus_cds_lo[chrom] = [iv[0] for iv in merged]

        return transcript_cds, gene_cds, gene_cds_end, gene_all_cds_starts, locus_cds_starts, locus_cds_intervals, locus_cds_lo

    @staticmethod
    def _position_in_any_cds_range(chrom, strand, position, locus_cds_intervals, locus_cds_lo):
        """
        Check if a position falls within ANY gene's CDS range at this chromosome.
        Uses merged intervals + bisect for O(log n) lookup.
        """
        lo_arr = locus_cds_lo.get(chrom)
        if not lo_arr:
            return False
        idx = bisect.bisect_right(lo_arr, position) - 1
        if idx < 0:
            return False
        return locus_cds_intervals[chrom][idx][0] <= position <= locus_cds_intervals[chrom][idx][1]

    def process_gtf_files(self):
        transcript_cds, gene_cds, gene_cds_end, gene_all_cds_starts, locus_cds_starts, locus_cds_intervals, locus_cds_lo = self._build_cds_start_lookup()

        priority_order = ['psORF', 'uoORF', 'uaoORF', 'doORF', 'daoORF', 'iORF', 'dORF', 'daORF', 'uORF', 'uaORF', 'lncRNA', 'riORF', 'aORF', 'eORF', 'udORF']

        def get_priority(annot):
            return priority_order.index(annot) if annot in priority_order else float('inf')

        gene_data = defaultdict(lambda: ('UA', 'Unknown', 'Unknown'))

        PSEUDO_BIOTYPES = {
            'processed_pseudogene', 'unprocessed_pseudogene', 'translated_unprocessed_pseudogene',
            'translated_processed_pseudogene', 'transcribed_processed_pseudogene',
            'transcribed_unprocessed_pseudogene', 'unitary_pseudogene', 'polymorphic_pseudogene'
        }

        NONCODING_BIOTYPES = {'lncRNA', 'lincRNA', 'antisense', 'sense_intronic', 'sense_overlapping'}

        with open(self.intersect_file, 'r') as file:
            for line in file:
                if not line.startswith('chr'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 18:
                    continue
                if parts[2] != 'CDS':
                    continue

                # Extract smORF info (columns 0-8)
                attrs = parse_attributes(parts[8])
                gene_id = attrs.get('gene_id', 'Unknown')

                # Extract the TRUE ORF start/end from the gene_id string
                # parts[3]/parts[4] are EXON coordinates (include UTR), not ORF coordinates
                smorf_strand = parts[6]
                coord_match = COORD_RE.search(gene_id)
                if coord_match:
                    orf_start = int(coord_match.group(2))
                    orf_end = int(coord_match.group(3))
                else:
                    orf_start = int(parts[3])
                    orf_end = int(parts[4])

                # Biological ORF start (strand-aware)
                orf_bio_start = orf_start if smorf_strand == '+' else orf_end

                # Extract ENSEMBL info (columns 9-17)
                gene_info = parts[17] if len(parts) > 17 else parts[8]
                gene_attrs = parse_attributes(gene_info)

                ensembl_gene_id = gene_attrs.get('gene_id', 'Unknown')
                gene_name = gene_attrs.get('gene_name', 'Unnamed')
                gene_biotype = gene_attrs.get('gene_biotype', 'Unnamed')
                transcript_biotype = gene_attrs.get('transcript_biotype', 'Unknown')
                feature_type = parts[11]

                # Chromosome for locus-level lookups
                smorf_chrom = parts[0]

                # === Annotation Logic ===
                if feature_type == 'three_prime_utr':
                    # Step 1: Check if smORF shares start codon with ANY CDS at this locus → iORF
                    #         (locus-level: catches readthrough genes, overlapping genes, etc.)
                    # Step 2: Check if smORF starts after the most-downstream CDS stop → dORF (confident)
                    #         Otherwise → daORF (downstream ambiguous, isoform artifact)
                    if (smorf_chrom, orf_bio_start) in locus_cds_starts:
                        annotation = 'iORF'
                    else:
                        is_confident = False
                        if ensembl_gene_id and ensembl_gene_id in gene_cds_end:
                            cds_end_pos, cds_strand = gene_cds_end[ensembl_gene_id]
                            # smORF starts after the most downstream CDS stop
                            if smorf_strand == '+' and orf_start > cds_end_pos:
                                is_confident = True
                            elif smorf_strand == '-' and orf_end < cds_end_pos:
                                is_confident = True

                        if is_confident:
                            annotation = 'dORF'
                        else:
                            # Before assigning daORF, check for readthrough gene artifacts:
                            # If the matched gene's CDS does NOT contain orf_bio_start,
                            # check whether another gene's CDS at this locus does.
                            # (e.g. CHKB-CPT1B readthrough has narrow CDS, but component
                            #  gene CHKB's CDS contains the smORF start → iORF)
                            matched_gene_contains_start = False
                            if ensembl_gene_id in gene_cds and ensembl_gene_id in gene_cds_end:
                                g_up = gene_cds[ensembl_gene_id][0]
                                g_down = gene_cds_end[ensembl_gene_id][0]
                                g_lo, g_hi = min(g_up, g_down), max(g_up, g_down)
                                if g_lo <= orf_bio_start <= g_hi:
                                    matched_gene_contains_start = True

                            if not matched_gene_contains_start and \
                               self._position_in_any_cds_range(smorf_chrom, smorf_strand, orf_bio_start, locus_cds_intervals, locus_cds_lo):
                                annotation = 'iORF'
                            else:
                                annotation = 'daORF'
                elif feature_type == 'five_prime_utr':
                    # Step 1: Check if smORF shares its start codon with ANY CDS at this locus → iORF
                    # Step 2: Compare ORF start against the GENE-LEVEL most-upstream CDS start.
                    #         If ORF starts before the earliest CDS start → uORF (confident)
                    #         Otherwise → uaORF (only upstream for some isoforms, ambiguous)
                    if (smorf_chrom, orf_bio_start) in locus_cds_starts:
                        annotation = 'iORF'
                    else:
                        is_confident = False
                        if ensembl_gene_id and ensembl_gene_id in gene_cds:
                            cds_start_pos, cds_strand = gene_cds[ensembl_gene_id]
                            if smorf_strand == '+' and orf_start < cds_start_pos:
                                is_confident = True
                            elif smorf_strand == '-' and orf_end > cds_start_pos:
                                is_confident = True

                        if is_confident:
                            annotation = 'uORF'
                        else:
                            # Readthrough gene check (same as three_prime_utr block)
                            matched_gene_contains_start = False
                            if ensembl_gene_id in gene_cds and ensembl_gene_id in gene_cds_end:
                                g_up = gene_cds[ensembl_gene_id][0]
                                g_down = gene_cds_end[ensembl_gene_id][0]
                                g_lo, g_hi = min(g_up, g_down), max(g_up, g_down)
                                if g_lo <= orf_bio_start <= g_hi:
                                    matched_gene_contains_start = True

                            if not matched_gene_contains_start and \
                               self._position_in_any_cds_range(smorf_chrom, smorf_strand, orf_bio_start, locus_cds_intervals, locus_cds_lo):
                                annotation = 'iORF'
                            else:
                                annotation = 'uaORF'
                elif feature_type == 'CDS':
                    # smORF overlaps an Ensembl CDS → at minimum iORF
                    # But also check if ORF extends upstream/downstream of
                    # the gene's CDS boundaries (e.g. smORF exon structure
                    # may not overlap a five_prime_utr feature directly)
                    # First check if start matches any CDS start or falls within any gene's CDS
                    if (smorf_chrom, orf_bio_start) in locus_cds_starts:
                        annotation = 'iORF'
                    elif self._position_in_any_cds_range(smorf_chrom, smorf_strand, orf_bio_start, locus_cds_intervals, locus_cds_lo):
                        annotation = 'iORF'
                    else:
                        annotation = 'iORF'
                        if ensembl_gene_id and ensembl_gene_id in gene_cds:
                            cds_start_pos, _ = gene_cds[ensembl_gene_id]
                            if (smorf_strand == '+' and orf_start < cds_start_pos) or \
                               (smorf_strand == '-' and orf_end > cds_start_pos):
                                annotation = 'uoORF'
                        if annotation == 'iORF' and ensembl_gene_id and ensembl_gene_id in gene_cds_end:
                            cds_end_pos, _ = gene_cds_end[ensembl_gene_id]
                            if (smorf_strand == '+' and orf_end > cds_end_pos) or \
                               (smorf_strand == '-' and orf_start < cds_end_pos):
                                annotation = 'doORF'
                elif feature_type == 'retrotransposed':
                    annotation = 'psORF'
                elif transcript_biotype in PSEUDO_BIOTYPES:
                    annotation = 'psORF'
                elif any(x in transcript_biotype for x in [
                    'non_stop_decay', 'nonsense_mediated_decay',
                    'ambiguous_orf', 'protein_coding_CDS_not_defined']):
                    annotation = 'aORF'
                elif 'retained_intron' in transcript_biotype:
                    annotation = 'riORF'
                elif gene_biotype in NONCODING_BIOTYPES:
                    annotation = 'lncRNA'
                elif feature_type == 'exon':
                    annotation = 'eORF'
                else:
                    annotation = 'UA'

                # === Conflict Resolution ===
                if gene_id not in gene_data:
                    gene_data[gene_id] = (annotation, gene_name, ensembl_gene_id)
                else:
                    existing_annotation, _, _ = gene_data[gene_id]

                    # Handle combined upstream (confident) + overlapping
                    if {existing_annotation, annotation} == {'uORF', 'iORF'}:
                        gene_data[gene_id] = ('uoORF', gene_name, ensembl_gene_id)

                    # Handle combined upstream ambiguous + overlapping
                    elif {existing_annotation, annotation} == {'uaORF', 'iORF'}:
                        gene_data[gene_id] = ('uaoORF', gene_name, ensembl_gene_id)

                    # Handle combined downstream ambiguous + overlapping
                    elif {existing_annotation, annotation} == {'daORF', 'iORF'}:
                        gene_data[gene_id] = ('daoORF', gene_name, ensembl_gene_id)

                    # Handle combined downstream + overlapping
                    elif {existing_annotation, annotation} == {'dORF', 'iORF'}:
                        gene_data[gene_id] = ('doORF', gene_name, ensembl_gene_id)

                    # Handle upstream + downstream (rare)
                    elif {existing_annotation, annotation} == {'uORF', 'dORF'}:
                        gene_data[gene_id] = ('udORF', gene_name, ensembl_gene_id)

                    # Handle upstream ambiguous + downstream (rare)
                    elif {existing_annotation, annotation} == {'uaORF', 'dORF'}:
                        gene_data[gene_id] = ('udORF', gene_name, ensembl_gene_id)

                    # Handle upstream + downstream ambiguous (rare)
                    elif {existing_annotation, annotation} == {'uORF', 'daORF'}:
                        gene_data[gene_id] = ('udORF', gene_name, ensembl_gene_id)

                    # Handle upstream ambiguous + downstream ambiguous (rare)
                    elif {existing_annotation, annotation} == {'uaORF', 'daORF'}:
                        gene_data[gene_id] = ('udORF', gene_name, ensembl_gene_id)

                    # dORF + daORF → keep dORF (confident wins)
                    elif {existing_annotation, annotation} == {'dORF', 'daORF'}:
                        gene_data[gene_id] = ('dORF', gene_name, ensembl_gene_id)

                    # uORF + uaORF → keep uORF (confident wins)
                    elif {existing_annotation, annotation} == {'uORF', 'uaORF'}:
                        gene_data[gene_id] = ('uORF', gene_name, ensembl_gene_id)

                    # Keep whichever has higher priority (lower index)
                    elif get_priority(annotation) < get_priority(existing_annotation):
                        gene_data[gene_id] = (annotation, gene_name, ensembl_gene_id)

        # === Add Intergenic Genes ===
        with open(self.non_intersect_file, 'r') as file:
            for line in file:
                if not line.startswith('chr'):
                    continue
                parts = line.strip().split('\t')
                attrs = parse_attributes(parts[8])
                gene_id = attrs.get('gene_id', 'Unknown')
                if gene_id not in gene_data:
                    gene_data[gene_id] = ('intergenic', 'None', 'None')

        # === Write to CSV with UCSC links ===
        with open(self.output_file, 'w', newline='') as out:
            writer = csv.writer(out)
            writer.writerow(['gene_id', 'annotation', 'gene_name', 'ensembl_gene_id', 'ucsc_link'])
            for gene_id, (annotation, gene_name, ensembl_gene_id) in gene_data.items():
                ucsc_link = make_ucsc_link(gene_id, db=self.genome)
                writer.writerow([gene_id, annotation, gene_name, ensembl_gene_id, ucsc_link])

        print(f"[✓] Output written to {self.output_file}")

        # === Print summary ===
        summary = defaultdict(int)
        for annot, _, _ in gene_data.values():
            summary[annot] += 1

        print("\n[📊] Annotation summary:")
        for annot in sorted(summary.keys()):
            print(f"  {annot:12} : {summary[annot]}")


def main():
    parser = argparse.ArgumentParser(description='Brain smORF Annotator - CSV output with UCSC links')
    parser.add_argument('--intersect_output', required=True, help='Path to bedtools intersect output')
    parser.add_argument('--non_intersect_output', required=True, help='Path to bedtools non-intersect output')
    parser.add_argument('--output_file', required=True, help='Path to output CSV file')
    parser.add_argument('--ensembl_gtf', required=True, help='Path to Ensembl GTF file')
    parser.add_argument('--genome', default='hg38', help='UCSC genome assembly (default: hg38)')
    args = parser.parse_args()

    annotator = smORFAnnotator(args)
    annotator.process_gtf_files()


if __name__ == '__main__':
    main()
