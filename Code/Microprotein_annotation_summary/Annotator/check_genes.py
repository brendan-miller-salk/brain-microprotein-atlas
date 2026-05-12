#!/usr/bin/env python3
import re

ATTR_RE = re.compile(r'(\S+) "([^"]+)"')
gene_cds_start = {}
gene_cds_end = {}
targets = {'ENSG00000140443', 'ENSG00000165506'}

gtf_path = '/Users/brendanmiller/Library/CloudStorage/Box-Box/post_shortstop_processing/annotation/ensemlb_hg38_filtered.gtf'

with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS':
            continue
        attrs = dict(ATTR_RE.findall(parts[8]))
        gid = attrs.get('gene_id')
        if gid not in targets:
            continue
        strand = parts[6]
        start, end = int(parts[3]), int(parts[4])
        cds_start_pos = start if strand == '+' else end
        cds_end_pos = end if strand == '+' else start

        if gid not in gene_cds_start:
            gene_cds_start[gid] = (cds_start_pos, strand)
            gene_cds_end[gid] = (cds_end_pos, strand)
        else:
            ex_s, _ = gene_cds_start[gid]
            ex_e, _ = gene_cds_end[gid]
            if strand == '+':
                gene_cds_start[gid] = (min(ex_s, start), strand)
                gene_cds_end[gid] = (max(ex_e, end), strand)
            else:
                gene_cds_start[gid] = (max(ex_s, end), strand)
                gene_cds_end[gid] = (min(ex_e, start), strand)

for gid in sorted(targets):
    print(f'{gid}:')
    print(f'  CDS start (most upstream): {gene_cds_start.get(gid)}')
    print(f'  CDS end (most downstream): {gene_cds_end.get(gid)}')

# Now check: for chr15:98960643-98961071 on + strand (IGF1R)
# orf_bio_start = 98960643
# Is it after cds_end? -> dORF?
# Does it fall in CDS range? -> why doORF?
print()
print("=== Case 1: IGF1R chr15:98960643-98961071 (+strand) ===")
gid = 'ENSG00000140443'
if gid in gene_cds_start:
    cs = gene_cds_start[gid][0]
    ce = gene_cds_end[gid][0]
    strand = gene_cds_start[gid][1]
    print(f"  Gene strand: {strand}")
    print(f"  Gene CDS start (upstream): {cs}")
    print(f"  Gene CDS end (downstream): {ce}")
    orf_start = 98960643
    orf_end = 98961071
    orf_bio_start = orf_start  # + strand
    print(f"  ORF bio start: {orf_bio_start}")
    print(f"  orf_start > cds_end? {orf_start > ce} -> dORF if true")
    print(f"  orf_bio_start in CDS range [{min(cs,ce)}, {max(cs,ce)}]? {min(cs,ce) <= orf_bio_start <= max(cs,ce)}")

print()
print("=== Case 2: DNAAF2 chr14:49633832-49634167 (-strand) ===")
gid = 'ENSG00000165506'
if gid in gene_cds_start:
    cs = gene_cds_start[gid][0]
    ce = gene_cds_end[gid][0]
    strand = gene_cds_start[gid][1]
    print(f"  Gene strand: {strand}")
    print(f"  Gene CDS start (upstream): {cs}")
    print(f"  Gene CDS end (downstream): {ce}")
    orf_start = 49633832
    orf_end = 49634167
    orf_bio_start = orf_end  # - strand
    print(f"  ORF bio start: {orf_bio_start}")
    print(f"  orf_end < cds_end? {orf_end < ce} -> dORF if true (- strand)")
    print(f"  orf_bio_start in CDS range [{min(cs,ce)}, {max(cs,ce)}]? {min(cs,ce) <= orf_bio_start <= max(cs,ce)}")
