# GTF_and_BED_files

Genome-browser and sequence files for the unreviewed brain microproteins — those absent
from SwissProt — from *"A microprotein atlas of the human frontal cortex in Alzheimer's
disease"* (Miller et al.).

## Two kinds of record

Every file here contains two record types.

| | What it is | Count |
|---|---|---|
| **Annotated-start** | The protein translated from the smORF's annotated start codon. | 4,814 |
| **Alternative-start** | A shorter protein translated from a downstream non-ATG start codon within the same smORF. | 867 |

Alternative-start proteoforms are additional products of smORFs already in the set, not
new loci. The 867 come from 681 parent smORFs.

**Telling them apart:**

| Format | Annotated-start | Alternative-start |
|---|---|---|
| GTF | column 2 is `GTF2FastaPatched` | column 2 is `AltInitiation` |
| FASTA | 4 header fields | 5th field is `ALT_INIT` |
| BED, IDs | — | name ends in `_alt_initiation_N` |

```bash
grep -v _alt_initiation_ Unreviewed_Brain_Microproteins_IDs.txt            # annotated-start only
awk -F'\t' '$2=="AltInitiation"' Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf
```

Where one smORF has two alternative-start entries (`_alt_initiation_1` and `_2`), they
are the same start site written twice — once with the native residue at position 1 and
once with methionine substituted. They share coordinates by design.

---

## Files

### Main set

| File | Size | Contents |
|---|---|---|
| `Unreviewed_Brain_Microproteins.fasta` | 762 K | Protein sequences. **5,681 records** = 4,814 annotated-start + 867 alternative-start. |
| `Unreviewed_Brain_Microproteins_IDs.txt` | 256 K | One identifier per line, matching the FASTA. 5,681 lines. |
| `Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed` | 522 K | CDS blocks, BED12. **4,553 records** = 3,686 annotated-start + 867 alternative-start. |
| `Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed` | 595 K | Every non-Swiss-Prot microprotein, BED12. **5,668 records** = the 4,553 above + 1,115 TrEMBL entries, colour-coded by type. See the note below. |
| `Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf` | 6.4 M | Gene features. **33,169 lines** = 32,002 annotated-start + 1,167 alternative-start CDS. |
| `Unreviewed_Brain_Microproteins_genomic_coordinates.txt` | 127 K | Unique `chrN:start-end` loci, one per line. 5,436 lines. |
| `Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv` | 507 K | Two columns: `genomic_coordinates`, `sequence`. 5,668 rows. |
| `Ensembl_and_Unreviewed_Brain_Microproteins.gtf` | 890 M | GENCODE v43 plus everything above, for loading a single annotation into a browser. ~1.98 M lines. Too large for GitHub — download from [Zenodo](https://doi.org/10.5281/zenodo.20045161), or rebuild locally by concatenating the GENCODE v43 GTF with `Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf`. |

`_All_Non_SwissProt.bed` carries a `track` line and colours its three record
types through `itemRgb` (column 9), so a browser separates them at a glance:

| Colour | itemRgb | Records | What it is |
|---|---|---|---|
| Blue | `31,78,121` | 3,686 | Annotated-start microprotein, real CDS blocks |
| Purple | `130,80,223` | 867 | Alternative-start proteoform, real CDS blocks |
| Amber | `217,119,6` | 1,115 | TrEMBL entry, single-interval CDS **envelope** |

```bash
awk -F'\t' '$9=="217,119,6"' Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed  # TrEMBL
awk -F'\t' '$9!="217,119,6" && /\t/' Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed  # real CDS
```

> **Amber blocks span introns.** TrEMBL microproteins are unreviewed UniProt
> entries carrying no transcript id, so nothing ties them to CDS blocks in our
> annotation. Each has one `chrN:start-end` interval, and that interval is a CDS
> *envelope* — start codon to stop codon with the introns left in — not a gene
> body: only 7 of the 1,115 equal an Ensembl gene span, while 536 equal some
> Ensembl transcript's CDS envelope exactly (and for 521 of those the transcript's
> CDS length is exactly 3x the microprotein length). Splicing is why an envelope
> runs a median ~15x the coding length. The block is drawn thick across the whole
> interval, so **colour is the only thing distinguishing it from real CDS** —
> strip column 9 and that information is gone.
>
> TrEMBL records are named by UniProt accession rather than by locus string.
> Strand comes from the Ensembl gene of the same symbol; 30 of the 1,115 had no
> match and carry `.`. A further 13 TrEMBL microproteins have no coordinates at
> all in the master table and appear only in the FASTA and ID list.


### `PROSIT_confidence_tiers/`

The same records split by PROSIT spectral-angle confidence — how closely the observed
MS/MS spectrum matches the predicted spectrum for the identifying peptide. Each tier has
a `.fasta`, `.bed`, `.gtf`, and an `_IDs_*.txt`.

| Tier | Records | annotated-start | alternative-start |
|---|---|---|---|
| `Strong` | 1,328 | 1,067 | 261 |
| `Moderate` | 1,847 | 1,433 | 414 |
| `Weak` | 456 | 343 | 113 |
| `Insufficient` | 203 | 158 | 45 |

Tiers are mutually exclusive, and they do **not** sum to the full set: records with no
PROSIT evidence appear in none of them (1,813 annotated-start, 34 alternative-start).

For an alternative-start proteoform the tier comes from its own supporting peptides —
those downstream of its start codon — not from the parent smORF.

> The `initiation_tier` attribute in the GTF (`1`, `2`, `3a`, `3b`, `3c`) is a different
> measure: confidence in the start site, not spectral match quality.

### `Ribo_ShortStop/`

Microproteins with Ribo-seq translation support that were **not** detected by DDA mass
spectrometry — the "Ribo-Seq Only" group in Figure 7, panel J.

| File | Records |
|---|---|
| `Unreviewed_Brain_Microproteins_Ribo_ShortStop.fasta` | 1,592 |
| `Unreviewed_Brain_Microproteins_IDs_Ribo_ShortStop.txt` | 1,592 |
| `Unreviewed_Brain_Microproteins_CDS_Ribo_ShortStop.bed` | 965 |
| `Unreviewed_Brain_Microproteins_Ribo_ShortStop.gtf` | 7,374 lines |

Contains no alternative-start records: every one of those was detected by mass
spectrometry, so none meet this set's criteria.

### `Alt_Proteoforms/`

The 867 alternative-start proteoforms on their own. The `.fasta`, `.gtf`, `.bed`, and
`_IDs.txt` here are identical to the alternative-start records already inside the main
files, extracted for convenience. One file is only available here:

| File | Contents |
|---|---|
| `Alt_Proteoforms_mapping.tsv` | Per-proteoform metadata, 20 columns: `alt_protein_id`, `parent_orf_id`, `parent_gene_id`, `gene_name`, `chromosome`, `strand`, `cds_start`, `cds_end`, `n_cds_blocks`, `aa_position`, `aa_length`, `codon`, `initiation_tier`, `cognate_status`, `context_strength`, `initiation_variant`, `includes_stop`, `prosit_grade`, `n_support_peptides`, `sequence`. |

---

## Formats

### FASTA

One sequence per line, unwrapped. Fields are `|`-delimited:
`id | gene_name | database | annotation_status`, then for alternative-start records
`ALT_INIT | codon= | init_tier=`. Files in `PROSIT_confidence_tiers/` add `PROSIT=<tier>`.

```
>ENST00000448905.6+chrX:150983407-150985604_F:1_P:0|HMGB3|Salk|MS
AAAAAAAAAAAAAAAAAPAQLSGPRRRERALPPPPPASTPHSLPAPSAFPPARRPGRKKQFSQDG
>ENST00000275227.9-chr6:132798497-132798556_F:1_P:1_alt_initiation_1|SLC18B1|Salk|MS|ALT_INIT|codon=GCG|init_tier=3a
MAARQSASGSLEVESKDVA
```

`codon=` is the genomic start triplet. Alternative-start sequences begin with `M` even
when that triplet is near-cognate, since initiator tRNA delivers methionine.

### BED

All `.bed` files are **BED12** — 12 columns, tab-separated, no header. Coordinates are
0-based half-open (`chromStart = GTF start − 1`).

```
chr22	50258300	50261070	ESPRESSO_chr22_19938_29-chr22:50258301-50261136_F:0_P:0_alt_initiation_1	0	-	50258300	50261070	0	2	1,242,	0,2528,
```

Columns: `chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd,
itemRgb, blockCount, blockSizes, blockStarts`. Score is always `0`;
`thickStart`/`thickEnd` equal `chromStart`/`chromEnd` because every block is coding.
Blocks ascend by genomic coordinate on both strands, per the BED spec.

**The blocks are the CDS exons, not the outer span.** 1,303 of the 4,553 records are
spliced, and across the whole file `chromStart`–`chromEnd` covers 11.7 Mb while the
blocks total 809 Kb — so ~93% of the enclosed range is intron. Tools that read blocks
(`bedtools getfasta -split`, IGV, the UCSC browser) extract only coding sequence;
tools that ignore them and use columns 2–3 alone will pull in introns.

### GTF

1-based inclusive coordinates.

```
chr11	AltInitiation	CDS	840407	840458	1000	+	0	gene_id "ESPRESSO_chr11_4460_1+chr11:840398-842483_F:2_P:2"; transcript_id "...F:2_P:2_alt_initiation_1"; protein_id "...F:2_P:2_alt_initiation_1"; gene_name "Intergenic"; start_codon_seq "CGA"; aa_position "4"; initiation_tier "3b"; cognate_status "non-near-cognate"; includes_stop "false"; prosit_grade "NA"; smorf_type "eORF";
```

Annotated-start records (`GTF2FastaPatched`) carry `transcript_id` and `gene_id` only,
with `.` in the frame column. Alternative-start records (`AltInitiation`) carry the
attributes above plus a computed CDS phase, and are **CDS features only** — no `exon` or
`transcript` lines. Tools that read `exon` features, such as RSEM, salmon, StringTie and
featureCounts, will therefore not see them. Tools that need a transcript-level feature,
such as VEP and SnpEff, will skip them.

---

## Why FASTA has more records than BED and GTF

1,128 of the 4,814 annotated-start microproteins come from TrEMBL rather than the smORF
discovery pipeline, so their identifiers are gene names and UniProt accessions (`ALDOC`,
`J3QKK1`) instead of coordinate-encoded IDs. No genomic locus can be derived for them,
so they appear in the FASTA and ID list but not in the BED or GTF. The same applies
within `Ribo_ShortStop/` (1,592 sequences, 965 intervals).
