# `browser_tracks/` — UCSC genome-browser tracks

Six tracks over the human frontal cortex (DLPFC), hg38: ribosome profiling,
ESPRESSO long-read isoform coverage, and smORF peptide / ORF-rules annotation.

## Loading them

UCSC Genome Browser → **My Data → Custom Tracks**, select **hg38**, and paste
the contents of [`brain_tracks.txt`](brain_tracks.txt), which sits in this
directory.

The track data streams from the public dataset
[`brmiller/brain-microprotein-atlas`](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas);
nothing needs to be downloaded first.

## The files

| File | Size | Track |
|------|-----:|-------|
| `adult_brain_riboseq_merged.fwd.bw` | 22 M | DLPFC Ribo (Fw) |
| `adult_brain_riboseq_merged.rev.bw` | 22 M | DLPFC Ribo (Rev) |
| `espresso_dlpfc_mean_cpm.bw` | 3.3 M | DLPFC ESPRESSO (CPM) |
| `espresso_dlpfc_isoforms.bb` | 3.1 M | DLPFC ESPRESSO Isoforms |
| `smorf_psm_peptides.bb` | 476 K | smORF PROSIT Peptides |
| `smorf_trx_priority_color_class.bb` | 628 K | smORF Rules |

Everything else in the directory is source or audit material for those six:

| File | What it is |
|------|------------|
| `smorf_trx_priority_color_class.bed` | source BED behind the smORF Rules track (5,668 features) |
| `smorf_trx_priority_color_class_assignments.tsv` | per-smORF rule call for the 5,621 smORFs with data — chosen transcript, decay class, disruption tier, and `nsd_gate_reason` |
| `smorf_psm_genomic.tsv` / `.bed` / `smorf_psm_peptides.gtf` | the peptide track's source at each stage (3,854 peptides over 3,154 smORFs) |
| `trembl_cds_phase_repair.tsv` | audits the 155 TrEMBL records whose CDS boundary was re-phased against the genome (see below) |
| `trembl_cds_rephase.tsv` | the 120 TrEMBL records whose CDS *frame* was corrected when the Rules track was built — a separate pass from the one above |

Colour schemes, scale choices, and per-tier feature counts are documented in the
comments of `brain_tracks.txt` itself.

## Re-hosting

UCSC fetches `bigDataUrl` byte-ranges **server-side** — genome.ucsc.edu, not
your machine, issues the request. So these files must sit at a public HTTP URL
that supports range requests; a local path or `file://` URL cannot work.

Hosts that work: Hugging Face, Zenodo, S3/GCS with public read, or any ordinary
web server. GitHub raw URLs do **not** — they do not honour range requests
reliably.

After moving the files, rewrite the six `bigDataUrl=` values in
`brain_tracks.txt` and check each one:

```bash
curl -s -o /dev/null -w '%{http_code}\n' -H 'Range: bytes=0-99' -L '<url>'
```

Want `206`. A `200` means ranges were ignored and UCSC will fail to load the
track. The `-L` matters — Hugging Face redirects to a CDN, and the `Range`
header has to survive the redirect, so a bare `curl -I` misleadingly reports
`200` on URLs that are actually fine.

## TrEMBL CDS repair (2026-08-28)

Two independent passes touched TrEMBL CDS boundaries. They are easy to confuse:
the first fixed *where* the CDS sat so peptides could be placed, the second fixed
the reading *frame* so UCSC could translate it.

### Pass 1 — thick boundary, for the peptide track

155 TrEMBL-derived smORFs carried peptides that could not be placed: their
annotated CDS did not encode the peptide in any position. The exon chains and
strands were right; only the thick (CDS) boundary was misplaced, which is what
you would expect from records that originate as blastp envelopes mapped onto a
transcript rather than as called ORFs. Two things identified it — for 132 of the
155 the CDS length exactly equalled the declared protein length, and translating
the full transcript showed 149 had one reading frame containing every one of
their peptides.

`transcript_rules/fix_trembl_cds_phase.py` (Box) re-phases those 149: it picks
the frame holding all the record's peptides and sets thick to a real ORF spanning
them (back to the nearest in-frame ATG, forward to the first stop). All 149 come
out as clean ORFs — length divisible by 3, no internal stop.

Effect on the peptide track: 3,670 -> 3,854 features over 3,154 smORFs. Eight
eligible peptides remain unplaced, all on smORFs with no CDS at all. Six records
were deliberately left alone (`A0A096LP52`, `A0A109PLW1`, `B4DQX9`, `F8VRZ4`,
`H3BPV0`, `J3QQY7`) because no single frame holds all their peptides — those need
the TrEMBL recovery revisited, not a boundary nudge.

**Every one of the 3,854 features was translated from hg38 and matches its
expected peptide exactly — 100%, exhaustive rather than sampled.**

### Pass 2 — reading frame, for the smORF Rules track

Building the Rules track re-checked every CDS frame against the genome and found
120 TrEMBL records that were out of frame: each had at least one internal stop
in its annotated frame, none ended in a stop codon, and all ran off the 3' end
of their transcript. Shifting the start by 1 nt (51 records) or 2 nt (69) cleared
them. The corrections are itemised in `trembl_cds_rephase.tsv`, one row per
record with the old and new CDS blocks.

After this pass all 5,621 coding features in the Rules track translate to a clean
ORF with no internal stop.
