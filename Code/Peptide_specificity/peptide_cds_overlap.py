"""Do a Salk smORF's MS peptides lie on genomic sequence that also encodes an annotated CDS?

Peptide specificity (make_peptide_flags.py) asks whether a peptide's *sequence*
occurs in another protein. A peptide can be unique by sequence and still be
read from the same nucleotides as an annotated CDS in another frame, e.g. the
part of a uoORF that overlaps the main ORF. This script answers the genomic
question, in any frame.

Peptides: the spliced genomic blocks of every detected tryptic peptide
  (Code/data/browser_tracks/smorf_psm_genomic.bed, BED12, name smorf_id|PEPTIDE).
Annotated CDS: the merged annotated (non-smORF) CDS on the smORF's own strand,
  within its span +/-150 nt (`annotated_cds` in Code/data/smorf_flank_structure.csv,
  written by Code/RP3_analysis/smorf_flank_psites.py from GENCODE v43).
A peptide overlaps when at least one of its nucleotides falls in that CDS.

Salk smORFs only. A TrEMBL entry's own CDS is itself an annotated Ensembl CDS,
so every TrEMBL peptide would overlap trivially; TrEMBL entries get no row.
Opposite-strand CDS are not considered.

usage:
    python peptide_cds_overlap.py [--peptides BED] [--structure CSV] [--outdir DIR]
"""
import argparse
import re
from pathlib import Path

import pandas as pd

repo = Path(__file__).resolve().parents[2]
salk_id_re = re.compile(r":\d+-\d+_F:")  # Salk gene_id; TrEMBL ids are UniProt accessions


def parse_blocks(s):
    if not isinstance(s, str) or not s:
        return []
    return [tuple(map(int, b.split("-"))) for b in s.split(";")]


def overlap_nt(blocks, ivs):
    return sum(max(0, min(b, e) - max(a, s)) for a, b in blocks for s, e in ivs)


def read_peptides(bed):
    rows = []
    with open(bed) as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            s = int(f[1])
            blocks = [(s + int(o), s + int(o) + int(z)) for z, o in
                      zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(","))]
            gene_id, peptide = f[3].split("|", 1)
            rows.append((gene_id, peptide, blocks))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--peptides", type=Path,
                   default=repo / "Code/data/browser_tracks/smorf_psm_genomic.bed")
    p.add_argument("--structure", type=Path, default=repo / "Code/data/smorf_flank_structure.csv")
    p.add_argument("--outdir", type=Path, default=repo / "Code/data")
    args = p.parse_args()

    st = pd.read_csv(args.structure, usecols=["gene_id", "annotated_cds"])
    ann = dict(zip(st.gene_id, st.annotated_cds.map(parse_blocks)))

    recs = []
    for gene_id, peptide, blocks in read_peptides(args.peptides):
        if not salk_id_re.search(gene_id) or gene_id not in ann:
            continue
        nt = sum(b - a for a, b in blocks)
        ov = overlap_nt(blocks, ann[gene_id])
        recs.append((gene_id, peptide, nt, ov))
    pep = pd.DataFrame(recs, columns=["gene_id", "peptide", "peptide_nt", "overlap_nt"])
    pep["overlaps"] = pep.overlap_nt.gt(0)
    pep["frac"] = pep.overlap_nt / pep.peptide_nt
    pep["label"] = (pep.peptide + "(" + pep.overlap_nt.astype(str) + "/"
                    + pep.peptide_nt.astype(str) + " nt)").where(pep.overlaps)

    out = pep.groupby("gene_id", sort=True).agg(
        n_peptides=("peptide", "size"),
        n_overlapping=("overlaps", "sum"),
        max_overlap_frac=("frac", "max"),
        overlapping_peptides=("label", lambda s: ";".join(s.dropna())),
    ).reset_index()
    out["max_overlap_frac"] = out.max_overlap_frac.round(3)
    # Not "None": pandas reads that back as NaN.
    out["status"] = "No peptides"
    out.loc[out.n_overlapping.gt(0), "status"] = "Some peptides"
    out.loc[out.n_overlapping.eq(out.n_peptides), "status"] = "All peptides"
    out = out[["gene_id", "status", "n_peptides", "n_overlapping", "max_overlap_frac",
               "overlapping_peptides"]]
    out.to_csv(args.outdir / "peptide_cds_overlap.csv", index=False)

    pd.DataFrame([
        ("gene_id", "smORF ID (master gene_id); Salk smORFs only"),
        ("status", "All / Some / No peptides: how many of the smORF's located peptides overlap an annotated CDS"),
        ("n_peptides", "Detected tryptic peptides placed on the genome"),
        ("n_overlapping", "Peptides with >=1 nt in an annotated same-strand CDS (any frame)"),
        ("max_overlap_frac", "Largest fraction of one peptide's nucleotides in annotated CDS"),
        ("overlapping_peptides", "PEPTIDE(overlap nt/peptide nt), ';'-separated"),
    ], columns=["column", "description"]).to_csv(
        args.outdir / "peptide_cds_overlap_data_dictionary.csv", index=False)

    print(f"{len(pep):,} peptides on {len(out):,} Salk smORFs; "
          f"{pep.overlaps.sum():,} peptides overlap annotated CDS")
    print(out.status.value_counts().to_string())


if __name__ == "__main__":
    main()
