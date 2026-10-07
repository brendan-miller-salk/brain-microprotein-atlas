"""Gene at the CDS locus of each TrEMBL microprotein, for its parent-gene name.

The master table names TrEMBL entries by their UniProt gene name, which can be
a placeholder (LOC..., the accession itself), outdated, or a paralogue of the
gene the CDS actually sits on (e.g. CALM1 for a CDS on CALM2). This picks the
gene at the locus instead, from the host genes listed for the entry in
smorf_transcript_flags.tsv (GENCODE v50), first matching rule wins:

  master name is a locus gene - the master name is one of the host genes: kept
  own gene                    - the host gene whose reference CDS is the
                                protein (CDS length = 3 x protein length, +3
                                with stop), or the Ensembl gene of the CDS the
                                entry is matched to in Psite_frame_by_sequence
                                .csv, when it has a symbol
  UniProt name                - the current UniProt gene name, when it is a
                                symbol (not LOC..., hCG_..., or an accession)
  named locus gene            - the only host gene with a symbol
  master name kept            - any other master name that is not the
                                accession itself (LOC..., hCG_...), when no
                                host gene has a symbol
  unnamed locus gene          - the host gene's Ensembl ID, when the master
                                name is the accession itself
Entries with no host genes in the flags keep the master name ("no locus data").

usage:
    python trembl_locus_genes.py [--master ZIP|CSV] [--flags TSV] [--uniprot TSV]
                                 [--psite-map CSV] [--psites TSV] [--outdir DIR]
"""
import argparse
import re
import zipfile
from pathlib import Path

import pandas as pd

repo = Path(__file__).resolve().parents[2]
UNNAMED = re.compile(r"(ENSG\d|LOC\d|hCG_|Unnamed$)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--master", type=Path, default=repo / "Code/data/microprotein_master.zip")
    p.add_argument("--flags", type=Path, default=repo / "Code/data/smorf_transcript_flags.tsv")
    p.add_argument("--uniprot", type=Path, default=repo / "Code/data/uniprotkb_proteome_UP000005640_2026_07_13.tsv")
    p.add_argument("--psite-map", type=Path, default=repo / "Results/RP3/Psite_frame_by_sequence.csv",
                   help="written by RP3_analysis/Psite_frame_mapping.py (tracked in git)")
    p.add_argument("--psites", type=Path, default=repo / "Code/data/psite_frame_counts_per_orf.tsv")
    p.add_argument("--outdir", type=Path, default=repo / "Code/data")
    return p.parse_args()


def is_symbol(g):
    return isinstance(g, str) and g != "" and not UNNAMED.match(g)


def read_master(path):
    cols = ["gene_id", "Database", "gene_name", "sequence"]
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".csv"))
            m = pd.read_csv(z.open(name), usecols=cols, low_memory=False)
    else:
        m = pd.read_csv(path, usecols=cols, low_memory=False)
    m = m[m.Database == "TrEMBL"].reset_index(drop=True)
    assert m.gene_id.is_unique
    return m


def locus_genes(path):
    """accession -> host genes [{id, name, ref_cds_len}] from the transcript flags."""
    f = pd.read_csv(path, sep="\t", usecols=["smorf_id", "smorf_record_type", "tx_source",
                                             "host_gene_id", "host_gene_name", "ref_cds_len"])
    f = f[(f.smorf_record_type == "TrEMBL_recovered") & (f.tx_source != "smORF_parent")]
    g = (f.groupby(["smorf_id", "host_gene_id"], sort=False)
         .agg(name=("host_gene_name", "first"), ref_cds_len=("ref_cds_len", "first")).reset_index())
    return {k: v.to_dict("records") for k, v in g.groupby("smorf_id", sort=False)}


def pick(accession, master_name, length, hosts, uniprot_name, psite_gene):
    names = [h["name"] for h in hosts]
    if not hosts:
        return master_name, "no locus data"
    if master_name in names:
        return master_name, "master name is a locus gene"
    own = {h["name"] for h in hosts if h["ref_cds_len"] in (3 * length, 3 * length + 3)}
    if psite_gene in names:
        own.add(psite_gene)
    own = sorted(n for n in own if is_symbol(n))
    named = sorted({n for n in names if is_symbol(n)})
    if len(own) == 1:
        return own[0], "own gene"
    if is_symbol(uniprot_name):
        return uniprot_name, "UniProt name"
    if len(named) == 1:
        return named[0], "named locus gene"
    assert not named, (master_name, named)  # several named genes and no tie-break: none today
    if master_name != accession:  # a LOC.../hCG_... name (a symbol would be a host gene or own gene)
        return master_name, "master name kept"
    return sorted(h["host_gene_id"] for h in hosts)[0].split(".")[0], "unnamed locus gene"


def main():
    args = parse_args()
    master = read_master(args.master)
    hosts = locus_genes(args.flags)
    u = pd.read_csv(args.uniprot, sep="\t", usecols=["Entry", "Gene Names"])
    uniprot = dict(zip(u.Entry, u["Gene Names"].fillna("").str.split().str[0]))
    pmap = pd.read_csv(args.psite_map, usecols=["sequence", "Psite_match", "Psite_transcript_id"])
    pmap = pmap[pmap.Psite_match == "coordinates"]
    psite_gene = pd.read_csv(args.psites, sep="\t", usecols=["transcript_id", "gene_name"]).set_index("transcript_id").gene_name
    own_gene = {s: psite_gene.get(sorted(ids.split(";"))[0]) for s, ids in zip(pmap.sequence, pmap.Psite_transcript_id)}

    rows = []
    for m in master.itertuples(index=False):
        seq = str(m.sequence).rstrip("*")
        hs = hosts.get(m.gene_id, [])
        locus, rule = pick(m.gene_id, m.gene_name, len(seq), hs, uniprot.get(m.gene_id), own_gene.get(seq))
        rows.append({"gene_id": m.gene_id, "master_name": m.gene_name, "locus_gene": locus, "rule": rule,
                     "host_genes": ";".join(h["name"] for h in hs), "uniprot_name": uniprot.get(m.gene_id)})
    out = pd.DataFrame(rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outdir / "trembl_locus_genes.csv", index=False)
    changed = out[out.locus_gene != out.master_name]
    print(f"done: {len(out)} TrEMBL entries; renamed {len(changed)}; rules {out.rule.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
