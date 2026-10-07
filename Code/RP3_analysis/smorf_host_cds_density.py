"""Ribo-seq density on each smORF relative to the main CDS of its host gene.

For every Salk and TrEMBL smORF on a protein-coding gene, compares ribosome
density on the smORF with density on the main CDS of its host gene, as a ratio
of RPKMs (P-sites per kb of CDS per million P-sites; no pseudocount), with
P-site counts from psite_frame_counts_per_orf.tsv. Two densities are written:
  frame 0    - each ORF's P-sites in its own reading frame (Psites_frame0_RPKM)
  all frames - all of each ORF's P-sites (Psites_total)
The call uses all frames for uORFs and dORFs that do not overlap the host CDS
(no host ribosome reads their codons), and frame 0 for every smORF that
overlaps or sits inside the host CDS, where all-frame counts would include the
host's own ribosomes.

smORF P-site row:
  Salk   - the GTF2FastaPatched row whose transcript_id is the master gene_id
  TrEMBL - the annotated CDS row matched by coordinates in
           Results/RP3/Psite_frame_by_sequence.csv; entries without one are
           skipped
Host CDS:
  1. the annotated CDS of the smORF's own transcript (the part of a Salk
     gene_id before "+chr"/"-chr"), when that transcript has one; else
  2. the MANE Select (else Ensembl canonical) CDS of a host gene listed for the
     smORF in smorf_transcript_flags.tsv, never the smORF's own row. Among
     several, prefer the gene named in the master, then one sharing an mRNA
     with the smORF, then the CDS overlapping (else nearest) the smORF, then
     the longest. A TrEMBL entry whose own CDS is some gene's reference CDS is
     hosted only by the gene the master names for it, else skipped.
Shared mRNA: the same transcript, or a compatible transcript in
smorf_transcript_flags.tsv (GENCODE v50, ENCODE4 or ESPRESSO long-read) that
contains both the smORF and the host start codon.

Call, first matching rule wins (P-site counts and ratio in the call's density):
  In-frame (not comparable) - in-frame isoforms (Iso, N-Iso, D-Iso) and TrEMBL
                              entries overlapping the host CDS: their P-sites
                              are largely the host's own
  Not on shared mRNA
  Low host coverage         - host CDS P-sites < MIN_HOST_PSITES
  ORF denser                - smORF P-sites >= MIN_ORF_PSITES, ratio > 1
  Low ORF coverage          - ratio > 1 on fewer smORF P-sites
  Host denser               - ratio <= 1
lncRNA, psORF and eORF smORFs are not compared (excluded by type) and get no
row. Densities are whole-ORF sums; P-site pile-up at start and stop codons is
not removed.

usage:
    python smorf_host_cds_density.py [--master ZIP|CSV] [--psites TSV] [--flags TSV]
                                     [--psite-map CSV] [--total-psites N] [--outdir DIR]
"""
import argparse
import re
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

repo = Path(__file__).resolve().parents[2]
MIN_HOST_PSITES = 20
MIN_ORF_PSITES = 10
UORF = {"uORF", "uaORF", "uoORF", "uaoORF"}
DORF = {"dORF", "daORF", "udORF", "daoORF", "doORF"}
ISO = {"Iso", "N-Iso", "D-Iso"}
SALK_TYPES = UORF | DORF | ISO | {"iORF"}
INFRAME_GROUPS = {"In-frame isoform", "TrEMBL overlapping CDS"}
ALL_FRAME_GROUPS = {"uORF", "dORF", "TrEMBL upstream", "TrEMBL downstream"}
COMPARABLE_CALLS = {"ORF denser", "Host denser", "Low ORF coverage"}
REF_LEVEL = {"MANE_Select": "MANE Select", "Ensembl_canonical": "Ensembl canonical"}
own_tx_re = re.compile(r"[+-]chr")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--master", type=Path, default=repo / "Code/data/microprotein_master.zip")
    p.add_argument("--psites", type=Path, default=repo / "Code/data/psite_frame_counts_per_orf.tsv")
    p.add_argument("--flags", type=Path, default=repo / "Code/data/smorf_transcript_flags.tsv",
                   help="one row per smORF x compatible transcript, with its host gene's reference CDS")
    p.add_argument("--psite-map", type=Path, default=repo / "Results/RP3/Psite_frame_by_sequence.csv",
                   help="written by Psite_frame_mapping.py; gives TrEMBL entries their P-site row")
    p.add_argument("--total-psites", type=int, default=116_287_730,
                   help="library P-sites behind psite_frame_counts_per_orf.tsv, for the all-frame RPKM")
    p.add_argument("--outdir", type=Path, default=repo / "Code/data")
    return p.parse_args()


def unversioned(tx):
    return tx.split(".")[0]


def read_master(path):
    cols = ["gene_id", "Database", "smorf_type", "gene_name", "sequence"]
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith(".csv"))
            m = pd.read_csv(z.open(name), usecols=cols, low_memory=False)
    else:
        m = pd.read_csv(path, usecols=cols, low_memory=False)
    m = m[m.Database.isin(["Salk", "TrEMBL"])].reset_index(drop=True)
    assert m.gene_id.is_unique
    m["subtype"] = np.where(m.Database == "TrEMBL", "TrEMBL", m.smorf_type)
    return m[(m.Database == "TrEMBL") | m.subtype.isin(SALK_TYPES)]


def host_candidates(path, annotated_u):
    """Host reference CDSs per smorf_id from the transcript flags, plus, per
    smorf_id, the transcripts (unversioned) whose own CDS is their gene's
    reference CDS (same length, wholly on the transcript)."""
    f = pd.read_csv(path, sep="\t", usecols=["smorf_id", "smorf_record_type", "tx_source", "transcript_id",
                                             "host_gene_name", "ref_tx", "ref_level", "ref_start_codon_in_tx",
                                             "ref_cds_len", "ref_cds_frac_in_tx", "tx_own_cds_len"])
    f = f[f.smorf_record_type.isin(["GTF2FastaPatched", "TrEMBL_recovered"])
          & (f.tx_source != "smORF_parent") & f.ref_tx.notna()].copy()
    f["ref_u"] = f.ref_tx.map(unversioned)
    f["start_in_tx"] = f.ref_start_codon_in_tx.astype(str).eq("True")
    is_ref = f[(f.ref_cds_frac_in_tx == 1) & (f.ref_cds_len == f.tx_own_cds_len)]
    ref_tx = {k: set(v.map(unversioned)) for k, v in is_ref.groupby("smorf_id").transcript_id}
    g = (f.groupby(["smorf_id", "ref_u"], sort=False)
         .agg(host_symbol=("host_gene_name", "first"), ref_level=("ref_level", "first"),
              shared=("start_in_tx", "any"))
         .reset_index())
    g = g[g.ref_u.isin(annotated_u)]
    return {k: v.to_dict("records") for k, v in g.groupby("smorf_id", sort=False)}, ref_tx


def position(s, e, hs, he, strand):
    """smORF CDS span (s, e) relative to the host CDS span (hs, he), strand-aware."""
    if s <= he and hs <= e:
        return "overlapping"
    return "upstream" if (e < hs) == (strand == "+") else "downstream"


def group_label(subtype, pos):
    if subtype == "TrEMBL":
        return "TrEMBL overlapping CDS" if pos == "overlapping" else f"TrEMBL {pos}"
    if subtype in ISO:
        return "In-frame isoform"
    if subtype == "iORF":
        return "iORF"
    family = "uORF" if subtype in UORF else "dORF"
    return f"{family} overlapping CDS" if pos == "overlapping" else family


def pick_host(cands, own, own_u, own_is_ref, gene_name, annotated):
    """Best host CDS row for one smORF, or (None, reason)."""
    if own_is_ref or any(c["ref_u"] in own_u for c in cands):
        # The smORF's own row is some gene's reference CDS: only its named gene can host it.
        cands = [c for c in cands if c["ref_u"] not in own_u and c["host_symbol"] == gene_name]
        if not cands:
            return None, "own row is the host CDS"
    if not cands:
        return None, "no host CDS in transcript flags"
    s, e = own.cds_start, own.cds_end

    def rank(c):
        h = annotated.loc[c["ref_u"]]
        overlap = max(0, min(e, h.cds_end) - max(s, h.cds_start) + 1)
        distance = max(0, h.cds_start - e, s - h.cds_end)
        return (c["host_symbol"] != gene_name, not c["shared"], -overlap, distance, -h.cds_len_nt, c["ref_u"])

    return min(cands, key=rank), None


def main():
    args = parse_args()
    master = read_master(args.master)

    ps = pd.read_csv(args.psites, sep="\t", low_memory=False)
    assert ps.transcript_id.is_unique
    smorf_rows = ps[ps.source == "GTF2FastaPatched"].set_index("transcript_id")
    annotated = ps[~ps.is_smorf & ~ps.transcript_id.str.endswith("_PAR_Y")].copy()
    assert annotated.source.isin(["HAVANA", "ENSEMBL"]).all()
    annotated["tx_u"] = annotated.transcript_id.map(unversioned)
    assert annotated.tx_u.is_unique
    par_u = set(ps.loc[ps.transcript_id.str.endswith("_PAR_Y"), "transcript_id"].map(unversioned))
    by_tx = annotated.set_index("transcript_id")
    by_u = annotated.set_index("tx_u")

    pmap = pd.read_csv(args.psite_map, usecols=["sequence", "Psite_match", "Psite_transcript_id"])
    pmap = pmap[pmap.Psite_match == "coordinates"]
    trembl_ids = {seq: sorted(ids.split(";")) for seq, ids in zip(pmap.sequence, pmap.Psite_transcript_id)}

    flags, ref_tx = host_candidates(args.flags, set(by_u.index))
    rpkm = lambda n, nt: n * 1e9 / (nt * args.total_psites)

    rows, skipped = [], Counter()
    for m in master.itertuples(index=False):
        seq = str(m.sequence).rstrip("*")
        if m.subtype == "TrEMBL":
            ids = trembl_ids.get(seq)
            if not ids:
                skipped["TrEMBL: no own P-site row"] += 1
                continue
            own = by_tx.loc[ids[0]]
            own_tx, own_u = ids[0], {unversioned(t) for t in ids}
        else:
            own = smorf_rows.loc[m.gene_id]
            own_tx, own_u = m.gene_id, set()
        assert own.n_codons == len(seq), m.gene_id

        same_tx = own_tx_re.split(m.gene_id, maxsplit=1)[0] if m.subtype != "TrEMBL" else None
        if same_tx in by_tx.index:
            host, host_tx = by_tx.loc[same_tx], same_tx
            choice, symbol, shared = "same transcript", host.gene_name, True
        else:
            own_is_ref = bool(own_u & ref_tx.get(m.gene_id, set()))
            c, why = pick_host(flags.get(m.gene_id, []), own, own_u, own_is_ref, m.gene_name, by_u)
            if c is None:
                skipped[f"{m.subtype if m.subtype == 'TrEMBL' else 'Salk'}: {why}"] += 1
                continue
            host = by_u.loc[c["ref_u"]]
            host_tx = host.transcript_id
            choice, symbol, shared = REF_LEVEL[c["ref_level"]], c["host_symbol"], bool(c["shared"])
        assert (own.chrom, own.strand) == (host.chrom, host.strand), m.gene_id
        assert unversioned(host_tx) not in par_u and unversioned(own_tx) not in par_u, m.gene_id

        pos = position(own.cds_start, own.cds_end, host.cds_start, host.cds_end, own.strand)
        rows.append({
            "gene_id": m.gene_id, "source": m.Database, "subtype": m.subtype,
            "group": group_label(m.subtype, pos), "position": pos, "shared_mrna": shared,
            "host_choice": choice, "host_symbol": symbol, "host_tx": host_tx, "orf_tx": own_tx,
            "orf_codons": int(own.n_codons),
            "orf_psites_f0": int(own.Psites_sum_frame0), "orf_rpkm_f0": float(own.Psites_frame0_RPKM),
            "orf_psites_all": int(own.Psites_total), "orf_rpkm_all": rpkm(own.Psites_total, own.cds_len_nt),
            "host_codons": int(host.n_codons),
            "host_psites_f0": int(host.Psites_sum_frame0), "host_rpkm_f0": float(host.Psites_frame0_RPKM),
            "host_psites_all": int(host.Psites_total), "host_rpkm_all": rpkm(host.Psites_total, host.cds_len_nt),
        })

    out = pd.DataFrame(rows)
    for d in ("f0", "all"):
        out[f"rpkm_ratio_{d}"] = out[f"orf_rpkm_{d}"] / out[f"host_rpkm_{d}"].where(out[f"host_rpkm_{d}"] > 0)
    use_all = out.group.isin(ALL_FRAME_GROUPS)
    out["density_basis"] = np.where(use_all, "all frames", "frame 0")
    ratio = out.rpkm_ratio_all.where(use_all, out.rpkm_ratio_f0)
    host_n = out.host_psites_all.where(use_all, out.host_psites_f0)
    orf_n = out.orf_psites_all.where(use_all, out.orf_psites_f0)
    out["rpkm_ratio"] = ratio
    out["log2_rpkm_ratio"] = np.log2(ratio.where(ratio > 0)).round(4)
    conds = [out.group.isin(INFRAME_GROUPS), ~out.shared_mrna, host_n < MIN_HOST_PSITES,
             (orf_n >= MIN_ORF_PSITES) & (ratio > 1), ratio > 1, ratio <= 1]
    calls = ["In-frame (not comparable)", "Not on shared mRNA", "Low host coverage",
             "ORF denser", "Low ORF coverage", "Host denser"]
    out["call"] = np.select(conds, calls, default="")
    assert (out.call != "").all()
    out["comparable"] = out.call.isin(COMPARABLE_CALLS)
    out["orf_denser"] = out.call.eq("ORF denser")
    assert (out.comparable == (out.shared_mrna & (host_n >= MIN_HOST_PSITES)
                               & ~out.group.isin(INFRAME_GROUPS))).all()
    assert (out.orf_denser == (out.comparable & (orf_n >= MIN_ORF_PSITES) & (ratio > 1))).all()

    args.outdir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outdir / "smorf_host_cds_density.csv", index=False)

    rpkm_def = "P-sites per kb of CDS per million P-sites"
    pd.DataFrame([
        ("gene_id", "smORF ID (master gene_id; UniProt accession for TrEMBL)"),
        ("source", "Salk (novel smORF) or TrEMBL (unreviewed UniProt microprotein)"),
        ("subtype", "Master smorf_type (TrEMBL for TrEMBL entries)"),
        ("group", "Comparison class: the master smorf_type family (uORF, dORF, iORF, In-frame isoform = "
                  "Iso/N-Iso/D-Iso), with uORF and dORF split by overlap with the host CDS (uORF / uORF "
                  "overlapping CDS, dORF / dORF overlapping CDS); TrEMBL entries by position (TrEMBL upstream / "
                  "overlapping CDS / downstream). Only TrEMBL labels follow position, so see position for the rest"),
        ("position", "smORF CDS relative to the host CDS on the genome, strand-aware: upstream, overlapping or downstream"),
        ("shared_mrna", "True when the smORF and the host start codon lie on one transcript: the same transcript, or a "
                        "compatible GENCODE v50, ENCODE4 or ESPRESSO long-read model in smorf_transcript_flags.tsv"),
        ("host_choice", "How the host CDS was chosen: same transcript (the smORF transcript's own annotated CDS), "
                        "or the host gene's MANE Select / Ensembl canonical CDS"),
        ("host_symbol", "Gene symbol of the host CDS"),
        ("host_tx", "Ensembl transcript of the host CDS (psite_frame_counts_per_orf.tsv transcript_id)"),
        ("orf_tx", "P-site row used for the smORF (psite_frame_counts_per_orf.tsv transcript_id)"),
        ("orf_codons", "smORF CDS codons, stop codon excluded (= protein length)"),
        ("orf_psites_f0", "Adult-brain P-sites in frame 0 of the smORF"),
        ("orf_rpkm_f0", f"smORF frame-0 {rpkm_def} (Psites_frame0_RPKM)"),
        ("orf_psites_all", "Adult-brain P-sites on the smORF CDS, all frames (Psites_total)"),
        ("orf_rpkm_all", f"smORF all-frame {rpkm_def}"),
        ("host_codons", "Host CDS codons, stop codon excluded"),
        ("host_psites_f0", "Adult-brain P-sites in frame 0 of the host CDS"),
        ("host_rpkm_f0", f"Host CDS frame-0 {rpkm_def} (Psites_frame0_RPKM)"),
        ("host_psites_all", "Adult-brain P-sites on the host CDS, all frames (Psites_total)"),
        ("host_rpkm_all", f"Host CDS all-frame {rpkm_def}"),
        ("rpkm_ratio_f0", "orf_rpkm_f0 / host_rpkm_f0, no pseudocount (empty when the host has no frame-0 P-sites)"),
        ("rpkm_ratio_all", "orf_rpkm_all / host_rpkm_all, no pseudocount (empty when the host has no P-sites)"),
        ("density_basis", "Density the call uses: all frames for groups " + ", ".join(sorted(ALL_FRAME_GROUPS)) +
                          " (no overlap with the host CDS); frame 0 otherwise, so host-frame ribosomes on shared "
                          "codons are not counted for the smORF"),
        ("rpkm_ratio", "The ratio the call uses (rpkm_ratio_all or rpkm_ratio_f0, per density_basis); whole-ORF "
                       "sums, P-site pile-up at start and stop codons not removed"),
        ("log2_rpkm_ratio", "log2 of rpkm_ratio (empty when rpkm_ratio is empty or 0)"),
        ("call", "First matching rule, with P-site counts in the density_basis: In-frame (not comparable) | "
                 f"Not on shared mRNA | Low host coverage (host P-sites < {MIN_HOST_PSITES}) | "
                 f"ORF denser (smORF P-sites >= {MIN_ORF_PSITES} and rpkm_ratio > 1) | "
                 "Low ORF coverage (rpkm_ratio > 1 on fewer P-sites) | Host denser (rpkm_ratio <= 1)"),
        ("comparable", "call is ORF denser, Host denser or Low ORF coverage"),
        ("orf_denser", "call is ORF denser"),
    ], columns=["column", "meaning"]).to_csv(args.outdir / "smorf_host_cds_data_dictionary.csv", index=False)

    print(f"done: {len(out)} smORFs {out.source.value_counts().to_dict()}; "
          f"skipped {dict(sorted(skipped.items()))}; calls {out.call.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
