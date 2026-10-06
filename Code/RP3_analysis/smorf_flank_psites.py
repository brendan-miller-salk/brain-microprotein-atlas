"""Ribo-seq P-sites in the windows flanking each smORF, with reading frame.

For every smORF, counts adult-brain P-sites (sum of the three frame tracks,
same strand) in the CDS and in 50/100/250/1000-nt windows immediately upstream
(5') and downstream (3'). Flanks walk along the spliced host transcript and are
clipped at its ends, so a window can be shorter than its nominal size.

smORFs:
  Salk   - novel smORFs; host transcript + CDS from the microprotein GTF
  TrEMBL - unreviewed UniProt microproteins; host transcript (BED12 blocks) and
           frame-repaired CDS (thick) from the smORF Rules track BED, whose
           names are "ACC|host"; records with an empty thick are skipped
Only gene_ids present in the master table are kept (this drops the
*_alt_initiation_N proteoforms). The stop codon always counts as CDS: where the
CDS coordinates stop short of it (CDS nt == 3 x protein length, as in the
smORF GTF), the first 3 nt of the 3' flank are moved into the CDS, unless the
transcript ends first.

Frame is relative to the smORF start codon in transcript coordinates (frame 0 =
in frame with the smORF). Frame percentages use "clean" positions only: not
overlapping an annotated GENCODE CDS on the same strand, and not among the 6 nt
just upstream of the start codon, where the initiation P-site peak spills over.

usage:
    python smorf_flank_psites.py [--total-psites N] [--bigwig-prefix PREFIX]
                                 [--full-gtf GTF] [--outdir DIR]
"""
import argparse
import bisect
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

repo = Path(__file__).resolve().parents[2]
WINDOWS = [50, 100, 250, 1000]
MAXW = max(WINDOWS)
START_SPILLOVER_NT = 6
tid_re = re.compile(r'transcript_id "([^"]+)"')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--smorf-gtf", type=Path,
                   default=repo / "GTF_and_BED_files/Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf")
    p.add_argument("--rules-bed", type=Path,
                   default=repo / "Code/data/browser_tracks/smorf_trx_priority_color_class.bed",
                   help="smORF Rules track BED12; source of TrEMBL host transcripts and CDS")
    p.add_argument("--full-gtf", type=Path,
                   default=repo / "GTF_and_BED_files/Ensembl_and_Unreviewed_Brain_Microproteins.gtf",
                   help="GENCODE + smORF GTF; non-smORF CDS lines mark annotated CDS")
    p.add_argument("--master", type=Path, default=repo / "Code/data/microprotein_master.zip")
    p.add_argument("--bigwig-prefix", type=Path,
                   default=repo / "Code/data/browser_tracks/adult_psite",
                   help="P-site bigwigs named <prefix>.f{0,1,2}.{fwd,rev}.bw, in CPM")
    p.add_argument("--total-psites", type=int, default=116_287_730,
                   help="total P-sites behind the bigwigs, to convert CPM back to counts")
    p.add_argument("--outdir", type=Path, default=repo / "Code/data")
    return p.parse_args()


def read_smorfs(smorf_gtf, rules_bed):
    recs = []  # (gene_id, source, chrom, strand, exons, cds), 0-based half-open
    tx = defaultdict(lambda: {"exon": [], "CDS": []})
    with open(smorf_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] not in ("exon", "CDS"):
                continue
            d = tx[tid_re.search(f[8]).group(1)]
            d["chrom"], d["strand"] = f[0], f[6]
            d[f[2]].append((int(f[3]) - 1, int(f[4])))
    for t, d in tx.items():
        recs.append((t, "Salk", d["chrom"], d["strand"], sorted(d["exon"]), sorted(d["CDS"])))
    with open(rules_bed) as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            f = line.rstrip("\n").split("\t")
            if "|" not in f[3] or f[3].startswith(("ENST", "ESPRESSO")):
                continue
            s, ts, te = int(f[1]), int(f[6]), int(f[7])
            if ts == te:
                continue
            exons = [(s + int(o), s + int(o) + int(z)) for z, o in
                     zip(f[10].rstrip(",").split(","), f[11].rstrip(",").split(","))]
            cds = [(max(a, ts), min(b, te)) for a, b in exons if min(b, te) > max(a, ts)]
            recs.append((f[3].split("|")[0], "TrEMBL", f[0], f[5], exons, cds))
    return recs


def read_annotated_cds(full_gtf):
    """Merged annotated (non-smORF) CDS intervals per (chrom, strand)."""
    ann = defaultdict(list)
    with open(full_gtf) as fh:
        for line in fh:
            if "\tCDS\t" not in line:
                continue
            f = line.split("\t", 8)
            if f[1] in ("GTF2FastaPatched", "AltInitiation"):
                continue
            ann[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    merged = {}
    for k, iv in ann.items():
        iv.sort()
        m = []
        for s, e in iv:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[k] = (np.array([a for a, _ in m]), np.array([b for _, b in m]))
    return merged


def overlapping(merged, chrom, strand, s, e):
    st, en = merged.get((chrom, strand), (np.array([], int), np.array([], int)))
    out, i = [], max(bisect.bisect_right(st, s) - 1, 0)
    while i < len(st) and st[i] < e:
        a, b = max(int(st[i]), s), min(int(en[i]), e)
        if a < b:
            out.append((a, b))
        i += 1
    return out


def tx_flank_segs(exons, cds, strand, maxw=MAXW):
    """Genomic segments of up to maxw nt of spliced transcript 5' and 3' of the CDS."""
    exons = sorted(exons)
    cs, ce = min(s for s, _ in cds), max(e for _, e in cds)
    left, right, need = [], [], maxw
    for s, e in reversed(exons):
        if s >= cs or need == 0:
            continue
        b = min(e, cs)
        a = max(s, b - need)
        left.insert(0, (a, b))
        need -= b - a
    need = maxw
    for s, e in exons:
        if e <= ce or need == 0:
            continue
        a = max(s, ce)
        b = min(e, a + need)
        right.append((a, b))
        need -= b - a
    return (left, right) if strand == "+" else (right, left)


def oriented(values_fn, mask_ivs, strand, segs):
    """Per-nt P-sites and annotated-CDS mask over segs, in 5'->3' transcript order."""
    vals, msk = [], []
    for s, e in segs:
        vals.append(values_fn(s, e))
        m = np.zeros(e - s, bool)
        for a, b in mask_ivs:
            if a < e and b > s:
                m[max(a, s) - s:min(b, e) - s] = True
        msk.append(m)
    if strand == "-":
        vals = [v[::-1] for v in vals][::-1]
        msk = [m[::-1] for m in msk][::-1]
    cat = lambda xs: np.concatenate(xs) if xs else np.zeros(0)
    return cat(vals), cat(msk).astype(bool)


def fmt_blocks(blocks):
    return ";".join(f"{a}-{b}" for a, b in blocks)


def main():
    args = parse_args()
    unit = 1e6 / args.total_psites  # CPM value of one P-site
    bws = {st: [pyBigWig.open(f"{args.bigwig_prefix}.f{k}.{st}.bw") for k in range(3)]
           for st in ("fwd", "rev")}

    with zipfile.ZipFile(args.master) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        meta = pd.read_csv(z.open(name), usecols=["gene_id", "Database", "sequence"],
                           low_memory=False)
    meta = meta[meta.Database.isin(["Salk", "TrEMBL"])].drop_duplicates("gene_id")
    aa_len = dict(zip(meta.gene_id, meta.sequence.astype(str).str.rstrip("*").str.len()))
    keep_ids = set(aa_len)

    recs = [r for r in read_smorfs(args.smorf_gtf, args.rules_bed) if r[0] in keep_ids]
    merged = read_annotated_cds(args.full_gtf)

    rows, struct = [], []
    for gid, src, c, st, exons, cds in recs:
        files = bws["fwd" if st == "+" else "rev"]
        values = lambda s, e: np.round(sum(np.nan_to_num(np.array(b.values(c, s, e)))
                                           for b in files) / unit) if e > s else np.zeros(0)
        stop_excluded = sum(b - a for a, b in cds) == 3 * aa_len[gid]
        up, dn = tx_flank_segs(exons, cds, st, maxw=MAXW + 3 * stop_excluded)
        span = [x for seg in (up, cds, dn) for x in seg]
        ann = overlapping(merged, c, st, min(a for a, _ in span), max(b for _, b in span))
        cds_v, _ = oriented(values, [], st, cds)
        uv, um = oriented(values, ann, st, up)
        dv, dm = oriented(values, ann, st, dn)
        stop_excluded = stop_excluded and len(dv) >= 3  # CDS running off the 3' end: no stop
        if stop_excluded:
            cds_v = np.concatenate([cds_v, dv[:3]])
            dv, dm = dv[3:], dm[3:]
        uv, um = uv[::-1], um[::-1]  # upstream: index 0 = nt adjacent to the start codon
        L = len(cds_v)
        cds_density = cds_v.sum() / L
        r = {"gene_id": gid, "source": src, "cds_len": L, "cds_psites": int(cds_v.sum()),
             "stop_moved_into_cds": stop_excluded}
        up_frame = (-(np.arange(len(uv)) + 1)) % 3
        dn_frame = (L + np.arange(len(dv))) % 3
        up_clean = ~um & (np.arange(len(uv)) >= START_SPILLOVER_NT)
        for side, v, m, fr, clean in (("up", uv, um, up_frame, up_clean),
                                      ("down", dv, dm, dn_frame, ~dm)):
            for w in WINDOWS:
                p = f"{side}{w}"
                n = min(len(v), w)
                ps = v[:w].sum()
                r[f"{p}_len"] = n
                r[f"{p}_psites"] = int(ps)
                r[f"{p}_density_over_cds"] = (ps / n) / cds_density if n and cds_density else np.nan
                r[f"{p}_annCDS_bp"] = int(m[:w].sum())
                fc = [v[:w][clean[:w] & (fr[:w] == k)].sum() for k in range(3)]
                tot = sum(fc)
                r[f"{p}_clean_n"] = int(tot)
                for k in range(3):
                    r[f"{p}_pct_f{k}"] = 100 * fc[k] / tot if tot else np.nan
        rows.append(r)
        struct.append({"gene_id": gid, "chrom": c, "strand": st,
                       "stop_moved_into_cds": stop_excluded,
                       "exons": fmt_blocks(exons), "cds": fmt_blocks(cds),
                       "annotated_cds": fmt_blocks(ann)})

    out = pd.DataFrame(rows)
    for col in out.columns:
        if col.endswith(("_density_over_cds", "_pct_f0", "_pct_f1", "_pct_f2")):
            out[col] = out[col].round(3)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.outdir / "smorf_flank_psites.csv", index=False)
    pd.DataFrame(struct).to_csv(args.outdir / "smorf_flank_structure.csv", index=False)

    w = "<side><window>, side = up|down, window = 50|100|250|1000 nt"
    pd.DataFrame([
        ("gene_id", "smORF ID (master gene_id; UniProt accession for TrEMBL)"),
        ("source", "Salk (novel smORF, microprotein GTF) or TrEMBL (smORF Rules track BED)"),
        ("cds_len", "CDS length in nt (spliced), including the stop codon"),
        ("stop_moved_into_cds", "True when the source CDS coordinates excluded the stop codon, so its 3 nt were taken from the 3' flank"),
        ("cds_psites", "Adult-brain P-sites in the CDS, all frames, same strand"),
        (f"{w}_len", "nt of host transcript available in the window (< window when the transcript ends first)"),
        (f"{w}_psites", "P-sites in the window, all frames"),
        (f"{w}_density_over_cds", "P-sites per nt in the window divided by P-sites per nt in the CDS"),
        (f"{w}_annCDS_bp", "nt of the window overlapping an annotated GENCODE CDS on the same strand"),
        (f"{w}_clean_n", f"P-sites used for the frame percentages: excludes annotated-CDS positions and the {START_SPILLOVER_NT} nt just upstream of the start codon"),
        (f"{w}_pct_f0/1/2", "% of clean P-sites in frame 0/1/2 relative to the smORF start (frame 0 = in frame; ~33/33/33 = no frame bias)"),
    ], columns=["column", "meaning"]).to_csv(args.outdir / "smorf_flank_data_dictionary.csv", index=False)

    print(f"done: {len(out)} smORFs {out.source.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
