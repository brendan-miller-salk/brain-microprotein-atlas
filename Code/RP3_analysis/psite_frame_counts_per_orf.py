"""P-site counts per ORF, split by codon position (frame 0/1/2 relative to the ORF start).

usage:
    python psite_frame_counts_per_orf.py --total-psites N [--gtf GTF] [--fasta FASTA]
                                         [--bigwig-prefix PREFIX] [--out OUT]
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import binomtest

repo = Path(__file__).resolve().parents[2]
keep_attrs = ["gene_id", "gene_name", "transcript_type", "protein_id", "parent_orf_id", "smorf_type"]
attr_re = re.compile(r'(\S+) "([^"]*)"')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--gtf", type=Path,
                   default=repo / "GTF_and_BED_files/Ensembl_and_Unreviewed_Brain_Microproteins.gtf",
                   help="GTF with CDS features")
    p.add_argument("--fasta", type=Path,
                   default=repo / "GTF_and_BED_files/Unreviewed_Brain_Microproteins.fasta",
                   help="smORF protein FASTA, header '>transcript_id|gene_name|...'")
    p.add_argument("--bigwig-prefix", type=Path,
                   default=repo / "Code/data/browser_tracks/adult_psite",
                   help="P-site bigwigs named <prefix>.f{0,1,2}.{fwd,rev}.bw, in CPM")
    p.add_argument("--total-psites", type=int, required=True,
                   help="total P-sites behind the bigwigs, to convert CPM back to counts")
    p.add_argument("--smorf-sources", nargs="+", default=["GTF2FastaPatched", "AltInitiation"],
                   help="GTF source values that mark smORF records")
    p.add_argument("--out", type=Path, default=repo / "Code/data/psite_frame_counts_per_orf.tsv")
    return p.parse_args()


def read_gtf_cds(path):
    orfs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "CDS":
                continue
            attrs = dict(attr_re.findall(f[8]))
            key = (f[1], attrs["transcript_id"])
            if key not in orfs:
                orfs[key] = {"chrom": f[0], "strand": f[6], "blocks": [],
                             "attrs": {k: attrs.get(k, "") for k in keep_attrs}}
            phase = 0 if f[7] == "." else int(f[7])
            orfs[key]["blocks"].append((int(f[3]) - 1, int(f[4]), phase))
    return orfs


def read_fasta(path):
    seqs = {}
    tid = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                fields = line[1:].split("|")
                tid = fields[0]
                seqs[tid] = [fields[1] if len(fields) > 1 else "", ""]
            elif tid:
                seqs[tid][1] += line
    return seqs


def get_psites(bws, chrom, start, end, total_psites):
    x = np.zeros(end - start)
    for bw in bws:
        if chrom not in bw.chroms():
            return x
        x += np.nan_to_num(bw.values(chrom, start, end, numpy=True))
    return np.round(x * total_psites / 1e6)


def count_frames(orf, bws, total_psites):
    minus = orf["strand"] == "-"
    bws = bws["rev"] if minus else bws["fwd"]
    blocks = sorted(set(orf["blocks"]), key=lambda b: b[0], reverse=minus)
    phase = blocks[0][2]

    x = []
    for start, end, _ in blocks:
        v = get_psites(bws, orf["chrom"], start, end, total_psites)
        x.append(v[::-1] if minus else v)
    x = np.concatenate(x)
    cds_len = len(x)
    x = x[phase:]
    n_codons = len(x) // 3
    codons = x[:n_codons * 3].reshape(n_codons, 3)

    res = {"n_blocks": len(blocks), "cds_len_nt": cds_len, "cds_phase": phase,
           "cds_len_mod3": cds_len % 3, "n_codons": n_codons}
    for i in range(3):
        res[f"Psites_sum_frame{i}"] = codons[:, i].sum()
    for i in range(3):
        res[f"Psites_coverage_frame{i}"] = (codons[:, i] > 0).mean() if n_codons else 0
    return res


def main():
    args = parse_args()
    smorf_sources = set(args.smorf_sources)

    orfs = read_gtf_cds(args.gtf)
    print(len(orfs), "CDS records")

    seqs = read_fasta(args.fasta) if args.fasta and args.fasta.exists() else {}
    for (source, tid), orf in orfs.items():
        if source in smorf_sources and tid in seqs:
            orf["attrs"]["gene_name"] = orf["attrs"]["gene_name"] or seqs[tid][0]
            orf["sequence"] = seqs[tid][1]

    prefix = str(args.bigwig_prefix)
    bws = {s: [pyBigWig.open(f"{prefix}.f{i}.{s}.bw") for i in range(3)] for s in ["fwd", "rev"]}

    rows = []
    for n, ((source, tid), orf) in enumerate(orfs.items()):
        if n % 10000 == 0:
            print(n)
        row = {"source": source, "transcript_id": tid, "is_smorf": source in smorf_sources}
        row.update(orf["attrs"])
        row["sequence"] = orf.get("sequence", "")
        row["chrom"] = orf["chrom"]
        row["strand"] = orf["strand"]
        row["cds_start"] = min(b[0] for b in orf["blocks"]) + 1
        row["cds_end"] = max(b[1] for b in orf["blocks"])
        row.update(count_frames(orf, bws, args.total_psites))
        rows.append(row)
    df = pd.DataFrame(rows)

    frames = ["Psites_sum_frame0", "Psites_sum_frame1", "Psites_sum_frame2"]
    df[frames] = df[frames].astype(int)
    df["Psites_total"] = df[frames].sum(axis=1)
    total = df["Psites_total"].where(df["Psites_total"] > 0)
    for i in range(3):
        df[f"Psites_pct_frame{i}"] = 100 * df[f"Psites_sum_frame{i}"] / total

    codons = df["n_codons"].where(df["n_codons"] > 0)
    df["Psites_per_codon"] = df["Psites_total"] / codons
    df["Psites_frame0_per_codon"] = df["Psites_sum_frame0"] / codons
    df["Psites_frame0_RPKM"] = df["Psites_sum_frame0"] * 1e9 / (df["cds_len_nt"] * args.total_psites)
    df["binom_p_frame0"] = [binomtest(int(k), int(t), 1 / 3, alternative="greater").pvalue if t > 0 else np.nan
                            for k, t in zip(df["Psites_sum_frame0"], df["Psites_total"])]

    cols = (["source", "transcript_id", "is_smorf"] + keep_attrs
            + ["sequence", "chrom", "strand", "cds_start", "cds_end",
               "n_blocks", "cds_len_nt", "cds_phase", "cds_len_mod3", "n_codons"]
            + frames + ["Psites_total"]
            + [f"Psites_pct_frame{i}" for i in range(3)]
            + [f"Psites_coverage_frame{i}" for i in range(3)]
            + ["Psites_per_codon", "Psites_frame0_per_codon", "Psites_frame0_RPKM", "binom_p_frame0"])
    df[cols].to_csv(args.out, sep="\t", index=False, float_format="%.4g")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
