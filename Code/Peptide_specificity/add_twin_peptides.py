#!/usr/bin/env python3
"""Add side-by-side peptide comparison columns to unreviewed_microproteins_peptide_flags.csv
(the output of make_peptide_flags.py):
   entry_peptides               all detected peptides of the entry, shown with flanking residues in the entry (^ = entry start, $ = entry end)
   matched_entry_peptide        the entry peptide(s) that have a match in another protein
   twin_peptide                 the matching stretch in the other protein, flanking residues shown (prev.PEPTIDE.next)
   twin_peptide_protein         UniProt accession|gene of the protein carrying the twin, and residue positions
   twin_peptide_trypsin         whether trypsin would release the twin from that protein
   differences                  position-by-position differences between entry peptide and twin (e.g. I5L, FP4-5PF)
Multiple matched peptides are separated by ' | ' in the same order across columns.

Usage: python add_twin_peptides.py --fasta up5640_iso.fasta.gz   # updates ../data/ flag table in place"""
import argparse, ast, bisect, gzip, itertools, os, re, zipfile
from multiprocessing import Pool
import numpy as np, pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
FLAGS_CSV = os.path.join(DATA_DIR, "unreviewed_microproteins_peptide_flags.csv")

# Parsed at module level: Pool workers (spawn on macOS) re-import this module and need FASTA.
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--master", default=os.path.join(DATA_DIR, "microprotein_master.zip"))
ap.add_argument("--fasta", required=True, help="UniProt UP000005640 FASTA incl. isoforms (.gz or plain)")
ap.add_argument("--flags", default=FLAGS_CSV, help="input flag table from make_peptide_flags.py")
ap.add_argument("--out", default=FLAGS_CSV, help="output path (default: overwrite --flags)")
ap.add_argument("--processes", type=int, default=8)
args = ap.parse_args()
MASTER, FASTA, FLAGS, OUT = args.master, args.fasta, args.flags, args.out

names, seqs, cur, buf = [], [], None, []
with (gzip.open if FASTA.endswith(".gz") else open)(FASTA, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if cur is not None: names.append(cur); seqs.append("".join(buf))
            h = line[1:].split(); acc = h[0].split("|")[1] if "|" in h[0] else h[0]
            gn = re.search(r"GN=(\S+)", line); cur = f"{acc}|{gn.group(1) if gn else acc}"; buf = []
        else: buf.append(line.strip())
    if cur is not None: names.append(cur); seqs.append("".join(buf))
ref = pd.DataFrame({"name": names, "sequence": seqs}).drop_duplicates("sequence")
refnames, refseqs = ref.name.tolist(), ref.sequence.tolist()
refIL = [s.replace("I", "L") for s in refseqs]
bigIL = "|".join(refIL)
offs = list(np.cumsum([0] + [len(s) + 1 for s in refseqs[:-1]]))

comp = {"G":(2,3,1,1,0),"A":(3,5,1,1,0),"S":(3,5,1,2,0),"P":(5,7,1,1,0),"V":(5,9,1,1,0),"T":(4,7,1,2,0),"C":(3,5,1,1,1),
        "L":(6,11,1,1,0),"I":(6,11,1,1,0),"N":(4,6,2,2,0),"D":(4,5,1,3,0),"Q":(5,8,2,2,0),"K":(6,12,2,1,0),"E":(5,7,1,3,0),
        "M":(5,9,1,1,1),"H":(6,7,3,1,0),"F":(9,9,1,1,0),"R":(6,12,4,1,0),"Y":(9,9,1,2,0),"W":(11,10,2,1,0)}
csum = lambda s: tuple(sum(comp[a][i] for a in s) for i in range(5))
aa = [a for a in comp if a != "I"]
groups = {}
for u in aa + ["".join(p) for p in itertools.product(aa, repeat=2)]:
    groups.setdefault(csum(u), []).append(u)
equiv = {u: [v for v in g if v != u] for g in groups.values() if len(g) > 1 for u in g}

def L(x):
    try: return ast.literal_eval(x)
    except Exception: return []

def find(q, ownIL):
    out = []
    for mt in re.finditer(re.escape(q), bigIL):
        k = bisect.bisect_right(offs, mt.start()) - 1
        if refIL[k] != ownIL: out.append((k, mt.start() - offs[k]))
    return out

def variants(q):
    yield q
    for i in range(len(q)):
        for w in (1, 2):
            for v in equiv.get(q[i:i+w], []): yield q[:i] + v + q[i+w:]
    for i, a in enumerate(q):
        if a == "N": yield q[:i] + "D" + q[i+1:]
        if a == "Q": yield q[:i] + "E" + q[i+1:]

def tryp(s, i, j):
    prev = s[i-1] if i > 0 else "^"; nxt = s[j] if j < len(s) else "$"
    n_ok = prev in "KR^"; c_ok = s[j-1] in "KR" or nxt == "$"
    return "tryptic" if n_ok and c_ok else ("semi-tryptic" if n_ok or c_ok else "non-tryptic")

def ctx(s, i, j):
    return f"{s[i-1] if i > 0 else '^'}.{s[i:j]}.{s[j] if j < len(s) else '$'}"

def diff(p, t):
    if len(p) != len(t):   # e.g. GG <-> N, AG <-> Q (two residues vs one, same atoms)
        import difflib
        ops = difflib.SequenceMatcher(None, p, t, autojunk=False).get_opcodes()
        return ",".join(f"{p[a:b]}{a+1}{'' if b-a == 1 else '-'+str(b)}{t[c:d]}" for tag, a, b, c, d in ops if tag != "equal") or "none"
    out, i = [], 0
    while i < len(p):
        if p[i] == t[i]: i += 1; continue
        j = i
        while j < len(p) and p[j] != t[j]: j += 1
        out.append(f"{p[i:j]}{i+1}{'' if j-i == 1 else '-'+str(j)}{t[i:j]}"); i = j
    return ",".join(out) if out else "none"

RANK = {"tryptic": 0, "semi-tryptic": 1, "non-tryptic": 2}
def process(args):
    gid, own, peps, search = args
    ownIL = own.replace("I", "L")
    allp = []
    for p in peps:
        i = own.find(p); allp.append(ctx(own, i, i + len(p)) if i >= 0 else p)
    m_e, m_t, m_p, m_s, m_d = [], [], [], [], []
    for p in (peps if search else []):
        for var in variants(p.replace("I", "L")):
            hits = find(var, ownIL)
            if not hits: continue
            best = min(hits, key=lambda h: (RANK[tryp(refseqs[h[0]], h[1], h[1] + len(var))], h[0]))
            k, i = best; j = i + len(var); s = refseqs[k]
            i0 = own.find(p)
            m_e.append(ctx(own, i0, i0 + len(p)) if i0 >= 0 else p)
            m_t.append(ctx(s, i, j)); m_p.append(f"{refnames[k]}:{i+1}-{j}"); m_s.append(tryp(s, i, j))
            m_d.append(diff(p, s[i:j]))
            break
    J = " | ".join
    return gid, " | ".join(allp), J(m_e), J(m_t), J(m_p), J(m_s), J(m_d)

if __name__ == "__main__":
    src = zipfile.ZipFile(MASTER).open("microprotein_master.csv") if MASTER.endswith(".zip") else MASTER
    m = pd.read_csv(src, low_memory=False).drop_duplicates("gene_id").set_index("gene_id")
    fl = pd.read_csv(FLAGS)
    todo = fl[fl.interpretation != "No MS evidence"]
    jobs = [(g, str(m.at[g, "sequence"]), L(m.at[g, "peptide_sequence"]), not str(u).startswith("unique"))
            for g, u in zip(todo.gene_id, todo.peptide_match_in_UniProt)]
    jobs.sort(key=lambda j: not j[3])
    with Pool(args.processes) as pool:
        res = pool.map(process, jobs, chunksize=20)
    cols = ["entry_peptides", "matched_entry_peptide", "twin_peptide", "twin_peptide_protein", "twin_peptide_trypsin", "differences"]
    add = pd.DataFrame(res, columns=["gene_id"] + cols).set_index("gene_id")
    out = fl.drop(columns=[c for c in cols if c in fl.columns]).join(add, on="gene_id")
    for c in cols: out[c] = out[c].fillna("")
    order = ["gene_id", "gene_name", "sequence", "interpretation", "entry_peptides", "matched_entry_peptide", "twin_peptide",
             "differences", "twin_peptide_protein", "twin_peptide_trypsin"]
    out = out[order + [c for c in out.columns if c not in order]]
    out.to_csv(OUT, index=False)
    print("done", len(out))
