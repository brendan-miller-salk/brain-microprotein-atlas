#!/usr/bin/env python3
"""
Peptide-evidence flags for the unreviewed (Salk + TrEMBL) entries of the
Brain Microprotein Atlas.

For every detected peptide of every unreviewed entry, ask whether the same
peptide (or a same-mass peptide) occurs in any OTHER human UniProt protein,
and if so whether trypsin would release it from that protein.

Inputs
  --master   microprotein_master.csv (or the .zip containing it);
             default ../data/microprotein_master.zip
  --fasta    UniProt human proteome UP000005640 incl. Swiss-Prot, TrEMBL and
             isoforms (gzipped FASTA). Download:
             curl -L -o up5640_iso.fasta.gz \
              "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640&format=fasta&includeIsoform=true&compressed=true"
Outputs (in --outdir, default ../data/)
  unreviewed_microproteins_peptide_flags.csv
  peptide_flags_data_dictionary.csv
  peptide_flags_summary.csv

Then run add_twin_peptides.py to add the side-by-side peptide/twin columns.

Requires: python>=3.9, pandas, numpy.   Runtime: ~1.5 h on one core.
"""
import argparse, ast, bisect, gzip, itertools, re, time, zipfile, os
import numpy as np
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--master", default=os.path.join(DATA_DIR, "microprotein_master.zip"))
ap.add_argument("--fasta", required=True, help="UniProt UP000005640 FASTA incl. isoforms (.gz or plain)")
ap.add_argument("--outdir", default=DATA_DIR)
ap.add_argument("--limit", type=int, default=None, help="test run on the first N entries")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
t0 = time.time()

# ---------------------------------------------------------------- inputs
if args.master.endswith(".zip"):
    m = pd.read_csv(zipfile.ZipFile(args.master).open("microprotein_master.csv"), low_memory=False)
else:
    m = pd.read_csv(args.master, low_memory=False)
U = m[m.Database.isin(["Salk", "TrEMBL"])].copy()
if args.limit:
    U = U.head(args.limit)

names, seqs, cur, buf = [], [], None, []
opener = gzip.open if args.fasta.endswith(".gz") else open
with opener(args.fasta, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if cur is not None:
                names.append(cur); seqs.append("".join(buf))
            h = line[1:].split()
            acc = h[0].split("|")[1] if "|" in h[0] else h[0]
            gn = re.search(r"GN=(\S+)", line)
            cur = gn.group(1) if gn else acc
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        names.append(cur); seqs.append("".join(buf))
ref = pd.DataFrame({"name": names, "sequence": seqs}).drop_duplicates("sequence")
refnames = ref.name.tolist()
refseqs = ref.sequence.tolist()
refIL = [s.replace("I", "L") for s in refseqs]
bigIL = "|".join(refIL)        # I and L merged (MS cannot tell them apart)
bigLT = "|".join(refseqs)      # letters as written
offs = list(np.cumsum([0] + [len(s) + 1 for s in refseqs[:-1]]))
print(f"reference sequences: {len(refseqs)}", flush=True)


def L(x):
    try:
        return ast.literal_eval(x)
    except Exception:
        return []


def find(q, big, own, own_is_IL):
    """indices of reference sequences containing q, excluding the entry itself"""
    out = []
    for mt in re.finditer(re.escape(q), big):
        k = bisect.bisect_right(offs, mt.start()) - 1
        if (refIL[k] if own_is_IL else refseqs[k]) == own:
            continue
        out.append(k)
    return out


# --------------------------------------------- same-mass substitution rules
# elemental composition (C,H,N,O,S) of residues; 1- and 2-residue units with
# identical composition are interchangeable without changing mass
comp = {"G": (2,3,1,1,0), "A": (3,5,1,1,0), "S": (3,5,1,2,0), "P": (5,7,1,1,0),
        "V": (5,9,1,1,0), "T": (4,7,1,2,0), "C": (3,5,1,1,1), "L": (6,11,1,1,0),
        "I": (6,11,1,1,0), "N": (4,6,2,2,0), "D": (4,5,1,3,0), "Q": (5,8,2,2,0),
        "K": (6,12,2,1,0), "E": (5,7,1,3,0), "M": (5,9,1,1,1), "H": (6,7,3,1,0),
        "F": (9,9,1,1,0), "R": (6,12,4,1,0), "Y": (9,9,1,2,0), "W": (11,10,2,1,0)}
csum = lambda s: tuple(sum(comp[a][i] for a in s) for i in range(5))
aa = [a for a in comp if a != "I"]
units = aa + ["".join(p) for p in itertools.product(aa, repeat=2)]
groups = {}
for u in units:
    groups.setdefault(csum(u), []).append(u)
equiv = {u: [v for v in groups[csum(u)] if v != u] for u in units if len(groups[csum(u)]) > 1}

mono = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
        "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293,
        "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
        "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931}


def n_distinguishing_ions(p, q):
    """number of b/y ions whose mass differs between sequences p and q"""
    def ions(s):
        return ([sum(mono[a] for a in s[:i]) for i in range(1, len(s))],
                [sum(mono[a] for a in s[i:]) for i in range(1, len(s))])
    bp, yp = ions(p); bq, yq = ions(q)
    f = lambda A, B: sum(all(abs(a - b) > 0.005 for b in B) for a in A)
    return f(bp, bq) + f(yp, yq)


def tryptic_status(var, ks):
    """would trypsin release `var` from any of the matching proteins?"""
    best = "non-tryptic"
    for k in ks:
        s = refseqs[k]
        for mt in re.finditer(re.escape(var), refIL[k]):
            i, j = mt.start(), mt.end()
            prev = s[i - 1] if i > 0 else "^"
            nxt = s[j] if j < len(s) else "$"
            n_ok = prev in "KR^"
            c_ok = s[j - 1] in "KR" or nxt == "$"
            st = "tryptic" if (n_ok and c_ok) else ("semi-tryptic" if (n_ok or c_ok) else "non-tryptic")
            if st == "tryptic":
                return st
            if st == "semi-tryptic":
                best = st
    return best


def classify_peptide(p, ownIL, own):
    """returns (match_type, matched_variant, reference indices)"""
    q = p.replace("I", "L")
    ks = find(q, bigIL, ownIL, True)
    if ks:
        kind = "identical" if find(p, bigLT, own, False) else "identical_except_IL"
        return kind, q, ks
    for i in range(len(q)):
        for w in (1, 2):
            u = q[i:i + w]
            for v in equiv.get(u, []):
                var = q[:i] + v + q[i + w:]
                ks = find(var, bigIL, ownIL, True)
                if ks:
                    return "same_mass_rearranged", var, ks
    for i, a in enumerate(q):
        for a2 in (["D"] if a == "N" else ["E"] if a == "Q" else []):
            var = q[:i] + a2 + q[i + 1:]
            ks = find(var, bigIL, ownIL, True)
            if ks:
                return "deamidation", var, ks
    return "unique", "", []


def n_sites(peps):
    """non-overlapping peptide locations (Met+/- and nested peptides collapsed)"""
    sites = []
    for p in sorted(set(peps), key=len, reverse=True):
        core = p[1:] if p.startswith("M") else p
        if not any(core in s or s in core for s in sites):
            sites.append(core)
    return len(sites)


severity = {"identical": 4, "identical_except_IL": 3, "same_mass_rearranged": 2, "deamidation": 1, "unique": 0}
MATCH_TEXT = {
    "identical": "identical sequence found in another protein",
    "identical_except_IL": "identical except I/L (same mass) in another protein",
    "same_mass_rearranged": "different sequence, same mass (residues rearranged) in another protein",
    "deamidation": "differs from another protein only by N/Q deamidation",
    "unique": "unique: no other UniProt protein contains these peptides"}
HOW_TEXT = {
    "tryptic": "normal trypsin cleavage (shared peptide)",
    "semi-tryptic": "needs one non-trypsin cut (fragment of that protein)",
    "non-tryptic": "needs two non-trypsin cuts (fragment of that protein)"}

# ---------------------------------------------------------------- main loop
rows = []
for n, (_, r) in enumerate(U.iterrows()):
    peps = L(r.peptide_sequence)
    base = dict(gene_id=r.gene_id, gene_name=r.gene_name, sequence=r.sequence)
    if not peps:
        rows.append({**base, "interpretation": "No MS evidence", "peptide_match_in_UniProt": "no MS peptides",
                     "peptides_with_this_match": "", "matching_protein": "", "matching_protein_is_same_gene": "",
                     "how_matching_protein_would_give_this_peptide": "", "n_fragment_ions_that_distinguish": np.nan,
                     "n_independent_peptide_sites": 0, "min_peptide_len": np.nan})
        continue
    own = str(r.sequence); ownIL = own.replace("I", "L")
    kinds, prots, stat, nd = [], set(), [], []
    for p in peps:
        kind, var, ks = classify_peptide(p, ownIL, own)
        kinds.append(kind)
        if ks:
            prots |= {refnames[k] for k in ks}
            stat.append(tryptic_status(var, ks))
            nd.append(n_distinguishing_ions(p.replace("I", "L"), var))
    matched = [k for k in kinds if k != "unique"]
    if not matched:
        match = "unique"
    elif set(matched) <= {"identical", "identical_except_IL"}:
        # sequence-level match: identical if any peptide is letter-identical
        match = "identical" if all(k == "identical" for k in matched) else (
            "identical_except_IL" if all(k == "identical_except_IL" for k in matched) else "identical")
    else:
        match = max(matched, key=lambda k: severity[k])
    some = bool(matched) and len(matched) < len(peps)
    ts = "tryptic" if "tryptic" in stat else ("semi-tryptic" if "semi-tryptic" in stat else ("non-tryptic" if stat else ""))
    plist = sorted(prots)[:8]
    same_gene = (str(r.gene_name) in plist) if plist else None

    if match == "unique":
        interp = "Supported: peptides unique to this entry"
    elif match in ("identical", "identical_except_IL"):
        if ts == "tryptic":
            interp = "Ambiguous: peptide shared with another protein"
        elif same_gene:
            interp = "Ambiguous: could be a fragment of the full-length protein from the same gene"
        else:
            interp = "Ambiguous: could be a fragment of a different protein"
    elif match == "same_mass_rearranged":
        interp = "Ambiguous: same-mass peptide in another protein; fragment ions may separate them"
    else:
        interp = "Mostly supported: differs only by a deamidation-like mass match"
    if some:
        interp += " (only some peptides; others unique)"

    rows.append({**base, "interpretation": interp, "peptide_match_in_UniProt": MATCH_TEXT[match],
                 "peptides_with_this_match": "" if match == "unique" else ("some" if some else "all"),
                 "matching_protein": ";".join(plist),
                 "matching_protein_is_same_gene": "" if same_gene is None else ("yes" if same_gene else "no"),
                 "how_matching_protein_would_give_this_peptide": HOW_TEXT.get(ts, ""),
                 "n_fragment_ions_that_distinguish": (min(nd) if (nd and match == "same_mass_rearranged") else np.nan),
                 "n_independent_peptide_sites": n_sites(peps),
                 "min_peptide_len": min(len(p) for p in peps)})
    if n % 500 == 0:
        print(f"{n}/{len(U)}  {round(time.time() - t0)} s", flush=True)

out = pd.DataFrame(rows)
out.to_csv(os.path.join(args.outdir, "unreviewed_microproteins_peptide_flags.csv"), index=False)
(out[out.interpretation != "No MS evidence"].groupby("interpretation").size()
    .rename("n_entries").sort_values(ascending=False).reset_index()
    .to_csv(os.path.join(args.outdir, "peptide_flags_summary.csv"), index=False))

# Full dictionary, including the twin columns that add_twin_peptides.py appends.
pd.DataFrame([
 ('interpretation',
  "One-line verdict for the entry's MS evidence"),
 ('entry_peptides',
  'All detected peptides of this entry, shown as previous.PEPTIDE.next residue within the entry (^ = entry start, $ = entry end)'),
 ('matched_entry_peptide',
  "The entry peptide(s) that also match another protein (' | ' separates multiple peptides; same order in the next columns)"),
 ('twin_peptide',
  'The matching stretch in the other protein, shown as previous.PEPTIDE.next residue in that protein'),
 ('differences',
  "Residue differences entry→twin, e.g. I3L (I at position 3 is L in the twin), GA4-5Q (GA at 4–5 is Q in the twin); 'none' = letter-identical"),
 ('twin_peptide_protein',
  'UniProt accession|gene of the protein carrying the twin, with residue positions'),
 ('twin_peptide_trypsin',
  'Whether trypsin releases the twin from that protein: tryptic, semi-tryptic (one non-trypsin cut needed), non-tryptic (two)'),
 ('peptide_match_in_UniProt',
  "Whether the detected peptides occur in any other UniProt human protein (Swiss-Prot + TrEMBL + isoforms, UP000005640). 'identical sequence' = letter-for-letter the same; 'identical except I/L' = same except isoleucine/leucine, which have identical mass and cannot be told apart by MS; 'different sequence, same mass' = same atoms, residues rearranged (e.g. QLLFPIVR vs QLLPFLVR); 'unique' = no match"),
 ('peptides_with_this_match',
  'all = every detected peptide has the match; some = only some (the rest are unique)'),
 ('matching_protein',
  'Gene name(s) of the other protein(s) containing the peptide'),
 ('matching_protein_is_same_gene',
  'yes = the other protein is from the same gene as this entry (e.g. the full-length protein)'),
 ('how_matching_protein_would_give_this_peptide',
  "Whether trypsin would release this peptide from the matching protein. The search allowed only trypsin cuts, so 'needs a non-trypsin cut' peptides were never considered for the matching protein and were credited to this entry, whose start/end creates the cut"),
 ('n_fragment_ions_that_distinguish',
  'Only for rearranged same-mass matches: number of b/y ions whose mass differs between the two sequences (0 = spectrum cannot separate them)'),
 ('n_independent_peptide_sites',
  'Number of non-overlapping peptide locations on the entry'),
 ('min_peptide_len',
  'Length of the shortest detected peptide'),
], columns=["column", "meaning"]).to_csv(os.path.join(args.outdir, "peptide_flags_data_dictionary.csv"), index=False)

print(f"done: {len(out)} entries in {round(time.time() - t0)} s")
print(out.interpretation.value_counts().to_string())
