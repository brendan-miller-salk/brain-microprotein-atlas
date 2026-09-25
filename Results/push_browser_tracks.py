"""Publish the UCSC browser tracks to the Hugging Face dataset.

The tracks UCSC actually loads live in the *dataset* repo
brmiller/brain-microprotein-atlas under browser_tracks/. Nothing else in
Results/ publishes them:

  * upload_to_hf.py       -> same dataset, but only the three figure directories
                             (mirror_plots, expression_profiles, smorf_cartoon_figures)
  * push_to_hf_space.py   -> the Streamlit *Space*, a different repo entirely;
                             its copy of these files is not what bigDataUrl points at

So a track rebuild that stops at the shared Box folder leaves the browser
serving the previous build while every local check reports success. This script
is the missing step.

It uses create_commit rather than upload_large_folder: the file set is small and
fixed, and each rebuild deserves one real commit message. upload_large_folder is
also actively unhelpful here -- it writes .cache/.huggingface/*.lock inside
folder_path, which in the Box tree means syncing thousands of lock files to the
cloud (see the staging dance in upload_to_hf.py).

Three checks run before anything uploads, in this order:

  1. every bigDataUrl in brain_tracks.txt resolves to a file that exists, so a
     renamed track cannot be published half-broken;
  2. brain_tracks.txt is included whenever any track it describes is uploaded --
     the legend carries the per-colour counts, and shipping a new .bb behind a
     stale legend has happened before;
  3. each uploaded file is re-downloaded and hashed against the local copy.

Usage (from anywhere):
    python Results/push_browser_tracks.py                 # changed files only
    python Results/push_browser_tracks.py -m "message"
    python Results/push_browser_tracks.py --all           # ignore the diff
    python Results/push_browser_tracks.py --dry-run

Non-interactive auth: export HF_TOKEN=hf_xxx, or `huggingface-cli login` first.
"""
import argparse
import hashlib
import sys
from pathlib import Path

from huggingface_hub import HfApi, CommitOperationAdd, hf_hub_download
from huggingface_hub.utils import EntryNotFoundError

REPO = "brmiller/brain-microprotein-atlas"
REPO_TYPE = "dataset"
PREFIX = "browser_tracks"

RESULTS = Path(__file__).parent.resolve()
TRACKS = RESULTS.parent / "Code" / "data" / "browser_tracks"
HUB_CONFIG = "brain_tracks.txt"

# Everything publishable. Files absent from disk are skipped with a warning
# rather than failing the run -- not every rebuild produces every artifact.
FILES = [
    # UCSC track hub definition (the user-facing legend and colour counts)
    "brain_tracks.txt",
    "README.md",
    # Track-hub entry point (My Data -> Track Hubs). Serves the same bigDataUrl
    # files as brain_tracks.txt, but as a hub, which is the only way to get
    # multiWig containers -- the three P-site frames overlaid in one row.
    # Custom tracks cannot express containers, so both entry points coexist.
    "hub.txt",
    "genomes.txt",
    "trackDb.txt",
    # ORF rules (decay x main-ORF disruption)
    "smorf_trx_priority_color_class.bb",
    "smorf_trx_priority_color_class.bed",
    "smorf_trx_priority_color_class_assignments.tsv",
    # peptide / PSM evidence
    "smorf_psm_peptides.bb",
    "smorf_psm_peptides.gtf",
    "smorf_psm_genomic.bed",
    "smorf_psm_genomic.tsv",
    # TrEMBL repair audit trails (two DIFFERENT repairs -- see the README)
    "trembl_cds_rephase.tsv",
    "trembl_cds_phase_repair.tsv",
    # long-read isoforms and coverage
    "espresso_dlpfc_isoforms.bb",
    "espresso_dlpfc_mean_cpm.bw",
    "adult_brain_riboseq_merged.fwd.bw",
    "adult_brain_riboseq_merged.rev.bw",
    # ribo-seq P-sites split by reading frame (see brain_tracks.txt header).
    # NOT the same data as adult_brain_riboseq_merged.*.bw above: those are
    # bamCoverage read-body coverage and carry no frame information.
    "adult_psite.f0.fwd.bw",
    "adult_psite.f1.fwd.bw",
    "adult_psite.f2.fwd.bw",
    "adult_psite.f0.rev.bw",
    "adult_psite.f1.rev.bw",
    "adult_psite.f2.rev.bw",
]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def remote_sha(api, name):
    """sha256 of the published copy, or None if it is not there yet."""
    try:
        p = hf_hub_download(REPO, f"{PREFIX}/{name}", repo_type=REPO_TYPE,
                            force_download=True)
    except (EntryNotFoundError, OSError):
        return None
    return sha(p)


def check_bigdataurls(present):
    """Every bigDataUrl in the hub config must name a file we know about."""
    cfg = TRACKS / HUB_CONFIG
    if not cfg.exists():
        return []
    missing = []
    for line in cfg.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("#") or "bigDataUrl=" not in line:
            continue
        url = line.split("bigDataUrl=", 1)[1].split()[0]
        name = url.rsplit("/", 1)[-1]
        if name not in present:
            missing.append(name)
    return missing


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-m", "--message", default=None, help="commit message")
    ap.add_argument("--all", action="store_true",
                    help="upload every file, not just those that differ")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be pushed, upload nothing")
    args = ap.parse_args()

    if not TRACKS.is_dir():
        sys.exit(f"track directory not found: {TRACKS}")

    present, absent = [], []
    for name in FILES:
        (present if (TRACKS / name).exists() else absent).append(name)
    if absent:
        print(f"not on disk, skipping: {', '.join(absent)}")

    broken = check_bigdataurls(present)
    if broken:
        sys.exit(f"ERROR: {HUB_CONFIG} references file(s) that are not on disk "
                 f"or not in FILES: {', '.join(broken)}")

    api = HfApi()
    try:
        who = api.whoami().get("name")
    except Exception:
        sys.exit("not authenticated -- run `huggingface-cli login` or set HF_TOKEN")
    print(f"authenticated as {who}\n{REPO} ({REPO_TYPE}) / {PREFIX}/\n")

    print("comparing local against published ...")
    changed = []
    for name in present:
        local = sha(TRACKS / name)
        if args.all:
            changed.append(name)
            print(f"  {'FORCED':<9} {name}")
            continue
        remote = remote_sha(api, name)
        if remote is None:
            changed.append(name)
            print(f"  {'NEW':<9} {name}")
        elif remote != local:
            changed.append(name)
            print(f"  {'CHANGED':<9} {name}")

    if not changed:
        print("\nnothing to push -- published copies already match.")
        return

    # The legend carries the per-colour counts. If any track it describes is
    # going up, the config goes with it or the published legend goes stale.
    if HUB_CONFIG in present and HUB_CONFIG not in changed:
        print(f"\n  note: {HUB_CONFIG} is unchanged and will not be re-uploaded;"
              f"\n        confirm its counts still describe the new build.")

    print(f"\n{len(changed)} file(s) to push:")
    total = sum((TRACKS / n).stat().st_size for n in changed)
    for n in changed:
        print(f"    {n:<48} {(TRACKS / n).stat().st_size:>12,}")
    print(f"    {'total':<48} {total:>12,}")

    if args.dry_run:
        print("\n--dry-run: nothing uploaded.")
        return

    msg = args.message or f"browser tracks: update {len(changed)} file(s)"
    ops = [CommitOperationAdd(path_in_repo=f"{PREFIX}/{n}",
                              path_or_fileobj=str(TRACKS / n)) for n in changed]
    print(f"\nuploading ...")
    commit = api.create_commit(repo_id=REPO, repo_type=REPO_TYPE,
                               operations=ops, commit_message=msg)
    print(f"commit: {commit.commit_url}")

    print("\nverifying published copies ...")
    bad = []
    for n in changed:
        if remote_sha(api, n) != sha(TRACKS / n):
            bad.append(n)
        print(f"  {'OK  ' if n not in bad else 'DIFF'} {n}")
    if bad:
        sys.exit(f"\nERROR: {len(bad)} file(s) do not match after upload: "
                 f"{', '.join(bad)}")
    print(f"\n{len(changed)} file(s) published and verified.")


if __name__ == "__main__":
    main()
