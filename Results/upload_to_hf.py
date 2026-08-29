"""Upload three figure directories to a Hugging Face dataset, preserving names.

Strategy: temporarily move the three directories into a staging parent so
``upload_large_folder`` picks up the parent layout, then move them back.

The staging parent lives OUTSIDE the Box-synced tree (default: ``~/Desktop``).
This matters because ``upload_large_folder`` writes resumable bookkeeping
(``.cache/.huggingface/`` with a ``.lock`` per file) *inside* ``folder_path`` —
if that sits in Box Drive, Box syncs thousands of ``*.lock`` files to the cloud.
Desktop and the Box tree are on the same filesystem, so the moves are instant
renames (no data copy). Override the location with ``HF_STAGE_DIR=/some/path``.
"""
import os
import shutil
from pathlib import Path

from huggingface_hub import HfApi

REPO = "brmiller/brain-microprotein-atlas"
RESULTS = Path(__file__).parent.resolve()
DIRS = ["mirror_plots", "expression_profiles", "smorf_cartoon_figures"]

# Staging parent OUTSIDE Box (keeps upload .lock/.cache off Box sync).
STAGE = Path(os.environ.get("HF_STAGE_DIR", Path.home() / "Desktop")).expanduser().resolve() / "hf_upload_stage"

# Guard: warn if the stage still ends up inside a synced CloudStorage tree.
if "CloudStorage" in STAGE.parts:
    print(f"WARNING: staging dir {STAGE} is inside a synced folder; "
          f"set HF_STAGE_DIR to a local path to avoid cloud-sync churn.")


def _move(src: Path, dst: Path):
    """Move src -> dst, falling back to a copy if across filesystems."""
    try:
        src.rename(dst)
    except OSError:
        shutil.move(str(src), str(dst))


api = HfApi()

# 1. Stage (move each figure dir into the out-of-Box staging parent).
STAGE.mkdir(parents=True, exist_ok=True)
for d in DIRS:
    src = RESULTS / d
    dst = STAGE / d
    if src.exists() and not dst.exists():
        print(f"Staging {d} -> {STAGE.name}/{d}")
        _move(src, dst)
    elif dst.exists():
        print(f"Already staged: {d}")
    else:
        print(f"WARNING: {d} not found at {src}")

# 2. Upload the staging parent (folder names preserved at repo root).
try:
    print(f"\nUploading {STAGE} -> {REPO} ...")
    api.upload_large_folder(
        repo_id=REPO,
        repo_type="dataset",
        folder_path=str(STAGE),
        num_workers=8,
    )
    print("\nUpload complete.")
finally:
    # 3. Always restore original locations (even if the upload failed).
    print("\nRestoring original locations...")
    for d in DIRS:
        src = STAGE / d
        dst = RESULTS / d
        if src.exists() and not dst.exists():
            _move(src, dst)
            print(f"Restored {d}")

    # 4. Remove the whole staging parent, including upload .cache/.lock files.
    stranded = [d for d in DIRS if (STAGE / d).exists()]
    if stranded:
        print(f"NOTE: {STAGE} still holds {stranded}; not removing.")
    elif STAGE.exists():
        shutil.rmtree(STAGE, ignore_errors=True)
        print(f"Cleaned up {STAGE}")