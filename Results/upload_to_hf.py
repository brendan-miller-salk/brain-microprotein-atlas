"""Upload three figure directories to Hugging Face dataset, preserving folder names.

Strategy: temporarily move the three directories into a staging parent so
upload_large_folder picks up the parent layout. Move them back when done.
"""
import shutil
from pathlib import Path
from huggingface_hub import HfApi

REPO = "brmiller/brain-microprotein-atlas"
RESULTS = Path(__file__).parent.resolve()
STAGE = RESULTS / "_hf_stage"
DIRS = ["mirror_plots", "expression_profiles", "smorf_cartoon_figures"]

api = HfApi()

# 1. Stage
STAGE.mkdir(exist_ok=True)
for d in DIRS:
    src = RESULTS / d
    dst = STAGE / d
    if src.exists() and not dst.exists():
        print(f"Staging {d} -> _hf_stage/{d}")
        src.rename(dst)
    elif dst.exists():
        print(f"Already staged: {d}")
    else:
        print(f"WARNING: {d} not found at {src}")

# 2. Upload the staging parent (folder names preserved at repo root)
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
    # 3. Always restore (even if upload failed)
    print("\nRestoring original locations...")
    for d in DIRS:
        src = STAGE / d
        dst = RESULTS / d
        if src.exists() and not dst.exists():
            src.rename(dst)
            print(f"Restored {d}")
    # Remove empty stage dir
    try:
        STAGE.rmdir()
    except OSError:
        print(f"NOTE: {STAGE} not empty; remove manually if desired.")