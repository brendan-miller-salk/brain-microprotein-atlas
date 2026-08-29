"""Deploy the microprotein dashboard to a public Hugging Face Streamlit Space.

This mirrors the app that already runs on Streamlit Community Cloud. It uploads
*exactly* the git-tracked file set (same content clean_and_push_repo.sh commits
to GitHub), so the large figure libraries (mirror_plots/, expression_profiles/,
smorf_cartoon_figures/) are NOT shipped — the dashboard streams those images
from the existing Hugging Face dataset brmiller/brain-microprotein-atlas.

Prerequisites:
    pip install huggingface_hub
    huggingface-cli login          # or set HF_TOKEN in the environment
    # `git ls-files` reads the index, so anything added since the last commit
    # must be staged first or it will not ship. launch_hf_dashboard.sh does this
    # for you (its step 2); if you run this script standalone, `git add -A` first.

Usage (from Results/):
    python push_to_hf_space.py

After the first deploy, in the Space UI go to Settings -> Variables and add a
plain variable  DASHBOARD_PUBLIC = 1  (NOT a secret) so the mirror is
password-free. Do NOT set DASHBOARD_PASSWORD_HASH on the Space.
"""
import subprocess
from pathlib import Path

from huggingface_hub import HfApi, create_repo

SPACE_REPO = "brmiller/brain-microprotein-atlas-app"
RESULTS = Path(__file__).resolve().parent
REPO_ROOT = RESULTS.parent  # the Github/ folder — becomes the Space root
APP_FILE = "Results/microproteins_dashboard.py"
# HF removed the native Streamlit SDK; Streamlit runs under the Docker SDK.
SPACE_SDK = "docker"
# Files needed by the Space that may not be git-tracked yet (the deploy runs
# before clean_and_push_repo.sh commits them). upload_folder walks the actual
# working tree, so listing them here ensures they ship regardless.
EXTRA_FILES = ["Dockerfile"]

# The GitHub README.md is kept free of YAML front matter (GitHub renders that
# block as an ugly metadata table). Hugging Face, however, requires the Space
# config in the README's front matter, so we prepend this block only on the
# copy pushed to the Space.
SPACE_README_YAML = """\
---
title: Brain Microprotein Atlas
emoji: 🧠
colorFrom: blue
colorTo: gray
sdk: docker
app_port: 8501
pinned: false
license: mit
---

"""


def tracked_files():
    """Return the list of git-tracked paths (relative to REPO_ROOT)."""
    out = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    files = [line.strip() for line in out.splitlines() if line.strip()]
    if not files:
        raise SystemExit(
            "No git-tracked files found. Run clean_and_push_repo.sh first so "
            "`git ls-files` reflects the committed set."
        )
    return files


def main():
    files = tracked_files()
    # Include extra (possibly untracked) files that exist on disk.
    for extra in EXTRA_FILES:
        if (REPO_ROOT / extra).exists() and extra not in files:
            files.append(extra)
    print(f"Deploying {len(files)} files from {REPO_ROOT}")
    print(f"  app_file = {APP_FILE}  (Docker SDK)")
    print(f"  (image libraries excluded; served from the HF dataset fallback)")

    api = HfApi()

    # 1. Create the Space if it does not exist (Docker SDK, public).
    create_repo(
        repo_id=SPACE_REPO,
        repo_type="space",
        space_sdk=SPACE_SDK,
        private=False,
        exist_ok=True,
    )

    # 2. Upload exactly the tracked set. LFS is handled automatically for the
    #    31 MB Code/data/microprotein_master.zip.
    print(f"\nUploading to space {SPACE_REPO} ...")
    api.upload_folder(
        repo_id=SPACE_REPO,
        repo_type="space",
        folder_path=str(REPO_ROOT),
        allow_patterns=files,
        commit_message="Deploy Brain Microprotein Atlas dashboard",
    )

    # 3. Streamlit on Spaces reads .streamlit/config.toml from the repo ROOT,
    #    but the tracked copy lives at Results/.streamlit/. Mirror it to the
    #    root so the dark theme applies (port stays 8501, the Spaces default).
    root_cfg = REPO_ROOT / "Results" / ".streamlit" / "config.toml"
    if root_cfg.exists():
        api.upload_file(
            path_or_fileobj=str(root_cfg),
            path_in_repo=".streamlit/config.toml",
            repo_id=SPACE_REPO,
            repo_type="space",
            commit_message="Add root Streamlit config (theme)",
        )

    # 4. The Space README needs YAML front matter for its config, but the
    #    GitHub README is kept YAML-free (GitHub renders it as a table).
    #    Upload a README that prepends the Space config to the GitHub content.
    readme = REPO_ROOT / "README.md"
    space_readme = SPACE_README_YAML + readme.read_text(encoding="utf-8")
    api.upload_file(
        path_or_fileobj=space_readme.encode("utf-8"),
        path_in_repo="README.md",
        repo_id=SPACE_REPO,
        repo_type="space",
        commit_message="Add Space config front matter to README",
    )

    print("\nUpload complete.")
    print(f"Space: https://huggingface.co/spaces/{SPACE_REPO}")
    print(
        "\nNEXT: In the Space -> Settings -> Variables, add\n"
        "    DASHBOARD_PUBLIC = 1   (plain variable, not a secret)\n"
        "so the mirror loads without the password gate."
    )


if __name__ == "__main__":
    main()
