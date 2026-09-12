#!/usr/bin/env bash
# Download and extract the large data assets that are too big for git.
#
# The repo ships the code and the master annotation table; the bulk data —
# the combined Ensembl+microprotein GTF, the three large CSVs, the RP3
# reference outputs, and the three figure libraries — lives on Zenodo. This
# script fetches those archives, verifies them, and extracts each one to the
# path the code expects.
#
# Usage:
#   ./download_zenodo_assets.sh                  # list the archives, download nothing
#   ./download_zenodo_assets.sh --all            # fetch everything (~5.1 GB)
#   ./download_zenodo_assets.sh large_data_tables large_reference_files
#   ./download_zenodo_assets.sh --check          # verify what is already cached
#   ./download_zenodo_assets.sh --cache DIR ...  # use a different cache location
#
# Re-running is safe and cheap: an archive that is already cached and passes
# its checksum is not downloaded again, so an interrupted multi-GB pull
# resumes at the next archive rather than starting over. Partial downloads are
# resumed byte-wise (curl -C -).
#
# Requires: curl, tar, python3, and one of shasum/md5/md5sum.

set -euo pipefail

# ── The published Zenodo record ─────────────────────────────────────────────
# This is the *version* record, pinned deliberately — NOT the concept DOI.
# zenodo.org/records/<concept-id> does not serve files, and a reader on an
# older clone should get the archives that clone's code expects rather than
# whatever is newest. Publishing a new version is a one-line edit here.
#
#   concept DOI (cite this):  10.5281/zenodo.20045161  -> always newest
#   version DOI (pinned):     10.5281/zenodo.21245767  -> v1.1, 2026-09-12
ZENODO_RECORD=21245767

BASE_URL="https://zenodo.org/records/${ZENODO_RECORD}/files"
API_URL="https://zenodo.org/api/records/${ZENODO_RECORD}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE_DIR="${ZENODO_CACHE_DIR:-$HOME/.cache/ad_microprotein_atlas}"

# ── Archive table ───────────────────────────────────────────────────────────
# name | extract-into (relative to repo root) | expected tar entries | size
#
# "extract-into" matters and is not uniform. Three archives store
# repo-relative paths internally (GTF_and_BED_files/…, Code/data/…) and so
# extract at the repo root; the figure archives are rooted at bare
# mirror_plots/ etc. and so extract into Results/; actin_micrographs is a flat
# bag of 25 PNGs with no directory of its own, so it MUST NOT be unpacked at
# the root — it goes into Code/Miscellaneous/, which is both where that
# module's README points --image-folder and what the existing .gitignore rule
# (Code/Miscellaneous/*.png) already covers.
#
# Entry counts are `tar -tz` lines: files PLUS directories. They are the real
# completeness check. Compressed size is not — v1.0.0 shipped two archives
# that looked plausible and were silently truncated, which is exactly what the
# checksum + entry-count verification here exists to catch on the reader side.
ARCHIVES=(
    "large_reference_files|.|44|76 MB"
    "large_data_tables|.|3|53 MB"
    "espresso_sequence_files|Code/Codon_context/espresso_sequences|4|158 MB"
    "figures_mirror_plots|Results|11379|2.3 GB"
    "figures_expression_profiles|Results|8655|1.8 GB"
    "figures_smorf_cartoons|Results|8689|371 MB"
    "actin_micrographs|Code/Miscellaneous|25|52 MB"
)

RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; YELLOW=$'\033[1;33m'
BLUE=$'\033[0;34m'; NC=$'\033[0m'

info()  { printf '%s\n' "$*"; }
ok()    { printf '%s✓%s %s\n' "$GREEN" "$NC" "$*"; }
warn()  { printf '%s!%s %s\n' "$YELLOW" "$NC" "$*"; }
fail()  { printf '%s✗%s %s\n' "$RED" "$NC" "$*" >&2; }

field() { printf '%s' "$1" | cut -d'|' -f"$2"; }

lookup() {
    # lookup <name> -> the matching ARCHIVES row, or empty
    local want="$1" row
    for row in "${ARCHIVES[@]}"; do
        [ "$(field "$row" 1)" = "$want" ] && { printf '%s' "$row"; return 0; }
    done
    return 1
}

# ── Checksums ───────────────────────────────────────────────────────────────
# Zenodo exposes MD5 in files[].checksum — note that despite what an earlier
# draft of the handoff plan said, it is MD5 and not SHA-256, so the repo's
# SHA256SUMS.txt cannot be compared against it directly. We read the expected
# digests from the live API rather than hardcoding them: the record is
# immutable once published, so the API is authoritative and a stale hardcoded
# table could only ever be wrong.
declare -a CK_NAMES=() CK_SUMS=() CK_SIZES=()

load_manifest() {
    local json
    info "==> Reading file manifest from Zenodo record ${ZENODO_RECORD}"
    if ! json="$(curl -sSfL --retry 3 --retry-delay 2 "$API_URL" 2>/dev/null)"; then
        fail "Could not reach the Zenodo API at $API_URL"
        info "    Check your connection, or download the archives by hand from"
        info "    https://doi.org/10.5281/zenodo.20045161"
        exit 1
    fi
    local parsed
    parsed="$(printf '%s' "$json" | python3 -c '
import sys, json
d = json.load(sys.stdin)
if not d.get("files"):
    sys.exit("no files in record")
for f in d["files"]:
    algo, _, digest = f.get("checksum", "").partition(":")
    print("%s\t%s\t%s\t%s" % (f["key"], algo, digest, f["size"]))')" || {
        fail "Could not parse the Zenodo record metadata."; exit 1; }

    local key algo digest size
    while IFS=$'\t' read -r key algo digest size; do
        [ -n "$key" ] || continue
        if [ "$algo" != "md5" ]; then
            warn "$key: unexpected checksum algorithm '$algo' — cannot verify"
            digest=""
        fi
        CK_NAMES+=("$key"); CK_SUMS+=("$digest"); CK_SIZES+=("$size")
    done <<< "$parsed"
    ok "manifest loaded (${#CK_NAMES[@]} files)"
}

expected_sum()  { local i; for i in "${!CK_NAMES[@]}"; do
    [ "${CK_NAMES[$i]}" = "$1" ] && { printf '%s' "${CK_SUMS[$i]}"; return; }; done; }
expected_size() { local i; for i in "${!CK_NAMES[@]}"; do
    [ "${CK_NAMES[$i]}" = "$1" ] && { printf '%s' "${CK_SIZES[$i]}"; return; }; done; }

# macOS ships `md5`; Linux ships `md5sum`. Both print the digest in a
# different column, hence the two branches.
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1;  then md5 -q "$1"
    else fail "Neither md5sum nor md5 found — cannot verify downloads."; exit 1; fi
}

# ── Per-archive operations ──────────────────────────────────────────────────

# Cached, complete, and checksum-clean?
cached_ok() {
    local name="$1" tarball="$CACHE_DIR/$1.tar.gz" want got
    [ -f "$tarball" ] || return 1
    want="$(expected_sum "$name.tar.gz")"
    [ -n "$want" ] || return 1
    got="$(md5_of "$tarball")"
    [ "$got" = "$want" ]
}

download() {
    local name="$1" url="$BASE_URL/$1.tar.gz?download=1"
    local tarball="$CACHE_DIR/$1.tar.gz" want got

    if cached_ok "$name"; then
        ok "$name — already downloaded and verified"
        return 0
    fi
    # A file that exists but failed the check is either a partial transfer or
    # corrupt in place. Resuming only helps the first case: if the file is
    # already the full published size, curl -C - considers it done and would
    # hand back the same bad bytes forever. So compare against the size Zenodo
    # reports and start over when there is nothing left to resume.
    if [ -f "$tarball" ]; then
        local have want_size
        have="$(wc -c < "$tarball" | tr -d ' ')"
        want_size="$(expected_size "$name.tar.gz")"
        if [ -n "$want_size" ] && [ "$have" -ge "$want_size" ]; then
            warn "$name — cached copy is complete but fails its checksum; re-downloading"
            rm -f "$tarball"
        else
            warn "$name — incomplete, resuming"
        fi
    fi

    info "==> Downloading $name.tar.gz"
    # -# gives a single-line progress bar instead of curl's multi-column
    # table, which scrolls unreadably in logs and CI.
    if ! curl -fL --progress-bar --retry 3 --retry-delay 2 -C - -o "$tarball" "$url"; then
        fail "$name — download failed"
        return 1
    fi

    want="$(expected_sum "$name.tar.gz")"
    if [ -z "$want" ]; then
        warn "$name — no checksum published; skipping verification"
        return 0
    fi
    got="$(md5_of "$tarball")"
    if [ "$got" != "$want" ]; then
        fail "$name — CHECKSUM MISMATCH"
        info "      expected md5 $want"
        info "      got      md5 $got"
        info "      The download is corrupt or truncated. The bad file is kept at"
        info "      $tarball"
        info "      so you can inspect it; delete it and re-run to try again."
        return 1
    fi
    ok "$name — checksum verified"
}

extract() {
    local name="$1" dest="$2" want_entries="$3"
    local tarball="$CACHE_DIR/$1.tar.gz" target got_entries

    target="$REPO_ROOT/$dest"
    mkdir -p "$target"

    # Verify the archive is whole BEFORE unpacking it over the working tree.
    # This is the reader-side equivalent of the check that caught the v1.0.0
    # truncation: a short entry count means a bad upload, and half-extracting
    # it would leave a confusing partial tree behind.
    info "==> Verifying $name"
    if ! got_entries="$(tar -tzf "$tarball" 2>/dev/null | wc -l | tr -d ' ')"; then
        fail "$name — archive is unreadable (corrupt gzip stream)"
        return 1
    fi
    if [ "$got_entries" != "$want_entries" ]; then
        fail "$name — entry count mismatch: expected $want_entries, found $got_entries"
        info "      Nothing was extracted. This usually means a bad or partial"
        info "      upload; report it rather than working around it."
        return 1
    fi
    ok "$name — $got_entries entries"

    info "==> Extracting $name -> $dest/"
    tar -xzf "$tarball" -C "$target"
    ok "$name — extracted"
}

check_only() {
    local row name want got tarball
    info "==> Checking cached archives in $CACHE_DIR"
    local any=0
    for row in "${ARCHIVES[@]}"; do
        name="$(field "$row" 1)"; tarball="$CACHE_DIR/$name.tar.gz"
        [ -f "$tarball" ] || continue
        any=1
        want="$(expected_sum "$name.tar.gz")"; got="$(md5_of "$tarball")"
        if [ "$got" = "$want" ]; then ok "$name — verified"
        else fail "$name — checksum mismatch (delete it and re-download)"; fi
    done
    [ "$any" -eq 1 ] || info "    (nothing cached yet)"
}

usage() {
    cat <<EOF
Download the Zenodo-hosted data assets for the Brain Microprotein Atlas.

Usage:
  $(basename "$0")                      list the archives and exit
  $(basename "$0") --all                download and extract everything
  $(basename "$0") NAME [NAME ...]      download and extract only these
  $(basename "$0") --check              verify already-cached archives
  $(basename "$0") --cache DIR ...      override the cache directory

Options:
  --all            fetch every archive (~5.1 GB download)
  --check          checksum what is already in the cache; download nothing
  --keep-cache     (default) leave tarballs in the cache after extracting
  --clean          delete each tarball once it extracts successfully
  --cache DIR      cache location (default: \$HOME/.cache/ad_microprotein_atlas,
                   override with \$ZENODO_CACHE_DIR)
  -h, --help       this message

Archives:
EOF
    local row
    for row in "${ARCHIVES[@]}"; do
        printf '  %-30s %8s   -> %s/\n' \
            "$(field "$row" 1)" "$(field "$row" 4)" "$(field "$row" 2)"
    done
    cat <<EOF

Data citation: https://doi.org/10.5281/zenodo.20045161  (concept DOI)
Pinned record: ${ZENODO_RECORD}  (v1.1)
EOF
}

# ── Argument parsing ────────────────────────────────────────────────────────
WANTED=(); DO_ALL=0; DO_CHECK=0; CLEAN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --all)        DO_ALL=1 ;;
        --check)      DO_CHECK=1 ;;
        --clean)      CLEAN=1 ;;
        --keep-cache) CLEAN=0 ;;
        --cache)      shift; [ $# -gt 0 ] || { fail "--cache needs a directory"; exit 1; }
                      CACHE_DIR="$1" ;;
        -h|--help)    usage; exit 0 ;;
        -*)           fail "Unknown option: $1"; echo; usage; exit 1 ;;
        *)            if lookup "$1" >/dev/null; then WANTED+=("$1")
                      else
                          fail "Unknown archive: $1"
                          info "    Run with --help to see the available names."
                          exit 1
                      fi ;;
    esac
    shift
done

mkdir -p "$CACHE_DIR"

if [ "$DO_CHECK" -eq 1 ]; then
    load_manifest; echo; check_only; exit 0
fi

if [ "$DO_ALL" -eq 0 ] && [ "${#WANTED[@]}" -eq 0 ]; then
    usage; exit 0
fi

if [ "$DO_ALL" -eq 1 ]; then
    WANTED=(); for row in "${ARCHIVES[@]}"; do WANTED+=("$(field "$row" 1)"); done
fi

# ── Run ─────────────────────────────────────────────────────────────────────
info "${BLUE}Brain Microprotein Atlas — Zenodo assets${NC}"
info "    record ${ZENODO_RECORD} · cache $CACHE_DIR"
echo
load_manifest
echo

failed=()
for name in "${WANTED[@]}"; do
    row="$(lookup "$name")"
    dest="$(field "$row" 2)"; entries="$(field "$row" 3)"
    echo "── $name ────────────────────────────────────────────"
    if ! download "$name";                 then failed+=("$name"); echo; continue; fi
    if ! extract "$name" "$dest" "$entries"; then failed+=("$name"); echo; continue; fi
    if [ "$CLEAN" -eq 1 ]; then
        rm -f "$CACHE_DIR/$name.tar.gz"; info "    cache cleared for $name"
    fi
    echo
done

if [ "${#failed[@]}" -gt 0 ]; then
    fail "Failed: ${failed[*]}"
    info "Re-run to retry — anything already verified will be skipped."
    exit 1
fi

ok "All requested assets are in place."
[ "$CLEAN" -eq 1 ] || info "Tarballs kept in $CACHE_DIR (--clean removes them)."
