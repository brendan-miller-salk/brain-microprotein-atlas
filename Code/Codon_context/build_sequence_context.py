#!/usr/bin/env python
"""
Build the sequence-context table both analysis steps of this pipeline read.

This is the *processing* half of this module and the only script here that touches the ESPRESSO
long-read sequence files. Those are ~1.5 GB and are not distributed with this repository (see
DATA_AVAILABILITY.md), so this script is documented rather than routinely run: its output,
`Code/data/codon_context_sequences.csv.gz`, is what ships, and `kozak_pipeline.py` /
`initiation_pipeline.py` read only that.

Three source files are scanned, each for something the others cannot give:

    <prefix>.split_nuc    each ORF's spliced CDS nucleotides
    <prefix>_ORFs.gtf     each ORF's exon blocks - a naive genomic slice would return intronic
                          sequence for a spliced ORF
    <prefix>.nuc          the assembled transcript, which is the only thing that makes the
                          upstream -1..-6 context available at all, since the ORF sequence
                          itself begins at +1

Only `FLANK_KEPT` (12) nt either side of each ORF is retained. That is twice the reach of any
scoring position in either step, so scoring against this table is identical to scoring
against the full transcript - not an approximation. Where the transcript provides less flank
than that, the short window is kept as-is and the unreachable positions score as '-', exactly as
they did before.

Run:
    python build_sequence_context.py --espresso-prefix /path/to/<prefix>

`--master` selects which ORFs to extract (default: the master table shipped with this repo).
Point it at the full annotation CSV to build the table for the complete discovery set.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from collections import defaultdict

import pandas as pd

import codon_context_common as common
from codon_context_common import FLANK_KEPT, ORF_ID_RE

log = logging.getLogger("sequence_context")


def _matching_key(name, wanted):
    """
    split_nuc / GTF ids may carry a trailing `_M` that the master CSV id lacks, and vice versa.

    An exact match always wins. The `_M` fallback only fires when the exact id is not wanted, so
    the ATG (`_M`) and stop-to-stop forms of the same ORF never collide even when both are in
    the cohort.
    """
    if name in wanted:
        return name
    alt = name[:-2] if name.endswith("_M") else name + "_M"
    return alt if alt in wanted else None


def read_orf_nucleotides(path, wanted):
    """ORF id -> spliced CDS nucleotides, from `.split_nuc`."""
    nt, header, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    nt[header] = "".join(buf).upper()
                header, buf = _matching_key(line[1:].strip(), wanted), []
            elif header:
                buf.append(line.strip())
    if header:
        nt[header] = "".join(buf).upper()
    log.info("      split_nuc: %d / %d ORF sequences", len(nt), len(wanted))
    return nt


def read_exon_blocks(path, wanted):
    """ORF id -> ([(lo, hi), ...], strand, chromosome), from `_ORFs.gtf`."""
    tid = re.compile(r'transcript_id "([^"]+)"')
    exons, strand, chrom = defaultdict(list), {}, {}
    with open(path) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            m = tid.search(p[8])
            key = _matching_key(m.group(1), wanted) if m else None
            if key:
                exons[key].append((int(p[3]), int(p[4])))
                strand[key], chrom[key] = p[6], p[0]
    log.info("      GTF: %d / %d exon structures", len(exons), len(wanted))
    return exons, strand, chrom


def read_transcript_flanks(path, wanted, nt):
    """
    ORF id -> {upstream, downstream, tx_offset, tx_copies}, from `.nuc`.

    The transcript is located by string search for the ORF's own CDS. `tx_copies` is carried
    through because an ORF occurring more than once in its transcript cannot be placed: the
    search cannot tell which copy the flanking window belongs to, and both steps treat that
    as a fatal QC flag rather than reporting a window that may be the wrong one.
    """
    flanks, header, buf = {}, None, []

    def flush(hdr, transcript):
        transcript = transcript.upper()
        # A `.nuc` header names every ORF carried by that transcript, whitespace separated.
        for raw in (hdr or "").split():
            key = _matching_key(raw, wanted)
            if key is None or key in flanks or key not in nt:
                continue
            i = transcript.find(nt[key])
            end = i + len(nt[key])
            flanks[key] = {
                "upstream": transcript[max(0, i - FLANK_KEPT):i] if i >= 0 else "",
                "downstream": transcript[end:end + FLANK_KEPT] if i >= 0 else "",
                "tx_offset": i, "tx_copies": transcript.count(nt[key])}

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(header, "".join(buf))
                header, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    flush(header, "".join(buf))
    log.info("      transcripts: %d / %d located", len(flanks), len(wanted))
    return flanks


def cohort_ids(master_csv):
    """Every unique smORF id in the master table, ATG and non-ATG alike."""
    ids = set()
    for chunk in pd.read_csv(master_csv, usecols=["gene_id"], chunksize=200_000,
                             low_memory=False):
        ids.update(s for s in chunk.gene_id.astype(str) if ORF_ID_RE.match(s))
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master", default=None,
                        help="master annotation CSV naming the ORFs to extract "
                             "(default: %s; env $%s)"
                             % (os.path.relpath(common.DEFAULT_MASTER, common.REPO_ROOT),
                                common.ENV_MASTER))
    parser.add_argument("--espresso-prefix", default=None,
                        help="path prefix of the ESPRESSO sequence files, i.e. the common part "
                             "before `.split_nuc` / `_ORFs.gtf` / `.nuc` (env $%s)"
                             % common.ENV_ESPRESSO)
    parser.add_argument("--out", default=None,
                        help="output table (default: %s; env $%s)"
                             % (os.path.relpath(common.DEFAULT_SEQUENCE_CONTEXT,
                                                common.REPO_ROOT), common.ENV_SEQUENCES))
    args = parser.parse_args()

    master = common.resolve(args.master, common.ENV_MASTER, common.DEFAULT_MASTER)
    prefix = common.resolve(args.espresso_prefix, common.ENV_ESPRESSO, None)
    out_path = common.resolve(args.out, common.ENV_SEQUENCES, common.DEFAULT_SEQUENCE_CONTEXT)
    if not prefix:
        raise SystemExit("--espresso-prefix (or $%s) is required: the ESPRESSO sequence files "
                         "are not distributed with this repository. See README.md."
                         % common.ENV_ESPRESSO)

    split_nuc, orfs_gtf, nuc = prefix + ".split_nuc", prefix + "_ORFs.gtf", prefix + ".nuc"
    for path in (master, split_nuc, orfs_gtf, nuc):
        if not os.path.exists(path):
            raise SystemExit("input not found: %s" % path)

    log_path = common.setup_logging("sequence_context",
                                    os.path.join(common.DEFAULT_RESULTS_ROOT, "logs"))
    log.info("=== sequence-context build ===")
    log.info("[1/3] reading %s", os.path.basename(master))
    wanted = cohort_ids(master)
    log.info("      %d smORF ids", len(wanted))

    log.info("[2/3] scanning ESPRESSO sequence files (~1.5 GB)")
    nt = read_orf_nucleotides(split_nuc, wanted)
    exons, strand, chrom = read_exon_blocks(orfs_gtf, wanted)
    flanks = read_transcript_flanks(nuc, wanted, nt)

    log.info("[3/3] writing %s", out_path)
    rows = []
    for orf in sorted(wanted):
        if orf not in nt or orf not in exons:
            # Both pipelines flag these as fatal QC. Recording them here as well would
            # duplicate that accounting with a row carrying no sequence.
            continue
        f = flanks.get(orf, {"upstream": "", "downstream": "", "tx_offset": -1, "tx_copies": 0})
        rows.append({"orf_id": orf, "chromosome": chrom.get(orf), "strand": strand.get(orf),
                     "exon_blocks": ";".join("%d-%d" % b for b in exons[orf]),
                     "orf_nt": nt[orf], "upstream_nt": f["upstream"],
                     "downstream_nt": f["downstream"], "tx_offset": f["tx_offset"],
                     "tx_copies": f["tx_copies"]})

    table = pd.DataFrame(rows, columns=common.SEQUENCE_CONTEXT_COLUMNS)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    table.to_csv(out_path, index=False)
    log.info("      %d / %d ORFs written (%.1f MB)", len(table), len(wanted),
             os.path.getsize(out_path) / 1e6)
    log.info("      %d ORFs could not be located in a transcript and have no flanking window",
             int((table.tx_offset < 0).sum()))
    log.info("done. log: %s", log_path)


if __name__ == "__main__":
    main()
