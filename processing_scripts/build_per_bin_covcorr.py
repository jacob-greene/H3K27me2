#!/usr/bin/env python3
"""Build data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz.

This is the minimal derivation of the one matrix `260801_nucleation_pub.ipynb` needs that
no other script in this repository produces. It is the data-assembly section (section 1)
of `260728_mtf2_correlation_nucleation_JG.ipynb` and nothing else — that notebook's
correlation analysis, candidate-nucleation-site analysis and figures are not reproduced
here and it remains in this directory, unmodified, as the record of them.

WHAT THE NUCLEATION NOTEBOOK ACTUALLY READS

    _cv = pd.read_csv(COVFULL, sep="\t", usecols=["bin_id", "MTF2", "SUZ12", "EZH2"])

i.e. four columns. Those four are always written. Each is the mean of a bigWig over the
1-kb bin, and `260801_nucleation_pub.ipynb` cell 4 re-checks that on a six-bin sample
against the same bigWigs, so the values must be the bigWig bin means and nothing else.

WHAT IS OPTIONAL

`260728_mtf2_correlation_nucleation_JG.ipynb` also reads `ME3`, `DNase`, `CGI` and
`MTF2_peak` from this file. They belong to that notebook's confound analysis, not to the
nucleation figures, and three of their four inputs have no public source (see README,
"Known gaps"). They are written when their inputs are present and skipped, loudly, when
they are not — so running this script on a clean clone gives you a file that runs the
nucleation notebook, and running it with the full archive in place gives you the file the
correlation notebook expects. Which columns were written, and which were skipped, is
printed at the end of the run.

The script will not overwrite an existing output unless you pass --force; the copy in the
Zenodo archive has every optional column, and a clean-clone run does not.

INPUTS
    data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz
        The bin universe: its first seven columns (domain_id, chr, start, end, bin_id,
        bin_index, tile_size) are copied through unchanged. Written by
        `replitag_nucleation_generate_matrices.sh`. The read-count matrix beside it is
        row-for-row identical in these columns (the rpkm matrix is derived from it by
        `awk`), so either may be used; the rpkm one is preferred here only because
        `260801_nucleation_pub.ipynb` already requires it.

    data/2024/MTF2_GSE164804/MTF2_shCT_hg38_coverage.bw          GEO GSE164804, lifted to hg38
    data/2024/MTF2_GSE164804/ENCODE_SUZ12_K562_ENCFF065KGU_fc_hg38.bw   ENCODE, fold-change over control
    data/2024/MTF2_GSE164804/ENCODE_EZH2_K562_ENCFF587SWK_fc_hg38.bw    ENCODE, fold-change over control

    Optional, for the four confound columns:
    data/2024/RepliTag/bws_bulk/K27me3_K_no_chrM.1-1000.bw       -> ME3
    data/2024/K562_annotations/EncodeRegDnaseUwK562Peak.bed      -> DNase
    data/2024/K562_annotations/processed_beds/hg38_CpG_Island.bed -> CGI
    data/2024/MTF2_GSE164804/MTF2_macs3_hg38_peaks.bed           -> MTF2_peak

OUTPUT
    data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz

Run from the repository root:  python3 processing_scripts/build_per_bin_covcorr.py
Requires: numpy, pandas, pyBigWig.
"""
import argparse
import bisect
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

P = Path(".")
RES = P / "data/results/replitag_sphase_nucleation"
OUT = P / "data/results/mtf2_nucleation"
MD = P / "data/2024/MTF2_GSE164804"

MATRIX = RES / "matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz"
META = ["domain_id", "chr", "start", "end", "bin_id", "bin_index", "tile_size"]

# read as the mean of the bigWig over each bin
REQUIRED_BW = {
    "MTF2": MD / "MTF2_shCT_hg38_coverage.bw",
    "SUZ12": MD / "ENCODE_SUZ12_K562_ENCFF065KGU_fc_hg38.bw",
    "EZH2": MD / "ENCODE_EZH2_K562_ENCFF587SWK_fc_hg38.bw",
}
OPTIONAL_BW = {
    "ME3": P / "data/2024/RepliTag/bws_bulk/K27me3_K_no_chrM.1-1000.bw",
}
# read as "does this bin touch any interval in the BED"
OPTIONAL_BED = {
    "DNase": P / "data/2024/K562_annotations/EncodeRegDnaseUwK562Peak.bed",
    "CGI": P / "data/2024/K562_annotations/processed_beds/hg38_CpG_Island.bed",
    "MTF2_peak": MD / "MTF2_macs3_hg38_peaks.bed",
}


def bw_mean(bw_path, coord):
    """Mean bigWig signal over each bin of `coord`, NaN where it cannot be measured.

    NaN, not 0, for: a chromosome the bigWig does not carry, a bin that falls off the end
    of its chromosome, and a bin the bigWig has no data over. The nucleation notebook
    distinguishes "not measurable" from "measured as zero" when it builds its site sets,
    so the two must not be conflated here.
    """
    bw = pyBigWig.open(str(bw_path))
    try:
        have = bw.chroms()
        chrom = coord["chr"].to_numpy()
        starts = coord["start"].to_numpy(np.int64)
        ends = coord["end"].to_numpy(np.int64)
        out = np.full(len(coord), np.nan)
        for ch in pd.unique(chrom):
            if ch not in have:
                continue
            clen = have[ch]
            for i in np.flatnonzero(chrom == ch):
                s, e = max(0, int(starts[i])), min(int(ends[i]), clen)
                if e <= s:
                    continue
                r = bw.stats(ch, s, e, type="mean")[0]
                if r is not None:
                    out[i] = r
        return out
    finally:
        bw.close()


def load_bed(path):
    """{chrom: (starts, prefix-max of ends)} for the `chr`-prefixed intervals of a BED."""
    by_chrom = {}
    opener = open
    if str(path).endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("track", "#", "browser")):
                continue
            f = line.split("\t")
            if len(f) >= 3 and f[0].startswith("chr"):
                by_chrom.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    out = {}
    for ch, ivs in by_chrom.items():
        arr = np.array(sorted(ivs), dtype=np.int64)
        # sorted by start, so "any interval starting at or before `be` ends after `bs`"
        # is a prefix query on the running maximum of the ends
        out[ch] = (arr[:, 0], np.maximum.accumulate(arr[:, 1]))
    return out


def overlaps(coord, bed):
    """Does each bin of `coord` overlap any interval in `bed`? (half-open coordinates)"""
    chrom = coord["chr"].to_numpy()
    starts = coord["start"].to_numpy(np.int64)
    ends = coord["end"].to_numpy(np.int64)
    flag = np.zeros(len(coord), bool)
    for ch in pd.unique(chrom):
        hit = bed.get(ch)
        if hit is None:
            continue
        a_start, a_endmax = hit
        rows = np.flatnonzero(chrom == ch)
        # bins whose start >= the last interval's end can be answered without a search
        for i in rows:
            n = bisect.bisect_right(a_start, int(ends[i]))
            if n > 0 and a_endmax[n - 1] > starts[i]:
                flag[i] = True
    return flag


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing per_bin_covcorr_full.tsv.gz")
    args = ap.parse_args()

    dest = OUT / "per_bin_covcorr_full.tsv.gz"
    if dest.exists() and not args.force:
        sys.exit(f"{dest} already exists (the Zenodo archive ships one). "
                 f"Pass --force to overwrite it.")
    if not MATRIX.exists():
        sys.exit(f"missing {MATRIX} — run processing_scripts/"
                 f"replitag_nucleation_generate_matrices.sh first")
    missing = [str(p) for p in REQUIRED_BW.values() if not p.exists()]
    if missing:
        sys.exit("missing required bigWig(s):\n  " + "\n  ".join(missing))

    coord = pd.read_csv(MATRIX, sep="\t", usecols=META)[META]
    assert (coord.tile_size == 1000).all(), "this script assumes the 1000 bp tiling"
    assert coord.bin_id.is_unique, "bin_id is not unique in the source matrix"
    print(f"{len(coord):,} bins from {MATRIX}")

    cov = coord.copy()
    written, skipped = [], []
    t0 = time.time()

    for name, path in {**REQUIRED_BW, **OPTIONAL_BW}.items():
        if not path.exists():
            skipped.append(f"{name} (no {path})")
            continue
        cov[name] = bw_mean(path, coord)
        written.append(name)
        print(f"  {name:<10} {time.time() - t0:6.1f}s  "
              f"{int(np.isfinite(cov[name]).sum()):,} bins measurable")

    for name, path in OPTIONAL_BED.items():
        if not path.exists():
            skipped.append(f"{name} (no {path})")
            continue
        cov[name] = overlaps(coord, load_bed(path))
        written.append(name)
        print(f"  {name:<10} {time.time() - t0:6.1f}s  {int(cov[name].sum()):,} bins overlap")

    need = ["bin_id", "MTF2", "SUZ12", "EZH2"]
    assert all(c in cov.columns for c in need), \
        f"260801_nucleation_pub.ipynb reads {need}; this frame is missing some"

    dest.parent.mkdir(parents=True, exist_ok=True)
    cov.to_csv(dest, sep="\t", index=False)

    print(f"\nwrote {dest}  {cov.shape[0]:,} x {cov.shape[1]}")
    print(f"  columns: {', '.join(cov.columns)}")
    if skipped:
        print("\n  NOT written, inputs absent — this file will run "
              "260801_nucleation_pub.ipynb but NOT the confound sections of "
              "260728_mtf2_correlation_nucleation_JG.ipynb:")
        for s in skipped:
            print(f"    {s}")


if __name__ == "__main__":
    main()
