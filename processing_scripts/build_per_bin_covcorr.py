#!/usr/bin/env python3
"""Build data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz.

`260801_nucleation_pub.ipynb` reads this file as

    _cv = pd.read_csv(COVFULL, sep="\t", usecols=["bin_id", "MTF2", "SUZ12", "EZH2"])

i.e. four columns. Each of the three signal columns is the mean of a bigWig over the 1-kb
bin, and cell 4 of that notebook re-checks it on a six-bin sample against the same bigWigs,
so the values must be the bigWig bin means and nothing else. This script writes those
columns, carrying the bin universe's seven coordinate columns through unchanged.

WHAT YOU NEED TO RUN IT

The three bigWigs are publicly obtainable (GEO and ENCODE, accessions below). The bin
universe is not: it comes from the Zenodo archive, or from
`replitag_nucleation_generate_matrices.sh` if you have run the upstream pipeline — which
needs duplicate-marked SAMs the authors produced with a lab-internal pipeline that is not
in this repository. **So this script does not run on a bare clone of the repository**; it
runs once you have unpacked the archive, and it is here so that you can rebuild this file
rather than having to trust the archived copy.

INPUTS
    data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz
        The bin universe: its first seven columns (domain_id, chr, start, end, bin_id,
        bin_index, tile_size) are copied through unchanged. Written by
        `replitag_nucleation_generate_matrices.sh`. The read-count matrix beside it is
        row-for-row identical in these columns — the rpkm matrix is derived from it by an
        `awk` that reprints $1..$7 verbatim, one output line per input line — so either may
        be used; the rpkm one is preferred here only because `260801_nucleation_pub.ipynb`
        already requires it.

    data/2024/MTF2_GSE164804/MTF2_shCT_hg38_coverage.bw          GEO GSE164804, lifted to hg38
    data/2024/MTF2_GSE164804/ENCODE_SUZ12_K562_ENCFF065KGU_fc_hg38.bw   ENCODE, fold-change over control
    data/2024/MTF2_GSE164804/ENCODE_EZH2_K562_ENCFF587SWK_fc_hg38.bw    ENCODE, fold-change over control

OUTPUT
    data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz

Run from the repository root:  python3 processing_scripts/build_per_bin_covcorr.py
Requires: numpy, pandas, pyBigWig.
"""
import argparse
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
SIGNAL_BW = {
    "MTF2": MD / "MTF2_shCT_hg38_coverage.bw",
    "SUZ12": MD / "ENCODE_SUZ12_K562_ENCFF065KGU_fc_hg38.bw",
    "EZH2": MD / "ENCODE_EZH2_K562_ENCFF587SWK_fc_hg38.bw",
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


def main():
    argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()

    if not MATRIX.exists():
        sys.exit(f"missing {MATRIX} — it comes from the Zenodo archive, or from "
                 f"processing_scripts/replitag_nucleation_generate_matrices.sh if you "
                 f"have run the upstream pipeline")
    missing = [str(p) for p in SIGNAL_BW.values() if not p.exists()]
    if missing:
        sys.exit("missing required bigWig(s):\n  " + "\n  ".join(missing))

    coord = pd.read_csv(MATRIX, sep="\t", usecols=META)[META]
    assert (coord.tile_size == 1000).all(), "this script assumes the 1000 bp tiling"
    assert coord.bin_id.is_unique, "bin_id is not unique in the source matrix"
    print(f"{len(coord):,} bins from {MATRIX}")

    cov = coord.copy()
    t0 = time.time()
    for name, path in SIGNAL_BW.items():
        cov[name] = bw_mean(path, coord)
        print(f"  {name:<6} {time.time() - t0:6.1f}s  "
              f"{int(np.isfinite(cov[name]).sum()):,} bins measurable")

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "per_bin_covcorr_full.tsv.gz"
    cov.to_csv(dest, sep="\t", index=False)
    print(f"\nwrote {dest}  {cov.shape[0]:,} x {cov.shape[1]}")
    print(f"  columns: {', '.join(cov.columns)}")


if __name__ == "__main__":
    main()
