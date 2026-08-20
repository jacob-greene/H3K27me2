# H3K27me2

Analysis code for the manuscript *"Histone H3K27 methylation states are sequentially
catalyzed in cycling cells"* by Greene, Ahmad, and Henikoff, 2026.

## Data availability

| Source | Accession | Contents |
|---|---|---|
| GEO | [`GSE327802`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327802) | Raw CUT&Tag sequencing data |
| Zenodo | [`10.5281/zenodo.19489947`](https://doi.org/10.5281/zenodo.19489947) | Processed data (the rest of the `data/` folder) |
| UCSC | [browser session](https://genome.ucsc.edu/s/jegreene/H3K27me2_manuscript) | Interactive tracks |

---

## Quick start

```bash
git clone https://github.com/jacob-greene/H3K27me2.git
cd H3K27me2

# Two notebooks run immediately, with no download at all:
jupyter lab K27me3_states_111epigenomes_pub.ipynb       # fetches its own data over HTTPS
jupyter lab 250214_growth_viability_cellcycle_pub.ipynb # reads only data/ as shipped

# For the other six, download the processed-data archive from Zenodo and unpack it so
# its contents land inside the data/ directory that ships with this repo. The result
# must look like "Expected data/ layout" below.
```

**Every path in this repository is relative to the repository root.** Start Jupyter from
the repository root, and run the shell/R scripts from the repository root
(`bash processing_scripts/RIFs_RepliTag.sh`).

Figures and figure-source tables are written **flat into `figures/`**. There are no
sub-directories, and no two cells write the same filename.

---

## Repository contents

```
H3K27me2/
├── *_pub.ipynb              8 publication notebooks (the figures)
├── processing_scripts/      23 scripts that build the processed data
├── data/                    input data — a few files ship here, the rest from Zenodo
├── figures/                 created by the notebooks; PDFs and figure-source TSVs land here
└── THIRD_PARTY_NOTICES.md   licence terms for the two redistributed supplement files
```

### What ships in `data/`

Everything else in `data/` comes from the Zenodo archive.

| File | Size | Provenance |
|---|---|---|
| `250214_K27me123ac_P2S25p_TAZEEDUTX_cell_counts.xlsx` | 17 KB | this study |
| `20260206_drug_8only_aurora/*.csv` (15 files) | 6.3 MB | this study — Aurora flow-cytometry exports |
| `Kc_merged_blacklist_hg38.bed` | 7.9 MB | this study — *Drosophila* spike-in blacklist, 286,386 intervals |
| `2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | 630 KB | this study — per-gene max RPKM + replication timing |
| `chr_lens.txt` | 12 KB | UCSC hg38 chromosome sizes |
| `2024/LAD_data/Dataset_S1_Promoter_SuRE_Classification.tsv` | 3.2 MB | **published supplement, not this study** |
| `2024/LAD_data/Dataset_S2_TRIP.tsv` | 1.1 MB | **published supplement, not this study** |

`Dataset_S1` / `Dataset_S2` are from Leemans et al., *Cell* 177:852–864 (2019), redistributed
unmodified under **CC BY-NC-ND 4.0** so that `260803_WT_timeseries_clean_chromHMM_pub.ipynb`
runs without a manual download. **Cite that paper, not this one, for those two files**, and
see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the terms.

---

## Notebooks

### These notebooks need no local data

- **`K27me3_states_111epigenomes_pub.ipynb`** needs **no `data/` files at all**. It
  downloads the Roadmap Epigenomics core-marks segmentations from S3 at run time. It does
  need `bedtools` on `PATH` (`conda install -c bioconda bedtools`; or set `BEDTOOLS_BIN` to
  prepend a specific directory). It *writes*
  `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed`, which one other notebook reads.
- **`250214_growth_viability_cellcycle_pub.ipynb`** reads only files that ship in this
  repository. No download required.

The other six need the Zenodo archive.

### Run order

Two notebooks are producers as well as consumers, so order matters:

1. **`K27me3_states_111epigenomes_pub.ipynb`** — writes
   `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed`.
2. **`250411_peakPrint_clean_pub.ipynb`** — writes
   `data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/*_top95regions.bed`.

Then, in any order: `250214_growth_viability_cellcycle_pub.ipynb`,
`250911_RIFs_FC_pub.ipynb`, `251203_drug_timeseries_pub.ipynb`,
`260801_nucleation_pub.ipynb`, `260803_WT_timeseries_clean_chromHMM_pub.ipynb`,
`bulk_RT_pub.ipynb`.

---

## Expected `data/` layout

```
data/
├── 2024/
│   ├── cCREs/intersected_annos/
│   │   ├── FRIPs.tsv
│   │   └── FC_out/*_featureCounts_bulk_counts.tsv, *.summary.tsv
│   ├── ENCODE_chromHMM/{FRIPs.tsv, OverlapEnrichments.tsv}
│   ├── ENCODE_states/K562_E123_15_coreMarks_domains.bed*
│   ├── K562_annotations/
│   │   ├── gencode.v44.basic.annotation.coding.gtf
│   │   └── processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}†
│   ├── LAD_data/
│   │   ├── gencode.v27.annotation.gtf.gz
│   │   └── *.saf
│   ├── MTF2_GSE164804/*.bw       MTF2_shCT + ENCFF587SWK / ENCFF065KGU fold-change
│   ├── RepliSeq/*.bigWig         UCSC-ENCODE UW Repli-seq K562 **hg19** —
│   │                             the INPUT to crossmaphg19tohg38.sh
│   ├── RepliTag/
│   │   ├── bws_bulk/*.bw
│   │   ├── filtered_sams_WT_1.0/, filtered_sams_drug_1.0/
│   │   ├── Matrices/{coverage_scores_*.tab, TAZEED/}
│   │   ├── refs/GOBP_CELL_CYCLE.v2025.1.Hs.json
│   │   └── RepliSeq/hg38_bws_crossmap/     (written by crossmaphg19tohg38.sh)
│   └── SH_all_data/
│       ├── sections4/, merged_bams4/, merged_bws4/
│       ├── peakPRINT_50/MERGED_downsampled_bins.tsv
│       └── downsample_peakPRINT_slope1_min300/
│           ├── reproducibility/FRiPs_merged_intervals_rep95.tsv
│           └── min350/{*.bed, *_meta_v*.tsv, top95/*}
└── results/
    ├── mtf2_nucleation/{nuc_cpgonly_spread_bins.tsv.gz, per_bin_covcorr_full.tsv.gz‡}
    └── replitag_sphase_nucleation/matrices/
            encode_coremarks.1000bp.{rpkm,read_count}.matrix.tsv.gz
```

`*` = written by `K27me3_states_111epigenomes_pub.ipynb`, which reads only Roadmap
Epigenomics S3 URLs and so has no blocked input. It is the one file here you can produce
from a bare clone, and it is not needed in the archive if you run the notebooks in order.

`†` = producible by `gtf2bed_TSS.sh` from the GENCODE v44 GTF, given `gtftools`.

`‡` = producible by `build_per_bin_covcorr.py`, but only once the rpkm matrix beside it is
in place.

`260803_WT_timeseries_clean_chromHMM_pub.ipynb` also *writes*
`data/2024/LAD_data/H3K27_TRIP.tsv`, which nothing reads; it is not required in the archive.

**External downloads** the cloner must fetch: GENCODE v27 and v44 GTFs; MSigDB
`GOBP_CELL_CYCLE.v2025.1.Hs.json`; UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; GEO
`GSE164804` (MTF2) and ENCODE `ENCFF587SWK` (EZH2) / `ENCFF065KGU` (SUZ12) fold-change
bigWigs.

---

## Known gaps

**Upstream of everything — read trimming and alignment.** `filter_sams*.sh` consume
duplicate-marked bowtie2 SAMs, which the authors produced with a lab-internal pipeline that
is not in this repository. `geo_to_sams.sh` can help process fastqs into sam files; the raw reads are in GEO `GSE327802`.

**Intermediate tables with no producer**, which must come from the Zenodo archive:

- `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz` — the 1-kb-bin analysis
  frame that all of `260801_nucleation_pub.ipynb` is built on (site class, super-domain,
  replication timing, per-fraction H3K27me3). Nothing here produces it.
- `data/2024/SH_all_data/peakPRINT_50/MERGED_downsampled_bins.tsv`,
  `.../downsample_peakPRINT_slope1_min300/reproducibility/FRiPs_merged_intervals_rep95.tsv`
- `data/2024/SH_all_data/sections4/*.txt` and `All_sams4.txt` (BAM path manifests)
- `data/2024/ENCODE_chromHMM/{FRIPs,OverlapEnrichments}.tsv` and
  `data/2024/cCREs/intersected_annos/FRIPs.tsv` (ChromHMM `OverlapEnrichment` output)
- `data/2024/RepliTag/bws_bulk/` — assembled by hand from the `bam2bed2bw.sh` bigWigs plus
  the lifted Repli-seq track; no script does the assembly

---



