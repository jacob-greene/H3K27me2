# H3K27me2

![H3K27 methylation states are sequentially catalyzed in cycling cells](Fig1.png)

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

conda env create -f envs/notebooks.yml && conda activate h3k27me2-notebooks

# Two notebooks run immediately, with no download at all:
jupyter lab K27me3_states_111epigenomes_pub.ipynb       # fetches its own data over HTTPS
jupyter lab 250214_growth_viability_cellcycle_pub.ipynb # reads only data/ as shipped

# For the other six, download the processed-data archive from Zenodo and unpack it so
# its contents land inside the data/ directory that ships with this repo. The result
# must look like "Expected data/ layout" below.
```

**Every path in this repo is relative to the repository root.** Start Jupyter from
the repository root, and run the processing scripts from the repository root
(`bash processing_scripts/run.sh <stage>`).

### Environments

| File | Covers |
|---|---|
| `envs/notebooks.yml` | the 8 publication notebooks |
| `envs/processing.yml` | `processing_scripts/`, except the two R scripts |
| `envs/r.yml` | `251113_DGE_RUVseq.R`, `251202_DGE_RUVseq_drug.R` |

Only `bowtie2` and `cutadapt` are pinned, to the versions the manuscript Methods name.
One dependency is not installable from conda: **`gtftools`**, needed by the `annotation`
stage, comes from <https://www.genemine.org/gtftools.php>.

Figures and figure-source tables are written **flat into `figures/`**. There are no
sub-directories, and no two cells write the same filename.

---

## Repository contents

```
H3K27me2/
├── *_pub.ipynb              8 publication notebooks (the figures)
├── processing_scripts/      24 scripts that build the processed data; run.sh drives them
├── envs/                    conda environment files
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

## Processing scripts

`processing_scripts/run.sh` drives them. Run it from the repository root:

```bash
conda env create -f envs/processing.yml && conda activate h3k27me2-processing

bash processing_scripts/run.sh list        # the stages, in pipeline order
bash processing_scripts/run.sh selftest    # ~1 min, needs no Zenodo data
bash processing_scripts/run.sh <stage>     # DRY_RUN=1 prints without running
```

`selftest` calls peaks on synthetic fragments; it exists to check that a clone is wired up
before you spend hours on real data. Every other stage needs the Zenodo archive unpacked
into `data/` first, and three carry caveats that `run.sh` repeats when you invoke them:
`align` was transcribed from the Methods and **has never been executed by the authors**;
`blacklist` reads four lab-internal files and its output already ships; `annotation` needs
`gtftools` (see Environments). The `dge` stage uses `envs/r.yml` instead.

`nucleation` submits a Slurm array and returns; run
`python3 processing_scripts/build_per_bin_covcorr.py` once those jobs finish.

---

## Expected `data/` layout

```
data/
├── 2024/
│   ├── cCREs/intersected_annos/
│   │   ├── FRIPs.tsv
│   │   └── FC_out/*_featureCounts_bulk_counts.tsv, *.summary.tsv
│   ├── ENCODE_chromHMM/{FRIPs.tsv, OverlapEnrichments.tsv}
│   ├── ENCODE_states/K562_E123_15_coreMarks_domains.bed
│   ├── K562_annotations/
│   │   ├── gencode.v44.basic.annotation.coding.gtf
│   │   └── processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}
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
    ├── mtf2_nucleation/{nuc_cpgonly_spread_bins.tsv.gz, per_bin_covcorr_full.tsv.gz}
    └── replitag_sphase_nucleation/matrices/
            encode_coremarks.1000bp.{rpkm,read_count}.matrix.tsv.gz
```

**External downloads** the cloner must fetch: GENCODE v27 and v44 GTFs; MSigDB
`GOBP_CELL_CYCLE.v2025.1.Hs.json`; UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; GEO
`GSE164804` (MTF2) and ENCODE `ENCFF587SWK` (EZH2) / `ENCFF065KGU` (SUZ12) fold-change
bigWigs.

---

## Known gaps

**Upstream of everything — read trimming and alignment.** `filter_sams*.sh` consume
duplicate-marked bowtie2 SAMs, which the authors produced with a lab-internal pipeline that
is not in this repository. `geo_to_sams.sh` can help process fastqs into sam files; the raw reads are in GEO `GSE327802`.

---

