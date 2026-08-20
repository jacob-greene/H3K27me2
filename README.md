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
├── paper.md                 methods excerpt / data-availability statement
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

### Which notebooks need no local data

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

`260803_WT_timeseries_clean_chromHMM_pub.ipynb` also writes
`data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv`, which three other notebooks read — but it
cannot reach that cell without the archive, so **that file ships in this repository** and
you do not need to run it first.

### Notebook inputs

| Notebook | Input | Source |
|---|---|---|
| **K27me3_states_111epigenomes** | Roadmap Epigenomics S3 bucket | downloaded at run time |
| **250214_growth_viability** | `data/250214_K27me123ac_P2S25p_TAZEEDUTX_cell_counts.xlsx` | **in repo** |
| | `data/20260206_drug_8only_aurora/*.csv` (15 flow exports) | **in repo** |
| **250411_peakPrint_clean** | `data/2024/SH_all_data/peakPRINT_50/MERGED_downsampled_bins.tsv` | Zenodo |
| | `.../min350/downsample_peakPRINT_slope1_min350_meta_v{1,2}.tsv`, `..._min300_meta_v3.tsv` | `downsample_peakPRINT.sh` |
| | `.../min350/*_reproducibility_95plus.bed`, `*_filtered.bed`, `K27me3_K_peak_coverage.bed` | `downsample_peakPRINT.sh` |
| | `.../reproducibility/FRiPs_merged_intervals_rep95.tsv` | Zenodo — **no producer** |
| | `data/2024/ENCODE_chromHMM/{FRIPs,OverlapEnrichments}.tsv` | Zenodo — **no producer** |
| | `data/2024/cCREs/intersected_annos/FRIPs.tsv` | Zenodo — **no producer** |
| **250911_RIFs_FC** | `data/2024/SH_all_data/sections4/*.txt` (BAM path lists) | Zenodo — **no producer** |
| | `data/2024/cCREs/intersected_annos/FC_out/*_bulk_counts.tsv`, `*.summary.tsv` | `Feature_counts_SAMs*.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | **in repo** |
| **251203_drug_timeseries** | `data/2024/RepliTag/Matrices/TAZEED/{mark}_drug_RUVr_series.tsv` | `251202_DGE_RUVseq_drug.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode.v44.basic.TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | `gtf2bed_TSS.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | **in repo** |
| **260801_nucleation** | `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz` | Zenodo — **no producer** |
| | `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz` | `build_per_bin_covcorr.py` |
| | `data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz` | `replitag_nucleation_generate_matrices.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | **in repo** |
| | `data/2024/MTF2_GSE164804/*.bw` | GEO `GSE164804`, ENCODE `ENCFF587SWK` / `ENCFF065KGU` |
| **260803_WT_timeseries** | `data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/{mark}/output.tsv` | `251113_DGE_RUVseq.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode...TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | `gtf2bed_TSS.sh` |
| | `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed` | 111epigenomes notebook |
| | `.../min350/top95/{K27me3,K27me2,K27me1,K27ac}_K_top95regions.bed` | peakPRINT notebook |
| | `data/2024/LAD_data/Dataset_S{1,2}_*.tsv` | **in repo** (published supplement) |
| | `data/2024/LAD_data/gencode.v27.annotation.gtf.gz` | GENCODE v27 |
| | `data/2024/RepliTag/refs/GOBP_CELL_CYCLE.v2025.1.Hs.json` | MSigDB |
| **bulk_RT** | `data/2024/RepliTag/Matrices/coverage_scores_bulk_10kb.tab` | `RIFs_RepliTag.sh` (10 kb block, currently commented out) |

The three `MTF2_GSE164804/*.bw` bigWigs are only needed by `build_per_bin_covcorr.py`; the
nucleation notebook re-checks its columns against them if `pyBigWig` and the files are both
present, and skips that check otherwise.

---

## Processing scripts

In `processing_scripts/`. **Run from the repository root**, except `gtf2bed_TSS.sh` as
noted. Most `module load` Fred Hutch Lmod modules (BEDTools, deepTools, SAMtools, Subread,
GNU parallel) and several submit Slurm jobs; substitute your own scheduler and module
system as needed.

### From GEO to SAMs

| Script | Produces | Requires |
|---|---|---|
| `geo_to_sams.sh` | `data/2024/RepliTag/geo_sams/<sample>.DTmarked.sam`, duplicate-marked, the input to `filter_sams.sh` | a two-column sample sheet (SRA run accession + sample name) and an hg38 bowtie2 index |

> **`geo_to_sams.sh` is provided for convenience and has never been run.** It is a
> transcription of the trimming and alignment methods in `paper.md` (cutadapt 4.4,
> Bowtie2 2.5.1) into runnable form, so that someone starting from GEO `GSE327802` can
> reach the input the rest of this pipeline expects. It is **not** the script the authors
> ran — the manuscript data were produced by a lab-internal pipeline that is not part of
> this repository. Check its output before trusting anything downstream of it. Every
> parameter the paper does not determine is marked in the source:
> `grep -n 'AMBIGUOUS (paper)' processing_scripts/geo_to_sams.sh`.

### Read filtering and coverage tracks

| Script | Produces | Requires |
|---|---|---|
| `blacklist_dm6.sh` | `data/Kc_merged_blacklist_hg38.bed` (**already shipped** — no need to re-run) | four *Drosophila*-on-hg38 BEDs, lab-internal and not redistributed |
| `filter_sams.sh` | `data/2024/RepliTag/filtered_sams_WT_1.0/` filtered SAMs, the combined featureCounts tables, the per-mark `_plusDUP_<mark>.tsv` split and the `.txt.summary.tsv` that `251113_DGE_RUVseq.R` reads | duplicate-marked bowtie2 SAMs (`geo_to_sams.sh`, or set `INPUT_DIR`), `data/Kc_merged_blacklist_hg38.bed`, `gencode.v44.basic.TSS-1000_TES_sorted.saf` (`gtf2bed_TSS.sh`) |
| `filter_sams_drug.sh` | same, into `filtered_sams_drug_1.0/` | as above |
| `filter_sams_drug_reseq.sh` | same, into `filtered_sams_drug_1.0/reseq/` | as above |
| `merge_sams2bam_slurm.sh` | `data/2024/SH_all_data/merged_bams4/` | `data/2024/SH_all_data/{All_sams4.txt,sections4/}` (Zenodo — **no producer**) |
| `bam2bed2bw.sh` | `filtered_sams_WT_1.0/*_{a,b}.bed`, merged bigWigs | filtered SAMs from `filter_sams.sh`; `data/chr_lens.txt`; `fraction_norm.csh` (in this directory) |
| `bam2bed2bw_drug.sh` | same for the drug series | as above |
| `fraction_norm.csh` | one normalized bigWig from one BED; called by `bam2bed2bw*.sh` | `bedtools`, `bedGraphToBigWig`, a chromosome-length file (`data/chr_lens.txt`) |

### Peak calling (peakPRINT)

| Script | Produces | Requires |
|---|---|---|
| `downsample_peakPRINT.sh` | `.../min350/*_merged_intervals.bed`, `*_reproducibility_95plus.bed`, the `_meta_v*.tsv` tables | bedGraphs from `bam2bed2bw*.sh`; `data/chr_lens.txt` |
| `bedtools_fingerprint_slope3.sh` | called by `downsample_peakPRINT.sh` | — |
| `compute_rank_threshold.py` | called by `bedtools_fingerprint_slope3.sh` | — |
| `Peak_heatmaps_v2.sh` | `figures/20250418_ALL_K_top95regions_10kb{,_profile}.pdf` | `top95/*_top95regions.bed` (peakPRINT notebook), `data/2024/SH_all_data/merged_bws4/` |

### Coverage matrices

| Script | Produces | Requires |
|---|---|---|
| `RIFs_RepliTag.sh` | `data/2024/RepliTag/Matrices/coverage_scores_gencode.v44.basic.TSS-1000_TES_sorted_bulk.tab`; and, if you uncomment the final block, `coverage_scores_bulk_10kb.tab` | `data/2024/RepliTag/bws_bulk/*.bw` (Zenodo — **no producer**), `gencode.v44.basic.TSS-1000_TES_sorted.bed` (`gtf2bed_TSS.sh`) |
| `Feature_counts_SAMs.sh` | `data/2024/cCREs/intersected_annos/FC_out/hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.tsv` | `sections4/` BAM lists (Zenodo), `gencode.v44.basic.TSS-1000_TES_sorted.saf` (`gtf2bed_TSS.sh`) |
| `Feature_counts_SAMs-allFeatures_k562.sh` | one `*_featureCounts_bulk_counts.tsv` + `.summary.tsv` per SAF, in the same directory | `data/2024/LAD_data/*.saf` (Zenodo), `sections4/` (Zenodo) |
| `replitag_nucleation_generate_matrices.sh` | `data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.{rpkm,read_count}.matrix.tsv.gz`, `config/samples.tsv` | `filtered_sams_WT_1.0/*_d2pcr_{a,b}.bed` (`bam2bed2bw.sh`), the hg38 Repli-seq bigWig (`crossmaphg19tohg38.sh`), `K562_E123_15_coreMarks_domains.bed` (111epigenomes notebook) |
| `run_replitag_nucleation_matrices_slurm.sh` | Slurm driver for the above | — |
| `build_per_bin_covcorr.py` | `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz` | `encode_coremarks.1000bp.rpkm.matrix.tsv.gz` (above) and the three `data/2024/MTF2_GSE164804/*.bw` bigWigs (GEO `GSE164804`, ENCODE `ENCFF587SWK` / `ENCFF065KGU`). Needs `numpy`, `pandas`, `pyBigWig`. |

`build_per_bin_covcorr.py` builds the one matrix `260801_nucleation_pub.ipynb` needs that
nothing else here produces — the four columns it reads (`bin_id`, `MTF2`, `SUZ12`, `EZH2`),
each a bigWig mean over a 1-kb bin. The three bigWigs are publicly obtainable; the bin
universe is not, so **this script does not run on a bare clone**. It runs once you have
unpacked the Zenodo archive, or once you have run
`replitag_nucleation_generate_matrices.sh` yourself.

```bash
python3 processing_scripts/build_per_bin_covcorr.py
```

### Differential expression

| Script | Produces | Requires |
|---|---|---|
| `251113_DGE_RUVseq.R` | `filtered_sams_WT_1.0/251113_DGE/<mark>/output.tsv` | the per-mark `_plusDUP_<mark>.tsv` tables and `_plusDUP.txt.summary.tsv`, both written by `filter_sams.sh` |
| `251202_DGE_RUVseq_drug.R` | `data/2024/RepliTag/Matrices/TAZEED/<mark>_drug_RUVr_series.tsv` | `_plusDUP.tsv` and `_plusDUP.txt.summary.tsv` for both the main and `reseq/` branches, written by `filter_sams_drug{,_reseq}.sh` |

### Reference preparation

| Script | Produces | Requires |
|---|---|---|
| `gtf2bed_TSS.sh` | `processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}` | `gtftools` on `PATH` and `gencode.v44.basic.annotation.coding.gtf` (GENCODE v44, protein-coding subset). **Uses relative paths — run it from `data/2024/K562_annotations/`, with `processed_beds/sorted/` already created.** |
| `crossmaphg19tohg38.sh` | `data/2024/RepliTag/RepliSeq/hg38_bws_crossmap/wgEncodeUwRepliSeqK562WaveSignalRep1_hg38.bigWig.bw` | UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; downloads the liftover chain itself |

---

## Expected `data/` layout

On top of what ships (see "What ships in `data/`" above), the Zenodo archive must supply:

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

Having a producer is **not** the same as being producible — the producer must also be
reachable. Everything not marked above must be in the archive.

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
is not in this repository. `geo_to_sams.sh` (above) is an untested reconstruction of it
from the parameters in `paper.md`; the raw reads are in GEO `GSE327802`.

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

## Open questions for the authors

Three things need an author's decision, and are recorded rather than guessed at.

- **Fig 3d–e cannot be regenerated here.** The cell in `260803_WT…ipynb` that writes
  `20260319_Fig3de_Native_vAutonomous_me3domain_me2_byRTbin.pdf` reads
  `per_gene_max_me2_cats`, which is assigned nowhere in any notebook — evidently defined in
  a cell dropped when the `_clean_pub` derivative was cut, so the cell raises `NameError`.
  It is left in place because it is the record of a manuscript panel. Restore the defining
  cell, or drop the panel.
- **Which *Drosophila* line the spike-in came from.** `paper.md` says **S2**; the four
  inputs to `blacklist_dm6.sh` are named for **Kc** and `BT_Dm`. The method and the shipped
  blacklist are unaffected either way, but manuscript and code disagree.
- **There is no `LICENSE`.** Without one, default copyright applies and readers have no
  permission to reuse the code. `THIRD_PARTY_NOTICES.md` covers only the two redistributed
  supplement files.

---

## Notes

- **If you build the Zenodo/supplementary archive, build it with `git archive` — not by
  zipping the working directory.** `.ipynb_checkpoints/` is gitignored and so is absent
  from any git-based export, but the checkpoint copies on disk are *pre-rewrite*: their
  source cells still carry lab-internal absolute paths. Zipping the working tree would
  publish exactly the paths the path rewrite removed.
- **Stored notebook outputs still contain lab-internal absolute paths**, of two kinds:
  `/fh/fast` paths inside rendered HTML tables and printed save confirmations, and
  `/loc/scratch/<slurm-jobid>/ipykernel_<pid>/` traceback frames inside `stderr` warning
  text. The *source* is clean in both cases; these are execution outputs. Clearing them
  changes the published record of what was run, so it is left as a deliberate decision
  rather than done silently.
