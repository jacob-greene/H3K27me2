# H3K27me2

Analysis code for the manuscript *"Histone H3K27 methylation states are sequentially
catalyzed in cycling cells"* by Greene, Ahmad, and Henikoff, 2026.

## Data availability

| Source | Accession | Contents |
|---|---|---|
| GEO | [`GSE327802`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327802) | Raw CUT&Tag sequencing data |
| Zenodo | [`10.5281/zenodo.19489947`](https://doi.org/10.5281/zenodo.19489947) | Processed data (`data/` folder contents) |
| UCSC | [browser session](https://genome.ucsc.edu/s/jegreene/H3K27me2_manuscript) | Interactive tracks |

---

## Quick start

```bash
git clone https://github.com/jacob-greene/H3K27me2.git
cd H3K27me2

# Download the processed-data archive from Zenodo and unpack it so that its
# contents land inside the data/ directory that ships with this repo.
# The result must look like the layout in "Expected data/ layout" below.
```

**Every path in this repository is relative to the repository root.** Start Jupyter from
the repository root, and run the shell/R scripts from the repository root
(`bash "processing scripts/RIFs_RepliTag.sh"`).

> **One exception, and it currently blocks five notebooks.** Cell 2 of
> `250411_peakPrint_clean_pub.ipynb`, `250911_RIFs_FC_pub.ipynb`,
> `251203_drug_timeseries_pub.ipynb`, `260803_WT_timeseries_clean_chromHMM_pub.ipynb`
> and `bulk_RT_pub.ipynb` does `sys.path.append(os.path.expanduser("~"))` and then
> `import color_palettes`. That module is **not in this repository** — see Known gaps.

---

## Repository contents

```
H3K27me2/
├── *_pub.ipynb              8 publication notebooks (the figures)
├── processing scripts/      19 scripts + 1 notebook that build the processed data
├── data/                    input data — mostly downloaded from Zenodo (see below)
├── figures/                 created by the notebooks; PDFs and figure-source TSVs land here
└── paper.md                 methods excerpt / data-availability statement
```

### Notebook run order

Three notebooks are **producers** as well as consumers, so order matters:

1. **`K27me3_states_111epigenomes_pub.ipynb`** — fully self-contained; downloads Roadmap
   Epigenomics core-marks segmentations from S3 at run time. Writes
   `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed`.
   *Requires `bedtools` on `PATH`; set `BEDTOOLS_BIN` to prepend a specific directory.*
2. **`250411_peakPrint_clean_pub.ipynb`** — writes
   `data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/min350/top95/*_top95regions.bed`.
3. **`260803_WT_timeseries_clean_chromHMM_pub.ipynb`** — writes
   `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv`, which three other notebooks read.
   **It cannot reach that cell on a clean clone**, so the file must come from the Zenodo
   archive regardless of run order — see "`JG_coding_genes_maxRPKM.tsv` has a producer it
   cannot run" under Known gaps.

Then, in any order:
`250214_growth_viability_cellcycle_pub.ipynb` (independent of all the above),
`250911_RIFs_FC_pub.ipynb`, `251203_drug_timeseries_pub.ipynb`,
`260801_nucleation_pub.ipynb`, `bulk_RT_pub.ipynb`.

---

## Notebook inputs

`250214_growth_viability_cellcycle_pub.ipynb` is the only notebook that runs on what
ships in this repository alone.

| Notebook | Input | Source |
|---|---|---|
| **250214_growth_viability** | `data/250214_K27me123ac_P2S25p_TAZEEDUTX_cell_counts.xlsx` | **in repo** |
| | `data/20260206_drug_8only_aurora/*.csv` (15 flow exports) | **in repo** |
| **K27me3_states_111epigenomes** | Roadmap Epigenomics S3 bucket | downloaded at run time |
| **250411_peakPrint_clean** | `data/2024/SH_all_data/peakPRINT_50/MERGED_downsampled_bins.tsv` | Zenodo |
| | `.../min350/downsample_peakPRINT_slope1_min350_meta_v{1,2}.tsv`, `..._min300_meta_v3.tsv` | `downsample_peakPRINT.sh` |
| | `.../min350/*_reproducibility_95plus.bed`, `*_filtered.bed`, `K27me3_K_peak_coverage.bed` | `downsample_peakPRINT.sh` |
| | `.../reproducibility/FRiPs_merged_intervals_rep95.tsv` | **blocked** — see below |
| | `data/2024/ENCODE_chromHMM/{FRIPs,OverlapEnrichments}.tsv` | **blocked** |
| | `data/2024/cCREs/intersected_annos/FRIPs.tsv` | **blocked** |
| **250911_RIFs_FC** | `data/2024/SH_all_data/sections4/*.txt` (BAM path lists) | **blocked** |
| | `data/2024/cCREs/intersected_annos/FC_out/*_bulk_counts.tsv`, `*.summary.tsv` | `Feature_counts_SAMs-allFeatures_k562.sh` |
| | `.../hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.tsv` | `Feature_counts_SAMs.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | Zenodo — WT notebook writes it but cannot run |
| **251203_drug_timeseries** | `data/2024/RepliTag/Matrices/TAZEED/{mark}_drug_RUVr_series.tsv` | `251202_DGE_RUVseq_drug.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode.v44.basic.TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | **blocked** |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | Zenodo — WT notebook writes it but cannot run |
| **260801_nucleation** | `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz` | **blocked** |
| | `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz` | `260728_mtf2_correlation_nucleation_JG.ipynb` |
| | `data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz` | `replitag_nucleation_generate_matrices.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | Zenodo — WT notebook writes it but cannot run |
| | `data/2024/MTF2_GSE164804/*.bw` | GEO `GSE164804`, ENCODE `ENCFF587SWK`/`ENCFF065KGU` — *optional*, skipped if `pyBigWig` is absent |
| **260803_WT_timeseries** | `data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/{mark}/output.tsv` | `251113_DGE_RUVseq.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode...TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | **blocked** |
| | `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed` | 111epigenomes notebook |
| | `.../min350/top95/{K27me3,K27me2,K27me1,K27ac}_K_top95regions.bed` | peakPRINT notebook |
| | `data/2024/LAD_data/Dataset_S1_Promoter_SuRE_Classification.tsv` | published supplement |
| | `data/2024/LAD_data/Dataset_S2_TRIP.tsv` | published supplement |
| | `data/2024/LAD_data/gencode.v27.annotation.gtf.gz` | GENCODE v27 |
| | `data/2024/RepliTag/refs/GOBP_CELL_CYCLE.v2025.1.Hs.json` | MSigDB |
| **bulk_RT** | `data/2024/RepliTag/Matrices/coverage_scores_bulk_10kb.tab` | `RIFs_RepliTag.sh` (10 kb block, currently commented out) |

---

## Processing scripts

In `processing scripts/`. Run from the repository root. Most `module load` Fred Hutch
Lmod modules (BEDTools, deepTools, SAMtools, Subread, GNU parallel, Subread) and several
submit Slurm jobs; substitute your own scheduler and module system as needed.

### Read filtering and coverage tracks

| Script | Produces | Requires |
|---|---|---|
| `blacklist_dm6.sh` | `data/Kc_merged_blacklist_hg38.bed` (**already shipped** — no need to re-run) | four *Drosophila*-on-hg38 BEDs under `/shared/ngs/...` (**blocked**) |
| `filter_sams.sh` | `data/2024/RepliTag/filtered_sams_WT_1.0/` filtered SAMs + `hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_{noDUP,plusDUP}.txt` | duplicate-marked bowtie2 SAMs (**blocked**), `data/Kc_merged_blacklist_hg38.bed`, `gencode.v44.basic.TSS-1000_TES_sorted.saf` (**blocked**) |
| `filter_sams_drug.sh` | same, into `filtered_sams_drug_1.0/` | as above |
| `filter_sams_drug_reseq.sh` | same, into `filtered_sams_drug_1.0/reseq/` | as above |
| `merge_sams2bam_slurm.sh` | `data/2024/SH_all_data/merged_bams4/` | `data/2024/SH_all_data/{All_sams4.txt,sections4/}` (**blocked**) |
| `bam2bed2bw.sh` | `filtered_sams_WT_1.0/*_{a,b}.bed`, merged bigWigs | `chr_lens.txt`, `fraction_norm.csh` (**blocked**) |
| `bam2bed2bw_drug.sh` | same for the drug series | as above |

### Peak calling (peakPRINT)

| Script | Produces | Requires |
|---|---|---|
| `downsample_peakPRINT.sh` | `.../min350/*_merged_intervals.bed`, `*_reproducibility_95plus.bed`, the `_meta_v*.tsv` tables | bedGraphs from `bam2bed2bw*.sh`; `chr_lens.txt` (**blocked**) |
| `bedtools_fingerprint_slope3.sh` | called by `downsample_peakPRINT.sh` | — |
| `compute_rank_threshold.py` | called by `bedtools_fingerprint_slope3.sh` | — |
| `Peak_heatmaps_v2.sh` | `figures/20250418_ALL_K_top95regions_10kb{,_profile}.pdf` | `top95/*_top95regions.bed` (peakPRINT notebook), `data/2024/SH_all_data/merged_bws4/` |

### Coverage matrices

| Script | Produces | Requires |
|---|---|---|
| `RIFs_RepliTag.sh` | `data/2024/RepliTag/Matrices/coverage_scores_gencode.v44.basic.TSS-1000_TES_sorted_bulk.tab`; and, if you uncomment the final block, `coverage_scores_bulk_10kb.tab` | `data/2024/RepliTag/bws_bulk/*.bw` (**blocked** — see below), `gencode.v44.basic.TSS-1000_TES_sorted.bed` (**blocked**) |
| `Feature_counts_SAMs.sh` | `hg38.gencodev44_HKP_genes-1000_featureCounts_bulk_counts.tsv` | `sections4/` BAM lists (**blocked**), `gencode.v44.basic.TSS-1000_TES_sorted.saf` (**blocked**) |
| `Feature_counts_SAMs-allFeatures_k562.sh` | one `*_featureCounts_bulk_counts.tsv` + `.summary.tsv` per SAF | `data/2024/LAD_data/*.saf` (**blocked**), `sections4/` (**blocked**) |
| `replitag_nucleation_generate_matrices.sh` | `data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz`, `config/samples.tsv` | `filtered_sams_WT_1.0/*_d2pcr_{a,b}.bed`, the hg38 Repli-seq bigWig, `K562_E123_15_coreMarks_domains.bed` |
| `run_replitag_nucleation_matrices_slurm.sh` | Slurm driver for the above | — |

### Nucleation coverage/correlation

| Script | Produces | Requires |
|---|---|---|
| `260728_mtf2_correlation_nucleation_JG.ipynb` | `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz` (section 1) | **Section 1 cannot run on a clean clone.** It reads `data/2024/RepliTag/bws_bulk/K27me3_K_no_chrM.1-1000.bw` (**blocked**, no producer), plus `encode_coremarks.1000bp.{rpkm,read_count}.matrix.tsv.gz` from `replitag_nucleation_generate_matrices.sh`, `data/2024/MTF2_GSE164804/{MTF2_shCT_hg38_coverage.bw, ENCODE_SUZ12_…ENCFF065KGU…bw, ENCODE_EZH2_…ENCFF587SWK…bw, MTF2_macs3_hg38_peaks.bed}`, `data/2024/K562_annotations/EncodeRegDnaseUwK562Peak.bed` and `.../processed_beds/hg38_CpG_Island.bed`. So `per_bin_covcorr_full.tsv.gz` **must ship in the archive** — having a producer and being producible are different things. Later sections additionally read `data/2024/JG_Hs_H1_K562_240110/`, `encode_coremarks_repliseq.tsv.gz`, `per_bin_input_mtf2mean0.tsv.gz` and `per_bin_specificity_tracks.tsv.gz` (all **blocked**; the last is a cache whose extract-from-bigWig fallback needs the two blocked directories, so it too must ship). |

### Differential expression

| Script | Produces | Requires |
|---|---|---|
| `251113_DGE_RUVseq.R` | `filtered_sams_WT_1.0/251113_DGE/<mark>/output.tsv` | `…_plusDUP_<mark>.**tsv**` per-mark tables and `…_plusDUP.txt.summary.**tsv**` — **both blocked**; `filter_sams.sh` writes one combined `.txt` plus featureCounts' own `.txt.summary`, so both the per-mark split *and* the `.txt.summary → .summary.tsv` rename are unproduced |
| `251202_DGE_RUVseq_drug.R` | `data/2024/RepliTag/Matrices/TAZEED/<mark>_drug_RUVr_series.tsv` | `filtered_sams_drug_1.0/…_plusDUP.**tsv**` and `…_plusDUP.txt.summary.**tsv**` for both the main and `reseq/` branches (**blocked** — `filter_sams_drug*.sh` write `.txt`, not `.tsv`), `data/2024/RepliTag/featureCounts_drug_dm6/...summary` (**blocked**) |

### Reference preparation

| Script | Produces | Requires |
|---|---|---|
| `crossmaphg19tohg38.sh` | `data/2024/RepliTag/RepliSeq/hg38_bws_crossmap/wgEncodeUwRepliSeqK562WaveSignalRep1_hg38.bigWig.bw` | UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; downloads the liftover chain itself |

---

## Expected `data/` layout

What ships in the repository:

```
data/
├── 250214_K27me123ac_P2S25p_TAZEEDUTX_cell_counts.xlsx
├── 20260206_drug_8only_aurora/           15 Aurora flow-cytometry CSV exports
└── Kc_merged_blacklist_hg38.bed          Drosophila spike-in blacklist (286,386 intervals)
```

What the Zenodo archive must add:

```
data/
├── 2024/
│   ├── cCREs/intersected_annos/
│   │   ├── FRIPs.tsv
│   │   └── FC_out/*_featureCounts_bulk_counts.tsv, *.summary.tsv
│   ├── ENCODE_chromHMM/{FRIPs.tsv, OverlapEnrichments.tsv}
│   ├── ENCODE_states/K562_E123_15_coreMarks_domains.bed*
│   ├── K562_annotations/
│   │   ├── gencode.v44.basic.annotation.gtf
│   │   ├── EncodeRegDnaseUwK562Peak.bed
│   │   ├── processed_beds/hg38_CpG_Island.bed
│   │   └── processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}
│   ├── LAD_data/
│   │   ├── Dataset_S1_Promoter_SuRE_Classification.tsv
│   │   ├── Dataset_S2_TRIP.tsv
│   │   ├── JG_coding_genes_maxRPKM.tsv   REQUIRED — its producer cannot run on a
│   │   │                                 clean clone; see Known gaps
│   │   ├── gencode.v27.annotation.gtf.gz
│   │   └── *.saf
│   ├── JG_Hs_H1_K562_240110/            (260728_mtf2_correlation… later sections, incl.
│   │                                     the cell-8 specificity fallback)
│   │                                    **blocked** — no producer
│   ├── MTF2_GSE164804/
│   │   ├── *.bw                        MTF2_shCT + ENCFF587SWK/ENCFF065KGU fold-change
│   │   ├── MTF2_macs3_hg38_peaks.bed   MACS3 peak calls — **blocked**, no producing script
│   │   └── ENCODE_{EZH2_K562_ENCFF080JPV,SUZ12_K562_ENCFF856HYC}_hg38.bed.gz
│   ├── RepliSeq/*.bigWig                 UCSC-ENCODE UW Repli-seq K562 **hg19** —
│   │                                     the INPUT to crossmaphg19tohg38.sh
│   ├── RepliTag/
│   │   ├── bws_bulk/*.bw
│   │   ├── featureCounts_drug_dm6/*.summary
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
    ├── mtf2_nucleation/{nuc_cpgonly_spread_bins.tsv.gz, per_bin_covcorr_full.tsv.gz,
    │                    per_bin_specificity_tracks.tsv.gz}
    └── replitag_sphase_nucleation/matrices/
            encode_coremarks.1000bp.{rpkm,read_count}.matrix.tsv.gz
```

`*` = written by a notebook in this repository, **and that notebook is verified to run on
a clean clone**; not required in the archive if you run the notebooks in the order given
above. There is exactly one such file, `K562_E123_15_coreMarks_domains.bed`: its producer
(`K27me3_states_111epigenomes_pub.ipynb` cell 9) imports no `color_palettes` and reads only
Roadmap Epigenomics S3 URLs, so it has no blocked input.

Having a producer is **not** sufficient for a star — the producer must also be reachable.
`JG_coding_genes_maxRPKM.tsv` carried a star until it was checked and does not: see below.

**Everything else in the block above must be in the archive.**

> **Packaging the code archive:** use `git archive`, not a zip of the working directory.
> See the last bullet under Notes — the untracked `.ipynb_checkpoints/` copies still contain
> the lab-internal absolute paths that the tracked notebooks no longer do.

---

## Known gaps

These inputs have no producing script in this repository and no public source. They must
be added to the Zenodo archive, or the corresponding script must be supplied.

**Upstream of everything — read trimming and alignment.** The `filter_sams*.sh` scripts
consume duplicate-marked bowtie2 SAMs at
`/shared/ngs/illumina/henikoff/{241106,250411,250910,250915}_bowtie2/JG_HsDm/marked/`.
No trimming or alignment script is included. The parameters are given in `paper.md`
(cutadapt 4.4, Bowtie2 2.5.1) and the raw reads are in GEO `GSE327802`.

**Lab-internal helpers**, referenced by `bam2bed2bw*.sh` and `downsample_peakPRINT.sh`:
`/shared/ngs/illumina/henikoff/databases/human/hg38/chr_lens.txt` and
`/shared/ngs/illumina/henikoff/bin/aligners/fraction_norm.csh`.

**Annotation derivations.** `gencode.v44.basic.TSS-1000_TES_sorted.bed` / `.saf` and the
`data/2024/LAD_data/*.saf` files are derived from GENCODE, but the derivation script is
not in this repository.

**`color_palettes.py` — a hard blocker, and the first thing a cloner will hit.**
Five of the eight notebooks import it at cell 2 from the user's home directory
(`sys.path.append(os.path.expanduser("~"))`). The module is not in this repository and
is not obtainable publicly, so those notebooks stop with `ModuleNotFoundError` at their
second cell — including `250411_peakPrint_clean_pub.ipynb`, which is **second in the
mandatory run order** and produces the `top95regions.bed` files two other notebooks need.
It defines the manuscript's figure colours, so it cannot be substituted without silently
changing published figures. It needs to be shipped with the repository.

**`JG_coding_genes_maxRPKM.tsv` has a producer it cannot run.** `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv`
is written by `260803_WT_timeseries_clean_chromHMM_pub.ipynb` at cell 26, and read by
`250911_RIFs_FC_pub.ipynb` (c16), `251203_drug_timeseries_pub.ipynb` (c12) and
`260801_nucleation_pub.ipynb` (c1). Cell 26 writes `rt_data`, which is built at cell 22 from
`bins_idx` — assigned at cell 12 from the **blocked**
`data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed`.
Two further blockers precede it: cell 2 imports `color_palettes` (above) and cell 10 reads
`data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/{mark}/output.tsv`, whose producer
`251113_DGE_RUVseq.R` is itself blocked. **So this file must ship in the Zenodo archive**; it
is not optional, and omitting it breaks three of the eight notebooks. (This is the same
distinction as `per_bin_covcorr_full.tsv.gz` in the merged notebook: having a producer and
being producible are different claims, and the archive manifest turns on the second.)

**Intermediate tables with no producer.**

- `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz`
- `data/2024/RepliTag/filtered_sams_WT_1.0/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP_<mark>.tsv`
  (the per-mark split of the combined `filter_sams.sh` table)
- `data/2024/RepliTag/featureCounts_drug_dm6/genes_v5_KABK_230607_featureCounts_counts_plusDUP.txt.summary`
- `data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/reproducibility/FRiPs_merged_intervals_rep95.tsv`
- `data/2024/SH_all_data/peakPRINT_50/MERGED_downsampled_bins.tsv`
- `data/2024/SH_all_data/sections4/*.txt` and `All_sams4.txt` (BAM path manifests)
- `data/2024/ENCODE_chromHMM/{FRIPs,OverlapEnrichments}.tsv` and
  `data/2024/cCREs/intersected_annos/FRIPs.tsv` (ChromHMM `OverlapEnrichment` output)
- `data/2024/RepliTag/filtered_sams_{WT,drug}_1.0/…_plusDUP.**tsv**` and
  `…_plusDUP.txt.summary.**tsv**` — `filter_sams*.sh` write `…_plusDUP.txt` and
  featureCounts' own `…_plusDUP.txt.summary`. Both DGE R scripts read `.tsv` versions of
  each. The `.txt → .tsv` conversion is unproduced on **both** the WT and drug branches.
- `data/2024/RepliTag/bws_bulk/` — assembled by hand from the `bam2bed2bw.sh` bigWigs plus
  the lifted Repli-seq track; no script does the assembly
- `data/2024/JG_Hs_H1_K562_240110/` — read by section 4 of
  `260728_mtf2_correlation_nucleation_JG.ipynb`
- `data/results/mtf2_nucleation/encode_coremarks_repliseq.tsv.gz` and
  `per_bin_input_mtf2mean0.tsv.gz` — read by that notebook's later sections, written
  nowhere in the corpus
- `data/results/mtf2_nucleation/per_bin_coverage_correlation.tsv.gz` — the optional
  coverage cache its section 1 reads if present; no producer, so the expensive bigWig
  path always runs
- `data/results/mtf2_nucleation/per_bin_specificity_tracks.tsv.gz` — the same shape, in the
  correction section (cell 8, lines 37–42): `if CACHE.exists()` read it, `else` extract from
  bigWigs. There is no producer outside that `else`, so on a clean clone the fallback always
  runs — and the fallback reads `data/2024/RepliTag/bws_bulk/{K27ac,K27me2}_K_no_chrM.1-1000.bw`
  and `data/2024/JG_Hs_H1_K562_240110/JG_Hs_K562_{K4me3,K9me3}_240110_25k_d2pcr_a.1-1000.bw`,
  **both directories blocked and listed above**. It must ship in the archive.
- `data/2024/MTF2_GSE164804/MTF2_macs3_hg38_peaks.bed` — MACS3 peak calls; no producing
  script, and not covered by the GEO/ENCODE bigWig downloads
- `data/2024/K562_annotations/EncodeRegDnaseUwK562Peak.bed` and
  `processed_beds/hg38_CpG_Island.bed`

**External downloads** the cloner must fetch: GENCODE v27 and v44 GTFs; MSigDB
`GOBP_CELL_CYCLE.v2025.1.Hs.json`; the published supplements `Dataset_S1` and
`Dataset_S2`; UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; GEO `GSE164804` (MTF2) and
ENCODE `ENCFF587SWK` (EZH2) / `ENCFF065KGU` (SUZ12) fold-change bigWigs, and
`ENCFF080JPV` (EZH2) / `ENCFF856HYC` (SUZ12) hg38 peak BEDs for the merged
`260728_mtf2_correlation_nucleation_JG.ipynb`.

---

## Notes

- **`FC_out` directory mismatch, both scripts.** `Feature_counts_SAMs-allFeatures_k562.sh`
  writes to `data/2024/LAD_data/FC_out/` and `Feature_counts_SAMs.sh` writes into
  `data/2024/SH_all_data/`, but `250911_RIFs_FC_pub.ipynb` globs both filename patterns
  inside `data/2024/cCREs/intersected_annos/FC_out/`. Filenames match exactly, so the
  links are certain; the outputs were evidently relocated after the run. The expected
  layout above sides with the notebook — **move both scripts' outputs there**, or repoint
  the scripts.
- **`H3K27_TRIP.tsv` is an output, not an input.** `260803_WT…` c063 *writes*
  `data/2024/LAD_data/H3K27_TRIP.tsv`; nothing reads it. It is not required in the archive.
- **`K562_IHEC_18_domains.bed` is not required.** Its only occurrence is a commented-out
  line in `260803_WT…` c017; the live line uses the E123 core-marks file.
- **Stored notebook outputs still contain 171 lab-internal absolute paths**, of two kinds.
  **82 `/fh/fast` paths** — 80 in `250411_peakPrint_clean_pub.ipynb` (cells 10/11/12,
  rendered HTML tables whose cell values are input filenames) and 2 in
  `251203_drug_timeseries_pub.ipynb` (cell 18, printed save confirmations naming a
  `setty_m` collaboration share). **89 `/loc/scratch/<slurm-jobid>/ipykernel_<pid>/<n>.py`
  traceback frames**, all inside `stderr` warning text, spread across five notebooks
  (`250411_peakPrint` 14, `250911_RIFs_FC` 6, `251203_drug` 38, `260803_WT` 30,
  `K27me3_states_111epigenomes` 1); between them they expose seven Slurm job IDs and nine
  kernel PIDs. Source is clean in both cases; these are execution *outputs*. Clearing them
  changes the published record of what was run, so it is left as a deliberate decision
  rather than done silently. Counts are over tracked files only, as literal
  `/fh/fast` and `/loc/scratch` substrings in each cell's `outputs` array.
- **If you build the Zenodo/supplementary archive, build it with `git archive` — not by
  zipping the working directory.** `.ipynb_checkpoints/` is gitignored (`.gitignore:1`) and
  so is absent from any git-based export, but the checkpoint copies on disk are *pre-rewrite*:
  their **source** cells still carry 44 lab-internal absolute paths (42 `/fh/fast`, 2
  `/home/jgreene3/…micromamba…`), on top of 82 in their stored outputs. Zipping the
  working tree would publish exactly the paths the path rewrite removed.
- **Two undefined names remain in `260803_WT_timeseries_clean_chromHMM_pub.ipynb`** and
  need an author decision, so they are left as-is rather than guessed at:
  `per_gene_max_me2_cats` (c040, used at L56 and L81 — c040's own header comment lists it
  under "Assumes these already exist", and the other three names in that list *are*
  defined) and `gl_df` (c053 L1). Neither is assigned anywhere in the repository. Both
  look like leftovers from cells dropped when the `_clean_pub` derivative was cut: either
  the defining cells should be restored, or the panels deleted. A related break in
  `250411_peakPrint_clean_pub.ipynb` (`ScalarMappable`/`Normalize` used unqualified at
  c026/c027/c029/c030) *was* fixed, since the missing imports were unambiguous.
- `260803_WT…` c074 writes `RNAPII_cellCycle_0.0001.csv` — a data product — into
  `figures/Fig1/`, following this project's convention of writing figure-source tables
  beside the figures. Nothing reads it, so the location is cosmetic.
- `processing scripts/260728_mtf2_correlation_nucleation_JG.ipynb` was copied
  byte-identical and then re-serialized with `ensure_ascii=False`, so a line-level `diff`
  against the lab-internal original shows many changed lines (unicode escapes became
  literal UTF-8 in the markdown cells) on top of the four intended path repoints. Only
  the four `source` edits are semantic: outputs, execution counts, cell ids, cell types
  and all metadata are unchanged.
- The `processing scripts/` directory name contains a space, so it must be quoted in
  shell commands.
