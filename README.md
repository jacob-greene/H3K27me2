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
(`bash "processing scripts/RIFs_RepliTag.sh"`). Nothing needs to be edited if you do.

---

## Repository contents

```
H3K27me2/
├── *_pub.ipynb              8 publication notebooks (the figures)
├── processing scripts/      19 scripts that build the processed data
├── data/                    input data — mostly downloaded from Zenodo (see below)
├── figures/                 created by the notebooks; PDFs land here
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
| | `data/20260206_drug_8only_aurora/*.csv` (16 flow exports) | **in repo** |
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
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | WT notebook (run it first) |
| **251203_drug_timeseries** | `data/2024/RepliTag/Matrices/TAZEED/{mark}_drug_RUVr_series.tsv` | `251202_DGE_RUVseq_drug.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode.v44.basic.TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | **blocked** |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | WT notebook |
| **260801_nucleation** | `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz` | **blocked** |
| | `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz` | **blocked** |
| | `data/results/replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz` | `replitag_nucleation_generate_matrices.sh` |
| | `data/2024/LAD_data/JG_coding_genes_maxRPKM.tsv` | WT notebook |
| | `data/2024/MTF2_GSE164804/*.bw` | GEO `GSE164804`, ENCODE `ENCFF587SWK`/`ENCFF065KGU` — *optional*, skipped if `pyBigWig` is absent |
| **260803_WT_timeseries** | `data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE/{mark}/output.tsv` | `251113_DGE_RUVseq.R` |
| | `data/2024/RepliTag/Matrices/coverage_scores_gencode...TSS-1000_TES_sorted_bulk.tab` | `RIFs_RepliTag.sh` |
| | `data/2024/K562_annotations/processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.bed` | **blocked** |
| | `data/2024/ENCODE_states/K562_E123_15_coreMarks_domains.bed` | 111epigenomes notebook |
| | `data/2024/ENCODE_states/K562_IHEC_18_domains.bed` | IHEC download — **needs a URL** |
| | `.../min350/top95/{K27me3,K27me2,K27me1,K27ac}_K_top95regions.bed` | peakPRINT notebook |
| | `data/2024/LAD_data/Dataset_S1_Promoter_SuRE_Classification.tsv` | published supplement |
| | `data/2024/LAD_data/Dataset_S2_TRIP.tsv` | published supplement |
| | `data/2024/LAD_data/H3K27_TRIP.tsv` | **blocked** |
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

### Differential expression

| Script | Produces | Requires |
|---|---|---|
| `251113_DGE_RUVseq.R` | `filtered_sams_WT_1.0/251113_DGE/<mark>/output.tsv` | `..._plusDUP_<mark>.tsv` per-mark tables (**blocked**) and `..._plusDUP.txt.summary.tsv` |
| `251202_DGE_RUVseq_drug.R` | `data/2024/RepliTag/Matrices/TAZEED/<mark>_drug_RUVr_series.tsv` | `filtered_sams_drug_1.0/`, `data/2024/RepliTag/featureCounts_drug_dm6/...summary` (**blocked**) |

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
├── 20260206_drug_8only_aurora/           16 Aurora flow-cytometry CSV exports
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
│   ├── ENCODE_states/{K562_IHEC_18_domains.bed, K562_E123_15_coreMarks_domains.bed*}
│   ├── K562_annotations/
│   │   ├── gencode.v44.basic.annotation.gtf
│   │   └── processed_beds/sorted/gencode.v44.basic.TSS-1000_TES_sorted.{bed,saf}
│   ├── LAD_data/
│   │   ├── Dataset_S1_Promoter_SuRE_Classification.tsv
│   │   ├── Dataset_S2_TRIP.tsv
│   │   ├── H3K27_TRIP.tsv
│   │   ├── JG_coding_genes_maxRPKM.tsv*
│   │   ├── gencode.v27.annotation.gtf.gz
│   │   └── *.saf
│   ├── MTF2_GSE164804/*.bw
│   ├── RepliTag/
│   │   ├── bws_bulk/*.bw
│   │   ├── featureCounts_drug_dm6/*.summary
│   │   ├── filtered_sams_WT_1.0/, filtered_sams_drug_1.0/
│   │   ├── Matrices/{coverage_scores_*.tab, TAZEED/}
│   │   ├── refs/GOBP_CELL_CYCLE.v2025.1.Hs.json
│   │   └── RepliSeq/hg38_bws_crossmap/
│   └── SH_all_data/
│       ├── sections4/, merged_bams4/, merged_bws4/
│       ├── peakPRINT_50/MERGED_downsampled_bins.tsv
│       └── downsample_peakPRINT_slope1_min300/
│           ├── reproducibility/FRiPs_merged_intervals_rep95.tsv
│           └── min350/{*.bed, *_meta_v*.tsv, top95/*}
└── results/
    ├── mtf2_nucleation/{nuc_cpgonly_spread_bins.tsv.gz, per_bin_covcorr_full.tsv.gz}
    └── replitag_sphase_nucleation/matrices/encode_coremarks.1000bp.rpkm.matrix.tsv.gz
```

`*` = written by a notebook in this repository; not required in the archive if you run
the notebooks in the order given above.

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

**Intermediate tables with no producer.**

- `data/results/mtf2_nucleation/nuc_cpgonly_spread_bins.tsv.gz`
- `data/results/mtf2_nucleation/per_bin_covcorr_full.tsv.gz`
- `data/2024/RepliTag/filtered_sams_WT_1.0/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP_<mark>.tsv`
  (the per-mark split of the combined `filter_sams.sh` table)
- `data/2024/RepliTag/featureCounts_drug_dm6/genes_v5_KABK_230607_featureCounts_counts_plusDUP.txt.summary`
- `data/2024/SH_all_data/downsample_peakPRINT_slope1_min300/reproducibility/FRiPs_merged_intervals_rep95.tsv`
- `data/2024/SH_all_data/peakPRINT_50/MERGED_downsampled_bins.tsv`
- `data/2024/SH_all_data/sections4/*.txt` and `All_sams4.txt` (BAM path manifests)
- `data/2024/ENCODE_chromHMM/{FRIPs,OverlapEnrichments}.tsv` and
  `data/2024/cCREs/intersected_annos/FRIPs.tsv` (ChromHMM `OverlapEnrichment` output)
- `data/2024/LAD_data/H3K27_TRIP.tsv`
- `data/2024/RepliTag/bws_bulk/` — assembled by hand from the `bam2bed2bw.sh` bigWigs plus
  the lifted Repli-seq track; no script does the assembly

**External downloads** the cloner must fetch: GENCODE v27 and v44 GTFs; MSigDB
`GOBP_CELL_CYCLE.v2025.1.Hs.json`; the published supplements `Dataset_S1` and
`Dataset_S2`; UCSC-ENCODE UW Repli-seq K562 hg19 bigWigs; GEO `GSE164804` (MTF2) and
ENCODE `ENCFF587SWK` (EZH2) / `ENCFF065KGU` (SUZ12); the IHEC K562 18-state segmentation.

---

## Notes

- `Feature_counts_SAMs-allFeatures_k562.sh` writes to `data/2024/LAD_data/FC_out/`, but
  `250911_RIFs_FC_pub.ipynb` reads from `data/2024/cCREs/intersected_annos/FC_out/`. The
  filenames match exactly; the outputs were evidently relocated after the run.
- The `processing scripts/` directory name contains a space, so it must be quoted in
  shell commands.
