#!/usr/bin/env Rscript

# 1) libraries
library(limma)
library(edgeR)
library(RUVSeq)
library(jsonlite)

# 2) marks
marks <- c("P2S25p", "K27me3", "K27me2", "K27me1", "K27ac")

# 3) paths
base_counts_path <- "data/2024/RepliTag/filtered_sams_WT_1.0/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP_%s.tsv"
base_output_dir  <- "data/2024/RepliTag/filtered_sams_WT_1.0/251113_DGE"

# 3b) complete library sizes from featureCounts summary
summary_file <- "data/2024/RepliTag/filtered_sams_WT_1.0/hg38.gencodev44_RepliTag_genes-1000_featureCounts_counts_plusDUP.txt.summary.tsv"

summary_tab <- read.delim(
  summary_file,
  sep         = "\t",
  check.names = FALSE,
  row.names   = 1
)

# total reads per column (assigned + all unassigned)
libsize_full <- colSums(summary_tab)

# make names match count_cols:
#   "/.../JG_HsDm_K562_P2S25p_241101_G1_d2pcr_a.DTmarked_filtered.sam"
#   -> "JG_HsDm_K562_P2S25p_241101_G1_d2pcr_a"
names(libsize_full) <- gsub(
  "\\.DTmarked_filtered\\.sam$",
  "",
  basename(colnames(summary_tab))
)

# 4) the 14 sample suffixes
suffixes <- c(
  "_241101_G1_d2pcr_a","_241101_G1_d2pcr_b",
  "_241101_G2_d2pcr_a","_241101_G2_d2pcr_b",
  "_241101_S1_d2pcr_a","_241101_S1_d2pcr_b",
  "_241101_S2_d2pcr_a","_241101_S2_d2pcr_b",
  "_241101_S3_d2pcr_a","_241101_S3_d2pcr_b",
  "_241101_S4_d2pcr_a","_241101_S4_d2pcr_b",
  "_241101_WT_d2pcr_a","_241101_WT_d2pcr_b"
)

# 5) phase metadata
phases <- c("G1","S1","S2","S3","S4","G2","SS")

sample_phases <- c(
  "G1","G1","G2","G2",
  "S1","S1","S2","S2",
  "S3","S3","S4","S4",
  "SS","SS"
)

# 6) loop over marks
for (mark in marks) {
  
  message("=== Processing ", mark, " ===")
  
  mark_label <- if (mark == "P2S25p") "RNAPII" else mark
  
  counts_file <- sprintf(base_counts_path, mark)
  output_dir  <- file.path(base_output_dir, mark)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  count_cols <- paste0("JG_HsDm_K562_", mark, suffixes)
  
  # read counts
  x <- read.delim(
    counts_file,
    sep          = "\t",
    check.names  = FALSE,
    colClasses   = "character",
    fileEncoding = "UTF-8-BOM"
  )
  colnames(x) <- trimws(gsub('^"|"$', "", colnames(x)))
  
  if (!"Geneid" %in% colnames(x)) {
    stop("Geneid column not found in counts file: ", counts_file)
  }
  rownames(x) <- x$Geneid
  
  x[, count_cols] <- lapply(x[, count_cols], as.numeric)
  counts <- as.matrix(x[, count_cols])
  
  rn_all <- rownames(counts)
  
  # get complete library sizes for these samples
  libsize_mark <- libsize_full[count_cols]
  if (any(is.na(libsize_mark))) {
    missing <- count_cols[is.na(libsize_mark)]
    stop("Missing library sizes for samples: ", paste(missing, collapse = ", "))
  }
  
  # design: one column per phase (degust style)
  phase_factor <- factor(sample_phases, levels = phases)
  design_raw <- model.matrix(~ 0 + phase_factor)
  colnames(design_raw) <- paste0(mark_label, "_", phases)
  design <- design_raw
  rownames(design) <- count_cols
  
  # DGEList with complete library sizes + RLE factors
  dge <- DGEList(counts = counts, lib.size = libsize_mark)
  dge <- calcNormFactors(dge, method = "RLE")
  
  # initial fit for RUVr
  y0   <- estimateDisp(dge, design)
  fit0 <- glmQLFit(y0, design)
  
  resids <- residuals(fit0)
  ctrls  <- rn_all
  ruv_r  <- RUVr(as.matrix(counts), ctrls, k = 1, res = resids)
  
  design_ruv <- cbind(design, ruv_r$W)
  colnames(design_ruv) <- c(
    colnames(design),
    paste0("W", seq_len(ncol(ruv_r$W)))
  )
  
  # RPKM summary
  gene_length <- as.numeric(x$Length)
  rpkm_all <- rpkm(
    dge,
    gene.length           = gene_length,
    normalized.lib.sizes  = TRUE
  )
  
  cond_means <- sapply(phases, function(ph) {
    idx <- sample_phases == ph
    rowMeans(rpkm_all[, idx, drop = FALSE])
  })
  time_phases <- phases[phases != "SS"]
  cond_means_noSS <- cond_means[, time_phases, drop = FALSE]
  max_mean_rpkm <- apply(cond_means_noSS, 1, max)
  min_mean_rpkm <- apply(cond_means_noSS, 1, min)
  
  high_expr <- max_mean_rpkm >= 0.5
  
  res_pairs <- list()
  res_all   <- NULL
  
  if (any(high_expr)) {
    
    dge_test <- dge[high_expr, , keep.lib.sizes = TRUE]
    
    y1   <- estimateDisp(dge_test, design_ruv)
    fit1 <- glmQLFit(y1, design_ruv)
    
    # 1) adjacent pairwise contrasts
    pair_list <- list(
      G2_vs_G1 = c("G2", "G1"),
      G1_vs_S1 = c("G1", "S1"),
      S1_vs_S2 = c("S1", "S2"),
      S2_vs_S3 = c("S2", "S3"),
      S3_vs_S4 = c("S3", "S4"),
      S4_vs_G2 = c("S4", "G2")
    )
    
    for (nm in names(pair_list)) {
      later   <- pair_list[[nm]][1]
      earlier <- pair_list[[nm]][2]
      
      v <- rep(0, ncol(design_ruv))
      names(v) <- colnames(design_ruv)
      
      col_later   <- paste0(mark_label, "_", later)
      col_earlier <- paste0(mark_label, "_", earlier)
      
      v[col_earlier] <- -1
      v[col_later]   <-  1
      
      qlf_pair <- glmQLFTest(fit1, contrast = v)
      res_pairs[[nm]] <- topTags(qlf_pair, n = Inf, sort.by = "none")$table
    }
    
    # 2) global time course: G1 vs all others
    baseline     <- "G1"
    other_phases <- phases[phases != baseline]
    
    cont_tc <- sapply(other_phases, function(ph) {
      v <- rep(0, ncol(design_ruv))
      names(v) <- colnames(design_ruv)
      v[paste0(mark_label, "_", baseline)] <- -1
      v[paste0(mark_label, "_", ph)]       <-  1
      v
    })
    colnames(cont_tc) <- paste0(other_phases, "_vs_", baseline)
    
    qlf_all <- glmQLFTest(fit1, contrast = cont_tc)
    res_all <- topTags(qlf_all, n = Inf, sort.by = "none")$table
  }
  
  # assemble output over all genes
  df_out <- data.frame(row.names = rn_all)
  
  pair_names <- c("G2_vs_G1", "G1_vs_S1", "S1_vs_S2",
                  "S2_vs_S3", "S3_vs_S4", "S4_vs_G2")
  
  for (nm in pair_names) {
    df_out[[paste0("logFC_", nm)]]  <- NA_real_
    df_out[[paste0("PValue_", nm)]] <- NA_real_
    df_out[[paste0("FDR_", nm)]]    <- NA_real_
  }
  
  tc_contrasts <- paste0(phases[phases != "G1"], "_vs_G1")
  for (nm in tc_contrasts) {
    df_out[[paste0("logFC_", nm)]] <- NA_real_
  }
  
  df_out[["AveExpr"]]      <- NA_real_
  df_out[["F_allTP"]]      <- NA_real_
  df_out[["PValue_allTP"]] <- NA_real_
  df_out[["FDR_allTP"]]    <- NA_real_
  
  if (!is.null(res_all)) {
    
    # high_expr genes are the only ones in res_pairs and res_all
    high_ids <- rownames(res_all)  # same subset as dge_test
    
    # fill pairwise adjacent stats
    for (nm in pair_names) {
      tab <- res_pairs[[nm]]
      rn  <- rownames(tab)
      df_out[rn, paste0("logFC_", nm)]  <- tab$logFC
      df_out[rn, paste0("PValue_", nm)] <- tab$PValue
      df_out[rn, paste0("FDR_", nm)]    <- tab$FDR
    }
    
    # global time course
    rn_tc <- rownames(res_all)
    
    lfc_mat <- as.matrix(res_all[, seq_along(tc_contrasts), drop = FALSE])
    colnames(lfc_mat) <- tc_contrasts
    
    for (j in seq_along(tc_contrasts)) {
      nm <- tc_contrasts[j]
      df_out[rn_tc, paste0("logFC_", nm)] <- lfc_mat[, j]
    }
    
    df_out[rn_tc, "AveExpr"]      <- res_all$logCPM
    df_out[rn_tc, "F_allTP"]      <- res_all$F
    df_out[rn_tc, "PValue_allTP"] <- res_all$PValue
    df_out[rn_tc, "FDR_allTP"]    <- res_all$FDR
  }
  
  # main P.Value and adj.P.Val from G2_vs_G1 contrast
  df_out[["P.Value"]]   <- df_out[["PValue_G2_vs_G1"]]
  df_out[["adj.P.Val"]] <- df_out[["FDR_G2_vs_G1"]]
  # time course F test lives in FDR_allTP and PValue_allTP
  
  # RPKM summaries for all genes
  df_out[["max_RPKM"]] <- max_mean_rpkm[rownames(df_out)]
  df_out[["min_RPKM"]] <- min_mean_rpkm[rownames(df_out)]
  
  # backend normalized CPM from RUVr
  backend_cpm <- cpm(
    ruv_r$normalizedCounts,
    log         = FALSE,
    prior.count = 0
  )
  
  for (j in seq_len(ncol(backend_cpm))) {
    samp <- colnames(backend_cpm)[j]
    colname <- paste0(samp, "_CPM")
    df_out[[colname]] <- backend_cpm[rownames(df_out), j]
  }
  
  # mean backend CPM per phase
  backend_cpm_means <- sapply(phases, function(ph) {
    idx <- sample_phases == ph
    rowMeans(backend_cpm[, idx, drop = FALSE])
  })
  
  for (ph in phases) {
    colname <- paste0("meanCPM_", ph)
    df_out[[colname]] <- backend_cpm_means[rownames(df_out), ph]
  }
  
  # annotate plus raw counts and write
  export_cols <- c("Geneid", "Length", count_cols)
  out2 <- cbind(df_out, x[, export_cols])
  
  write.table(
    out2,
    sep       = "\t",
    file      = file.path(output_dir, "output.tsv"),
    row.names = FALSE,
    na        = ""
  )
}
