# Primary raw-count differential-expression analysis.
#
# This script runs only when data/raw_counts.csv is present. The normalized
# abundance matrix in data/expression.csv is not a valid DESeq2 input.

if (!file.exists(files$raw_counts)) {
  note <- c(
    "DESeq2 analysis was skipped because data/raw_counts.csv is absent.",
    "Add the unnormalized integer count matrix with columns:",
    "Gene,A1,A2,A3,B1,B2,B3,C1,C2,C3",
    "and rerun source('run_all.R').",
    "Do not substitute data/expression.csv for raw counts."
  )
  writeLines(note, file.path(output_dir, "01_DESeq2_SKIPPED.txt"))
  message(note[1])
} else {
  require_package("DESeq2")

  count_data <- read_csv_checked(files$raw_counts)
  if (ncol(count_data) != 10L) {
    stop("raw_counts.csv must have one gene column and nine sample columns.")
  }

  gene_id <- as.character(count_data[[1]])
  counts <- as.matrix(count_data[, -1, drop = FALSE])
  suppressWarnings(storage.mode(counts) <- "numeric")
  rownames(counts) <- gene_id

  if (anyNA(counts) || any(!is.finite(counts))) {
    stop("raw_counts.csv contains missing or non-finite values.")
  }
  if (any(counts < 0) || any(counts != round(counts))) {
    stop("DESeq2 input must contain non-negative integer counts.")
  }
  if (anyDuplicated(gene_id)) stop("Duplicated gene IDs in raw_counts.csv.")

  expected_samples <- names(sample_groups)
  if (!setequal(colnames(counts), expected_samples)) {
    stop(
      "Count columns must be: ",
      paste(expected_samples, collapse = ", ")
    )
  }
  counts <- counts[, expected_samples, drop = FALSE]

  sample_data <- data.frame(
    row.names = expected_samples,
    group = factor(
      unname(sample_groups[expected_samples]),
      levels = c("Control", "F1", "F20")
    )
  )

  keep <- rowSums(counts >= 10) >= 3
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts[keep, , drop = FALSE]),
    colData = sample_data,
    design = ~ group
  )
  dds <- DESeq2::DESeq(dds)

  contrast_definitions <- list(
    F1_vs_Control = c("group", "F1", "Control"),
    F20_vs_Control = c("group", "F20", "Control"),
    F20_vs_F1 = c("group", "F20", "F1")
  )

  deseq_results <- list()
  for (contrast_name in names(contrast_definitions)) {
    result <- DESeq2::results(
      dds,
      contrast = contrast_definitions[[contrast_name]],
      alpha = thresholds$primary_fdr
    )
    table <- as.data.frame(result)
    table$Gene <- rownames(table)
    table$DEG <- (
      !is.na(table$padj) &
        table$padj <= thresholds$primary_fdr &
        abs(table$log2FoldChange) >= thresholds$primary_abs_log2fc
    )
    table$Direction <- ifelse(
      table$DEG & table$log2FoldChange > 0, "Up",
      ifelse(table$DEG & table$log2FoldChange < 0, "Down", "Not_DEG")
    )
    table <- table[
      ,
      c(
        "Gene", "baseMean", "log2FoldChange", "lfcSE", "stat",
        "pvalue", "padj", "DEG", "Direction"
      )
    ]
    rownames(table) <- NULL

    deseq_results[[contrast_name]] <- table
    write_result(
      table,
      file.path("DESeq2", paste0(contrast_name, "_all_genes.csv"))
    )
    write_result(
      table[table$DEG, , drop = FALSE],
      file.path("DESeq2", paste0(contrast_name, "_DEGs.csv"))
    )
  }

  saveRDS(
    deseq_results,
    file.path(intermediate_dir, "deseq2_results.rds")
  )
  saveRDS(dds, file.path(intermediate_dir, "deseq2_dataset.rds"))
}

