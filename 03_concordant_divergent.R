# Concordant and divergent transcriptomic responses.
#
# Linear contrasts:
#   C = [(F1 - Control) + (F20 - Control)] / 2
#   D = (F20 - Control) - (F1 - Control) = F20 - F1
#
# The reported concordant gene set was defined from directionally consistent
# pairwise responses, whereas the divergent set corresponds to the F20-F1
# comparison. When raw-count DESeq2 results are available, this script
# reconstructs those DEG sets. Otherwise, it validates the supplied reported
# gene sets and all contrast coefficients.

require_package("limma")

expression_analysis <- readRDS(
  file.path(intermediate_dir, "expression_pattern_analysis.rds")
)
log_expression_unscaled <- log2(expression_analysis$expression_matrix + 1)
design_cd <- expression_analysis$design

contrast_cd <- limma::makeContrasts(
  Concordant_Contrast = ((F1 - Control) + (F20 - Control)) / 2,
  Divergent_Contrast = (F20 - Control) - (F1 - Control),
  levels = design_cd
)

fit_cd <- limma::lmFit(log_expression_unscaled, design_cd)
fit_cd <- limma::contrasts.fit(fit_cd, contrast_cd)
fit_cd <- limma::eBayes(fit_cd)

coefficient_table <- data.frame(
  Gene = rownames(log_expression_unscaled),
  Concordant_Contrast =
    fit_cd$coefficients[, "Concordant_Contrast"],
  Divergent_Contrast =
    fit_cd$coefficients[, "Divergent_Contrast"],
  F1_vs_Control =
    rowMeans(log_expression_unscaled[, expression_analysis$group == "F1"]) -
    rowMeans(
      log_expression_unscaled[
        , expression_analysis$group == "Control", drop = FALSE
      ]
    ),
  F20_vs_Control =
    rowMeans(log_expression_unscaled[, expression_analysis$group == "F20"]) -
    rowMeans(
      log_expression_unscaled[
        , expression_analysis$group == "Control", drop = FALSE
      ]
    ),
  stringsAsFactors = FALSE
)
write_result(
  coefficient_table,
  file.path(
    "concordant_divergent",
    "linear_contrast_coefficients_all_genes.csv"
  )
)

reported_cd <- read_csv_checked(
  files$reported_concordant_divergent,
  c(
    "Concordant_Contrast", "Divergent_Contrast",
    "Hot_vs_Anc", "Cold_vs_Anc", "Classification"
  )
)
names(reported_cd)[1] <- "Gene"

reported_cd$Formula_C_error <- reported_cd$Concordant_Contrast -
  (reported_cd$Hot_vs_Anc + reported_cd$Cold_vs_Anc) / 2
reported_cd$Formula_D_error <- reported_cd$Divergent_Contrast -
  (reported_cd$Hot_vs_Anc - reported_cd$Cold_vs_Anc)

validation <- merge(
  reported_cd,
  coefficient_table,
  by = "Gene",
  suffixes = c("_reported", "_recalculated"),
  all.x = TRUE,
  sort = FALSE
)
validation$Recalculated_C_error <- (
  validation$Concordant_Contrast_reported -
    validation$Concordant_Contrast_recalculated
)
validation$Recalculated_D_error <- (
  validation$Divergent_Contrast_reported -
    validation$Divergent_Contrast_recalculated
)
write_result(
  validation,
  file.path(
    "concordant_divergent",
    "reported_contrast_validation.csv"
  )
)

# The reported CSV stores rounded decimal coefficients; 1e-7 accommodates
# only that text-export rounding and is far below the displayed precision.
tolerance <- 1e-7
error_columns <- c(
  "Formula_C_error", "Formula_D_error",
  "Recalculated_C_error", "Recalculated_D_error"
)
maximum_errors <- vapply(
  validation[, error_columns, drop = FALSE],
  function(x) max(abs(x), na.rm = TRUE),
  numeric(1)
)
if (any(maximum_errors > tolerance)) {
  stop(
    "Concordant/divergent coefficient validation failed: ",
    paste(names(maximum_errors), signif(maximum_errors, 4), collapse = "; ")
  )
}

reported_counts <- table(
  factor(
    reported_cd$Classification,
    levels = names(expected_concordant_divergent_counts)
  )
)
if (!identical(
  as.integer(reported_counts),
  as.integer(expected_concordant_divergent_counts)
)) {
  stop("Reported concordant/divergent counts are inconsistent.")
}

write_result(
  data.frame(
    Response = names(expected_concordant_divergent_counts),
    Count = as.integer(reported_counts)
  ),
  file.path(
    "concordant_divergent",
    "reported_response_counts.csv"
  )
)

deseq2_rds <- file.path(intermediate_dir, "deseq2_results.rds")
if (file.exists(deseq2_rds)) {
  deseq2_results <- readRDS(deseq2_rds)
  f1 <- deseq2_results[["F1_vs_Control"]]
  f20 <- deseq2_results[["F20_vs_Control"]]
  f20_f1 <- deseq2_results[["F20_vs_F1"]]

  pairwise <- merge(
    f1[, c("Gene", "log2FoldChange", "padj", "DEG")],
    f20[, c("Gene", "log2FoldChange", "padj", "DEG")],
    by = "Gene",
    suffixes = c("_F1_Control", "_F20_Control")
  )
  pairwise$Concordant <- (
    pairwise$DEG_F1_Control &
      pairwise$DEG_F20_Control &
      sign(pairwise$log2FoldChange_F1_Control) ==
        sign(pairwise$log2FoldChange_F20_Control)
  )

  divergent <- f20_f1[, c(
    "Gene", "log2FoldChange", "padj", "DEG"
  )]
  names(divergent) <- c(
    "Gene", "log2FoldChange_F20_F1", "padj_F20_F1", "Divergent"
  )

  derived_sets <- merge(pairwise, divergent, by = "Gene", all = TRUE)
  write_result(
    derived_sets,
    file.path(
      "concordant_divergent",
      "DESeq2_derived_concordant_divergent_sets.csv"
    )
  )
} else {
  message(
    "Raw-count DESeq2 output unavailable: using the reported gene sets ",
    "for downstream concordant/divergent summaries."
  )
}

saveRDS(
  list(
    coefficients = coefficient_table,
    reported = reported_cd,
    validation = validation
  ),
  file.path(intermediate_dir, "concordant_divergent_analysis.rds")
)

message("Concordant/divergent formulas and coefficients validated.")
