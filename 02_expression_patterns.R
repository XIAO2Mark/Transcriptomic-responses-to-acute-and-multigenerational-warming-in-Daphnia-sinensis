# Expression-pattern classification from normalized abundance values.
#
# The historical trajectory analysis used:
#   - expression > 0 in at least three samples;
#   - global scaling to an identical sample mean;
#   - log2(expression + 1);
#   - limma with empirical-Bayes moderation;
#   - moderated P.Value < 0.05 for the reported pattern assignments.
#
# Adjusted P values are exported but are not substituted for the historical
# threshold because doing so would not reproduce the reported gene sets.

require_package("limma")

expression_data <- read_csv_checked(files$expression)
if (ncol(expression_data) != 10L) {
  stop("expression.csv must have one gene column and nine sample columns.")
}

gene_id <- as.character(expression_data[[1]])
expression_matrix_all <- as.matrix(
  expression_data[, -1, drop = FALSE]
)
suppressWarnings(storage.mode(expression_matrix_all) <- "numeric")
rownames(expression_matrix_all) <- gene_id

if (anyDuplicated(gene_id)) stop("Duplicated gene IDs in expression.csv.")
if (anyNA(expression_matrix_all) || any(!is.finite(expression_matrix_all))) {
  stop("expression.csv contains missing or non-finite values.")
}
if (any(expression_matrix_all < 0)) {
  stop("expression.csv contains negative abundance values.")
}

expected_samples <- names(sample_groups)
if (!setequal(colnames(expression_matrix_all), expected_samples)) {
  stop(
    "Expression columns must be: ",
    paste(expected_samples, collapse = ", ")
  )
}
expression_matrix_all <- expression_matrix_all[
  , expected_samples, drop = FALSE
]

expressed <- rowSums(expression_matrix_all > 0) >= 3
expression_matrix <- expression_matrix_all[expressed, , drop = FALSE]

target_mean <- mean(colMeans(expression_matrix))
normalized_matrix <- sweep(
  expression_matrix,
  2,
  colMeans(expression_matrix),
  FUN = "/"
) * target_mean
log_expression <- log2(normalized_matrix + 1)

group <- factor(
  unname(sample_groups[colnames(log_expression)]),
  levels = c("Control", "F1", "F20")
)
design <- stats::model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(log_expression)

contrast_matrix <- limma::makeContrasts(
  F1_vs_Control = F1 - Control,
  F20_vs_Control = F20 - Control,
  F20_vs_F1 = F20 - F1,
  levels = design
)

run_pattern_limma <- function(matrix) {
  base_fit <- limma::lmFit(matrix, design)
  contrast_fit <- limma::contrasts.fit(base_fit, contrast_matrix)
  contrast_fit <- limma::eBayes(contrast_fit)

  result_list <- vector("list", ncol(contrast_matrix))
  names(result_list) <- colnames(contrast_matrix)

  for (index in seq_len(ncol(contrast_matrix))) {
    table <- limma::topTable(
      contrast_fit,
      coef = index,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )
    table$Gene <- rownames(table)
    table <- table[
      ,
      c(
        "Gene", "logFC", "AveExpr", "t", "P.Value",
        "adj.P.Val", "B"
      )
    ]
    rownames(table) <- NULL
    result_list[[index]] <- table
  }

  list(
    base_fit = base_fit,
    contrast_fit = contrast_fit,
    results = result_list
  )
}

classify_patterns <- function(result_list) {
  f1 <- result_list[["F1_vs_Control"]]
  f20 <- result_list[["F20_vs_Control"]]
  f20_f1 <- result_list[["F20_vs_F1"]]

  stopifnot(
    identical(f1$Gene, f20$Gene),
    identical(f1$Gene, f20_f1$Gene)
  )

  significant_f1 <- f1$P.Value < thresholds$pattern_p
  significant_f20 <- f20$P.Value < thresholds$pattern_p
  significant_f20_f1 <- f20_f1$P.Value < thresholds$pattern_p

  same_direction <- sign(f1$logFC) == sign(f20$logFC)
  f20_stronger <- abs(f20$logFC) > abs(f1$logFC)

  pattern <- rep("Other", nrow(f1))

  # The order below is explicit and the rules are mutually interpretable.
  pattern[!significant_f1 & significant_f20] <- "Novel Adaptation"
  pattern[significant_f1 & !significant_f20] <- "Reverse"
  pattern[
    significant_f1 & significant_f20 &
      same_direction & f20_stronger
  ] <- "Reinforcement"
  pattern[
    significant_f1 & significant_f20 &
      same_direction & !f20_stronger & !significant_f20_f1
  ] <- "Assimilation"

  data.frame(
    Gene = f1$Gene,
    LFC_F1_Ctrl = f1$logFC,
    P_F1_Ctrl = f1$P.Value,
    FDR_F1_Ctrl = f1$adj.P.Val,
    LFC_F20_Ctrl = f20$logFC,
    P_F20_Ctrl = f20$P.Value,
    FDR_F20_Ctrl = f20$adj.P.Val,
    LFC_F20_F1 = f20_f1$logFC,
    P_F20_F1 = f20_f1$P.Value,
    FDR_F20_F1 = f20_f1$adj.P.Val,
    Pattern_internal = pattern,
    Pattern = display_pattern(pattern),
    stringsAsFactors = FALSE
  )
}

expression_fit <- run_pattern_limma(log_expression)
pattern_assignments <- classify_patterns(expression_fit$results)

for (contrast_name in names(expression_fit$results)) {
  write_result(
    expression_fit$results[[contrast_name]],
    file.path(
      "expression_patterns",
      paste0("pairwise_", contrast_name, ".csv")
    )
  )
}

write_result(
  pattern_assignments,
  file.path(
    "expression_patterns",
    "expression_pattern_assignments.csv"
  )
)

pattern_counts_internal <- table(
  factor(
    pattern_assignments$Pattern_internal,
    levels = names(expected_pattern_counts)
  )
)
if (!identical(
  as.integer(pattern_counts_internal),
  as.integer(expected_pattern_counts)
)) {
  stop(
    "Reported expression-pattern counts were not reproduced. Obtained: ",
    paste(
      names(expected_pattern_counts),
      as.integer(pattern_counts_internal),
      collapse = "; "
    )
  )
}

pattern_count_table <- data.frame(
  Pattern_internal = names(expected_pattern_counts),
  Pattern = display_pattern(names(expected_pattern_counts)),
  Count = as.integer(pattern_counts_internal),
  stringsAsFactors = FALSE
)
write_result(
  pattern_count_table,
  file.path("expression_patterns", "expression_pattern_counts.csv")
)

# Independent gene-by-gene comparison with the saved reported assignments.
reported_patterns <- read_csv_checked(
  files$reported_patterns,
  c("Gene", "Pattern")
)
reported_patterns$Pattern_reported <- display_pattern(
  reported_patterns$Pattern
)
reproduction_check <- merge(
  pattern_assignments[, c("Gene", "Pattern")],
  reported_patterns[, c("Gene", "Pattern_reported")],
  by = "Gene",
  all = TRUE,
  sort = FALSE
)
reproduction_check$Match <- (
  !is.na(reproduction_check$Pattern) &
    !is.na(reproduction_check$Pattern_reported) &
    reproduction_check$Pattern == reproduction_check$Pattern_reported
)
write_result(
  reproduction_check,
  file.path("expression_patterns", "gene_level_reproduction_check.csv")
)

if (!all(reproduction_check$Match)) {
  stop(
    "Gene-level expression-pattern reproduction failed for ",
    sum(!reproduction_check$Match), " genes."
  )
}

saveRDS(
  list(
    expression_matrix = expression_matrix,
    normalized_matrix = normalized_matrix,
    log_expression = log_expression,
    group = group,
    design = design,
    contrast_matrix = contrast_matrix,
    fit = expression_fit,
    assignments = pattern_assignments
  ),
  file.path(intermediate_dir, "expression_pattern_analysis.rds")
)

message("Expression-pattern assignments reproduced exactly.")

