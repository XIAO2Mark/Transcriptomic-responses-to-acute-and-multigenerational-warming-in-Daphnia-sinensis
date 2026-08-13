# Sensitivity analysis excluding the upper 10% of genes with the greatest
# pooled within-group residual variance.

require_package("limma")

analysis <- readRDS(
  file.path(intermediate_dir, "expression_pattern_analysis.rds")
)
log_expression <- analysis$log_expression
original_fit <- analysis$fit
original_assignments <- analysis$assignments

residual_variance <- original_fit$base_fit$sigma^2
variance_cutoff <- as.numeric(
  stats::quantile(
    residual_variance,
    probs = 1 - thresholds$high_variance_fraction,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
)
high_variance <- residual_variance > variance_cutoff
variance_eligible <- !high_variance

variance_table <- data.frame(
  Gene = rownames(log_expression),
  Within_group_residual_variance = residual_variance,
  High_variance_top_10_percent = high_variance,
  Original_pattern_internal = original_assignments$Pattern_internal,
  Original_pattern = original_assignments$Pattern,
  stringsAsFactors = FALSE
)
variance_table <- variance_table[
  order(variance_table$Within_group_residual_variance, decreasing = TRUE),
]
write_result(
  variance_table,
  file.path(
    "variance_sensitivity",
    "within_group_residual_variance_all_genes.csv"
  )
)
write_result(
  variance_table[
    variance_table$High_variance_top_10_percent, , drop = FALSE
  ],
  file.path(
    "variance_sensitivity",
    "high_variance_genes_removed.csv"
  )
)

filtered_log_expression <- log_expression[
  variance_eligible, , drop = FALSE
]
filtered_fit <- run_pattern_limma(filtered_log_expression)
filtered_assignments <- classify_patterns(filtered_fit$results)

write_result(
  filtered_assignments,
  file.path(
    "variance_sensitivity",
    "variance_filtered_pattern_assignments.csv"
  )
)

comparison <- merge(
  original_assignments[
    , c("Gene", "Pattern_internal", "Pattern"), drop = FALSE
  ],
  filtered_assignments[
    , c("Gene", "Pattern_internal", "Pattern"), drop = FALSE
  ],
  by = "Gene",
  all.x = TRUE,
  suffixes = c("_original", "_filtered"),
  sort = FALSE
)

summary_rows <- lapply(
  names(expected_pattern_counts),
  function(pattern_name) {
    original_scope <- comparison$Pattern_internal_original == pattern_name
    eligible_scope <- original_scope &
      !is.na(comparison$Pattern_internal_filtered)
    retained_scope <- eligible_scope &
      comparison$Pattern_internal_filtered == pattern_name

    data.frame(
      Pattern_internal = pattern_name,
      Pattern = display_pattern(pattern_name),
      Original_n = sum(original_scope),
      High_variance_genes_removed =
        sum(original_scope) - sum(eligible_scope),
      Variance_eligible_n = sum(eligible_scope),
      Retained_same_pattern_n = sum(retained_scope),
      Retention_percent = 100 * sum(retained_scope) / sum(eligible_scope),
      Reclassified_n = sum(
        filtered_assignments$Pattern_internal == pattern_name
      ),
      stringsAsFactors = FALSE
    )
  }
)
sensitivity_summary <- do.call(rbind, summary_rows)
write_result(
  sensitivity_summary,
  file.path(
    "variance_sensitivity",
    "Table_S17_pattern_robustness.csv"
  )
)

expected_s17 <- data.frame(
  Pattern_internal = names(expected_pattern_counts),
  Original_n = c(2438L, 839L, 439L, 260L, 8877L),
  High_variance_genes_removed = c(206L, 113L, 21L, 28L, 918L),
  Variance_eligible_n = c(2232L, 726L, 418L, 232L, 7959L),
  Retained_same_pattern_n = c(2114L, 667L, 413L, 228L, 7506L),
  Reclassified_n = c(2462L, 769L, 519L, 270L, 7547L)
)

check <- merge(
  sensitivity_summary,
  expected_s17,
  by = "Pattern_internal",
  suffixes = c("_observed", "_expected"),
  sort = FALSE
)
integer_metrics <- c(
  "Original_n",
  "High_variance_genes_removed",
  "Variance_eligible_n",
  "Retained_same_pattern_n",
  "Reclassified_n"
)
for (metric in integer_metrics) {
  observed <- check[[paste0(metric, "_observed")]]
  expected <- check[[paste0(metric, "_expected")]]
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop("Table S17 reproduction failed for metric: ", metric)
  }
}

transition_matrix <- table(
  factor(
    comparison$Pattern_internal_original,
    levels = names(expected_pattern_counts)
  ),
  factor(
    comparison$Pattern_internal_filtered,
    levels = names(expected_pattern_counts)
  )
)
write_result(
  as.data.frame.matrix(transition_matrix),
  file.path(
    "variance_sensitivity",
    "pattern_transition_matrix.csv"
  ),
  row.names = TRUE
)

saveRDS(
  list(
    variance_table = variance_table,
    filtered_fit = filtered_fit,
    filtered_assignments = filtered_assignments,
    comparison = comparison,
    summary = sensitivity_summary
  ),
  file.path(intermediate_dir, "variance_sensitivity_analysis.rds")
)

message("Table S17 sensitivity results reproduced exactly.")

