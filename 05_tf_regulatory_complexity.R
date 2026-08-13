# Reproduce the statistical source data for Fig. 5A-E.
#
# TF load is the number of distinct high-confidence motif clusters detected in
# the promoter of each gene after promoter QC. It is a predicted regulatory
# metric, not a direct measurement of TF binding.

tf_metrics <- read_csv_checked(
  files$tf_gene_metrics,
  c(
    "Gene", "Pattern_display", "Direction", "F0_mean",
    "gc_fraction_3kb", "predicted_TF_motif_load",
    "genomic_block_100kb"
  )
)
tf_metrics$Pattern_display <- display_pattern(
  tf_metrics$Pattern_display
)
tf_metrics$Direction <- factor(
  tf_metrics$Direction,
  levels = c("Down", "Up")
)
tf_metrics$Pattern_display <- factor(
  tf_metrics$Pattern_display,
  levels = pattern_levels_display
)
tf_metrics$predicted_TF_motif_load <- as.numeric(
  tf_metrics$predicted_TF_motif_load
)
tf_metrics$Has_predicted_TF <- (
  tf_metrics$predicted_TF_motif_load > 0
)

if (anyNA(tf_metrics$Pattern_display) || anyNA(tf_metrics$Direction)) {
  stop("Unexpected pattern or direction labels in TF gene-level data.")
}
if (anyDuplicated(tf_metrics$Gene)) {
  stop("TF gene-level data must contain exactly one row per gene.")
}

# --------------------------------------------------------------------------
# Fig. 5A: predicted TF coverage by expression direction within each pattern.
# --------------------------------------------------------------------------
panel_a_rows <- list()
panel_a_tests <- data.frame(
  Pattern = character(),
  fisher_p_up_vs_down = numeric(),
  stringsAsFactors = FALSE
)

for (pattern_name in pattern_levels_display) {
  scope <- tf_metrics[
    tf_metrics$Pattern_display == pattern_name, , drop = FALSE
  ]
  contingency <- table(
    factor(scope$Direction, levels = c("Down", "Up")),
    factor(scope$Has_predicted_TF, levels = c(FALSE, TRUE))
  )
  p_value <- stats::fisher.test(contingency)$p.value

  panel_a_tests <- rbind(
    panel_a_tests,
    data.frame(
      Pattern = pattern_name,
      fisher_p_up_vs_down = p_value
    )
  )

  for (direction_name in c("Down", "Up")) {
    direction_scope <- scope[
      scope$Direction == direction_name, , drop = FALSE
    ]
    panel_a_rows[[paste(pattern_name, direction_name)]] <- data.frame(
      Pattern = pattern_name,
      Direction = direction_name,
      n_genes = nrow(direction_scope),
      n_with_predicted_TF_input =
        sum(direction_scope$Has_predicted_TF),
      percent_with_predicted_TF_input =
        100 * mean(direction_scope$Has_predicted_TF),
      mean_TF_load =
        mean(direction_scope$predicted_TF_motif_load),
      median_TF_load =
        stats::median(direction_scope$predicted_TF_motif_load),
      stringsAsFactors = FALSE
    )
  }
}
panel_a <- do.call(rbind, panel_a_rows)
rownames(panel_a) <- NULL
panel_a_tests$fisher_q_up_vs_down <- stats::p.adjust(
  panel_a_tests$fisher_p_up_vs_down,
  method = "BH"
)
panel_a <- merge(
  panel_a,
  panel_a_tests,
  by = "Pattern",
  all.x = TRUE,
  sort = FALSE
)
panel_a$Pattern <- factor(
  panel_a$Pattern,
  levels = pattern_levels_display
)
panel_a <- panel_a[
  order(panel_a$Pattern, panel_a$Direction), , drop = FALSE
]
panel_a$Pattern <- as.character(panel_a$Pattern)

# --------------------------------------------------------------------------
# Fig. 5B/C: overall TF coverage and mean TF load by pattern.
# --------------------------------------------------------------------------
other_scope <- tf_metrics[
  tf_metrics$Pattern_display == "Other genes", , drop = FALSE
]
panel_bc_rows <- list()

for (pattern_name in pattern_levels_display) {
  scope <- tf_metrics[
    tf_metrics$Pattern_display == pattern_name, , drop = FALSE
  ]

  if (pattern_name == "Other genes") {
    fisher_p <- NA_real_
    fisher_or <- NA_real_
    wilcoxon_p <- NA_real_
  } else {
    a <- sum(scope$Has_predicted_TF)
    b <- sum(!scope$Has_predicted_TF)
    c_value <- sum(other_scope$Has_predicted_TF)
    d <- sum(!other_scope$Has_predicted_TF)
    fisher_p <- stats::fisher.test(
      matrix(c(a, b, c_value, d), nrow = 2, byrow = TRUE)
    )$p.value
    fisher_or <- raw_odds_ratio(a, b, c_value, d)
    wilcoxon_p <- stats::wilcox.test(
      scope$predicted_TF_motif_load,
      other_scope$predicted_TF_motif_load,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )$p.value
  }

  panel_bc_rows[[pattern_name]] <- data.frame(
    Pattern = pattern_name,
    n_genes = nrow(scope),
    n_with_predicted_TF_input = sum(scope$Has_predicted_TF),
    percent_with_predicted_TF_input =
      100 * mean(scope$Has_predicted_TF),
    fisher_odds_ratio_vs_Other = fisher_or,
    fisher_p_value_vs_Other = fisher_p,
    mean_TF_load = mean(scope$predicted_TF_motif_load),
    median_TF_load = stats::median(scope$predicted_TF_motif_load),
    maximum_TF_load = max(scope$predicted_TF_motif_load),
    wilcoxon_p_value_load_vs_Other = wilcoxon_p,
    stringsAsFactors = FALSE
  )
}
panel_bc <- do.call(rbind, panel_bc_rows)
rownames(panel_bc) <- NULL

non_reference <- panel_bc$Pattern != "Other genes"
panel_bc$fisher_q_value_vs_Other <- NA_real_
panel_bc$fisher_q_value_vs_Other[non_reference] <- stats::p.adjust(
  panel_bc$fisher_p_value_vs_Other[non_reference],
  method = "BH"
)
panel_bc$wilcoxon_q_value_load_vs_Other <- NA_real_
panel_bc$wilcoxon_q_value_load_vs_Other[non_reference] <-
  stats::p.adjust(
    panel_bc$wilcoxon_p_value_load_vs_Other[non_reference],
    method = "BH"
  )

panel_b <- panel_bc[
  ,
  c(
    "Pattern", "n_genes", "n_with_predicted_TF_input",
    "percent_with_predicted_TF_input",
    "fisher_odds_ratio_vs_Other",
    "fisher_p_value_vs_Other",
    "fisher_q_value_vs_Other"
  )
]
panel_c <- panel_bc[
  ,
  c(
    "Pattern", "n_genes", "mean_TF_load", "median_TF_load",
    "maximum_TF_load", "wilcoxon_p_value_load_vs_Other",
    "wilcoxon_q_value_load_vs_Other"
  )
]

# --------------------------------------------------------------------------
# Fig. 5D: one-versus-rest logistic regression.
# --------------------------------------------------------------------------
model_data <- tf_metrics
model_data$baseline_z <- population_z(model_data$F0_mean)
model_data$tf_load_z <- population_z(
  model_data$predicted_TF_motif_load
)
model_data$promoter_gc_z <- population_z(
  model_data$gc_fraction_3kb
)

model_rows <- lapply(
  pattern_levels_display,
  function(pattern_name) {
    fit_clustered_membership_model(
      data = model_data,
      outcome = model_data$Pattern_display == pattern_name,
      cluster_column = "genomic_block_100kb",
      pattern_name = pattern_name
    )
  }
)
panel_d <- do.call(rbind, model_rows)
panel_d$q_value_all_10_effects <- stats::p.adjust(
  panel_d$p_value,
  method = "BH"
)

# --------------------------------------------------------------------------
# Fig. 5E: TF-encoding genes and TF-family composition.
# --------------------------------------------------------------------------
tf_annotation <- read_csv_checked(
  files$tf_annotation,
  c("Gene ID", "Family")
)
tf_annotation <- unique(
  tf_annotation[, c("Gene ID", "Family"), drop = FALSE]
)
names(tf_annotation)[1] <- "Gene"

# Panel E describes genes that themselves encode TFs. It therefore uses the
# full expression-pattern gene universe and is not restricted by promoter-QC
# eligibility, which applies only to the predicted promoter-load analyses in
# panels A-D.
expression_patterns <- readRDS(
  file.path(intermediate_dir, "expression_pattern_analysis.rds")
)$assignments
pattern_direction <- expression_patterns[
  ,
  c("Gene", "Pattern", "LFC_F20_Ctrl"),
  drop = FALSE
]
names(pattern_direction)[2] <- "Pattern_display"
pattern_direction$Pattern_display <- display_pattern(
  pattern_direction$Pattern_display
)
pattern_direction$Direction <- ifelse(
  pattern_direction$LFC_F20_Ctrl >= 0,
  "Up",
  "Down"
)
pattern_direction$LFC_F20_Ctrl <- NULL
tf_family <- merge(
  tf_annotation,
  pattern_direction,
  by = "Gene",
  all = FALSE
)
tf_family <- tf_family[
  tf_family$Pattern_display %in%
    c("Novel Adaptation", "Reverse"), ,
  drop = FALSE
]

panel_e <- stats::aggregate(
  tf_family$Gene,
  by = list(
    Pattern = tf_family$Pattern_display,
    Family = tf_family$Family,
    Direction = as.character(tf_family$Direction)
  ),
  FUN = function(x) length(unique(x))
)
names(panel_e)[4] <- "n_TF_encoding_genes"
panel_e <- panel_e[
  order(panel_e$Pattern, panel_e$Family, panel_e$Direction),
  ,
  drop = FALSE
]

# Write Fig. 5 source tables.
write_result(panel_a, file.path("TF_analysis", "Figure5A_source_data.csv"))
write_result(panel_b, file.path("TF_analysis", "Figure5B_source_data.csv"))
write_result(panel_c, file.path("TF_analysis", "Figure5C_source_data.csv"))
write_result(panel_d, file.path("TF_analysis", "Figure5D_source_data.csv"))
write_result(panel_e, file.path("TF_analysis", "Figure5E_source_data.csv"))

# Chen-inspired supplementary summaries across every expression category.
chen_rows <- lapply(
  pattern_levels_display,
  function(pattern_name) {
    scope <- tf_metrics[
      tf_metrics$Pattern_display == pattern_name, , drop = FALSE
    ]
    data.frame(
      Pattern = pattern_name,
      n_genes = nrow(scope),
      mean_TF_load = mean(scope$predicted_TF_motif_load),
      median_TF_load = stats::median(scope$predicted_TF_motif_load),
      fraction_with_no_predicted_TF =
        mean(scope$predicted_TF_motif_load == 0),
      mean_ancestral_expression = mean(scope$F0_mean),
      median_ancestral_expression = stats::median(scope$F0_mean),
      stringsAsFactors = FALSE
    )
  }
)
chen_summary <- do.call(rbind, chen_rows)
write_result(
  chen_summary,
  file.path(
    "TF_analysis",
    "supplementary_all_pattern_TF_expression_summary.csv"
  )
)

spearman <- stats::cor.test(
  tf_metrics$predicted_TF_motif_load,
  tf_metrics$F0_mean,
  method = "spearman",
  exact = FALSE
)
zero_vs_positive <- stats::wilcox.test(
  tf_metrics$F0_mean[tf_metrics$predicted_TF_motif_load == 0],
  tf_metrics$F0_mean[tf_metrics$predicted_TF_motif_load > 0],
  exact = FALSE,
  correct = FALSE
)
association_summary <- data.frame(
  Analysis = c(
    "Spearman: TF load vs ancestral expression",
    "Wilcoxon: ancestral expression, zero vs positive TF load"
  ),
  Estimate = c(unname(spearman$estimate), NA_real_),
  P_value = c(spearman$p.value, zero_vs_positive$p.value),
  stringsAsFactors = FALSE
)
write_result(
  association_summary,
  file.path(
    "TF_analysis",
    "supplementary_TF_load_expression_associations.csv"
  )
)

saveRDS(
  list(
    gene_level = tf_metrics,
    panel_a = panel_a,
    panel_b = panel_b,
    panel_c = panel_c,
    panel_d = panel_d,
    panel_e = panel_e,
    chen_summary = chen_summary,
    association_summary = association_summary
  ),
  file.path(intermediate_dir, "tf_regulatory_complexity_analysis.rds")
)

message("Fig. 5A-E statistical source data generated.")
