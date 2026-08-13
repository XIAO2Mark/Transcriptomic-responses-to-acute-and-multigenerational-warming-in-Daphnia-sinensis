require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Required package '", package, "' is not installed.\n",
      "See README.md for installation requirements."
    )
  }
}

assert_columns <- function(data, required, label = deparse(substitute(data))) {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      label, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

read_csv_checked <- function(path, required = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(required)) assert_columns(data, required, basename(path))
  data
}

write_result <- function(data, filename, row.names = FALSE) {
  path <- file.path(output_dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(data, path, row.names = row.names, na = "")
  message("Wrote: ", path)
  invisible(path)
}

population_z <- function(x) {
  x <- as.numeric(x)
  center <- mean(x, na.rm = TRUE)
  scale_value <- sqrt(mean((x - center)^2, na.rm = TRUE))
  if (!is.finite(scale_value) || scale_value == 0) {
    stop("Cannot standardize a constant or non-finite variable.")
  }
  (x - center) / scale_value
}

significance_label <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.001, "***",
      ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
    )
  )
}

raw_odds_ratio <- function(a, b, c, d) {
  if (b * c == 0) return(NA_real_)
  (a * d) / (b * c)
}

cluster_vcov_glm <- function(fit, cluster) {
  x <- stats::model.matrix(fit)
  y <- stats::model.response(stats::model.frame(fit))
  mu <- stats::fitted(fit)
  w <- as.numeric(mu * (1 - mu))

  if (length(cluster) != nrow(x)) {
    stop("Cluster vector length does not match the fitted model.")
  }
  cluster <- factor(cluster)

  bread <- solve(crossprod(x, x * w))
  score_rows <- x * as.numeric(y - mu)
  score_by_cluster <- rowsum(score_rows, cluster, reorder = FALSE)
  meat <- crossprod(score_by_cluster)

  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))

  bread %*% meat %*% bread * correction
}

fit_clustered_membership_model <- function(
    data,
    outcome,
    cluster_column,
    pattern_name) {
  model_data <- data
  model_data$Target <- as.integer(outcome)

  fit <- stats::glm(
    Target ~ baseline_z + tf_load_z + promoter_gc_z,
    data = model_data,
    family = stats::binomial()
  )

  vcov_clustered <- cluster_vcov_glm(
    fit,
    model_data[[cluster_column]]
  )
  coefficient_names <- c("baseline_z", "tf_load_z")

  rows <- lapply(coefficient_names, function(term) {
    beta <- stats::coef(fit)[term]
    standard_error <- sqrt(vcov_clustered[term, term])
    z_value <- beta / standard_error
    p_value <- 2 * stats::pnorm(-abs(z_value))

    data.frame(
      Pattern = pattern_name,
      Predictor = ifelse(
        term == "baseline_z", "Ancestral expression", "TF load"
      ),
      n_genes = nrow(model_data),
      n_pattern = sum(model_data$Target),
      covariance = "100-kb genomic-block clustered",
      model_adjustment = "ancestral expression + promoter GC",
      beta = unname(beta),
      standard_error = unname(standard_error),
      odds_ratio = exp(beta),
      ci_low = exp(beta - 1.96 * standard_error),
      ci_high = exp(beta + 1.96 * standard_error),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

