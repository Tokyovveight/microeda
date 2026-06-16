#' Build a tidy QC summary from a phyloseq object or count matrix
#'
#' `microeda_qc()` returns tidy per-sample, per-feature, per-rank, and
#' metadata-completeness tables derived directly from the extracted count,
#' metadata, and taxonomy data. It does not run the rule engine, so for
#' diagnostics and recommendations use [microeda_check()].
#'
#' @inheritParams microeda_check
#' @param min_prevalence Features with sample prevalence strictly below this
#'   threshold are reported but flagged as `above_prevalence_threshold = FALSE`
#'   in the per-feature table. The per-sample `n_features_above_prevalence`
#'   count is restricted to features at or above the threshold.
#'
#' @return A `microeda_qc` object with the following tidy data frames and
#'   compact summaries:
#'   \describe{
#'     \item{`per_sample`}{One row per sample with `sample_id`, `library_size`,
#'       `zero_fraction`, `n_features_detected`, `n_features_above_prevalence`.}
#'     \item{`per_feature`}{One row per feature with `feature_id`, `total_reads`,
#'       `prevalence`, `n_samples_detected`, `above_prevalence_threshold`.}
#'     \item{`library_size_summary`}{A named list with compact library-size
#'       diagnostics, including total reads, quartiles, zero-library samples,
#'       and imbalance ratios.}
#'     \item{`sparsity_summary`}{A named list with matrix-wide zero fraction,
#'       zero-library samples, zero-abundance features, and median sample/feature
#'       zero fractions.}
#'     \item{`prevalence_summary`}{A named list with compact prevalence and
#'       filtering diagnostics, including features above and below the
#'       `min_prevalence` threshold.}
#'     \item{`qc_flags`}{A data frame of conservative QC flags. Empty when no
#'       simple library-size, sparsity, or prevalence flags are triggered.}
#'     \item{`qc_observations`}{A data frame of concise human-readable QC
#'       observations derived from the compact summaries and QC flags.}
#'     \item{`per_rank`}{One row per taxonomy rank with `rank`, `n_assigned`,
#'       `n_unique`, `missing_fraction`. `NULL` when no taxonomy is supplied.}
#'     \item{`metadata_completeness`}{One row per metadata column with
#'       `column`, `missing_fraction`, `n_unique`, `is_constant`, `is_group`.
#'       `NULL` when no metadata is supplied.}
#'   }
#' @export
microeda_qc <- function(x,
                        metadata = NULL,
                        taxonomy = NULL,
                        group = NULL,
                        taxa_are_rows = TRUE,
                        min_prevalence = 0.05) {
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 ||
      is.na(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
    stop("`min_prevalence` must be a single number in [0, 1].", call. = FALSE)
  }

  extracted <- microeda_extract(x, metadata, taxonomy, taxa_are_rows)
  counts <- extracted$counts

  per_sample <- .qc_per_sample(counts, min_prevalence = min_prevalence)
  per_feature <- .qc_per_feature(counts, min_prevalence = min_prevalence)
  library_size_summary <- .qc_library_size_summary(counts)
  sparsity_summary <- .qc_sparsity_summary(counts)
  prevalence_summary <- .qc_prevalence_summary(
    per_feature = per_feature,
    min_prevalence = min_prevalence,
    n_samples = nrow(counts)
  )
  qc_flags <- .qc_flags(
    library_size_summary,
    sparsity_summary,
    prevalence_summary
  )
  per_rank <- .qc_per_rank(extracted$taxonomy)
  meta_complete <- .qc_metadata_completeness(extracted$metadata, group = group)
  qc_observations <- .qc_observations(
    library_size_summary = library_size_summary,
    sparsity_summary = sparsity_summary,
    prevalence_summary = prevalence_summary,
    qc_flags = qc_flags,
    per_rank = per_rank,
    metadata_completeness = meta_complete
  )

  structure(
    list(
      per_sample = per_sample,
      per_feature = per_feature,
      library_size_summary = library_size_summary,
      sparsity_summary = sparsity_summary,
      prevalence_summary = prevalence_summary,
      qc_flags = qc_flags,
      qc_observations = qc_observations,
      per_rank = per_rank,
      metadata_completeness = meta_complete,
      params = list(
        min_prevalence = min_prevalence,
        group = group
      ),
      call = match.call()
    ),
    class = "microeda_qc"
  )
}

.qc_observations <- function(library_size_summary,
                             sparsity_summary,
                             prevalence_summary,
                             qc_flags,
                             per_rank,
                             metadata_completeness) {
  observations <- list(
    .qc_observation(
      observation_id = "input_dimensions",
      category = "input",
      severity = "info",
      message = paste0(
        "Input contains ",
        library_size_summary$n_samples,
        " sample(s) and ",
        prevalence_summary$n_features,
        " feature(s)."
      )
    ),
    .qc_observation(
      observation_id = "total_reads",
      category = "library_size",
      severity = "info",
      message = paste0(
        "Total reads across all samples: ",
        library_size_summary$total_reads,
        "."
      )
    ),
    .qc_observation(
      observation_id = "overall_zero_fraction",
      category = "sparsity",
      severity = "info",
      message = paste0(
        "Overall zero fraction: ",
        .qc_format_percent(sparsity_summary$overall_zero_fraction),
        "."
      )
    ),
    .qc_observation(
      observation_id = "features_above_prevalence",
      category = "prevalence",
      severity = "info",
      message = paste0(
        prevalence_summary$n_features_above_threshold,
        " feature(s) are at or above the min_prevalence threshold (",
        prevalence_summary$min_prevalence_threshold,
        ")."
      )
    )
  )

  if (is.null(metadata_completeness)) {
    observations[[length(observations) + 1]] <- .qc_observation(
      observation_id = "metadata_absent",
      category = "metadata",
      severity = "info",
      message = "Metadata table is absent."
    )
  } else {
    observations[[length(observations) + 1]] <- .qc_observation(
      observation_id = "metadata_columns",
      category = "metadata",
      severity = "info",
      message = paste0(
        nrow(metadata_completeness),
        " metadata column(s) are available for completeness checks."
      )
    )
  }

  if (is.null(per_rank)) {
    observations[[length(observations) + 1]] <- .qc_observation(
      observation_id = "taxonomy_absent",
      category = "taxonomy",
      severity = "info",
      message = "Taxonomy table is absent."
    )
  } else {
    observations[[length(observations) + 1]] <- .qc_observation(
      observation_id = "taxonomy_ranks",
      category = "taxonomy",
      severity = "info",
      message = paste0(
        nrow(per_rank),
        " taxonomy rank(s) are available for completeness checks."
      )
    )
  }

  if (nrow(qc_flags) > 0) {
    flag_observations <- lapply(seq_len(nrow(qc_flags)), function(i) {
      flag <- qc_flags[i, , drop = FALSE]
      .qc_observation(
        observation_id = paste0("flag_", flag$flag_id),
        category = .qc_flag_category(flag$flag_id),
        severity = flag$severity,
        message = flag$message
      )
    })
    observations <- c(observations, flag_observations)
  }

  out <- do.call(rbind, observations)
  rownames(out) <- NULL
  out[!duplicated(out$observation_id), , drop = FALSE]
}

.qc_observation <- function(observation_id, category, severity, message) {
  data.frame(
    observation_id = observation_id,
    category = category,
    severity = severity,
    message = message,
    stringsAsFactors = FALSE
  )
}

.qc_flag_category <- function(flag_id) {
  if (flag_id %in% c("zero_library_samples", "library_size_imbalance")) {
    return("library_size")
  }

  if (flag_id %in% c("zero_abundance_features", "high_sparsity")) {
    return("sparsity")
  }

  if (flag_id %in% c(
    "many_features_below_prevalence",
    "many_single_sample_features"
  )) {
    return("prevalence")
  }

  "input"
}

.qc_format_percent <- function(x) {
  paste0(format(round(100 * x, 2), trim = TRUE, scientific = FALSE), "%")
}

.qc_prevalence_summary <- function(per_feature, min_prevalence, n_samples) {
  prevalence <- per_feature$prevalence
  n_features <- nrow(per_feature)
  above_threshold <- per_feature$above_prevalence_threshold
  n_features_above_threshold <- sum(above_threshold)
  n_features_below_threshold <- sum(!above_threshold)
  detected_in_one_sample <- per_feature$n_samples_detected == 1

  list(
    n_features = n_features,
    n_samples = n_samples,
    min_prevalence_threshold = min_prevalence,
    n_features_above_threshold = unname(n_features_above_threshold),
    n_features_below_threshold = unname(n_features_below_threshold),
    fraction_features_above_threshold = .qc_safe_ratio(
      n_features_above_threshold,
      n_features
    ),
    fraction_features_below_threshold = .qc_safe_ratio(
      n_features_below_threshold,
      n_features
    ),
    min_prevalence = unname(min(prevalence)),
    q1_prevalence = unname(stats::quantile(prevalence, 0.25, names = FALSE)),
    median_prevalence = unname(stats::median(prevalence)),
    mean_prevalence = unname(mean(prevalence)),
    q3_prevalence = unname(stats::quantile(prevalence, 0.75, names = FALSE)),
    max_prevalence = unname(max(prevalence)),
    n_features_detected_in_all_samples = unname(sum(
      per_feature$n_samples_detected == n_samples
    )),
    n_features_detected_in_one_sample = unname(sum(detected_in_one_sample)),
    fraction_features_detected_in_one_sample = .qc_safe_ratio(
      sum(detected_in_one_sample),
      n_features
    )
  )
}

.qc_library_size_summary <- function(counts) {
  library_sizes <- rowSums(counts)
  positive_library_sizes <- library_sizes[library_sizes > 0]
  max_library_size <- unname(max(library_sizes))
  median_library_size <- unname(stats::median(library_sizes))

  if (length(positive_library_sizes) == 0) {
    nonzero_min <- NA_real_
  } else {
    nonzero_min <- unname(min(positive_library_sizes))
  }

  list(
    n_samples = nrow(counts),
    total_reads = unname(sum(library_sizes)),
    min = unname(min(library_sizes)),
    q1 = unname(stats::quantile(library_sizes, 0.25, names = FALSE)),
    median = median_library_size,
    mean = unname(mean(library_sizes)),
    q3 = unname(stats::quantile(library_sizes, 0.75, names = FALSE)),
    max = max_library_size,
    zero_library_samples = unname(sum(library_sizes == 0)),
    nonzero_min = nonzero_min,
    max_to_median_ratio = .qc_safe_ratio(max_library_size, median_library_size),
    max_to_min_nonzero_ratio = .qc_safe_ratio(max_library_size, nonzero_min)
  )
}

.qc_sparsity_summary <- function(counts) {
  library_sizes <- rowSums(counts)
  feature_totals <- colSums(counts)
  zero_library_samples <- sum(library_sizes == 0)
  zero_abundance_features <- sum(feature_totals == 0)

  list(
    overall_zero_fraction = unname(mean(counts == 0)),
    zero_library_samples = unname(zero_library_samples),
    zero_library_sample_fraction = .qc_safe_ratio(zero_library_samples, nrow(counts)),
    zero_abundance_features = unname(zero_abundance_features),
    zero_abundance_feature_fraction = .qc_safe_ratio(
      zero_abundance_features,
      ncol(counts)
    ),
    median_sample_zero_fraction = unname(stats::median(rowMeans(counts == 0))),
    median_feature_zero_fraction = unname(stats::median(colMeans(counts == 0)))
  )
}

.qc_flags <- function(library_size_summary,
                      sparsity_summary,
                      prevalence_summary) {
  flags <- list()

  if (library_size_summary$zero_library_samples > 0) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "zero_library_samples",
      severity = "warning",
      message = paste0(
        library_size_summary$zero_library_samples,
        " sample(s) have zero total counts."
      )
    )
  }

  if (!is.na(library_size_summary$max_to_median_ratio) &&
      library_size_summary$max_to_median_ratio >= 10) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "library_size_imbalance",
      severity = "warning",
      message = "Maximum library size is at least 10 times the median."
    )
  }

  if (sparsity_summary$zero_abundance_features > 0) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "zero_abundance_features",
      severity = "warning",
      message = paste0(
        sparsity_summary$zero_abundance_features,
        " feature(s) have zero total abundance."
      )
    )
  }

  if (sparsity_summary$overall_zero_fraction >= 0.7) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "high_sparsity",
      severity = "warning",
      message = "At least 70% of count matrix entries are zero."
    )
  }

  if (prevalence_summary$fraction_features_below_threshold >= 0.5) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "many_features_below_prevalence",
      severity = "warning",
      message = paste0(
        "At least 50% of features are below the min_prevalence threshold (",
        prevalence_summary$min_prevalence_threshold,
        ")."
      )
    )
  }

  if (prevalence_summary$fraction_features_detected_in_one_sample >= 0.3) {
    flags[[length(flags) + 1]] <- .qc_flag(
      flag_id = "many_single_sample_features",
      severity = "warning",
      message = "At least 30% of features are detected in only one sample."
    )
  }

  if (length(flags) == 0) {
    return(data.frame(
      flag_id = character(),
      severity = character(),
      message = character(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, flags)
}

.qc_flag <- function(flag_id, severity, message) {
  data.frame(
    flag_id = flag_id,
    severity = severity,
    message = message,
    stringsAsFactors = FALSE
  )
}

.qc_safe_ratio <- function(numerator, denominator) {
  if (length(denominator) == 0 || is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}

.qc_per_sample <- function(counts, min_prevalence) {
  feature_prevalence <- colMeans(counts > 0)
  keep_features <- feature_prevalence >= min_prevalence
  counts_above <- counts[, keep_features, drop = FALSE]

  data.frame(
    sample_id = rownames(counts),
    library_size = unname(rowSums(counts)),
    zero_fraction = unname(rowMeans(counts == 0)),
    n_features_detected = unname(rowSums(counts > 0)),
    n_features_above_prevalence = unname(rowSums(counts_above > 0)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.qc_per_feature <- function(counts, min_prevalence) {
  prevalence <- colMeans(counts > 0)
  data.frame(
    feature_id = colnames(counts),
    total_reads = unname(colSums(counts)),
    prevalence = unname(prevalence),
    n_samples_detected = unname(colSums(counts > 0)),
    above_prevalence_threshold = unname(prevalence >= min_prevalence),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.qc_per_rank <- function(taxonomy) {
  if (is.null(taxonomy)) {
    return(NULL)
  }

  missing_terms <- c(
    "", "na", "n/a", "unknown", "unclassified", "uncultured", "unassigned"
  )

  is_missing_col <- function(column) {
    values <- tolower(trimws(as.character(column)))
    is.na(column) | values %in% missing_terms
  }

  ranks <- colnames(taxonomy)
  n_assigned <- vapply(taxonomy, function(col) {
    sum(!is_missing_col(col))
  }, integer(1))

  n_unique <- vapply(taxonomy, function(col) {
    missing <- is_missing_col(col)
    length(unique(trimws(as.character(col[!missing]))))
  }, integer(1))

  missing_fraction <- vapply(taxonomy, function(col) {
    mean(is_missing_col(col))
  }, numeric(1))

  data.frame(
    rank = ranks,
    n_assigned = n_assigned,
    n_unique = n_unique,
    missing_fraction = missing_fraction,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.qc_metadata_completeness <- function(metadata, group = NULL) {
  if (is.null(metadata)) {
    return(NULL)
  }

  missing_fraction <- vapply(metadata, function(col) {
    mean(is.na(col) | trimws(as.character(col)) == "")
  }, numeric(1))

  n_unique <- vapply(metadata, function(col) {
    length(unique(col[!(is.na(col) | trimws(as.character(col)) == "")]))
  }, integer(1))

  is_constant <- n_unique <= 1

  is_group <- if (is.null(group)) {
    rep(FALSE, ncol(metadata))
  } else {
    colnames(metadata) == group
  }

  data.frame(
    column = colnames(metadata),
    missing_fraction = missing_fraction,
    n_unique = n_unique,
    is_constant = is_constant,
    is_group = is_group,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
