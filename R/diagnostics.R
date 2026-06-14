diagnose_microbiome_data <- function(counts,
                                     metadata = NULL,
                                     taxonomy = NULL,
                                     group = NULL,
                                     min_prevalence = 0.05,
                                     feature_read_n = 50) {
  if (any(!is.finite(counts))) {
    stop("`counts` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(counts < 0)) {
    stop("`counts` cannot contain negative values.", call. = FALSE)
  }

  library_sizes <- rowSums(counts)
  feature_reads <- colSums(counts)
  feature_prevalence <- colMeans(counts > 0)
  sample_zero_fraction <- rowMeans(counts == 0)
  positive_library_sizes <- library_sizes[library_sizes > 0]

  list(
    n_samples = nrow(counts),
    n_features = ncol(counts),
    taxa_sample_ratio = safe_ratio(ncol(counts), nrow(counts)),
    count_type = diagnose_count_type(counts, library_sizes),
    library_size = summarize_library_sizes(library_sizes, positive_library_sizes),
    reads = summarize_reads(
      library_sizes = library_sizes,
      feature_reads = feature_reads,
      feature_read_n = feature_read_n
    ),
    sparsity = summarize_sparsity(
      counts = counts,
      feature_prevalence = feature_prevalence,
      sample_zero_fraction = sample_zero_fraction,
      min_prevalence = min_prevalence
    ),
    metadata = summarize_metadata(metadata, group = group),
    taxonomy = summarize_taxonomy(taxonomy)
  )
}

diagnose_count_type <- function(counts, library_sizes) {
  integerish <- all(abs(counts - round(counts)) < sqrt(.Machine$double.eps))
  row_sum_1 <- all(abs(library_sizes - 1) < 1e-6)
  row_sum_100 <- all(abs(library_sizes - 100) < 1e-4)

  if (row_sum_1) {
    input_scale <- "relative_1"
  } else if (row_sum_100) {
    input_scale <- "relative_100"
  } else if (integerish) {
    input_scale <- "integer_counts"
  } else {
    input_scale <- "non_integer"
  }

  list(
    integerish = integerish,
    input_scale = input_scale,
    looks_relative = input_scale %in% c("relative_1", "relative_100"),
    has_zero_sum_samples = any(library_sizes == 0)
  )
}

summarize_library_sizes <- function(library_sizes, positive_library_sizes) {
  if (length(positive_library_sizes) == 0) {
    imbalance_ratio <- NA_real_
  } else {
    imbalance_ratio <- safe_ratio(
      max(positive_library_sizes),
      min(positive_library_sizes)
    )
  }

  list(
    min = unname(min(library_sizes)),
    q1 = unname(stats::quantile(library_sizes, 0.25, names = FALSE)),
    median = unname(stats::median(library_sizes)),
    mean = unname(mean(library_sizes)),
    q3 = unname(stats::quantile(library_sizes, 0.75, names = FALSE)),
    max = unname(max(library_sizes)),
    sd = unname(stats::sd(library_sizes)),
    cv = safe_ratio(stats::sd(library_sizes), mean(library_sizes)),
    imbalance_ratio = imbalance_ratio
  )
}

summarize_reads <- function(library_sizes, feature_reads, feature_read_n) {
  if (feature_read_n < 1) {
    stop("`feature_read_n` must be at least 1.", call. = FALSE)
  }

  n_features <- length(feature_reads)
  n_used <- min(feature_read_n, n_features)
  first_n_reads <- sum(feature_reads[seq_len(n_used)])
  top_n_reads <- sum(sort(feature_reads, decreasing = TRUE)[seq_len(n_used)])
  total_reads <- sum(feature_reads)

  min_sample_index <- which.min(library_sizes)
  max_sample_index <- which.max(library_sizes)

  list(
    total_reads = unname(total_reads),
    feature_read_n = n_used,
    first_n_reads = unname(first_n_reads),
    first_n_fraction = safe_ratio(first_n_reads, total_reads),
    top_n_reads = unname(top_n_reads),
    top_n_fraction = safe_ratio(top_n_reads, total_reads),
    min_sample = list(
      sample_id = names(library_sizes)[min_sample_index],
      sample_index = min_sample_index,
      reads = unname(library_sizes[min_sample_index])
    ),
    max_sample = list(
      sample_id = names(library_sizes)[max_sample_index],
      sample_index = max_sample_index,
      reads = unname(library_sizes[max_sample_index])
    )
  )
}

summarize_sparsity <- function(counts,
                               feature_prevalence,
                               sample_zero_fraction,
                               min_prevalence) {
  list(
    zero_fraction = unname(mean(counts == 0)),
    sample_zero_fraction_median = unname(stats::median(sample_zero_fraction)),
    feature_prevalence_median = unname(stats::median(feature_prevalence)),
    low_prevalence_fraction = unname(mean(feature_prevalence < min_prevalence)),
    min_prevalence = min_prevalence
  )
}

summarize_metadata <- function(metadata, group = NULL) {
  if (is.null(metadata)) {
    return(list(
      provided = FALSE,
      n_columns = 0,
      missing_fraction_by_column = NULL,
      constant_columns = character(),
      group = NULL
    ))
  }

  missing_fraction <- vapply(metadata, function(column) {
    mean(is.na(column) | trimws(as.character(column)) == "")
  }, numeric(1))

  constant_columns <- names(Filter(isTRUE, lapply(metadata, function(column) {
    length(unique(column[!(is.na(column) | trimws(as.character(column)) == "")])) <= 1
  })))

  group_summary <- NULL
  if (!is.null(group)) {
    if (!group %in% colnames(metadata)) {
      stop("`group` is not a column in `metadata`.", call. = FALSE)
    }

    group_values <- metadata[[group]]
    group_sizes <- sort(table(group_values, useNA = "ifany"))
    group_summary <- list(
      column = group,
      n_groups = length(group_sizes),
      sizes = group_sizes,
      min_n = unname(min(group_sizes)),
      max_n = unname(max(group_sizes))
    )
  }

  list(
    provided = TRUE,
    n_columns = ncol(metadata),
    missing_fraction_by_column = missing_fraction,
    constant_columns = constant_columns,
    group = group_summary
  )
}

summarize_taxonomy <- function(taxonomy) {
  if (is.null(taxonomy)) {
    return(list(
      provided = FALSE,
      missing_fraction_by_rank = NULL,
      assigned_features_by_rank = NULL,
      unique_taxa_by_rank = NULL
    ))
  }

  missing_terms <- c(
    "",
    "na",
    "n/a",
    "unknown",
    "unclassified",
    "uncultured",
    "unassigned"
  )

  missing_by_cell <- lapply(taxonomy, function(column) {
    values <- tolower(trimws(as.character(column)))
    is.na(column) | values %in% missing_terms
  })

  missing_fraction <- vapply(missing_by_cell, mean, numeric(1))

  assigned_features <- vapply(missing_by_cell, function(is_missing) {
    sum(!is_missing)
  }, integer(1))

  unique_taxa <- vapply(names(taxonomy), function(rank_name) {
    column <- taxonomy[[rank_name]]
    is_missing <- missing_by_cell[[rank_name]]
    length(unique(trimws(as.character(column[!is_missing]))))
  }, numeric(1))

  list(
    provided = TRUE,
    missing_fraction_by_rank = missing_fraction,
    assigned_features_by_rank = assigned_features,
    unique_taxa_by_rank = unique_taxa
  )
}

safe_ratio <- function(numerator, denominator) {
  if (length(denominator) == 0 || is.na(denominator) || denominator == 0) {
    return(NA_real_)
  }

  numerator / denominator
}
