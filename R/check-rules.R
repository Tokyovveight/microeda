check_rule_context <- function(diagnostics) {
  validate_check_diagnostics(diagnostics)

  metadata <- diagnostics$metadata
  taxonomy <- diagnostics$taxonomy
  group <- metadata$group

  summary <- data.frame(
    n_samples = diagnostics$n_samples,
    n_features = diagnostics$n_features,
    input_scale = diagnostics$count_type$input_scale,
    zero_fraction = diagnostics$sparsity$zero_fraction,
    library_size_imbalance_ratio = diagnostics$library_size$imbalance_ratio,
    has_group = !is.null(group),
    min_group_n = check_group_min_n(group),
    has_metadata = isTRUE(metadata$provided),
    max_metadata_missing_fraction = check_max_metadata_missing_fraction(metadata),
    has_taxonomy = isTRUE(taxonomy$provided),
    max_taxonomy_missing_fraction = check_max_taxonomy_missing_fraction(taxonomy),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    context = check_rule_context_rows(diagnostics)
  )
}

validate_check_diagnostics <- function(diagnostics) {
  required_fields <- c(
    "n_samples",
    "n_features",
    "count_type",
    "library_size",
    "sparsity",
    "metadata",
    "taxonomy"
  )

  if (!is.list(diagnostics) ||
      any(!required_fields %in% names(diagnostics))) {
    stop(
      "`diagnostics` must be a diagnostics object from microeda_check().",
      call. = FALSE
    )
  }

  invisible(diagnostics)
}

check_rule_context_rows <- function(diagnostics) {
  metadata <- diagnostics$metadata
  taxonomy <- diagnostics$taxonomy
  group <- metadata$group

  rows <- list(
    check_rule_context_row(
      context_id = "count_scale",
      topic = "input",
      metric = "input_scale",
      value = diagnostics$count_type$input_scale,
      label = "Input count scale classification."
    ),
    check_rule_context_row(
      context_id = "sample_count",
      topic = "input",
      metric = "n_samples",
      value = diagnostics$n_samples,
      numeric_value = diagnostics$n_samples,
      label = "Number of samples."
    ),
    check_rule_context_row(
      context_id = "feature_count",
      topic = "input",
      metric = "n_features",
      value = diagnostics$n_features,
      numeric_value = diagnostics$n_features,
      label = "Number of features."
    ),
    check_rule_context_row(
      context_id = "zero_fraction",
      topic = "sparsity",
      metric = "zero_fraction",
      value = diagnostics$sparsity$zero_fraction,
      numeric_value = diagnostics$sparsity$zero_fraction,
      label = "Fraction of count table entries equal to zero."
    ),
    check_rule_context_row(
      context_id = "library_size_imbalance",
      topic = "library_size",
      metric = "imbalance_ratio",
      value = diagnostics$library_size$imbalance_ratio,
      numeric_value = diagnostics$library_size$imbalance_ratio,
      label = "Maximum positive library size divided by minimum positive library size."
    ),
    check_rule_context_row(
      context_id = "group_available",
      topic = "group_design",
      metric = "has_group",
      value = !is.null(group),
      label = "Whether group diagnostics are available."
    ),
    check_rule_context_row(
      context_id = "metadata_available",
      topic = "metadata",
      metric = "has_metadata",
      value = isTRUE(metadata$provided),
      label = "Whether sample metadata was supplied."
    ),
    check_rule_context_row(
      context_id = "taxonomy_available",
      topic = "taxonomy",
      metric = "has_taxonomy",
      value = isTRUE(taxonomy$provided),
      label = "Whether taxonomy data was supplied."
    )
  )

  if (!is.null(group)) {
    rows[[length(rows) + 1L]] <- check_rule_context_row(
      context_id = "minimum_group_size",
      topic = "group_design",
      metric = "min_group_n",
      value = group$min_n,
      numeric_value = group$min_n,
      label = "Minimum observed group size."
    )
  }

  if (isTRUE(metadata$provided)) {
    max_missing <- check_max_metadata_missing_fraction(metadata)
    rows[[length(rows) + 1L]] <- check_rule_context_row(
      context_id = "metadata_missingness",
      topic = "metadata",
      metric = "max_missing_fraction",
      value = max_missing,
      numeric_value = max_missing,
      label = "Maximum missing fraction across metadata columns."
    )
  }

  if (isTRUE(taxonomy$provided)) {
    max_missing <- check_max_taxonomy_missing_fraction(taxonomy)
    rows[[length(rows) + 1L]] <- check_rule_context_row(
      context_id = "taxonomy_missingness",
      topic = "taxonomy",
      metric = "max_missing_fraction",
      value = max_missing,
      numeric_value = max_missing,
      label = "Maximum missing fraction across taxonomy ranks."
    )
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

check_rule_context_row <- function(context_id,
                                   topic,
                                   metric,
                                   value,
                                   label,
                                   numeric_value = NA_real_) {
  data.frame(
    context_id = context_id,
    topic = topic,
    metric = metric,
    value = as.character(value),
    numeric_value = numeric_value,
    label = label,
    stringsAsFactors = FALSE
  )
}

check_group_min_n <- function(group) {
  if (is.null(group)) {
    return(NA_real_)
  }

  group$min_n
}

check_max_metadata_missing_fraction <- function(metadata) {
  if (!isTRUE(metadata$provided) ||
      length(metadata$missing_fraction_by_column) == 0) {
    return(NA_real_)
  }

  unname(max(metadata$missing_fraction_by_column))
}

check_max_taxonomy_missing_fraction <- function(taxonomy) {
  if (!isTRUE(taxonomy$provided) ||
      length(taxonomy$missing_fraction_by_rank) == 0) {
    return(NA_real_)
  }

  unname(max(taxonomy$missing_fraction_by_rank))
}
