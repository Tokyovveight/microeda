da_prepare_context <- function(x,
                               metadata = NULL,
                               taxonomy = NULL,
                               group,
                               contrast,
                               methods = c("aldex2", "ancombc2", "deseq2"),
                               tax_rank = NULL,
                               prevalence_filter = NULL,
                               min_count = NULL,
                               p_adjust_method = NULL,
                               taxa_are_rows = TRUE) {
  if (missing(group)) {
    stop("`group` is required.", call. = FALSE)
  }
  if (missing(contrast)) {
    stop("`contrast` is required.", call. = FALSE)
  }

  extracted <- microeda_extract(
    x = x,
    metadata = metadata,
    taxonomy = taxonomy,
    taxa_are_rows = taxa_are_rows
  )
  counts <- extracted$counts
  metadata <- extracted$metadata
  taxonomy <- extracted$taxonomy

  da_validate_counts(counts)
  group <- da_validate_group(metadata, group)
  group_values <- metadata[[group]]
  names(group_values) <- rownames(metadata)
  contrast <- da_validate_contrast(
    contrast = contrast,
    group_values = group_values,
    group = group
  )
  methods <- da_validate_methods(methods)
  tax_rank <- da_validate_tax_rank(tax_rank, taxonomy)
  filters <- da_validate_filters(
    prevalence_filter = prevalence_filter,
    min_count = min_count
  )
  p_adjust_method <- da_validate_p_adjust_method(p_adjust_method)
  contrast_label <- paste0(contrast[1], "_vs_", contrast[2])
  sample_ids <- rownames(counts)
  feature_ids <- colnames(counts)

  structure(
    list(
      counts = counts,
      metadata = metadata,
      taxonomy = taxonomy,
      group = group,
      contrast = contrast,
      contrast_label = contrast_label,
      methods = methods,
      tax_rank = tax_rank,
      filters = filters,
      p_adjust_method = p_adjust_method,
      feature_ids = feature_ids,
      sample_ids = sample_ids,
      group_values = group_values,
      caveats = da_context_caveats(
        counts = counts,
        group_values = group_values,
        methods = methods,
        tax_rank = tax_rank,
        taxonomy = taxonomy
      ),
      params = list(
        methods = methods,
        tax_rank = tax_rank,
        filters = filters,
        p_adjust_method = p_adjust_method,
        taxa_are_rows = taxa_are_rows
      ),
      call = match.call()
    ),
    class = "microeda_da_context"
  )
}

da_validate_methods <- function(methods) {
  supported_methods <- da_supported_methods()
  if (!is.character(methods) || length(methods) < 1 ||
      any(is.na(methods)) || any(!nzchar(methods))) {
    stop(
      "`methods` must contain one or more supported DA method IDs.",
      call. = FALSE
    )
  }

  if (anyDuplicated(methods)) {
    stop("`methods` cannot contain duplicate values.", call. = FALSE)
  }

  unknown_methods <- setdiff(methods, supported_methods)
  if (length(unknown_methods) > 0) {
    stop(
      "Unknown DA method(s): ",
      paste(unknown_methods, collapse = ", "),
      ". Supported methods: ",
      paste(supported_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  methods
}

da_validate_contrast <- function(contrast, group_values = NULL, group = "group") {
  if (!is.character(contrast) || length(contrast) != 2 ||
      any(is.na(contrast)) || any(!nzchar(contrast))) {
    stop(
      "`contrast` must be a length-2 character vector of group levels.",
      call. = FALSE
    )
  }

  if (identical(contrast[1], contrast[2])) {
    stop("`contrast` must contain two different group levels.", call. = FALSE)
  }

  if (!is.null(group_values)) {
    available_levels <- unique(as.character(group_values))
    missing_levels <- setdiff(contrast, available_levels)
    if (length(missing_levels) > 0) {
      stop(
        "`contrast` level(s) not found in `",
        group,
        "`: ",
        paste(missing_levels, collapse = ", "),
        ". Available levels: ",
        paste(available_levels, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  contrast
}

da_empty_result <- function() {
  data.frame(
    feature_id = character(),
    taxon_label = character(),
    rank = character(),
    method = character(),
    contrast = character(),
    group1 = character(),
    group2 = character(),
    effect = numeric(),
    effect_type = character(),
    log_fold_change = numeric(),
    statistic = numeric(),
    standard_error = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    p_value = numeric(),
    p_adjusted = numeric(),
    p_adjust_method = character(),
    p_adjust_scope = character(),
    significance = character(),
    direction = character(),
    method_note = character(),
    stringsAsFactors = FALSE
  )
}

da_standard_result <- function(...) {
  values <- list(...)
  out <- da_empty_result()
  if (length(values) == 0) {
    return(out)
  }

  if (is.null(names(values)) || any(!nzchar(names(values)))) {
    stop("`da_standard_result()` inputs must be named.", call. = FALSE)
  }

  unknown_columns <- setdiff(names(values), names(out))
  if (length(unknown_columns) > 0) {
    stop(
      "Unknown standardized DA result column(s): ",
      paste(unknown_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  n <- max(vapply(values, length, integer(1)))
  if (n == 0) {
    return(out)
  }

  columns <- lapply(names(out), function(column) {
    if (column %in% names(values)) {
      return(da_recycle_result_value(values[[column]], n, column))
    }

    if (is.numeric(out[[column]])) {
      return(rep(NA_real_, n))
    }
    rep(NA_character_, n)
  })
  names(columns) <- names(out)

  data.frame(columns, stringsAsFactors = FALSE, check.names = FALSE)
}

da_method_notes <- function(methods = da_supported_methods()) {
  methods <- da_validate_methods(methods)
  rows <- lapply(methods, da_method_note)
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_supported_methods <- function() {
  c("aldex2", "ancombc2", "deseq2")
}

da_validate_counts <- function(counts) {
  if (any(!is.finite(counts))) {
    stop("`counts` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(counts < 0)) {
    stop("`counts` cannot contain negative values.", call. = FALSE)
  }

  invisible(counts)
}

da_validate_group <- function(metadata, group) {
  if (is.null(metadata)) {
    stop("`metadata` is required for DA group validation.", call. = FALSE)
  }

  if (!is.character(group) || length(group) != 1 ||
      is.na(group) || !nzchar(group)) {
    stop("`group` must be a single non-missing character string.", call. = FALSE)
  }

  if (!group %in% colnames(metadata)) {
    stop("`group` is not a column in `metadata`.", call. = FALSE)
  }

  group_values <- metadata[[group]]
  group_labels <- as.character(group_values)
  if (any(is.na(group_labels)) || any(!nzchar(group_labels))) {
    stop("Group labels must be present for all samples.", call. = FALSE)
  }

  if (length(unique(group_labels)) < 2) {
    stop("`group` must contain at least two levels.", call. = FALSE)
  }

  group
}

da_validate_tax_rank <- function(tax_rank, taxonomy) {
  if (is.null(tax_rank)) {
    return(NULL)
  }

  if (!is.character(tax_rank) || length(tax_rank) != 1 ||
      is.na(tax_rank) || !nzchar(tax_rank)) {
    stop("`tax_rank` must be NULL or a single character string.", call. = FALSE)
  }

  if (!is.null(taxonomy) && !tax_rank %in% colnames(taxonomy)) {
    stop("`tax_rank` is not a column in `taxonomy`.", call. = FALSE)
  }

  tax_rank
}

da_validate_filters <- function(prevalence_filter = NULL, min_count = NULL) {
  if (!is.null(prevalence_filter) &&
      (!is.numeric(prevalence_filter) || length(prevalence_filter) != 1 ||
       is.na(prevalence_filter) || !is.finite(prevalence_filter) ||
       prevalence_filter < 0 || prevalence_filter > 1)) {
    stop("`prevalence_filter` must be NULL or a number in [0, 1].", call. = FALSE)
  }

  if (!is.null(min_count) &&
      (!is.numeric(min_count) || length(min_count) != 1 ||
       is.na(min_count) || !is.finite(min_count) ||
       min_count < 0 || min_count != floor(min_count))) {
    stop("`min_count` must be NULL or a non-negative whole number.", call. = FALSE)
  }

  list(
    prevalence_filter = prevalence_filter,
    min_count = if (is.null(min_count)) NULL else as.integer(min_count),
    applied = FALSE
  )
}

da_validate_p_adjust_method <- function(p_adjust_method) {
  if (is.null(p_adjust_method)) {
    return(NULL)
  }

  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1 ||
      is.na(p_adjust_method) || !p_adjust_method %in% stats::p.adjust.methods) {
    stop(
      "`p_adjust_method` must be NULL or one of stats::p.adjust.methods.",
      call. = FALSE
    )
  }

  p_adjust_method
}

da_context_caveats <- function(counts, group_values, methods, tax_rank, taxonomy) {
  rows <- list(da_caveat(
    method = NA_character_,
    caveat_id = "method_native_p_adjustment",
    topic = "differential_abundance",
    severity = "info",
    message = paste(
      "By default, DA backends should use method-native p-value adjustment;",
      "microeda should not globally re-adjust backend outputs."
    )
  ))

  group_sizes <- table(as.character(group_values))
  if (min(group_sizes) < 5) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = NA_character_,
      caveat_id = "small_group_size",
      topic = "group_design",
      severity = "warning",
      message = paste(
        "At least one contrast group has fewer than five samples;",
        "DA results may be unstable."
      )
    )
  }

  if (mean(counts == 0) >= 0.7) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = NA_character_,
      caveat_id = "high_sparsity",
      topic = "sparsity",
      severity = "warning",
      message = "The count table is highly sparse; DA method behavior may be sensitive to filtering and zero handling."
    )
  }

  if (!all(abs(counts - round(counts)) < sqrt(.Machine$double.eps))) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = NA_character_,
      caveat_id = "non_integer_counts",
      topic = "input",
      severity = "warning",
      message = "Counts are not integer-like; count-based DA methods may not be appropriate."
    )
  }

  if (!is.null(tax_rank) && is.null(taxonomy)) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = NA_character_,
      caveat_id = "taxonomy_unavailable",
      topic = "taxonomy",
      severity = "warning",
      message = "`tax_rank` was requested but no taxonomy table is available."
    )
  }

  rows <- c(rows, split(da_method_notes(methods), seq_along(methods)))
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_method_note <- function(method) {
  if (identical(method, "aldex2")) {
    return(da_caveat(
      method = method,
      caveat_id = "aldex2_compositional_note",
      topic = "differential_abundance",
      severity = "info",
      message = paste(
        "ALDEx2 is treated as a compositional-aware method;",
        "interpret results with its Monte Carlo Dirichlet assumptions."
      )
    ))
  }

  if (identical(method, "ancombc2")) {
    return(da_caveat(
      method = method,
      caveat_id = "ancombc2_compositional_note",
      topic = "differential_abundance",
      severity = "info",
      message = paste(
        "ANCOM-BC2 is treated as a compositional-aware method;",
        "interpret results with its model and bias-correction assumptions."
      )
    ))
  }

  da_caveat(
    method = method,
    caveat_id = "deseq2_sensitivity_note",
    topic = "differential_abundance",
    severity = "info",
    message = paste(
      "DESeq2 is treated as a comparison/sensitivity method;",
      "microbiome compositionality and sparsity can violate its assumptions."
    )
  )
}

da_caveat <- function(method, caveat_id, topic, severity, message) {
  data.frame(
    method = method,
    caveat_id = caveat_id,
    topic = topic,
    severity = severity,
    message = message,
    stringsAsFactors = FALSE
  )
}

da_recycle_result_value <- function(value, n, column) {
  if (length(value) == n) {
    return(value)
  }

  if (length(value) == 1) {
    return(rep(value, n))
  }

  stop(
    "`",
    column,
    "` must have length 1 or match the number of result rows.",
    call. = FALSE
  )
}
