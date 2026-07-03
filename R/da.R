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
  contrast_plan <- da_validate_contrast(
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
  contrast_label <- da_contrast_label(contrast_plan)
  sample_ids <- rownames(counts)
  feature_ids <- colnames(counts)

  structure(
    list(
      counts = counts,
      metadata = metadata,
      taxonomy = taxonomy,
      group = group,
      contrast = contrast,
      contrast_plan = contrast_plan,
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
        contrast_plan = contrast_plan,
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
  if (!is.character(contrast) || length(contrast) < 1 ||
      any(is.na(contrast)) || any(!nzchar(contrast))) {
    stop(
      "`contrast` must be a length-2 character vector of group levels or \"pairwise\".",
      call. = FALSE
    )
  }

  if (length(contrast) == 1) {
    if (identical(contrast, "pairwise")) {
      return(da_pairwise_contrast_plan(group_values = group_values, group = group))
    }

    stop(
      "`contrast` must be a length-2 character vector of group levels or \"pairwise\".",
      call. = FALSE
    )
  }

  if (length(contrast) != 2) {
    stop(
      "`contrast` must be a length-2 character vector of group levels or \"pairwise\".",
      call. = FALSE
    )
  }

  if (identical(contrast[1], contrast[2])) {
    stop("`contrast` must contain two different group levels.", call. = FALSE)
  }

  available_levels <- NULL
  if (!is.null(group_values)) {
    available_levels <- da_group_levels(group_values)
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

  data.frame(
    contrast = paste0(contrast[1], "_vs_", contrast[2]),
    group1 = contrast[1],
    group2 = contrast[2],
    contrast_type = "explicit",
    stringsAsFactors = FALSE
  )
}

da_pairwise_contrast_plan <- function(group_values, group = "group") {
  if (is.null(group_values)) {
    stop(
      "`contrast = \"pairwise\"` requires group labels.",
      call. = FALSE
    )
  }

  levels <- da_group_levels(group_values)
  if (length(levels) < 2) {
    stop(
      "`contrast = \"pairwise\"` requires at least two levels in `",
      group,
      "`.",
      call. = FALSE
    )
  }

  pairs <- utils::combn(levels, 2)
  data.frame(
    contrast = paste0(pairs[1, ], "_vs_", pairs[2, ]),
    group1 = pairs[1, ],
    group2 = pairs[2, ],
    contrast_type = "pairwise",
    stringsAsFactors = FALSE
  )
}

da_group_levels <- function(group_values) {
  group_labels <- as.character(group_values)
  group_labels <- group_labels[!is.na(group_labels) & nzchar(group_labels)]
  unique(group_labels)
}

da_contrast_label <- function(contrast_plan) {
  if (nrow(contrast_plan) == 1 &&
      identical(contrast_plan$contrast_type, "explicit")) {
    return(contrast_plan$contrast)
  }

  "pairwise"
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

da_run_aldex2 <- function(context,
                          mc.samples = 128,
                          denom = "all",
                          paired.test = FALSE) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }

  if (!"aldex2" %in% context$methods) {
    stop("`context$methods` must include \"aldex2\".", call. = FALSE)
  }

  if (nrow(context$contrast_plan) != 1 ||
      !identical(context$contrast_plan$contrast_type, "explicit")) {
    stop(
      "ALDEx2 pairwise contrast execution is not implemented yet.",
      call. = FALSE
    )
  }

  if (!da_optional_package_available("ALDEx2")) {
    stop(
      "`da_run_aldex2()` requires the optional package `ALDEx2`.",
      call. = FALSE
    )
  }

  params <- da_validate_aldex2_params(
    mc.samples = mc.samples,
    denom = denom,
    paired.test = paired.test
  )

  da_run_aldex2_contrast(
    context = context,
    contrast_row = context$contrast_plan[1, , drop = FALSE],
    params = params
  )
}

da_run_aldex2_contrast <- function(context, contrast_row, params) {
  group_values <- as.character(context$group_values)
  sample_ids <- names(context$group_values)
  group1 <- contrast_row$group1
  group2 <- contrast_row$group2
  keep <- group_values %in% c(group1, group2)
  selected_samples <- sample_ids[keep]
  conditions <- group_values[keep]

  if (!all(c(group1, group2) %in% conditions)) {
    stop("Both contrast groups must have samples for ALDEx2.", call. = FALSE)
  }

  conditions <- factor(conditions, levels = c(group1, group2))
  reads <- t(context$counts[selected_samples, , drop = FALSE])
  warnings <- character()
  messages <- character()
  clr <- NULL
  ttest <- NULL
  effect <- NULL

  stdout <- utils::capture.output({
    withCallingHandlers(
      {
        clr <- ALDEx2::aldex.clr(
          reads,
          as.character(conditions),
          mc.samples = params$mc.samples,
          denom = params$denom,
          verbose = FALSE
        )
        ttest <- ALDEx2::aldex.ttest(
          clr,
          paired.test = params$paired.test,
          verbose = FALSE
        )
        effect <- ALDEx2::aldex.effect(
          clr,
          verbose = FALSE,
          paired.test = params$paired.test
        )
      },
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  })
  messages <- c(messages, stdout[nzchar(stdout)])

  ttest <- as.data.frame(ttest, stringsAsFactors = FALSE, check.names = FALSE)
  effect <- as.data.frame(effect, stringsAsFactors = FALSE, check.names = FALSE)
  combined <- da_combine_aldex2_tables(ttest = ttest, effect = effect)
  results <- da_standardize_aldex2_result(
    combined = combined,
    context = context,
    contrast_row = contrast_row
  )

  raw_output <- list(
    clr = clr,
    ttest = ttest,
    effect = effect,
    combined = combined,
    conditions = as.character(conditions),
    contrast_row = contrast_row,
    input_orientation = "feature_by_sample",
    transposed_from_context = TRUE,
    params = params,
    warnings = warnings,
    messages = messages
  )

  da_backend_result(
    method = "aldex2",
    results = results,
    raw_output = raw_output,
    notes = da_method_notes("aldex2"),
    params = params
  )
}

da_standardize_aldex2_result <- function(combined, context, contrast_row) {
  feature_ids <- rownames(combined)
  p_adjusted <- da_column_or_na(combined, "we.eBH")
  method_note <- da_method_note("aldex2")$message

  da_standard_result(
    feature_id = feature_ids,
    taxon_label = da_taxon_labels(context, feature_ids),
    rank = da_result_rank(context),
    method = "aldex2",
    contrast = contrast_row$contrast,
    group1 = contrast_row$group1,
    group2 = contrast_row$group2,
    effect = da_column_or_na(combined, "effect"),
    effect_type = "aldex2_effect",
    log_fold_change = NA_real_,
    statistic = NA_real_,
    standard_error = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = da_column_or_na(combined, "we.ep"),
    p_adjusted = p_adjusted,
    p_adjust_method = "aldex2_native_BH",
    p_adjust_scope = "method_contrast",
    significance = da_p_significance(p_adjusted),
    direction = NA_character_,
    method_note = method_note
  )
}

da_backend_result <- function(method,
                              results = da_empty_result(),
                              raw_output = NULL,
                              notes = NULL,
                              params = list()) {
  method <- da_validate_backend_method(method)
  if (is.null(notes)) {
    notes <- da_method_notes(method)
  }

  out <- structure(
    list(
      method = method,
      results = results,
      raw_output = raw_output,
      notes = notes,
      params = params
    ),
    class = "microeda_da_backend_result"
  )

  da_validate_backend_result(out)
}

da_validate_backend_result <- function(x) {
  if (!inherits(x, "microeda_da_backend_result")) {
    stop("`x` must be a microeda_da_backend_result object.", call. = FALSE)
  }

  expected_fields <- c("method", "results", "raw_output", "notes", "params")
  missing_fields <- setdiff(expected_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "`x` is missing backend result field(s): ",
      paste(missing_fields, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  da_validate_backend_method(x$method)
  da_validate_standard_result_table(x$results)
  da_validate_note_table(x$notes)
  if (!is.list(x$params)) {
    stop("`x$params` must be a list.", call. = FALSE)
  }

  invisible(x)
}

da_standardize_backend_result <- function(x) {
  da_validate_backend_result(x)
  x$results
}

da_combine_method_results <- function(backend_results) {
  if (!is.list(backend_results) || length(backend_results) < 1) {
    stop("`backend_results` must be a non-empty list.", call. = FALSE)
  }

  backend_results <- lapply(backend_results, function(result) {
    da_validate_backend_result(result)
    result
  })
  methods <- unname(vapply(
    backend_results,
    function(result) result$method,
    character(1)
  ))
  if (anyDuplicated(methods)) {
    stop("`backend_results` cannot contain duplicate methods.", call. = FALSE)
  }

  results <- do.call(
    rbind,
    lapply(backend_results, da_standardize_backend_result)
  )
  row.names(results) <- NULL

  names(backend_results) <- methods
  raw_outputs <- lapply(backend_results, function(result) result$raw_output)
  names(raw_outputs) <- methods

  caveats <- da_combine_backend_notes(backend_results)

  list(
    results = results,
    method_results = backend_results,
    raw_outputs = raw_outputs,
    caveats = caveats,
    methods = methods
  )
}

da_build_result_object <- function(context, backend_results) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }

  combined <- da_combine_method_results(backend_results)
  if (!identical(combined$methods, context$methods)) {
    stop(
      "`backend_results` methods must match `context$methods` in order.",
      call. = FALSE
    )
  }

  caveats <- da_deduplicate_caveats(rbind(context$caveats, combined$caveats))

  structure(
    list(
      results = combined$results,
      method_results = combined$method_results,
      raw_outputs = combined$raw_outputs,
      methods = combined$methods,
      group = context$group,
      contrast = context$contrast,
      contrast_plan = context$contrast_plan,
      contrast_label = context$contrast_label,
      tax_rank = context$tax_rank,
      feature_metadata = da_feature_metadata(context),
      filters = context$filters,
      caveats = caveats,
      params = list(
        context = context$params,
        backend = lapply(combined$method_results, function(result) result$params)
      ),
      call = match.call()
    ),
    class = "microeda_da"
  )
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

da_validate_backend_method <- function(method) {
  if (!is.character(method) || length(method) != 1 ||
      is.na(method) || !nzchar(method)) {
    stop("`method` must be a single supported DA method ID.", call. = FALSE)
  }

  da_validate_methods(method)
}

da_validate_standard_result_table <- function(results) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame.", call. = FALSE)
  }

  expected_columns <- names(da_empty_result())
  if (!identical(names(results), expected_columns)) {
    stop(
      "`results` must have exactly the standardized DA result columns.",
      call. = FALSE
    )
  }

  invisible(results)
}

da_validate_note_table <- function(notes) {
  if (is.null(notes)) {
    return(invisible(notes))
  }

  if (!is.data.frame(notes)) {
    stop("`notes` must be NULL or a data frame.", call. = FALSE)
  }

  expected_columns <- names(da_empty_caveats())
  if (!identical(names(notes), expected_columns)) {
    stop(
      "`notes` must have columns: ",
      paste(expected_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(notes)
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

da_validate_aldex2_params <- function(mc.samples, denom, paired.test) {
  if (!is.numeric(mc.samples) || length(mc.samples) != 1 ||
      is.na(mc.samples) || !is.finite(mc.samples) ||
      mc.samples < 1 || mc.samples != floor(mc.samples)) {
    stop("`mc.samples` must be a positive whole number.", call. = FALSE)
  }

  if (!is.character(denom) || length(denom) != 1 ||
      is.na(denom) || !nzchar(denom)) {
    stop("`denom` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.logical(paired.test) || length(paired.test) != 1 ||
      is.na(paired.test)) {
    stop("`paired.test` must be TRUE or FALSE.", call. = FALSE)
  }

  list(
    mc.samples = as.integer(mc.samples),
    denom = denom,
    paired.test = paired.test
  )
}

da_optional_package_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
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

da_empty_caveats <- function() {
  data.frame(
    method = character(),
    caveat_id = character(),
    topic = character(),
    severity = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

da_deduplicate_caveats <- function(caveats) {
  if (is.null(caveats)) {
    return(da_empty_caveats())
  }

  da_validate_note_table(caveats)
  if (nrow(caveats) == 0) {
    return(caveats)
  }

  out <- caveats[!duplicated(caveats[names(da_empty_caveats())]), , drop = FALSE]
  row.names(out) <- NULL
  out
}

da_combine_backend_notes <- function(backend_results) {
  rows <- lapply(backend_results, function(result) {
    notes <- result$notes
    if (is.null(notes) || nrow(notes) == 0) {
      return(da_empty_caveats())
    }

    missing_method <- is.na(notes$method) | !nzchar(notes$method)
    notes$method[missing_method] <- result$method
    notes
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_feature_metadata <- function(context) {
  out <- data.frame(
    feature_id = context$feature_ids,
    stringsAsFactors = FALSE
  )

  if (!is.null(context$taxonomy)) {
    out <- cbind(
      out,
      as.data.frame(context$taxonomy, stringsAsFactors = FALSE)
    )
  }

  row.names(out) <- NULL
  out
}

da_combine_aldex2_tables <- function(ttest, effect) {
  if (is.null(rownames(ttest)) || is.null(rownames(effect)) ||
      any(!nzchar(rownames(ttest))) || any(!nzchar(rownames(effect)))) {
    stop("ALDEx2 output must include feature row names.", call. = FALSE)
  }

  missing_effect <- setdiff(rownames(ttest), rownames(effect))
  if (length(missing_effect) > 0) {
    stop("ALDEx2 t-test and effect outputs have different features.", call. = FALSE)
  }

  effect <- effect[rownames(ttest), , drop = FALSE]
  out <- data.frame(ttest, effect, check.names = FALSE)
  row.names(out) <- rownames(ttest)
  out
}

da_column_or_na <- function(data, column) {
  if (column %in% names(data)) {
    return(as.numeric(data[[column]]))
  }

  rep(NA_real_, nrow(data))
}

da_taxon_labels <- function(context, feature_ids) {
  if (is.null(context$tax_rank) || is.null(context$taxonomy)) {
    return(rep(NA_character_, length(feature_ids)))
  }

  taxonomy <- as.data.frame(context$taxonomy, stringsAsFactors = FALSE)
  taxon_values <- rep(NA_character_, length(feature_ids))
  matched <- match(feature_ids, rownames(taxonomy))
  found <- !is.na(matched)
  taxon_values[found] <- as.character(taxonomy[matched[found], context$tax_rank])
  taxon_values
}

da_result_rank <- function(context) {
  if (is.null(context$tax_rank)) {
    return(NA_character_)
  }

  context$tax_rank
}

da_p_significance <- function(p) {
  out <- rep("ns", length(p))
  out[is.na(p)] <- NA_character_
  out[!is.na(p) & p <= 0.05] <- "*"
  out[!is.na(p) & p <= 0.01] <- "**"
  out[!is.na(p) & p <= 0.001] <- "***"
  out
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
