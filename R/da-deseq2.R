da_validate_deseq2_public_arguments <- function(mc_samples_supplied,
                                                 denom_supplied,
                                                 paired.test,
                                                 pair_id) {
  if (!is.logical(paired.test) || length(paired.test) != 1 ||
      is.na(paired.test)) {
    stop("`paired.test` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(paired.test)) {
    stop(
      "Paired/repeated DESeq2 execution is not implemented; ",
      "`paired.test = TRUE` can currently be used only when ALDEx2 is the ",
      "sole requested method.",
      call. = FALSE
    )
  }
  if (!is.null(pair_id)) {
    stop(
      "Paired/repeated DESeq2 execution is not implemented; `pair_id` can ",
      "currently be supplied only when ALDEx2 is the sole requested method.",
      call. = FALSE
    )
  }
  if (isTRUE(mc_samples_supplied)) {
    stop(
      "`mc.samples` is an ALDEx2-specific argument and cannot be supplied ",
      "when DESeq2 is requested without ALDEx2.",
      call. = FALSE
    )
  }
  if (isTRUE(denom_supplied)) {
    stop(
      "`denom` is an ALDEx2-specific argument and cannot be supplied when ",
      "DESeq2 is requested without ALDEx2.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

da_run_deseq2 <- function(context, contrast_row = NULL, params) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }
  if (!"deseq2" %in% context$methods) {
    stop("`context$methods` must include \"deseq2\".", call. = FALSE)
  }
  if (!da_optional_package_available("DESeq2")) {
    stop(
      "`da_run_deseq2()` requires the optional package `DESeq2`. ",
      "Install it with `BiocManager::install(\"DESeq2\")`.",
      call. = FALSE
    )
  }

  params <- da_validate_deseq2_params(params)
  if (!is.null(contrast_row)) {
    return(da_run_deseq2_contrast(
      context = context,
      contrast_row = contrast_row,
      params = params
    ))
  }

  contrast_plan <- context$contrast_plan
  if (nrow(contrast_plan) == 1 &&
      identical(contrast_plan$contrast_type, "explicit")) {
    return(da_run_deseq2_contrast(
      context = context,
      contrast_row = contrast_plan[1, , drop = FALSE],
      params = params
    ))
  }
  if (!all(contrast_plan$contrast_type == "pairwise")) {
    stop(
      "DESeq2 supports explicit or pairwise contrast plans only.",
      call. = FALSE
    )
  }

  contrast_results <- lapply(seq_len(nrow(contrast_plan)), function(i) {
    contrast_row <- contrast_plan[i, , drop = FALSE]
    tryCatch(
      da_run_deseq2_contrast(
        context = context,
        contrast_row = contrast_row,
        params = params
      ),
      error = function(e) {
        stop(
          "DESeq2 pairwise execution failed for contrast `",
          contrast_row$contrast,
          "`: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  })
  contrast_labels <- contrast_plan$contrast
  names(contrast_results) <- contrast_labels

  results <- do.call(
    rbind,
    lapply(contrast_results, function(result) result$results)
  )
  row.names(results) <- NULL

  raw_output <- lapply(contrast_results, function(result) result$raw_output)
  names(raw_output) <- contrast_labels

  note_rows <- lapply(seq_along(contrast_results), function(i) {
    da_deseq2_pairwise_notes(
      notes = contrast_results[[i]]$notes,
      contrast = contrast_labels[[i]]
    )
  })
  notes <- da_deduplicate_caveats(do.call(rbind, note_rows))
  contrast_params <- lapply(contrast_results, function(result) result$params)
  names(contrast_params) <- contrast_labels

  da_backend_result(
    method = "deseq2",
    results = results,
    raw_output = raw_output,
    notes = notes,
    params = list(
      sf_type = params$sf_type,
      fit_type = params$fit_type,
      independent_filtering = params$independent_filtering,
      alpha = params$alpha,
      p_adjust_method = "BH",
      contrast_plan = contrast_plan,
      contrasts = contrast_params
    )
  )
}

da_run_deseq2_contrast <- function(context, contrast_row, params) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }
  if (!"deseq2" %in% context$methods) {
    stop("`context$methods` must include \"deseq2\".", call. = FALSE)
  }
  if (!da_optional_package_available("DESeq2")) {
    stop(
      "`da_run_deseq2_contrast()` requires the optional package `DESeq2`. ",
      "Install it with `BiocManager::install(\"DESeq2\")`.",
      call. = FALSE
    )
  }

  params <- da_validate_deseq2_params(params)
  contrast_row <- da_validate_deseq2_contrast_row(contrast_row)
  input <- da_prepare_deseq2_input(context, contrast_row)
  deseq_params <- da_deseq2_fit_params(params)
  results_params <- da_deseq2_results_params(input, params)

  warnings <- character()
  messages <- character()
  dds <- NULL
  native_results <- NULL
  backend_error <- NULL
  stdout <- utils::capture.output({
    tryCatch(
      withCallingHandlers(
        {
          dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = input$counts,
            colData = input$metadata,
            design = input$design
          )
          dds <- do.call(
            DESeq2::DESeq,
            c(list(object = dds), deseq_params)
          )
          native_results <- do.call(
            DESeq2::results,
            c(list(object = dds), results_params)
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
      ),
      error = function(e) {
        backend_error <<- conditionMessage(e)
      }
    )
  })
  messages <- unique(c(messages, stdout[nzchar(trimws(stdout))]))
  warnings <- unique(warnings)

  if (!is.null(backend_error)) {
    stop(
      "DESeq2 failed for contrast `",
      contrast_row$contrast,
      "`: ",
      backend_error,
      call. = FALSE
    )
  }
  if (is.null(dds) || is.null(native_results)) {
    stop(
      "DESeq2 returned an incomplete result for contrast `",
      contrast_row$contrast,
      "`.",
      call. = FALSE
    )
  }

  results_table <- as.data.frame(native_results)
  results <- da_standardize_deseq2_result(
    results_table = results_table,
    context = context,
    contrast_row = contrast_row,
    feature_order = input$feature_order
  )

  actual_fit_type <- da_deseq2_actual_fit_type(dds, params$fit_type)
  dds_mcols <- as.data.frame(da_deseq2_mcols(dds))
  results_metadata <- da_deseq2_metadata(native_results)
  results_mcols <- as.data.frame(da_deseq2_mcols(native_results))
  assay_names <- da_deseq2_assay_names(dds)
  cooks <- if ("cooks" %in% assay_names) {
    da_deseq2_assay(dds, "cooks")
  } else {
    NULL
  }
  cooks_cutoff <- da_deseq2_cooks_cutoff(dds)
  replacement_diagnostics <- da_deseq2_replacement_diagnostics(
    dds_mcols = dds_mcols,
    assay_names = assay_names
  )
  feature_diagnostics <- da_deseq2_feature_diagnostics(
    results_table = results_table,
    dds_mcols = dds_mcols,
    cooks_cutoff = cooks_cutoff,
    independent_filtering = params$independent_filtering
  )
  notes <- da_deseq2_notes(
    warnings = warnings,
    contrast = contrast_row$contrast,
    requested_fit_type = params$fit_type,
    actual_fit_type = actual_fit_type,
    feature_diagnostics = feature_diagnostics,
    replacement_diagnostics = replacement_diagnostics
  )

  backend_params <- list(
    sf_type = params$sf_type,
    requested_fit_type = params$fit_type,
    actual_fit_type = actual_fit_type,
    independent_filtering = params$independent_filtering,
    alpha = params$alpha,
    p_adjust_method = "BH",
    deseq = deseq_params,
    results = results_params
  )

  raw_output <- list(
    native_dds = dds,
    native_results = native_results,
    results_table = results_table,
    results_metadata = results_metadata,
    results_mcols = results_mcols,
    dds_mcols = dds_mcols,
    results_names = DESeq2::resultsNames(dds),
    requested_contrast = results_params$contrast,
    reference_level = contrast_row$group1,
    factor_levels = c(contrast_row$group1, contrast_row$group2),
    original_conditions = input$original_conditions,
    backend_conditions = input$backend_conditions,
    sample_order = input$sample_order,
    feature_order = input$feature_order,
    size_factors = DESeq2::sizeFactors(dds),
    normalization_factors = DESeq2::normalizationFactors(dds),
    dispersions = DESeq2::dispersions(dds),
    design = DESeq2::design(dds),
    backend_group_column = input$backend_group_column,
    backend_metadata = input$metadata,
    contrast_row = contrast_row,
    input_orientation = "feature_by_sample",
    transposed_from_context = TRUE,
    requested_fit_type = params$fit_type,
    actual_fit_type = actual_fit_type,
    deseq_params = deseq_params,
    results_params = results_params,
    cooks = cooks,
    cooks_cutoff = cooks_cutoff,
    replacement_diagnostics = replacement_diagnostics,
    feature_diagnostics = feature_diagnostics,
    warnings = warnings,
    messages = messages,
    params = backend_params,
    package_version = da_deseq2_package_version()
  )

  da_backend_result(
    method = "deseq2",
    results = results,
    raw_output = raw_output,
    notes = notes,
    params = backend_params
  )
}

da_validate_deseq2_contrast_row <- function(contrast_row) {
  required <- c("contrast", "group1", "group2", "contrast_type")
  if (!is.data.frame(contrast_row) || nrow(contrast_row) != 1 ||
      !all(required %in% names(contrast_row))) {
    stop(
      "`contrast_row` must contain one DESeq2 contrast.",
      call. = FALSE
    )
  }
  if (!as.character(contrast_row$contrast_type) %in%
      c("explicit", "pairwise")) {
    stop(
      "`contrast_row$contrast_type` must be \"explicit\" or \"pairwise\".",
      call. = FALSE
    )
  }

  values <- unlist(contrast_row[1, c("contrast", "group1", "group2")])
  if (any(is.na(values)) || any(!nzchar(values)) ||
      identical(as.character(contrast_row$group1),
                as.character(contrast_row$group2))) {
    stop("`contrast_row` contains invalid group labels.", call. = FALSE)
  }

  contrast_row[1, required, drop = FALSE]
}

da_validate_deseq2_params <- function(params) {
  if (!is.list(params)) {
    stop("DESeq2 `params` must be a list.", call. = FALSE)
  }
  expected <- c("sf_type", "fit_type", "independent_filtering", "alpha")
  unknown <- setdiff(names(params), expected)
  missing <- setdiff(expected, names(params))
  if (length(unknown) > 0) {
    stop(
      "Unknown DESeq2 parameter(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (length(missing) > 0) {
    stop(
      "Missing DESeq2 parameter(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  sf_types <- c("ratio", "poscounts", "iterate")
  if (!is.character(params$sf_type) || length(params$sf_type) != 1 ||
      is.na(params$sf_type) || !params$sf_type %in% sf_types) {
    stop(
      "`deseq2_sf_type` must be one of: ",
      paste(sf_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  fit_types <- c("parametric", "local", "mean", "glmGamPoi")
  if (!is.character(params$fit_type) || length(params$fit_type) != 1 ||
      is.na(params$fit_type) || !params$fit_type %in% fit_types) {
    stop(
      "`deseq2_fit_type` must be one of: ",
      paste(fit_types, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (identical(params$fit_type, "glmGamPoi") &&
      !da_optional_package_available("glmGamPoi")) {
    stop(
      "`deseq2_fit_type = \"glmGamPoi\"` requires the optional package ",
      "`glmGamPoi`.",
      call. = FALSE
    )
  }

  if (!is.logical(params$independent_filtering) ||
      length(params$independent_filtering) != 1 ||
      is.na(params$independent_filtering)) {
    stop(
      "`deseq2_independent_filtering` must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  if (!is.numeric(params$alpha) || length(params$alpha) != 1 ||
      is.na(params$alpha) || !is.finite(params$alpha) ||
      params$alpha <= 0 || params$alpha >= 1) {
    stop(
      "`deseq2_alpha` must be a single finite number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  list(
    sf_type = params$sf_type,
    fit_type = params$fit_type,
    independent_filtering = params$independent_filtering,
    alpha = as.numeric(params$alpha)
  )
}

da_prepare_deseq2_input <- function(context, contrast_row) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }
  contrast_row <- da_validate_deseq2_contrast_row(contrast_row)

  sample_ids <- context$sample_ids
  feature_ids <- context$feature_ids
  da_validate_deseq2_ids(sample_ids, "sample")
  da_validate_deseq2_ids(feature_ids, "feature")

  if (!identical(rownames(context$counts), sample_ids)) {
    stop(
      "DESeq2 sample IDs are not aligned with the context count table.",
      call. = FALSE
    )
  }
  if (is.null(context$metadata) ||
      !identical(rownames(context$metadata), sample_ids)) {
    stop(
      "DESeq2 sample metadata must be uniquely aligned to count-table ",
      "sample IDs.",
      call. = FALSE
    )
  }
  if (!identical(names(context$group_values), sample_ids)) {
    stop(
      "DESeq2 group labels are not aligned with count-table sample IDs.",
      call. = FALSE
    )
  }

  group1 <- as.character(contrast_row$group1)
  group2 <- as.character(contrast_row$group2)
  group_values <- as.character(context$group_values)
  keep <- group_values %in% c(group1, group2)
  selected_samples <- sample_ids[keep]
  conditions <- group_values[keep]
  if (!all(c(group1, group2) %in% conditions)) {
    stop(
      "Both contrast groups must have samples for DESeq2.",
      call. = FALSE
    )
  }

  selected_counts <- context$counts[selected_samples, , drop = FALSE]
  da_validate_deseq2_counts(selected_counts)

  group_column <- da_deseq2_internal_name(
    "microeda_deseq2_group",
    colnames(context$metadata)
  )
  backend_conditions <- factor(conditions, levels = c(group1, group2))
  backend_metadata <- data.frame(row.names = selected_samples)
  backend_metadata[[group_column]] <- backend_conditions
  design <- stats::reformulate(group_column)

  list(
    counts = t(selected_counts),
    metadata = backend_metadata,
    design = design,
    original_conditions = stats::setNames(conditions, selected_samples),
    backend_conditions = stats::setNames(
      backend_conditions,
      selected_samples
    ),
    sample_order = selected_samples,
    feature_order = feature_ids,
    backend_group_column = group_column
  )
}

da_validate_deseq2_ids <- function(ids, type) {
  if (!is.character(ids) || length(ids) < 1 ||
      any(is.na(ids)) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(
      "DESeq2 requires unique, non-missing ",
      type,
      " IDs.",
      call. = FALSE
    )
  }

  invisible(ids)
}

da_validate_deseq2_counts <- function(counts) {
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("DESeq2 requires a numeric count matrix.", call. = FALSE)
  }
  if (any(!is.finite(counts))) {
    stop("DESeq2 counts must contain only finite values.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("DESeq2 counts cannot contain negative values.", call. = FALSE)
  }
  integer_like <- abs(counts - round(counts)) < sqrt(.Machine$double.eps)
  if (!all(integer_like)) {
    stop(
      "DESeq2 requires integer-like counts; microeda does not round counts.",
      call. = FALSE
    )
  }

  library_sizes <- rowSums(counts)
  if (any(library_sizes <= 0)) {
    stop(
      "DESeq2 requires positive library sizes for every selected sample; ",
      "zero-library sample(s): ",
      paste(rownames(counts)[library_sizes <= 0], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(counts)
}

da_deseq2_internal_name <- function(base, used_names) {
  candidate <- base
  suffix <- 0L
  while (candidate %in% used_names) {
    suffix <- suffix + 1L
    candidate <- paste0(base, "_", suffix)
  }
  candidate
}

da_deseq2_fit_params <- function(params) {
  list(
    test = "Wald",
    fitType = params$fit_type,
    sfType = params$sf_type,
    betaPrior = FALSE,
    quiet = TRUE,
    minReplicatesForReplace = Inf,
    useT = FALSE,
    parallel = FALSE
  )
}

da_deseq2_results_params <- function(input, params) {
  factor_levels <- levels(input$metadata[[input$backend_group_column]])
  list(
    contrast = c(
      input$backend_group_column,
      factor_levels[[2]],
      factor_levels[[1]]
    ),
    lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    independentFiltering = params$independent_filtering,
    alpha = params$alpha,
    pAdjustMethod = "BH"
  )
}

da_standardize_deseq2_result <- function(results_table,
                                         context,
                                         contrast_row,
                                         feature_order) {
  if (!is.data.frame(results_table)) {
    stop("DESeq2 results must be coercible to a data frame.", call. = FALSE)
  }
  contrast_row <- da_validate_deseq2_contrast_row(contrast_row)
  required <- c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  missing <- setdiff(required, names(results_table))
  if (length(missing) > 0) {
    stop(
      "DESeq2 results are missing required column(s): ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  feature_ids <- rownames(results_table)
  da_validate_deseq2_ids(feature_ids, "result feature")
  if (!identical(feature_ids, feature_order)) {
    stop(
      "DESeq2 result feature IDs or order differ from the backend input.",
      call. = FALSE
    )
  }

  numeric_values <- lapply(required, function(column) {
    values <- results_table[[column]]
    if (!is.numeric(values) || length(values) != nrow(results_table)) {
      stop(
        "DESeq2 result column `",
        column,
        "` must be numeric and match result rows.",
        call. = FALSE
      )
    }
    as.numeric(values)
  })
  names(numeric_values) <- required
  p_adjusted <- numeric_values$padj
  method_note <- da_method_note("deseq2")$message

  da_standard_result(
    feature_id = feature_ids,
    taxon_label = da_taxon_labels(context, feature_ids),
    rank = da_result_rank(context),
    method = "deseq2",
    contrast = contrast_row$contrast,
    group1 = contrast_row$group1,
    group2 = contrast_row$group2,
    effect = numeric_values$log2FoldChange,
    effect_type = "deseq2_log2_fold_change_group2_vs_group1",
    log_fold_change = numeric_values$log2FoldChange,
    statistic = numeric_values$stat,
    standard_error = numeric_values$lfcSE,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = numeric_values$pvalue,
    p_adjusted = p_adjusted,
    p_adjust_method = "BH",
    p_adjust_scope = "method_contrast",
    significance = da_p_significance(p_adjusted),
    direction = NA_character_,
    method_note = method_note
  )
}

da_deseq2_actual_fit_type <- function(dds, requested_fit_type) {
  fit_type <- attr(DESeq2::dispersionFunction(dds), "fitType")
  if (is.null(fit_type) || length(fit_type) != 1 ||
      is.na(fit_type) || !nzchar(fit_type)) {
    return(requested_fit_type)
  }

  as.character(fit_type)
}

da_deseq2_cooks_cutoff <- function(dds) {
  model_matrix <- attr(dds, "dispModelMatrix")
  if (is.null(model_matrix)) {
    return(NA_real_)
  }
  m <- nrow(model_matrix)
  p <- ncol(model_matrix)
  if (m <= p) {
    return(NA_real_)
  }

  stats::qf(0.99, p, m - p)
}

da_deseq2_feature_diagnostics <- function(results_table,
                                          dds_mcols,
                                          cooks_cutoff,
                                          independent_filtering) {
  feature_ids <- rownames(results_table)
  n <- length(feature_ids)
  all_zero <- if ("allZero" %in% names(dds_mcols)) {
    as.logical(dds_mcols$allZero)
  } else {
    rep(FALSE, n)
  }
  all_zero[is.na(all_zero)] <- FALSE
  max_cooks <- if ("maxCooks" %in% names(dds_mcols)) {
    as.numeric(dds_mcols$maxCooks)
  } else {
    rep(NA_real_, n)
  }

  p_value_missing <- is.na(results_table$pvalue)
  p_adjusted_missing <- is.na(results_table$padj)
  cooks_outlier <- p_value_missing & !all_zero &
    is.finite(max_cooks) & is.finite(cooks_cutoff) &
    max_cooks > cooks_cutoff
  independently_filtered <- isTRUE(independent_filtering) &
    !p_value_missing & p_adjusted_missing
  other_na <- (p_value_missing | p_adjusted_missing) &
    !all_zero & !cooks_outlier & !independently_filtered

  data.frame(
    feature_id = feature_ids,
    all_zero = all_zero,
    cooks_outlier = cooks_outlier,
    independently_filtered = independently_filtered,
    other_na = other_na,
    p_value_missing = p_value_missing,
    p_adjusted_missing = p_adjusted_missing,
    max_cooks = max_cooks,
    stringsAsFactors = FALSE
  )
}

da_deseq2_replacement_diagnostics <- function(dds_mcols, assay_names) {
  replace_flags <- if ("replace" %in% names(dds_mcols)) {
    as.logical(dds_mcols$replace)
  } else {
    rep(FALSE, nrow(dds_mcols))
  }
  replace_flags[is.na(replace_flags)] <- FALSE

  list(
    min_replicates_for_replace = Inf,
    replace_counts_present = "replaceCounts" %in% assay_names,
    replace_cooks_present = "replaceCooks" %in% assay_names,
    replacement_flag_count = sum(replace_flags),
    replacement_flags = replace_flags
  )
}

da_deseq2_notes <- function(warnings,
                            contrast,
                            requested_fit_type,
                            actual_fit_type,
                            feature_diagnostics,
                            replacement_diagnostics) {
  rows <- list(da_method_note("deseq2"))

  warnings <- unique(warnings[nzchar(warnings)])
  if (length(warnings) > 0) {
    warning_rows <- lapply(seq_along(warnings), function(i) {
      da_caveat(
        method = "deseq2",
        caveat_id = paste0("deseq2_backend_warning_", i),
        topic = "backend",
        severity = "warning",
        message = paste0("DESeq2 reported: ", warnings[[i]])
      )
    })
    rows <- c(rows, warning_rows)
  }

  if (!identical(requested_fit_type, actual_fit_type)) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = "deseq2",
      caveat_id = "deseq2_dispersion_fit_fallback",
      topic = "model_fit",
      severity = "warning",
      message = paste0(
        "DESeq2 used dispersion fit `",
        actual_fit_type,
        "` instead of requested `",
        requested_fit_type,
        "` for contrast ",
        contrast,
        "."
      )
    )
  }

  diagnostic_specs <- list(
    all_zero = c(
      "deseq2_all_zero_features",
      "all-zero feature(s) were retained with native missing test results"
    ),
    cooks_outlier = c(
      "deseq2_cooks_outliers",
      "feature(s) had native p-values suppressed by Cook's diagnostics"
    ),
    independently_filtered = c(
      "deseq2_independently_filtered",
      paste(
        "feature(s) retained native p-values but have missing adjusted",
        "p-values after independent filtering"
      )
    ),
    other_na = c(
      "deseq2_other_na",
      "feature(s) had other native missing p-value or adjusted-p diagnostics"
    )
  )
  for (diagnostic in names(diagnostic_specs)) {
    count <- sum(feature_diagnostics[[diagnostic]], na.rm = TRUE)
    if (count == 0) {
      next
    }
    spec <- diagnostic_specs[[diagnostic]]
    rows[[length(rows) + 1L]] <- da_caveat(
      method = "deseq2",
      caveat_id = spec[[1]],
      topic = "backend_diagnostics",
      severity = if (identical(diagnostic, "independently_filtered")) {
        "info"
      } else {
        "warning"
      },
      message = paste0(
        "Contrast ",
        contrast,
        ": ",
        count,
        " ",
        spec[[2]],
        "."
      )
    )
  }

  replacement_present <- isTRUE(
    replacement_diagnostics$replace_counts_present
  ) || replacement_diagnostics$replacement_flag_count > 0
  if (replacement_present) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = "deseq2",
      caveat_id = "deseq2_unexpected_count_replacement",
      topic = "backend_diagnostics",
      severity = "warning",
      message = paste0(
        "DESeq2 unexpectedly recorded count replacement for contrast ",
        contrast,
        " despite `minReplicatesForReplace = Inf`."
      )
    )
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_deseq2_pairwise_notes <- function(notes, contrast) {
  if (is.null(notes) || nrow(notes) == 0) {
    return(da_empty_caveats())
  }

  method_note <- notes$caveat_id == "deseq2_sensitivity_note"
  contrast_specific <- !method_note
  notes$caveat_id[contrast_specific] <- paste0(
    notes$caveat_id[contrast_specific],
    "_",
    contrast
  )
  notes
}

da_deseq2_mcols <- function(x) {
  getExportedValue("S4Vectors", "mcols")(x)
}

da_deseq2_metadata <- function(x) {
  getExportedValue("S4Vectors", "metadata")(x)
}

da_deseq2_assay_names <- function(x) {
  getExportedValue("SummarizedExperiment", "assayNames")(x)
}

da_deseq2_assay <- function(x, name) {
  getExportedValue("SummarizedExperiment", "assay")(x, name)
}

da_deseq2_package_version <- function() {
  as.character(utils::packageVersion("DESeq2"))
}
