da_validate_ancombc2_public_arguments <- function(mc_samples_supplied,
                                                   denom_supplied,
                                                   paired.test,
                                                   pair_id) {
  if (!is.logical(paired.test) || length(paired.test) != 1 ||
      is.na(paired.test)) {
    stop("`paired.test` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(paired.test)) {
    stop(
      "`paired.test = TRUE` is only supported by the ALDEx2 backend.",
      call. = FALSE
    )
  }
  if (!is.null(pair_id)) {
    stop(
      "`pair_id` is only supported by the paired ALDEx2 backend.",
      call. = FALSE
    )
  }
  if (isTRUE(mc_samples_supplied)) {
    stop(
      "`mc.samples` is an ALDEx2-specific argument and cannot be supplied ",
      "with `methods = \"ancombc2\"`.",
      call. = FALSE
    )
  }
  if (isTRUE(denom_supplied)) {
    stop(
      "`denom` is an ALDEx2-specific argument and cannot be supplied with ",
      "`methods = \"ancombc2\"`.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

da_run_ancombc2 <- function(context, contrast_row, params) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }
  if (!"ancombc2" %in% context$methods) {
    stop("`context$methods` must include \"ancombc2\".", call. = FALSE)
  }
  contrast_row <- da_validate_ancombc2_contrast_row(contrast_row)
  if (!da_optional_package_available("ANCOMBC")) {
    stop(
      "`da_run_ancombc2()` requires the optional package `ANCOMBC`. ",
      "Install it with `BiocManager::install(\"ANCOMBC\")`.",
      call. = FALSE
    )
  }

  params <- da_validate_ancombc2_params(params)
  input <- da_prepare_ancombc2_input(context, contrast_row)
  actual_params <- da_ancombc2_call_params(
    group_column = input$backend_group_column,
    p_adj_method = params$p_adj_method
  )

  warnings <- character()
  messages <- character()
  native_result <- NULL
  backend_error <- NULL
  stdout <- utils::capture.output({
    tryCatch(
      withCallingHandlers(
        {
          native_result <- da_call_ancombc2(
            data = input$counts,
            meta_data = input$metadata,
            params = actual_params
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
  messages <- unique(c(messages, stdout[nzchar(stdout)]))
  warnings <- unique(warnings)

  if (!is.null(backend_error)) {
    stop(
      "ANCOM-BC2 failed for contrast `",
      contrast_row$contrast,
      "`: ",
      backend_error,
      call. = FALSE
    )
  }
  if (!is.list(native_result) || is.null(native_result$res)) {
    stop(
      "ANCOM-BC2 returned an invalid result for contrast `",
      contrast_row$contrast,
      "`: `res` is missing.",
      call. = FALSE
    )
  }

  results <- da_standardize_ancombc2_result(
    res = native_result$res,
    context = context,
    contrast_row = contrast_row,
    coefficient = input$coefficient,
    p_adj_method = params$p_adj_method
  )
  exclusions <- da_ancombc2_feature_exclusions(
    input_features = input$feature_order,
    output_features = results$feature_id
  )
  notes <- da_ancombc2_notes(
    warnings = warnings,
    excluded_features = exclusions$excluded_features,
    contrast = contrast_row$contrast
  )

  raw_output <- list(
    native_result = native_result,
    res = native_result$res,
    res_global = native_result[["res_global"]],
    res_pair = native_result[["res_pair"]],
    res_dunn = native_result[["res_dunn"]],
    res_trend = native_result[["res_trend"]],
    zero_ind = native_result[["zero_ind"]],
    ss_tab = native_result[["ss_tab"]],
    coefficient = input$coefficient,
    reference_level = contrast_row$group1,
    factor_levels = c(contrast_row$group1, contrast_row$group2),
    original_conditions = input$original_conditions,
    backend_conditions = input$backend_conditions,
    sample_order = input$sample_order,
    feature_order = input$feature_order,
    backend_excluded_features = exclusions$excluded_features,
    backend_excluded_feature_count = length(exclusions$excluded_features),
    contrast_row = contrast_row,
    input_orientation = "feature_by_sample",
    transposed_from_context = TRUE,
    backend_group_column = input$backend_group_column,
    inert_metadata_column = input$inert_metadata_column,
    backend_metadata = input$metadata,
    params = actual_params,
    warnings = warnings,
    messages = messages,
    package_version = da_ancombc2_package_version()
  )

  da_backend_result(
    method = "ancombc2",
    results = results,
    raw_output = raw_output,
    notes = notes,
    params = actual_params
  )
}

da_validate_ancombc2_contrast_row <- function(contrast_row) {
  required <- c("contrast", "group1", "group2", "contrast_type")
  if (!is.data.frame(contrast_row) || nrow(contrast_row) != 1 ||
      !all(required %in% names(contrast_row))) {
    stop(
      "`contrast_row` must contain one explicit ANCOM-BC2 contrast.",
      call. = FALSE
    )
  }
  if (!identical(as.character(contrast_row$contrast_type), "explicit")) {
    stop(
      "ANCOM-BC2 currently supports exactly one explicit contrast; ",
      "pairwise execution is not implemented.",
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

da_validate_ancombc2_params <- function(params) {
  if (!is.list(params)) {
    stop("ANCOM-BC2 `params` must be a list.", call. = FALSE)
  }
  unknown <- setdiff(names(params), "p_adj_method")
  if (length(unknown) > 0) {
    stop(
      "Unknown ANCOM-BC2 parameter(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  p_adj_method <- params$p_adj_method
  if (!is.character(p_adj_method) || length(p_adj_method) != 1 ||
      is.na(p_adj_method) ||
      !p_adj_method %in% stats::p.adjust.methods) {
    stop(
      "`ancombc2_p_adj_method` must be one of stats::p.adjust.methods.",
      call. = FALSE
    )
  }

  list(p_adj_method = p_adj_method)
}

da_prepare_ancombc2_input <- function(context, contrast_row) {
  if (!inherits(context, "microeda_da_context")) {
    stop("`context` must be a microeda_da_context object.", call. = FALSE)
  }
  contrast_row <- da_validate_ancombc2_contrast_row(contrast_row)

  sample_ids <- context$sample_ids
  feature_ids <- context$feature_ids
  da_validate_ancombc2_ids(sample_ids, "sample")
  da_validate_ancombc2_ids(feature_ids, "feature")

  if (!identical(rownames(context$counts), sample_ids)) {
    stop(
      "ANCOM-BC2 sample IDs are not aligned with the context count table.",
      call. = FALSE
    )
  }
  if (is.null(context$metadata) ||
      !identical(rownames(context$metadata), sample_ids)) {
    stop(
      "ANCOM-BC2 sample metadata must be uniquely aligned to count-table ",
      "sample IDs.",
      call. = FALSE
    )
  }
  if (!identical(names(context$group_values), sample_ids)) {
    stop(
      "ANCOM-BC2 group labels are not aligned with count-table sample IDs.",
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
      "Both explicit contrast groups must have samples for ANCOM-BC2.",
      call. = FALSE
    )
  }

  selected_counts <- context$counts[selected_samples, , drop = FALSE]
  da_validate_ancombc2_counts(selected_counts)

  used_names <- colnames(context$metadata)
  group_column <- da_ancombc2_internal_name(
    "microeda_ancombc2_group",
    used_names
  )
  inert_column <- da_ancombc2_internal_name(
    "microeda_ancombc2_metadata",
    c(used_names, group_column)
  )
  backend_conditions <- factor(conditions, levels = c(group1, group2))
  backend_metadata <- data.frame(row.names = selected_samples)
  backend_metadata[[group_column]] <- backend_conditions
  backend_metadata[[inert_column]] <- rep("microeda_internal", length(conditions))

  model_columns <- colnames(stats::model.matrix(
    stats::reformulate(group_column),
    data = backend_metadata
  ))
  coefficient <- setdiff(model_columns, "(Intercept)")
  if (length(coefficient) != 1) {
    stop(
      "ANCOM-BC2 expected exactly one non-intercept group coefficient.",
      call. = FALSE
    )
  }

  list(
    counts = t(selected_counts),
    metadata = backend_metadata,
    original_conditions = stats::setNames(conditions, selected_samples),
    backend_conditions = stats::setNames(
      backend_conditions,
      selected_samples
    ),
    sample_order = selected_samples,
    feature_order = feature_ids,
    coefficient = unname(coefficient),
    backend_group_column = group_column,
    inert_metadata_column = inert_column
  )
}

da_validate_ancombc2_ids <- function(ids, type) {
  if (!is.character(ids) || length(ids) < 1 ||
      any(is.na(ids)) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop(
      "ANCOM-BC2 requires unique, non-missing ",
      type,
      " IDs.",
      call. = FALSE
    )
  }

  invisible(ids)
}

da_validate_ancombc2_counts <- function(counts) {
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("ANCOM-BC2 requires a numeric count matrix.", call. = FALSE)
  }
  if (any(!is.finite(counts))) {
    stop("ANCOM-BC2 counts must contain only finite values.", call. = FALSE)
  }
  if (any(counts < 0)) {
    stop("ANCOM-BC2 counts cannot contain negative values.", call. = FALSE)
  }
  integer_like <- abs(counts - round(counts)) < sqrt(.Machine$double.eps)
  if (!all(integer_like)) {
    stop(
      "ANCOM-BC2 requires integer-like counts; microeda does not round counts.",
      call. = FALSE
    )
  }

  library_sizes <- rowSums(counts)
  if (any(library_sizes <= 0)) {
    stop(
      "ANCOM-BC2 requires positive library sizes for every selected sample; ",
      "zero-library sample(s): ",
      paste(rownames(counts)[library_sizes <= 0], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  feature_totals <- colSums(counts)
  if (any(feature_totals == 0)) {
    zero_features <- colnames(counts)[feature_totals == 0]
    stop(
      "ANCOM-BC2 input contains ",
      length(zero_features),
      " all-zero feature(s) in the selected contrast: ",
      paste(zero_features, collapse = ", "),
      ". microeda does not remove features silently.",
      call. = FALSE
    )
  }

  invisible(counts)
}

da_ancombc2_internal_name <- function(base, used_names) {
  candidate <- base
  suffix <- 0L
  while (candidate %in% used_names) {
    suffix <- suffix + 1L
    candidate <- paste0(base, "_", suffix)
  }
  candidate
}

da_ancombc2_call_params <- function(group_column, p_adj_method) {
  list(
    taxa_are_rows = TRUE,
    rank = NULL,
    tax_level = NULL,
    aggregate_data = NULL,
    fix_formula = group_column,
    rand_formula = NULL,
    p_adj_method = p_adj_method,
    pseudo = 0,
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    s0_perc = 0.05,
    group = NULL,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE,
    iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
    em_control = list(tol = 1e-05, max_iter = 100)
  )
}

da_call_ancombc2 <- function(data, meta_data, params) {
  do.call(
    ANCOMBC::ancombc2,
    c(
      list(
        data = data,
        meta_data = meta_data
      ),
      params
    )
  )
}

da_ancombc2_package_version <- function() {
  as.character(utils::packageVersion("ANCOMBC"))
}

da_standardize_ancombc2_result <- function(res,
                                           context,
                                           contrast_row,
                                           coefficient,
                                           p_adj_method) {
  if (!is.data.frame(res)) {
    stop("ANCOM-BC2 `res` must be a data frame.", call. = FALSE)
  }
  contrast_row <- da_validate_ancombc2_contrast_row(contrast_row)
  if (!is.character(coefficient) || length(coefficient) != 1 ||
      is.na(coefficient) || !nzchar(coefficient)) {
    stop("ANCOM-BC2 coefficient name is invalid.", call. = FALSE)
  }

  native_columns <- c(
    lfc = paste0("lfc_", coefficient),
    se = paste0("se_", coefficient),
    statistic = paste0("W_", coefficient),
    p = paste0("p_", coefficient),
    q = paste0("q_", coefficient)
  )
  required <- c("taxon", unname(native_columns))
  missing_columns <- setdiff(required, names(res))
  if (length(missing_columns) > 0) {
    stop(
      "ANCOM-BC2 `res` is missing required coefficient column(s): ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  feature_ids <- as.character(res$taxon)
  if (any(is.na(feature_ids)) || any(!nzchar(feature_ids)) ||
      anyDuplicated(feature_ids)) {
    stop(
      "ANCOM-BC2 `res$taxon` must contain unique, non-missing feature IDs.",
      call. = FALSE
    )
  }
  unknown_features <- setdiff(feature_ids, context$feature_ids)
  if (length(unknown_features) > 0) {
    stop(
      "ANCOM-BC2 returned feature ID(s) absent from its input: ",
      paste(unknown_features, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  numeric_values <- lapply(native_columns, function(column) {
    values <- res[[column]]
    if (!is.numeric(values) || length(values) != nrow(res)) {
      stop(
        "ANCOM-BC2 coefficient column `",
        column,
        "` must be numeric and match `res` rows.",
        call. = FALSE
      )
    }
    as.numeric(values)
  })
  lfc <- numeric_values$lfc
  p_adjusted <- numeric_values$q
  method_note <- da_method_note("ancombc2")$message

  da_standard_result(
    feature_id = feature_ids,
    taxon_label = da_taxon_labels(context, feature_ids),
    rank = da_result_rank(context),
    method = "ancombc2",
    contrast = contrast_row$contrast,
    group1 = contrast_row$group1,
    group2 = contrast_row$group2,
    effect = lfc,
    effect_type = "ancombc2_log_fold_change_group2_vs_group1",
    log_fold_change = lfc,
    statistic = numeric_values$statistic,
    standard_error = numeric_values$se,
    ci_low = NA_real_,
    ci_high = NA_real_,
    p_value = numeric_values$p,
    p_adjusted = p_adjusted,
    p_adjust_method = p_adj_method,
    p_adjust_scope = "method_contrast",
    significance = da_p_significance(p_adjusted),
    direction = NA_character_,
    method_note = method_note
  )
}

da_ancombc2_feature_exclusions <- function(input_features, output_features) {
  unknown <- setdiff(output_features, input_features)
  if (length(unknown) > 0) {
    stop(
      "ANCOM-BC2 returned feature ID(s) absent from its input: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(excluded_features = setdiff(input_features, output_features))
}

da_ancombc2_notes <- function(warnings,
                              excluded_features,
                              contrast) {
  rows <- list(da_method_note("ancombc2"))

  warnings <- unique(warnings[nzchar(warnings)])
  if (length(warnings) > 0) {
    warning_rows <- lapply(seq_along(warnings), function(i) {
      da_caveat(
        method = "ancombc2",
        caveat_id = paste0("ancombc2_backend_warning_", i),
        topic = "backend",
        severity = "warning",
        message = paste0("ANCOM-BC2 reported: ", warnings[[i]])
      )
    })
    rows <- c(rows, warning_rows)
  }

  if (length(excluded_features) > 0) {
    rows[[length(rows) + 1L]] <- da_caveat(
      method = "ancombc2",
      caveat_id = "ancombc2_backend_feature_exclusion",
      topic = "backend",
      severity = "warning",
      message = paste0(
        "ANCOM-BC2 returned ",
        length(excluded_features),
        " fewer feature(s) than supplied for contrast ",
        contrast,
        "; standardized results include only native `res` rows."
      )
    )
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}
