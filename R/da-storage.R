da_validate_progress <- function(progress) {
  if (!is.logical(progress) || length(progress) != 1 || is.na(progress)) {
    stop("`progress` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  progress
}

da_empty_timings <- function() {
  data.frame(
    method = character(),
    contrast = character(),
    elapsed_seconds = numeric(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

da_execution_state <- function(progress, total_methods) {
  state <- new.env(parent = emptyenv())
  state$progress <- da_validate_progress(progress)
  state$total_methods <- as.integer(total_methods)
  state$method_index <- NA_integer_
  state$timings <- da_empty_timings()
  state
}

da_execution_context <- function(context, state, method_index) {
  state$method_index <- as.integer(method_index)
  context$execution <- state
  context
}

da_progress_method <- function(context, method) {
  state <- context$execution
  if (!is.environment(state) || !isTRUE(state$progress)) {
    return(invisible(NULL))
  }

  message(
    "[microeda DA] Method ",
    state$method_index,
    "/",
    state$total_methods,
    ": ",
    method
  )
  invisible(NULL)
}

da_execute_contrast <- function(context, method, contrast_row, runner) {
  if (!is.function(runner)) {
    stop("`runner` must be a function.", call. = FALSE)
  }

  contrast <- as.character(contrast_row$contrast[[1]])
  contrast_index <- match(contrast, context$contrast_plan$contrast)
  contrast_total <- nrow(context$contrast_plan)
  state <- context$execution
  show_progress <- is.environment(state) && isTRUE(state$progress)

  if (show_progress) {
    message(
      "[microeda DA]   Contrast ",
      contrast_index,
      "/",
      contrast_total,
      ": ",
      contrast,
      " - started"
    )
  }

  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    runner(),
    error = function(e) {
      elapsed <- max(0, proc.time()[["elapsed"]] - started)
      if (show_progress) {
        message(
          "[microeda DA]   Contrast ",
          contrast_index,
          "/",
          contrast_total,
          ": ",
          contrast,
          " - failed after ",
          format(round(elapsed, 1), nsmall = 1, trim = TRUE),
          " s"
        )
      }
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  elapsed <- max(0, proc.time()[["elapsed"]] - started)

  if (is.environment(state)) {
    state$timings <- rbind(
      state$timings,
      data.frame(
        method = method,
        contrast = contrast,
        elapsed_seconds = as.numeric(elapsed),
        status = "success",
        stringsAsFactors = FALSE
      )
    )
    row.names(state$timings) <- NULL
  }

  if (show_progress) {
    message(
      "[microeda DA]   Contrast ",
      contrast_index,
      "/",
      contrast_total,
      ": ",
      contrast,
      " - finished in ",
      format(round(elapsed, 1), nsmall = 1, trim = TRUE),
      " s"
    )
  }

  result
}

da_apply_raw_storage <- function(backend_result,
                                 raw_storage,
                                 contrast_plan,
                                 timings = da_empty_timings()) {
  da_validate_backend_result(backend_result)
  raw_storage <- match.arg(raw_storage, c("full", "compact", "none"))

  if (identical(raw_storage, "full")) {
    return(backend_result)
  }
  if (identical(raw_storage, "none")) {
    backend_result["raw_output"] <- list(NULL)
    return(backend_result)
  }

  backend_result$raw_output <- da_compact_backend_raw(
    method = backend_result$method,
    raw = backend_result$raw_output,
    contrast_plan = contrast_plan,
    timings = timings
  )
  backend_result
}

da_compact_backend_raw <- function(method, raw, contrast_plan, timings) {
  is_pairwise <- nrow(contrast_plan) > 1 ||
    all(contrast_plan$contrast_type == "pairwise")

  if (!is_pairwise) {
    return(da_compact_contrast_raw(
      method = method,
      raw = raw,
      timing = timings[
        timings$contrast == contrast_plan$contrast[[1]],
        ,
        drop = FALSE
      ]
    ))
  }

  contrast_labels <- contrast_plan$contrast
  if (identical(method, "aldex2")) {
    compact_contrasts <- lapply(contrast_labels, function(contrast) {
      da_compact_contrast_raw(
        method = method,
        raw = raw$contrasts[[contrast]],
        timing = timings[timings$contrast == contrast, , drop = FALSE]
      )
    })
    names(compact_contrasts) <- contrast_labels
    return(list(
      raw_storage = "compact",
      contrasts = compact_contrasts,
      contrast_plan = contrast_plan,
      params = raw$params,
      input_orientation = raw$input_orientation,
      transposed_from_context = raw$transposed_from_context,
      timings = timings
    ))
  }

  compact <- lapply(contrast_labels, function(contrast) {
    da_compact_contrast_raw(
      method = method,
      raw = raw[[contrast]],
      timing = timings[timings$contrast == contrast, , drop = FALSE]
    )
  })
  names(compact) <- contrast_labels
  compact
}

da_compact_contrast_raw <- function(method, raw, timing = da_empty_timings()) {
  if (!is.list(raw)) {
    stop(
      "Cannot create compact raw output for method `",
      method,
      "` because its raw output is not a list.",
      call. = FALSE
    )
  }

  keep <- switch(
    method,
    aldex2 = c(
      "ttest", "effect", "combined", "conditions", "backend_conditions",
      "condition_mapping", "sample_order", "feature_order", "pair_id",
      "pairing", "contrast_row", "input_orientation",
      "transposed_from_context", "params", "warnings", "messages",
      "package_version"
    ),
    ancombc2 = c(
      "res", "zero_ind", "ss_tab", "coefficient", "reference_level",
      "factor_levels", "original_conditions", "backend_conditions",
      "sample_order", "feature_order", "backend_excluded_features",
      "backend_excluded_feature_count", "contrast_row", "input_orientation",
      "transposed_from_context", "backend_group_column",
      "inert_metadata_column", "params", "warnings", "messages",
      "package_version"
    ),
    deseq2 = c(
      "results_table", "results_metadata", "results_mcols", "results_names",
      "requested_contrast", "reference_level", "factor_levels",
      "original_conditions", "backend_conditions", "sample_order",
      "feature_order", "size_factors", "dispersions", "design",
      "backend_group_column", "contrast_row", "input_orientation",
      "transposed_from_context", "requested_fit_type", "actual_fit_type",
      "deseq_params", "results_params", "cooks_cutoff",
      "replacement_diagnostics", "feature_diagnostics", "warnings",
      "messages", "params", "package_version"
    ),
    character()
  )
  keep <- intersect(keep, names(raw))
  out <- raw[keep]

  if (identical(method, "aldex2")) {
    if (!"feature_order" %in% names(out) && is.data.frame(out$combined)) {
      out$feature_order <- rownames(out$combined)
    }
    if (!"package_version" %in% names(out) &&
        da_optional_package_available("ALDEx2")) {
      out$package_version <- as.character(utils::packageVersion("ALDEx2"))
    }
  }

  out <- c(
    list(raw_storage = "compact"),
    out,
    list(timing = timing)
  )
  out
}

da_raw_storage_mode <- function(x) {
  mode <- x$params$raw_storage
  if (is.character(mode) && length(mode) == 1 && !is.na(mode) &&
      mode %in% c("full", "compact", "none")) {
    return(mode)
  }

  "full"
}
