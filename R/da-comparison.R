#' Extract one method- and contrast-specific DA raw output
#'
#' `as_da_raw_output()` provides a stable access path across the different raw
#' output nesting used by ALDEx2, ANCOM-BC2, and DESeq2. The returned object is
#' not transformed; its structure remains specific to the selected backend.
#'
#' @param x A `microeda_da` object.
#' @param method Optional single method ID. It may be omitted only when `x`
#'   contains one method.
#' @param contrast Optional single contrast label. It may be omitted only when
#'   `x` contains one contrast.
#'
#' @return The selected method-specific raw backend object.
#'
#' @seealso [as_da_results()], [as_da_comparison()]
#'
#' @examples
#' \dontrun{
#' raw <- as_da_raw_output(da, method = "deseq2", contrast = "A_vs_B")
#' }
#'
#' @export
as_da_raw_output <- function(x, method = NULL, contrast = NULL) {
  da_comparison_validate_x(x)
  method <- da_raw_output_method(x, method)
  contrast <- da_raw_output_contrast(x, contrast)

  if (is.null(x$raw_outputs) || !is.list(x$raw_outputs) ||
      !method %in% names(x$raw_outputs)) {
    stop(
      "Raw output for method `",
      method,
      "` is missing from `x$raw_outputs`.",
      call. = FALSE
    )
  }

  raw <- x$raw_outputs[[method]]
  contrast_row <- x$contrast_plan[
    x$contrast_plan$contrast == contrast,
    ,
    drop = FALSE
  ]
  if (nrow(contrast_row) != 1) {
    stop(
      "The contrast plan does not contain exactly one row for `",
      contrast,
      "`.",
      call. = FALSE
    )
  }

  if (identical(as.character(contrast_row$contrast_type), "explicit")) {
    return(raw)
  }

  if (identical(method, "aldex2")) {
    if (!is.list(raw) || !is.list(raw$contrasts) ||
        !contrast %in% names(raw$contrasts)) {
      stop(
        "Malformed ALDEx2 pairwise raw output: contrast `",
        contrast,
        "` is missing from `x$raw_outputs$aldex2$contrasts`.",
        call. = FALSE
      )
    }
    return(raw$contrasts[[contrast]])
  }

  if (!is.list(raw) || !contrast %in% names(raw)) {
    stop(
      "Malformed ",
      method,
      " pairwise raw output: contrast `",
      contrast,
      "` is missing.",
      call. = FALSE
    )
  }

  raw[[contrast]]
}

#' Align standardized DA results side by side
#'
#' `as_da_comparison()` creates one row per feature and contrast and aligns
#' selected method-specific standardized outputs in prefixed column blocks.
#' It is a descriptive view, not a consensus analysis or method ranking.
#'
#' ALDEx2 effects retain the native ALDEx2 scale. ANCOM-BC2 effects are native
#' natural-log estimates for `group2 - group1`. DESeq2 effects are native
#' unshrunk log2 fold changes for `group2 / group1`. Signs share the
#' group2-relative orientation, but effect magnitudes are not directly
#' comparable across methods.
#'
#' @details
#' Identifier columns are `feature_id`, `taxon_label`, `rank`, `contrast`,
#' `group1`, and `group2`. Each selected method adds prefixed `effect`,
#' `effect_type`, `p_value`, `p_adjusted`, `significance`, `effect_sign`,
#' `tested`, and `p_adjusted_le_alpha` columns.
#'
#' `tested` is `TRUE` when the method returned the feature with a non-missing
#' native p-value, `FALSE` when the row exists but its native p-value is
#' missing, and `NA` when the backend did not return that feature for the
#' contrast. Missing values are not converted to non-significant results.
#'
#' @param x A `microeda_da` object.
#' @param contrast Optional contrast labels. `NULL` selects all contrasts.
#' @param methods Optional method IDs. `NULL` selects all methods.
#' @param alpha Strictly positive adjusted p-value threshold below one. It is
#'   used only for each method's descriptive `p_adjusted_le_alpha` field.
#' @param features Optional feature IDs to retain.
#'
#' @return A base data frame with identifier columns followed by one prefixed
#'   column block per selected method. Native p-values, adjusted p-values,
#'   significance labels, and effects are not recalculated.
#'
#' @seealso [as_da_results()], [as_da_summary()], [as_da_raw_output()],
#'   [microeda_da_comparison_report()], [write_da_comparison()]
#'
#' @examples
#' \dontrun{
#' comparison <- as_da_comparison(da, contrast = "A_vs_B")
#' }
#'
#' @export
as_da_comparison <- function(x,
                             contrast = NULL,
                             methods = NULL,
                             alpha = 0.05,
                             features = NULL) {
  da_comparison_validate_x(x)
  selected_contrasts <- da_comparison_select_contrasts(x, contrast)
  selected_methods <- da_comparison_select_methods(x, methods)
  alpha <- da_comparison_validate_alpha(alpha)
  selected_features <- da_comparison_select_features(x, features)

  results <- as_da_results(x)
  results <- results[
    results$method %in% selected_methods &
      results$contrast %in% selected_contrasts,
    ,
    drop = FALSE
  ]
  feature_order <- da_comparison_feature_order(x)
  rows <- list()

  for (contrast_label in selected_contrasts) {
    contrast_results <- results[
      results$contrast == contrast_label,
      ,
      drop = FALSE
    ]
    universe <- unique(contrast_results$feature_id)
    universe <- feature_order[feature_order %in% universe]
    if (!is.null(selected_features)) {
      universe <- universe[universe %in% selected_features]
    }

    contrast_row <- x$contrast_plan[
      x$contrast_plan$contrast == contrast_label,
      ,
      drop = FALSE
    ]
    for (feature_id in universe) {
      rows[[length(rows) + 1L]] <- da_comparison_row(
        results = contrast_results,
        methods = selected_methods,
        contrast_row = contrast_row,
        feature_id = feature_id,
        alpha = alpha
      )
    }
  }

  if (length(rows) == 0) {
    out <- da_empty_comparison(selected_methods)
  } else {
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
  }

  attr(out, "alpha") <- alpha
  attr(out, "methods") <- selected_methods
  attr(out, "contrasts") <- selected_contrasts
  out
}

#' Create a readable side-by-side DA comparison report
#'
#' `microeda_da_comparison_report()` reports method-specific results next to
#' one another without combining p-values, harmonizing effects, ranking
#' methods, or creating consensus calls.
#'
#' @inheritParams as_da_comparison
#' @param max_features Maximum number of displayed feature rows per contrast.
#'   Use `Inf` to show all selected rows.
#'
#' @return A single character string suitable for `cat()`.
#'
#' @details
#' Unless `features` is supplied, feature rows are displayed when at least one
#' selected method has a finite native adjusted p-value no greater than
#' `alpha`. This controls display only and is not a consensus rule. Truncation
#' by `max_features` preserves the original feature order.
#'
#' @seealso [as_da_comparison()], [as_da_raw_output()],
#'   [microeda_da_report()]
#'
#' @examples
#' \dontrun{
#' cat(microeda_da_comparison_report(
#'   da,
#'   contrast = "A_vs_B",
#'   max_features = 15
#' ))
#' }
#'
#' @export
microeda_da_comparison_report <- function(x,
                                          contrast = NULL,
                                          methods = NULL,
                                          alpha = 0.05,
                                          features = NULL,
                                          max_features = 20) {
  comparison <- as_da_comparison(
    x = x,
    contrast = contrast,
    methods = methods,
    alpha = alpha,
    features = features
  )
  selected_methods <- attr(comparison, "methods")
  selected_contrasts <- attr(comparison, "contrasts")
  max_features <- da_comparison_validate_max_features(max_features)
  results <- as_da_results(x)

  lines <- c(
    "microeda differential representation comparison report",
    "",
    paste0("Methods: ", paste(selected_methods, collapse = ", ")),
    paste0("Contrasts: ", paste(selected_contrasts, collapse = ", ")),
    paste0("Display alpha: ", da_report_number(alpha, digits = 3)),
    "",
    paste(
      "Effect scales differ: ALDEx2 uses its native effect, ANCOM-BC2 uses",
      "a natural-log group2-minus-group1 estimate, and DESeq2 uses an",
      "unshrunk log2 group2/group1 fold change."
    ),
    paste(
      "Effect magnitudes must not be compared directly across methods;",
      "matching signs remain model-specific descriptive sensitivity patterns."
    )
  )

  for (contrast_label in selected_contrasts) {
    contrast_table <- comparison[
      comparison$contrast == contrast_label,
      ,
      drop = FALSE
    ]
    contrast_features <- contrast_table$feature_id
    lines <- c(
      lines,
      "",
      "============================================================",
      paste0("Contrast: ", contrast_label),
      "Method summaries:"
    )

    for (method in selected_methods) {
      lines <- c(
        lines,
        da_comparison_method_report_lines(
          x = x,
          results = results,
          method = method,
          contrast = contrast_label,
          features = contrast_features,
          alpha = alpha
        )
      )
    }

    if (is.null(features)) {
      display <- da_comparison_display_rows(
        contrast_table,
        methods = selected_methods,
        alpha = alpha
      )
      lines <- c(
        lines,
        paste0(
          "Display criterion: at least one selected method has a finite ",
          "native adjusted p-value <= ",
          da_report_number(alpha, digits = 3),
          ". This is a display filter, not a consensus call."
        )
      )
    } else {
      display <- contrast_table
      lines <- c(
        lines,
        "Display criterion: explicitly requested feature IDs."
      )
    }

    if (nrow(display) == 0) {
      lines <- c(
        lines,
        "No features meet the display criterion for this contrast."
      )
      next
    }

    shown_n <- if (is.infinite(max_features)) {
      nrow(display)
    } else {
      min(nrow(display), max_features)
    }
    hidden_n <- nrow(display) - shown_n
    display <- display[seq_len(shown_n), , drop = FALSE]
    lines <- c(
      lines,
      "Side-by-side standardized values:",
      da_report_table_lines(
        da_comparison_report_table(display, selected_methods),
        digits = 3
      )
    )
    if (hidden_n > 0) {
      lines <- c(
        lines,
        paste0(
          hidden_n,
          " feature row(s) hidden by `max_features`; original feature ",
          "order was preserved."
        )
      )
    }
  }

  selected_caveats <- x$caveats[
    is.na(x$caveats$method) | !nzchar(x$caveats$method) |
      x$caveats$method %in% selected_methods,
    ,
    drop = FALSE
  ]
  lines <- c(
    lines,
    "",
    "Object caveats:",
    da_report_caveat_lines(selected_caveats),
    "",
    "Interpretation caveats:",
    "- Native p-value adjustment and filtering differ by method.",
    "- Missing values are not evidence of non-significance.",
    "- Effect magnitudes are not directly comparable across methods.",
    paste(
      "- Agreement among methods is a descriptive sensitivity pattern,",
      "not proof of a biological discovery."
    ),
    paste(
      "- DESeq2 is a count-model sensitivity/comparison method;",
      "size-factor normalization is not a compositional correction."
    )
  )

  paste(lines, collapse = "\n")
}

#' Write a side-by-side DA comparison table to CSV
#'
#' `write_da_comparison()` writes only the table returned by
#' [as_da_comparison()]. It does not serialize `raw_outputs`,
#' `method_results`, or native backend objects.
#'
#' @inheritParams as_da_comparison
#' @param file Single non-empty character path for the output CSV.
#' @param ... Named CSV formatting arguments. Supported values are `quote`,
#'   `na`, `eol`, `qmethod`, and `fileEncoding`. The exported table, file path,
#'   and `row.names = FALSE` cannot be overridden.
#'
#' @return The output file path, invisibly.
#'
#' @seealso [write_da_results()], [as_da_comparison()]
#'
#' @examples
#' \dontrun{
#' write_da_comparison(da, "da_comparison.csv", contrast = "A_vs_B")
#' }
#'
#' @export
write_da_comparison <- function(x,
                                file,
                                contrast = NULL,
                                methods = NULL,
                                alpha = 0.05,
                                features = NULL,
                                ...) {
  da_comparison_validate_file(file)
  dots <- da_comparison_validate_csv_arguments(list(...))
  comparison <- as_da_comparison(
    x = x,
    contrast = contrast,
    methods = methods,
    alpha = alpha,
    features = features
  )

  do.call(
    utils::write.csv,
    c(
      list(
        x = comparison,
        file = file,
        row.names = FALSE
      ),
      dots
    )
  )

  invisible(file)
}

da_comparison_validate_x <- function(x) {
  if (!inherits(x, "microeda_da")) {
    stop("`x` must be a microeda_da object.", call. = FALSE)
  }
  if (!is.character(x$methods) || length(x$methods) < 1 ||
      any(is.na(x$methods)) || any(!nzchar(x$methods)) ||
      anyDuplicated(x$methods)) {
    stop("`x$methods` must contain unique method IDs.", call. = FALSE)
  }
  required_plan <- c("contrast", "group1", "group2", "contrast_type")
  if (!is.data.frame(x$contrast_plan) || nrow(x$contrast_plan) < 1 ||
      !all(required_plan %in% names(x$contrast_plan)) ||
      anyDuplicated(x$contrast_plan$contrast)) {
    stop("`x$contrast_plan` is malformed.", call. = FALSE)
  }

  invisible(x)
}

da_raw_output_method <- function(x, method) {
  available <- x$methods
  if (is.null(method)) {
    if (length(available) == 1) {
      return(available[[1]])
    }
    stop(
      "`method` must be specified for a multi-method object. Available ",
      "methods: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  da_comparison_validate_scalar_selection(method, "method")
  if (!method %in% available) {
    stop(
      "Unknown `method` `",
      method,
      "`. Available methods: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  method
}

da_raw_output_contrast <- function(x, contrast) {
  available <- as.character(x$contrast_plan$contrast)
  if (is.null(contrast)) {
    if (length(available) == 1) {
      return(available[[1]])
    }
    stop(
      "`contrast` must be specified for an object with multiple contrasts. ",
      "Available contrasts: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  da_comparison_validate_scalar_selection(contrast, "contrast")
  if (!contrast %in% available) {
    stop(
      "Unknown `contrast` `",
      contrast,
      "`. Available contrasts: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  contrast
}

da_comparison_validate_scalar_selection <- function(x, name) {
  if (!is.character(x) || length(x) != 1 ||
      is.na(x) || !nzchar(x)) {
    stop(
      "`",
      name,
      "` must be a single non-empty character value.",
      call. = FALSE
    )
  }

  invisible(x)
}

da_comparison_select_methods <- function(x, methods) {
  available <- x$methods
  if (is.null(methods)) {
    return(available)
  }
  da_comparison_validate_vector_selection(methods, "methods")
  unknown <- setdiff(methods, available)
  if (length(unknown) > 0) {
    stop(
      "Unknown `methods`: ",
      paste(unknown, collapse = ", "),
      ". Available methods: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  available[available %in% methods]
}

da_comparison_select_contrasts <- function(x, contrast) {
  available <- as.character(x$contrast_plan$contrast)
  if (is.null(contrast)) {
    return(available)
  }
  da_comparison_validate_vector_selection(contrast, "contrast")
  unknown <- setdiff(contrast, available)
  if (length(unknown) > 0) {
    stop(
      "Unknown `contrast`: ",
      paste(unknown, collapse = ", "),
      ". Available contrasts: ",
      paste(available, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  available[available %in% contrast]
}

da_comparison_select_features <- function(x, features) {
  if (is.null(features)) {
    return(NULL)
  }
  da_comparison_validate_vector_selection(features, "features")
  available <- da_comparison_feature_order(x)
  unknown <- setdiff(features, available)
  if (length(unknown) > 0) {
    stop(
      "Unknown `features`: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  available[available %in% features]
}

da_comparison_validate_vector_selection <- function(x, name) {
  if (!is.character(x) || length(x) < 1 ||
      any(is.na(x)) || any(!nzchar(x))) {
    stop(
      "`",
      name,
      "` must be a non-empty character vector without missing values.",
      call. = FALSE
    )
  }
  if (anyDuplicated(x)) {
    stop("`", name, "` cannot contain duplicates.", call. = FALSE)
  }

  invisible(x)
}

da_comparison_validate_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop(
      "`alpha` must be a single finite number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  as.numeric(alpha)
}

da_comparison_validate_max_features <- function(max_features) {
  if (!is.numeric(max_features) || length(max_features) != 1 ||
      is.na(max_features) || max_features <= 0 ||
      (!is.infinite(max_features) &&
       (!is.finite(max_features) || max_features != floor(max_features)))) {
    stop(
      "`max_features` must be a positive integer-like number or Inf.",
      call. = FALSE
    )
  }

  if (is.infinite(max_features)) Inf else as.integer(max_features)
}

da_comparison_feature_order <- function(x) {
  metadata_ids <- character()
  if (is.data.frame(x$feature_metadata) &&
      "feature_id" %in% names(x$feature_metadata)) {
    metadata_ids <- as.character(x$feature_metadata$feature_id)
  }
  result_ids <- as.character(as_da_results(x)$feature_id)
  ids <- unique(c(metadata_ids, result_ids))
  ids[!is.na(ids) & nzchar(ids)]
}

da_comparison_row <- function(results,
                              methods,
                              contrast_row,
                              feature_id,
                              alpha) {
  feature_rows <- results[results$feature_id == feature_id, , drop = FALSE]
  out <- data.frame(
    feature_id = feature_id,
    taxon_label = da_comparison_metadata_value(
      feature_rows,
      "taxon_label",
      feature_id
    ),
    rank = da_comparison_metadata_value(feature_rows, "rank", feature_id),
    contrast = as.character(contrast_row$contrast),
    group1 = as.character(contrast_row$group1),
    group2 = as.character(contrast_row$group2),
    stringsAsFactors = FALSE
  )

  for (method in methods) {
    method_row <- feature_rows[
      feature_rows$method == method,
      ,
      drop = FALSE
    ]
    if (nrow(method_row) > 1) {
      stop(
        "Standardized results contain duplicate rows for method `",
        method,
        "`, contrast `",
        contrast_row$contrast,
        "`, and feature `",
        feature_id,
        "`.",
        call. = FALSE
      )
    }
    out <- cbind(
      out,
      da_comparison_method_block(method_row, method, alpha)
    )
  }

  out
}

da_comparison_metadata_value <- function(rows, column, feature_id) {
  values <- unique(as.character(
    rows[[column]][!is.na(rows[[column]]) & nzchar(rows[[column]])]
  ))
  if (length(values) > 1) {
    stop(
      "Conflicting `",
      column,
      "` values for feature `",
      feature_id,
      "`.",
      call. = FALSE
    )
  }
  if (length(values) == 0) NA_character_ else values[[1]]
}

da_comparison_method_block <- function(row, method, alpha) {
  if (nrow(row) == 0) {
    values <- list(
      effect = NA_real_,
      effect_type = NA_character_,
      p_value = NA_real_,
      p_adjusted = NA_real_,
      significance = NA_character_,
      effect_sign = NA_character_,
      tested = NA,
      p_adjusted_le_alpha = NA
    )
  } else {
    effect <- row$effect[[1]]
    p_value <- row$p_value[[1]]
    p_adjusted <- row$p_adjusted[[1]]
    values <- list(
      effect = effect,
      effect_type = row$effect_type[[1]],
      p_value = p_value,
      p_adjusted = p_adjusted,
      significance = row$significance[[1]],
      effect_sign = da_comparison_effect_sign(effect),
      tested = !is.na(p_value),
      p_adjusted_le_alpha = if (is.na(p_adjusted)) {
        NA
      } else {
        is.finite(p_adjusted) && p_adjusted <= alpha
      }
    )
  }

  names(values) <- paste0(method, "_", names(values))
  as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
}

da_comparison_effect_sign <- function(effect) {
  if (!is.numeric(effect) || length(effect) != 1 ||
      is.na(effect) || !is.finite(effect)) {
    return(NA_character_)
  }
  if (effect > 0) {
    return("positive")
  }
  if (effect < 0) {
    return("negative")
  }

  "zero"
}

da_empty_comparison <- function(methods) {
  out <- data.frame(
    feature_id = character(),
    taxon_label = character(),
    rank = character(),
    contrast = character(),
    group1 = character(),
    group2 = character(),
    stringsAsFactors = FALSE
  )
  for (method in methods) {
    out[[paste0(method, "_effect")]] <- numeric()
    out[[paste0(method, "_effect_type")]] <- character()
    out[[paste0(method, "_p_value")]] <- numeric()
    out[[paste0(method, "_p_adjusted")]] <- numeric()
    out[[paste0(method, "_significance")]] <- character()
    out[[paste0(method, "_effect_sign")]] <- character()
    out[[paste0(method, "_tested")]] <- logical()
    out[[paste0(method, "_p_adjusted_le_alpha")]] <- logical()
  }

  out
}

da_comparison_method_report_lines <- function(x,
                                              results,
                                              method,
                                              contrast,
                                              features,
                                              alpha) {
  rows <- results[
    results$method == method &
      results$contrast == contrast &
      results$feature_id %in% features,
    ,
    drop = FALSE
  ]
  p_value_present <- !is.na(rows$p_value)
  p_adjusted_present <- !is.na(rows$p_adjusted)
  native_na <- is.na(rows$p_value) | is.na(rows$p_adjusted)
  adjustment <- paste0(
    da_summary_collapse_unique(rows$p_adjust_method),
    " (scope: ",
    da_summary_collapse_unique(rows$p_adjust_scope),
    ")"
  )
  effect_type <- da_summary_collapse_unique(rows$effect_type)
  method_note <- da_summary_collapse_unique(rows$method_note)

  lines <- c(
    "-----------------------------------------",
    paste0("Method: ", method),
    paste0("Total rows: ", nrow(rows)),
    paste0("Features with p_value: ", sum(p_value_present)),
    paste0("Features with p_adjusted: ", sum(p_adjusted_present)),
    paste0(
      "Features with p_adjusted <= ",
      da_report_number(alpha, digits = 3),
      ": ",
      sum(p_adjusted_present & rows$p_adjusted <= alpha)
    ),
    paste0("Features with native NA: ", sum(native_na)),
    paste0("Native adjustment: ", adjustment),
    paste0("Effect type: ", effect_type),
    paste0("Interpretation: ", method_note)
  )

  c(
    lines,
    da_comparison_backend_diagnostic_lines(
      x = x,
      method = method,
      contrast = contrast,
      features = features
    )
  )
}

da_comparison_backend_diagnostic_lines <- function(x,
                                                   method,
                                                   contrast,
                                                   features) {
  raw <- as_da_raw_output(x, method = method, contrast = contrast)
  if (identical(method, "deseq2") &&
      is.data.frame(raw$feature_diagnostics)) {
    diagnostics <- raw$feature_diagnostics
    if ("feature_id" %in% names(diagnostics)) {
      diagnostics <- diagnostics[
        diagnostics$feature_id %in% features,
        ,
        drop = FALSE
      ]
    }
    labels <- c(
      all_zero = "all-zero",
      cooks_outlier = "Cook's-filtered",
      independently_filtered = "independently filtered",
      other_na = "other NA/model issue"
    )
    available <- names(labels)[names(labels) %in% names(diagnostics)]
    return(paste0(
      "DESeq2 diagnostics: ",
      paste(
        paste0(
          labels[available],
          "=",
          vapply(
            diagnostics[available],
            function(value) sum(value, na.rm = TRUE),
            integer(1)
          )
        ),
        collapse = "; "
      )
    ))
  }

  if (identical(method, "ancombc2") &&
      is.numeric(raw$backend_excluded_feature_count) &&
      length(raw$backend_excluded_feature_count) == 1) {
    return(paste0(
      "ANCOM-BC2 backend-excluded features: ",
      raw$backend_excluded_feature_count,
      "."
    ))
  }

  character()
}

da_comparison_display_rows <- function(data, methods, alpha) {
  flags <- lapply(methods, function(method) {
    values <- data[[paste0(method, "_p_adjusted")]]
    !is.na(values) & is.finite(values) & values <= alpha
  })
  keep <- Reduce(`|`, flags)
  data[keep, , drop = FALSE]
}

da_comparison_report_table <- function(data, methods) {
  columns <- c("feature_id", "taxon_label")
  for (method in methods) {
    columns <- c(
      columns,
      paste0(
        method,
        c("_effect", "_effect_sign", "_p_adjusted", "_tested")
      )
    )
  }
  data[, columns, drop = FALSE]
}

da_comparison_validate_file <- function(file) {
  if (!is.character(file) || length(file) != 1 ||
      is.na(file) || !nzchar(file)) {
    stop("`file` must be a single non-empty character path.", call. = FALSE)
  }

  invisible(file)
}

da_comparison_validate_csv_arguments <- function(dots) {
  if (length(dots) == 0) {
    return(dots)
  }
  dot_names <- names(dots)
  if (is.null(dot_names) || any(!nzchar(dot_names)) ||
      anyDuplicated(dot_names)) {
    stop("CSV arguments in `...` must be uniquely named.", call. = FALSE)
  }
  allowed <- c("quote", "na", "eol", "qmethod", "fileEncoding")
  unknown <- setdiff(dot_names, allowed)
  if (length(unknown) > 0) {
    stop(
      "Unsupported CSV argument(s): ",
      paste(unknown, collapse = ", "),
      ". Supported arguments: ",
      paste(allowed, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  dots
}
