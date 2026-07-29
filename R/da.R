#' Run exploratory differential representation with ALDEx2
#'
#' `microeda_da()` runs a cautious differential representation workflow for
#' microbiome count data. In this version, the only implemented backend is
#' ALDEx2. Results are standardized for downstream use while preserving the raw
#' backend output in the returned object.
#'
#' This helper is exploratory: it uses method-native p-value adjustment, does
#' not globally re-adjust backend outputs, does not rank methods, and does not
#' claim ground-truth differentially abundant taxa.
#'
#' @param x A `phyloseq` object, count matrix, or count data frame.
#' @param metadata Optional sample metadata when `x` is a count table.
#' @param taxonomy Optional taxonomy table.
#' @param group Metadata column containing group labels.
#' @param contrast Either a length-2 character vector of group levels or
#'   `"pairwise"` to run all pairwise group comparisons.
#' @param methods Differential representation methods to run. Only
#'   `"aldex2"` is implemented in this version.
#' @param tax_rank Optional taxonomy rank used for `taxon_label`.
#' @param prevalence_filter Optional filter value recorded for future use.
#'   This wrapper does not apply filtering.
#' @param min_count Optional count filter recorded for future use. This wrapper
#'   does not apply filtering.
#' @param p_adjust_method Optional future backend override request. `NULL`
#'   keeps backend-native/default p-value adjustment.
#' @param taxa_are_rows Logical; for count tables, whether taxa/features are
#'   rows.
#' @param mc.samples Number of Monte Carlo samples passed to ALDEx2.
#' @param denom Denominator strategy passed to `ALDEx2::aldex.clr()`.
#' @param paired.test Logical passed to ALDEx2 test/effect helpers.
#' @param pair_id Optional metadata column identifying matched sample pairs.
#'   Required when `paired.test = TRUE` and rejected otherwise. Within every
#'   analyzed contrast, each pair must contain exactly one `group1` sample and
#'   one `group2` sample.
#'
#' @details
#' For ALDEx2 results, `effect` retains the native ALDEx2 effect size. microeda
#' supplies deterministic internal condition labels so the between-group
#' difference is `group2 - group1`: positive values indicate higher CLR
#' abundance in `group2`, and negative values indicate higher CLR abundance in
#' `group1`.
#'
#' For paired analyses, samples are matched by `pair_id` separately within
#' every contrast. Pairs are ordered deterministically by their character
#' labels, with the `group1` and `group2` sample blocks passed to ALDEx2 in the
#' same pair order. Incomplete or ambiguous pairs cause an error; samples are
#' never removed silently.
#'
#' @return A `microeda_da` object containing standardized results, method
#'   results, preserved raw backend outputs, caveats, and parameters.
#'
#' @examples
#' counts <- matrix(
#'   c(40, 8, 20, 2,
#'     38, 9, 22, 3,
#'     42, 7, 19, 2,
#'     12, 35, 18, 4,
#'     10, 37, 17, 5,
#'     11, 33, 20, 4),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", seq_len(6))
#' colnames(counts) <- paste0("ASV", seq_len(4))
#' metadata <- data.frame(
#'   group = c("A", "A", "A", "B", "B", "B"),
#'   row.names = rownames(counts)
#' )
#'
#' if (requireNamespace("ALDEx2", quietly = TRUE)) {
#'   da <- microeda_da(
#'     counts,
#'     metadata = metadata,
#'     group = "group",
#'     contrast = c("A", "B"),
#'     taxa_are_rows = FALSE,
#'     mc.samples = 16
#'   )
#'   as_da_results(da)
#' }
#'
#' @export
microeda_da <- function(x,
                        metadata = NULL,
                        taxonomy = NULL,
                        group,
                        contrast = "pairwise",
                        methods = "aldex2",
                        tax_rank = NULL,
                        prevalence_filter = NULL,
                        min_count = NULL,
                        p_adjust_method = NULL,
                        taxa_are_rows = TRUE,
                        mc.samples = 128,
                        denom = "all",
                        paired.test = FALSE,
                        pair_id = NULL) {
  methods <- da_validate_methods(methods)
  unsupported <- setdiff(methods, "aldex2")
  if (length(unsupported) > 0) {
    stop(
      "Only methods = \"aldex2\" is implemented in this version of ",
      "microeda_da(). ANCOM-BC2 and DESeq2 are planned for later optional ",
      "backends.",
      call. = FALSE
    )
  }
  pair_id <- da_validate_pair_id_argument(
    pair_id = pair_id,
    paired.test = paired.test
  )

  context <- da_prepare_context(
    x = x,
    metadata = metadata,
    taxonomy = taxonomy,
    group = group,
    contrast = contrast,
    methods = methods,
    tax_rank = tax_rank,
    prevalence_filter = prevalence_filter,
    min_count = min_count,
    p_adjust_method = p_adjust_method,
    taxa_are_rows = taxa_are_rows,
    pair_id = pair_id
  )
  backend <- da_run_aldex2(
    context,
    mc.samples = mc.samples,
    denom = denom,
    paired.test = paired.test
  )
  da_build_result_object(context, list(aldex2 = backend))
}

#' Extract standardized differential representation results
#'
#' `as_da_results()` returns the standardized result table from a
#' `microeda_da` object. Raw backend outputs are preserved separately in
#' `x$raw_outputs`.
#'
#' @param x A `microeda_da` object.
#'
#' @return A data frame with the standardized differential representation
#'   schema.
#'
#' @export
as_da_results <- function(x) {
  if (!inherits(x, "microeda_da")) {
    stop("`x` must be a microeda_da object.", call. = FALSE)
  }

  x$results
}

#' Summarize standardized differential representation results by contrast
#'
#' `as_da_summary()` returns one compact row per method and contrast from the
#' standardized table returned by `as_da_results(x)`. Raw backend outputs are
#' not inspected or serialized by this helper; they remain available in
#' `x$raw_outputs`.
#'
#' Summary rows describe exploratory differential representation outputs, not
#' confirmed biological discoveries.
#'
#' @param x A `microeda_da` object.
#' @param alpha Adjusted p-value threshold used for summary counts.
#'
#' @return A base `data.frame` with one row per method and contrast.
#'
#' @examples
#' counts <- matrix(
#'   c(40, 8, 20, 2,
#'     38, 9, 22, 3,
#'     42, 7, 19, 2,
#'     12, 35, 18, 4,
#'     10, 37, 17, 5,
#'     11, 33, 20, 4),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", seq_len(6))
#' colnames(counts) <- paste0("ASV", seq_len(4))
#' metadata <- data.frame(
#'   group = c("A", "A", "A", "B", "B", "B"),
#'   row.names = rownames(counts)
#' )
#'
#' if (requireNamespace("ALDEx2", quietly = TRUE)) {
#'   da <- microeda_da(
#'     counts,
#'     metadata = metadata,
#'     group = "group",
#'     contrast = c("A", "B"),
#'     taxa_are_rows = FALSE,
#'     mc.samples = 16
#'   )
#'   as_da_summary(da)
#' }
#'
#' @export
as_da_summary <- function(x, alpha = 0.05) {
  if (!inherits(x, "microeda_da")) {
    stop("`x` must be a microeda_da object.", call. = FALSE)
  }
  alpha <- da_validate_report_alpha(alpha)

  results <- as_da_results(x)
  da_summary_from_results(results, alpha = alpha)
}

#' Create a compact differential representation text report
#'
#' `microeda_da_report()` formats a `microeda_da` object as a plain-text,
#' console-friendly report. It summarizes contrasts, adjusted p-value counts,
#' selected standardized rows, caveats, and where to find standardized and raw
#' backend outputs.
#'
#' The report is descriptive and exploratory. It does not change differential
#' representation calculations, globally re-adjust p-values, rank methods, or
#' claim confirmed biological discoveries.
#'
#' @param x A `microeda_da` object.
#' @param top_n Number of top standardized rows to show after sorting by
#'   `p_adjusted` ascending with missing values last. Use `0` to omit top rows.
#' @param alpha Adjusted p-value threshold used for per-contrast counts.
#' @param digits Number of digits used for numeric report values.
#'
#' @return A single character string suitable for `cat()`.
#'
#' @examples
#' counts <- matrix(
#'   c(40, 8, 20, 2,
#'     38, 9, 22, 3,
#'     42, 7, 19, 2,
#'     12, 35, 18, 4,
#'     10, 37, 17, 5,
#'     11, 33, 20, 4),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", seq_len(6))
#' colnames(counts) <- paste0("ASV", seq_len(4))
#' metadata <- data.frame(
#'   group = c("A", "A", "A", "B", "B", "B"),
#'   row.names = rownames(counts)
#' )
#'
#' if (requireNamespace("ALDEx2", quietly = TRUE)) {
#'   da <- microeda_da(
#'     counts,
#'     metadata = metadata,
#'     group = "group",
#'     contrast = c("A", "B"),
#'     taxa_are_rows = FALSE,
#'     mc.samples = 16
#'   )
#'   cat(microeda_da_report(da, top_n = 5))
#' }
#'
#' @export
microeda_da_report <- function(x, top_n = 10, alpha = 0.05, digits = 3) {
  if (!inherits(x, "microeda_da")) {
    stop("`x` must be a microeda_da object.", call. = FALSE)
  }

  top_n <- da_validate_report_integer(top_n, "top_n")
  alpha <- da_validate_report_alpha(alpha)
  digits <- da_validate_report_integer(digits, "digits")

  results <- as_da_results(x)
  contrast_values <- unique(results$contrast)
  contrast_values <- contrast_values[!is.na(contrast_values) & nzchar(contrast_values)]
  p_adjust_methods <- unique(results$p_adjust_method)
  p_adjust_methods <- da_report_value(p_adjust_methods)

  lines <- c(
    "microeda differential representation report",
    "",
    paste0("Methods: ", paste(x$methods, collapse = ", ")),
    paste0("Group: ", x$group)
  )

  if (all(x$contrast_plan$contrast_type == "pairwise")) {
    lines <- c(lines, paste0("Pairwise contrasts: ", nrow(x$contrast_plan)))
  } else {
    lines <- c(lines, paste0("Contrast: ", x$contrast_label))
  }

  lines <- c(
    lines,
    paste0("Result rows: ", nrow(results)),
    paste0("Unique features: ", length(unique(results$feature_id))),
    paste0("Contrasts: ", length(contrast_values)),
    "",
    "P-value adjustment:",
    paste0("p_adjust_method: ", paste(p_adjust_methods, collapse = ", ")),
    paste(
      "Backend-native/method-specific adjustment is used;",
      "microeda does not globally re-adjust backend outputs."
    ),
    "",
    "Per-contrast summary:"
  )

  lines <- c(
    lines,
    da_report_table_lines(
      da_report_contrast_summary(results, alpha = alpha),
      digits = digits
    ),
    "",
    "Top standardized rows by adjusted p-value:"
  )

  if (top_n == 0L) {
    lines <- c(lines, "No top rows requested.")
  } else {
    top_rows <- da_report_top_rows(results, top_n = top_n)
    lines <- c(
      lines,
      da_report_table_lines(
        top_rows[
          ,
          c(
            "method",
            "contrast",
            "feature_id",
            "taxon_label",
            "effect",
            "p_value",
            "p_adjusted",
            "significance"
          ),
          drop = FALSE
        ],
        digits = digits
      )
    )
  }

  lines <- c(
    lines,
    "",
    "Caveats:",
    da_report_caveat_lines(x$caveats),
    "",
    "These rows are exploratory method outputs, not confirmed biological discoveries.",
    "",
    "Raw output:",
    "Raw backend outputs are in x$raw_outputs.",
    "Standardized table is available with as_da_results(x)."
  )

  paste(lines, collapse = "\n")
}

#' Write standardized differential representation results to CSV
#'
#' `write_da_results()` writes only the standardized result table returned by
#' `as_da_results(x)`. Raw backend outputs are not exported or serialized by
#' this helper; they remain available in `x$raw_outputs`.
#'
#' Exported rows are exploratory method outputs, not confirmed biological
#' discoveries.
#'
#' @param x A `microeda_da` object.
#' @param file Single non-empty character path for the output CSV file.
#' @param na String used for missing values in the CSV.
#' @param quote Logical; passed to `utils::write.csv()`.
#'
#' @return The output file path, invisibly.
#'
#' @examples
#' counts <- matrix(
#'   c(40, 8, 20, 2,
#'     38, 9, 22, 3,
#'     42, 7, 19, 2,
#'     12, 35, 18, 4,
#'     10, 37, 17, 5,
#'     11, 33, 20, 4),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", seq_len(6))
#' colnames(counts) <- paste0("ASV", seq_len(4))
#' metadata <- data.frame(
#'   group = c("A", "A", "A", "B", "B", "B"),
#'   row.names = rownames(counts)
#' )
#'
#' if (requireNamespace("ALDEx2", quietly = TRUE)) {
#'   da <- microeda_da(
#'     counts,
#'     metadata = metadata,
#'     group = "group",
#'     contrast = c("A", "B"),
#'     taxa_are_rows = FALSE,
#'     mc.samples = 16
#'   )
#'   write_da_results(da, tempfile(fileext = ".csv"))
#' }
#'
#' @export
write_da_results <- function(x, file, na = "", quote = TRUE) {
  if (!inherits(x, "microeda_da")) {
    stop("`x` must be a microeda_da object.", call. = FALSE)
  }
  if (!is.character(file) || length(file) != 1 ||
      is.na(file) || !nzchar(file)) {
    stop("`file` must be a single non-empty character path.", call. = FALSE)
  }
  if (!is.character(na) || length(na) != 1 || is.na(na)) {
    stop("`na` must be a single character value.", call. = FALSE)
  }
  if (!is.logical(quote) || length(quote) != 1 || is.na(quote)) {
    stop("`quote` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  utils::write.csv(
    as_da_results(x),
    file = file,
    row.names = FALSE,
    na = na,
    quote = quote
  )

  invisible(file)
}

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
                               taxa_are_rows = TRUE,
                               pair_id = NULL) {
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
  pair_id <- da_validate_pair_id_column(metadata, pair_id)
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
      pair_id = pair_id,
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
        contrast_plan = contrast_plan,
        methods = methods,
        tax_rank = tax_rank,
        taxonomy = taxonomy
      ),
      params = list(
        methods = methods,
        contrast_plan = contrast_plan,
        pair_id = pair_id,
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

  if (!da_optional_package_available("ALDEx2")) {
    stop(
      "`da_run_aldex2()` requires the optional package `ALDEx2`.",
      call. = FALSE
    )
  }

  params <- da_validate_aldex2_params(
    mc.samples = mc.samples,
    denom = denom,
    paired.test = paired.test,
    pair_id = context$pair_id
  )

  if (nrow(context$contrast_plan) == 1 &&
      identical(context$contrast_plan$contrast_type, "explicit")) {
    return(da_run_aldex2_contrast(
      context = context,
      contrast_row = context$contrast_plan[1, , drop = FALSE],
      params = params
    ))
  }

  if (!all(context$contrast_plan$contrast_type == "pairwise")) {
    stop(
      "ALDEx2 supports explicit or pairwise contrast plans only.",
      call. = FALSE
    )
  }

  contrast_results <- lapply(seq_len(nrow(context$contrast_plan)), function(i) {
    da_run_aldex2_contrast(
      context = context,
      contrast_row = context$contrast_plan[i, , drop = FALSE],
      params = params
    )
  })
  contrast_labels <- context$contrast_plan$contrast
  names(contrast_results) <- contrast_labels

  results <- do.call(
    rbind,
    lapply(contrast_results, function(result) result$results)
  )
  row.names(results) <- NULL

  raw_contrasts <- lapply(contrast_results, function(result) result$raw_output)
  names(raw_contrasts) <- contrast_labels

  raw_output <- list(
    contrasts = raw_contrasts,
    contrast_plan = context$contrast_plan,
    params = params,
    input_orientation = "feature_by_sample",
    transposed_from_context = TRUE
  )

  da_backend_result(
    method = "aldex2",
    results = results,
    raw_output = raw_output,
    notes = da_method_notes("aldex2"),
    params = params
  )
}

da_run_aldex2_contrast <- function(context, contrast_row, params) {
  sample_plan <- da_aldex2_sample_plan(
    context = context,
    contrast_row = contrast_row,
    params = params
  )
  selected_samples <- sample_plan$sample_order
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
          sample_plan$backend_conditions,
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
    conditions = sample_plan$conditions,
    backend_conditions = sample_plan$backend_conditions,
    condition_mapping = sample_plan$condition_mapping,
    sample_order = sample_plan$sample_order,
    pair_id = params$pair_id,
    pairing = sample_plan$pairing,
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

da_aldex2_sample_plan <- function(context, contrast_row, params) {
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

  pairing <- NULL
  if (isTRUE(params$paired.test)) {
    pairing <- da_aldex2_pairing_plan(
      context = context,
      contrast_row = contrast_row,
      pair_id = params$pair_id
    )
    selected_samples <- c(
      pairing$group1_sample_id,
      pairing$group2_sample_id
    )
    conditions <- c(
      rep(group1, nrow(pairing)),
      rep(group2, nrow(pairing))
    )
  }

  backend_conditions <- ifelse(
    conditions == group1,
    "microeda_group1",
    "microeda_group2"
  )

  list(
    sample_order = selected_samples,
    conditions = conditions,
    backend_conditions = backend_conditions,
    condition_mapping = data.frame(
      backend_condition = c("microeda_group1", "microeda_group2"),
      group = c(group1, group2),
      stringsAsFactors = FALSE
    ),
    pairing = pairing
  )
}

da_aldex2_pairing_plan <- function(context, contrast_row, pair_id) {
  group_values <- as.character(context$group_values)
  sample_ids <- names(context$group_values)
  group1 <- contrast_row$group1
  group2 <- contrast_row$group2
  contrast_label <- contrast_row$contrast
  keep <- group_values %in% c(group1, group2)
  selected_samples <- sample_ids[keep]
  selected_groups <- group_values[keep]
  pair_values <- context$metadata[selected_samples, pair_id, drop = TRUE]
  pair_labels <- as.character(pair_values)
  missing_pair <- is.na(pair_values) |
    is.na(pair_labels) |
    !nzchar(trimws(pair_labels))

  if (any(missing_pair)) {
    stop(
      "`pair_id` column `",
      pair_id,
      "` contains missing or empty values for contrast `",
      contrast_label,
      "` in sample(s): ",
      paste(selected_samples[missing_pair], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  pair_order <- sort(unique(pair_labels), method = "radix")
  pair_counts <- lapply(pair_order, function(pair_label) {
    in_pair <- pair_labels == pair_label
    c(
      group1 = sum(selected_groups[in_pair] == group1),
      group2 = sum(selected_groups[in_pair] == group2),
      total = sum(in_pair)
    )
  })
  invalid <- vapply(pair_counts, function(counts) {
    counts[["group1"]] != 1L ||
      counts[["group2"]] != 1L ||
      counts[["total"]] != 2L
  }, logical(1))

  if (any(invalid)) {
    invalid_details <- vapply(which(invalid), function(i) {
      counts <- pair_counts[[i]]
      paste0(
        pair_order[i],
        " (",
        group1,
        "=",
        counts[["group1"]],
        ", ",
        group2,
        "=",
        counts[["group2"]],
        ", total=",
        counts[["total"]],
        ")"
      )
    }, character(1))
    stop(
      "Invalid paired samples for contrast `",
      contrast_label,
      "` in `",
      pair_id,
      "`: each pair must contain exactly one `",
      group1,
      "` sample and one `",
      group2,
      "` sample. Invalid pair(s): ",
      paste(invalid_details, collapse = "; "),
      ".",
      call. = FALSE
    )
  }

  group1_samples <- vapply(pair_order, function(pair_label) {
    selected_samples[
      pair_labels == pair_label & selected_groups == group1
    ]
  }, character(1))
  group2_samples <- vapply(pair_order, function(pair_label) {
    selected_samples[
      pair_labels == pair_label & selected_groups == group2
    ]
  }, character(1))

  data.frame(
    pair_id = pair_order,
    group1_sample_id = unname(group1_samples),
    group2_sample_id = unname(group2_samples),
    stringsAsFactors = FALSE
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

da_validate_pair_id_argument <- function(pair_id, paired.test) {
  if (!is.logical(paired.test) || length(paired.test) != 1 ||
      is.na(paired.test)) {
    stop("`paired.test` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(paired.test)) {
    if (is.null(pair_id)) {
      stop(
        "`pair_id` must name a metadata column when `paired.test = TRUE`.",
        call. = FALSE
      )
    }
  } else if (!is.null(pair_id)) {
    stop(
      "`pair_id` can only be supplied when `paired.test = TRUE`.",
      call. = FALSE
    )
  }

  if (!is.null(pair_id) &&
      (!is.character(pair_id) || length(pair_id) != 1 ||
       is.na(pair_id) || !nzchar(pair_id))) {
    stop("`pair_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  pair_id
}

da_validate_pair_id_column <- function(metadata, pair_id) {
  if (is.null(pair_id)) {
    return(NULL)
  }

  if (!is.character(pair_id) || length(pair_id) != 1 ||
      is.na(pair_id) || !nzchar(pair_id)) {
    stop("`pair_id` must be NULL or a single non-empty character string.", call. = FALSE)
  }

  if (!pair_id %in% colnames(metadata)) {
    stop(
      "`pair_id` column `",
      pair_id,
      "` is not present in `metadata`.",
      call. = FALSE
    )
  }

  pair_id
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

da_validate_aldex2_params <- function(mc.samples,
                                      denom,
                                      paired.test,
                                      pair_id = NULL) {
  if (!is.numeric(mc.samples) || length(mc.samples) != 1 ||
      is.na(mc.samples) || !is.finite(mc.samples) ||
      mc.samples < 1 || mc.samples != floor(mc.samples)) {
    stop("`mc.samples` must be a positive whole number.", call. = FALSE)
  }

  if (!is.character(denom) || length(denom) != 1 ||
      is.na(denom) || !nzchar(denom)) {
    stop("`denom` must be a single non-empty character string.", call. = FALSE)
  }

  pair_id <- da_validate_pair_id_argument(pair_id, paired.test)

  list(
    mc.samples = as.integer(mc.samples),
    denom = denom,
    paired.test = paired.test,
    pair_id = pair_id
  )
}

da_validate_report_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 ||
      is.na(x) || !is.finite(x) || x < 0 || x != floor(x)) {
    stop(
      "`",
      name,
      "` must be a single non-negative integer-like number.",
      call. = FALSE
    )
  }

  as.integer(x)
}

da_validate_report_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      is.na(alpha) || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single numeric value between 0 and 1.", call. = FALSE)
  }

  alpha
}

da_optional_package_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

da_context_caveats <- function(counts,
                               group_values,
                               contrast_plan,
                               methods,
                               tax_rank,
                               taxonomy) {
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

  small_group_caveats <- da_small_group_caveats(
    group_values = group_values,
    contrast_plan = contrast_plan
  )
  if (nrow(small_group_caveats) > 0) {
    rows <- c(
      rows,
      split(small_group_caveats, seq_len(nrow(small_group_caveats)))
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

da_small_group_caveats <- function(group_values, contrast_plan) {
  rows <- lapply(seq_len(nrow(contrast_plan)), function(i) {
    contrast_row <- contrast_plan[i, , drop = FALSE]
    group1 <- contrast_row$group1
    group2 <- contrast_row$group2
    group_labels <- as.character(group_values)
    group1_n <- sum(group_labels == group1)
    group2_n <- sum(group_labels == group2)

    if (min(group1_n, group2_n) >= 5) {
      return(NULL)
    }

    da_caveat(
      method = NA_character_,
      caveat_id = "small_group_size",
      topic = "group_design",
      severity = "warning",
      message = paste0(
        "Contrast ",
        contrast_row$contrast,
        " has a group with fewer than five samples (",
        group1,
        "=",
        group1_n,
        "; ",
        group2,
        "=",
        group2_n,
        "); DA results may be unstable."
      )
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(da_empty_caveats())
  }

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
        "interpret results with its Monte Carlo Dirichlet assumptions.",
        "Its native effect is oriented as group2 minus group1:",
        "positive values favor group2 and negative values favor group1."
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

da_report_contrast_summary <- function(results, alpha) {
  out <- data.frame(
    contrast = character(),
    tested_features = integer(),
    p_adjusted_le_alpha = integer(),
    min_p_adjusted = numeric(),
    stringsAsFactors = FALSE
  )

  if (nrow(results) == 0) {
    return(out)
  }

  contrast_values <- unique(results$contrast)
  rows <- lapply(contrast_values, function(contrast) {
    rows <- results[results$contrast == contrast, , drop = FALSE]
    p_adjusted <- rows$p_adjusted
    p_adjusted_present <- !is.na(p_adjusted)
    min_p_adjusted <- if (any(p_adjusted_present)) {
      min(p_adjusted[p_adjusted_present])
    } else {
      NA_real_
    }

    data.frame(
      contrast = contrast,
      tested_features = length(unique(rows$feature_id)),
      p_adjusted_le_alpha = sum(p_adjusted_present & p_adjusted <= alpha),
      min_p_adjusted = min_p_adjusted,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_summary_from_results <- function(results, alpha) {
  out <- da_empty_summary()
  if (nrow(results) == 0) {
    return(out)
  }

  keys <- unique(results[c("method", "contrast")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    method <- keys$method[[i]]
    contrast <- keys$contrast[[i]]
    group_rows <- results[
      results$method == method & results$contrast == contrast,
      ,
      drop = FALSE
    ]
    da_summary_row(group_rows, alpha = alpha)
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

da_empty_summary <- function() {
  data.frame(
    method = character(),
    contrast = character(),
    group1 = character(),
    group2 = character(),
    tested_features = integer(),
    n_p_value = integer(),
    n_p_adjusted = integer(),
    n_p_adjusted_le_alpha = integer(),
    min_p_adjusted = numeric(),
    top_feature_id = character(),
    top_taxon_label = character(),
    top_effect = numeric(),
    top_p_adjusted = numeric(),
    p_adjust_method = character(),
    p_adjust_scope = character(),
    stringsAsFactors = FALSE
  )
}

da_summary_row <- function(rows, alpha) {
  p_value_present <- !is.na(rows$p_value)
  p_adjusted_present <- !is.na(rows$p_adjusted)
  top_index <- da_summary_top_index(rows$p_adjusted)

  if (is.na(top_index)) {
    top_feature_id <- NA_character_
    top_taxon_label <- NA_character_
    top_effect <- NA_real_
    top_p_adjusted <- NA_real_
  } else {
    top_feature_id <- rows$feature_id[[top_index]]
    top_taxon_label <- rows$taxon_label[[top_index]]
    top_effect <- rows$effect[[top_index]]
    top_p_adjusted <- rows$p_adjusted[[top_index]]
  }

  data.frame(
    method = da_summary_collapse_unique(rows$method),
    contrast = da_summary_collapse_unique(rows$contrast),
    group1 = da_summary_collapse_unique(rows$group1),
    group2 = da_summary_collapse_unique(rows$group2),
    tested_features = length(unique(rows$feature_id)),
    n_p_value = sum(p_value_present),
    n_p_adjusted = sum(p_adjusted_present),
    n_p_adjusted_le_alpha = sum(p_adjusted_present & rows$p_adjusted <= alpha),
    min_p_adjusted = if (any(p_adjusted_present)) {
      min(rows$p_adjusted[p_adjusted_present])
    } else {
      NA_real_
    },
    top_feature_id = top_feature_id,
    top_taxon_label = top_taxon_label,
    top_effect = top_effect,
    top_p_adjusted = top_p_adjusted,
    p_adjust_method = da_summary_collapse_unique(rows$p_adjust_method),
    p_adjust_scope = da_summary_collapse_unique(rows$p_adjust_scope),
    stringsAsFactors = FALSE
  )
}

da_summary_top_index <- function(p_adjusted) {
  present <- which(!is.na(p_adjusted))
  if (length(present) == 0) {
    return(NA_integer_)
  }

  present[which.min(p_adjusted[present])]
}

da_summary_collapse_unique <- function(x) {
  values <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
  if (length(values) == 0) {
    return(NA_character_)
  }

  paste(values, collapse = "; ")
}

da_report_top_rows <- function(results, top_n) {
  if (nrow(results) == 0) {
    return(results[0, , drop = FALSE])
  }

  order_index <- order(is.na(results$p_adjusted), results$p_adjusted)
  sorted <- results[order_index, , drop = FALSE]
  sorted[seq_len(min(top_n, nrow(sorted))), , drop = FALSE]
}

da_report_table_lines <- function(data, digits) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  headers <- names(data)
  if (length(headers) == 0) {
    return("No columns.")
  }

  formatted <- lapply(data, function(column) {
    if (is.numeric(column)) {
      return(da_report_number(column, digits = digits))
    }

    da_report_value(column)
  })

  widths <- vapply(seq_along(headers), function(i) {
    max(nchar(c(headers[i], formatted[[i]])), na.rm = TRUE)
  }, integer(1))

  format_row <- function(values) {
    paste(
      vapply(seq_along(values), function(i) {
        paste0(values[[i]], strrep(" ", widths[[i]] - nchar(values[[i]])))
      }, character(1)),
      collapse = "  "
    )
  }

  lines <- format_row(headers)
  if (nrow(data) == 0) {
    return(c(lines, "No rows available."))
  }

  row_lines <- vapply(seq_len(nrow(data)), function(i) {
    format_row(vapply(formatted, `[[`, character(1), i))
  }, character(1))

  c(lines, row_lines)
}

da_report_number <- function(x, digits) {
  out <- format(round(x, digits), trim = TRUE, scientific = FALSE)
  out[is.na(x)] <- "NA"
  out
}

da_report_value <- function(x) {
  if (length(x) == 0) {
    return("NA")
  }

  out <- as.character(x)
  out[is.na(out) | !nzchar(out)] <- "NA"
  out
}

da_report_caveat_lines <- function(caveats) {
  if (is.null(caveats) || nrow(caveats) == 0) {
    return("None.")
  }

  vapply(seq_len(nrow(caveats)), function(i) {
    row <- caveats[i, , drop = FALSE]
    method <- if (is.na(row$method) || !nzchar(row$method)) {
      "input"
    } else {
      row$method
    }
    paste0(
      "- [",
      row$severity,
      "] ",
      method,
      "/",
      row$caveat_id,
      ": ",
      row$message
    )
  }, character(1))
}
