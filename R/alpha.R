#' Compute alpha diversity indices
#'
#' `microeda_alpha()` calculates per-sample alpha diversity metrics from raw
#' count-like microbiome tables. It reports classic indices and their Hill
#' number interpretation so Shannon and Simpson can be read as effective taxa.
#'
#' @inheritParams microeda_check
#'
#' @return A `microeda_alpha` object with per-sample indices and optional
#'   group summaries.
#' @export
microeda_alpha <- function(x,
                           metadata = NULL,
                           group = NULL,
                           taxa_are_rows = TRUE) {
  extracted <- microeda_extract(
    x = x,
    metadata = metadata,
    taxonomy = NULL,
    taxa_are_rows = taxa_are_rows
  )

  counts <- extracted$counts
  metadata <- extracted$metadata
  library_sizes <- rowSums(counts)
  count_type <- diagnose_count_type(counts, library_sizes)
  count_based_ok <- isTRUE(count_type$integerish) &&
    !isTRUE(count_type$looks_relative)

  alpha <- calculate_alpha_indices(
    counts = counts,
    count_based_ok = count_based_ok
  )

  group_values <- NULL
  if (!is.null(group)) {
    if (is.null(metadata)) {
      stop("`metadata` is required when `group` is supplied.", call. = FALSE)
    }
    if (!group %in% colnames(metadata)) {
      stop("`group` is not a column in `metadata`.", call. = FALSE)
    }
    group_values <- metadata[[group]]
    alpha <- insert_group_column(alpha, group = group, values = group_values)
  }

  notes <- alpha_notes(count_type, count_based_ok)

  structure(
    list(
      indices = alpha,
      group_summary = summarize_alpha_groups(alpha, group = group),
      notes = notes,
      group = group,
      count_type = count_type,
      call = match.call()
    ),
    class = "microeda_alpha"
  )
}

#' Extract per-sample alpha diversity indices
#'
#' @param x A `microeda_alpha` object.
#'
#' @return A data frame with one row per sample.
#' @export
as_alpha_table <- function(x) {
  if (!inherits(x, "microeda_alpha")) {
    stop("`x` must be a microeda_alpha object.", call. = FALSE)
  }

  x$indices
}

#' Extract alpha diversity group summaries
#'
#' @param x A `microeda_alpha` object.
#'
#' @return A data frame with summary statistics by group and index.
#' @export
as_alpha_summary <- function(x) {
  if (!inherits(x, "microeda_alpha")) {
    stop("`x` must be a microeda_alpha object.", call. = FALSE)
  }

  x$group_summary
}

#' Plot an alpha diversity metric
#'
#' `microeda_alpha_plot()` draws a minimal base R barplot for one numeric alpha
#' diversity metric across samples.
#'
#' @param x A `microeda_alpha` object.
#' @param metric Optional metric column to plot. If `NULL`, the first available
#'   alpha diversity metric is used.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return The value returned by [graphics::barplot()], invisibly.
#' @examples
#' counts <- matrix(c(10, 0, 0, 5, 20, 0, 1, 0), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("S1", "S2")
#' colnames(counts) <- paste0("ASV", 1:4)
#'
#' alpha <- microeda_alpha(counts, taxa_are_rows = FALSE)
#' microeda_alpha_plot(alpha)
#' microeda_alpha_plot(alpha, metric = "shannon")
#' @export
microeda_alpha_plot <- function(x, metric = NULL, ...) {
  if (!inherits(x, "microeda_alpha")) {
    stop("`x` must be a microeda_alpha object.", call. = FALSE)
  }

  alpha_table <- x$indices
  metric_names <- alpha_numeric_metric_names(alpha_table, group = x$group)

  if (is.null(metric)) {
    metric <- default_alpha_plot_metric(metric_names)
  } else if (!is.character(metric) || length(metric) != 1 ||
             is.na(metric) || !nzchar(metric)) {
    stop("`metric` must be a single non-missing character string.", call. = FALSE)
  }

  if (!metric %in% metric_names) {
    stop(
      "`metric` must name a numeric alpha metric. Available metrics: ",
      paste(metric_names, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  heights <- alpha_table[[metric]]
  names(heights) <- alpha_table$sample_id

  dots <- list(...)
  if (is.null(dots$names.arg)) {
    dots$names.arg <- names(heights)
  }
  if (is.null(dots$xlab)) {
    dots$xlab <- "Sample"
  }
  if (is.null(dots$ylab)) {
    dots$ylab <- metric
  }
  if (is.null(dots$main)) {
    dots$main <- paste("Alpha diversity:", metric)
  }

  invisible(do.call(graphics::barplot, c(list(height = heights), dots)))
}

#' Compare alpha diversity indices across groups
#'
#' `microeda_alpha_compare()` runs exploratory group comparisons for alpha
#' diversity indices. It uses Kruskal-Wallis tests as omnibus checks and, by
#' default, pairwise Wilcoxon rank-sum tests with p-value correction.
#'
#' @param x A `microeda_alpha` object, `phyloseq` object, matrix, or data frame.
#' @param metadata Optional sample metadata when `x` is not a
#'   `microeda_alpha` object.
#' @param group Grouping column. Optional when `x` is a `microeda_alpha` object
#'   that already contains a group column.
#' @param taxa_are_rows For matrix/data frame inputs, whether rows are taxa.
#' @param indices Alpha indices to compare.
#' @param p_adjust_method P-value adjustment method passed to `p.adjust()`.
#' @param pairwise Whether to run pairwise Wilcoxon tests.
#'
#' @return A `microeda_alpha_compare` object.
#' @export
microeda_alpha_compare <- function(x,
                                   metadata = NULL,
                                   group = NULL,
                                   taxa_are_rows = TRUE,
                                   indices = c(
                                     "observed",
                                     "chao1",
                                     "hill_q1",
                                     "hill_q2",
                                     "shannon",
                                     "inverse_simpson",
                                     "pielou_evenness",
                                     "goods_coverage"
                                   ),
                                   p_adjust_method = "BH",
                                   pairwise = TRUE) {
  alpha <- coerce_to_alpha(
    x = x,
    metadata = metadata,
    group = group,
    taxa_are_rows = taxa_are_rows
  )

  alpha_table <- alpha$indices
  group <- resolve_alpha_group(alpha_table, group = group, stored_group = alpha$group)
  indices <- validate_alpha_indices(alpha_table, indices)

  tests <- do.call(rbind, lapply(indices, function(index) {
    run_alpha_omnibus(alpha_table, group = group, index = index)
  }))
  rownames(tests) <- NULL
  tests$p_value_adjusted <- stats::p.adjust(tests$p_value, method = p_adjust_method)
  tests$p_adjust_method <- p_adjust_method

  pairwise_tests <- data.frame()
  if (isTRUE(pairwise)) {
    pairwise_tests <- do.call(rbind, lapply(indices, function(index) {
      run_alpha_pairwise(
        alpha_table = alpha_table,
        group = group,
        index = index,
        p_adjust_method = p_adjust_method
      )
    }))
    rownames(pairwise_tests) <- NULL
  }

  diagnostics <- alpha_compare_diagnostics(
    alpha_table = alpha_table,
    group = group,
    p_adjust_method = p_adjust_method
  )

  notes <- alpha_compare_notes(
    alpha = alpha,
    tests = tests,
    diagnostics = diagnostics
  )

  structure(
    list(
      tests = tests,
      pairwise = pairwise_tests,
      group_summary = alpha$group_summary,
      diagnostics = diagnostics,
      notes = notes,
      alpha = alpha,
      group = group,
      call = match.call()
    ),
    class = "microeda_alpha_compare"
  )
}

#' Extract alpha diversity omnibus tests
#'
#' @param x A `microeda_alpha_compare` object.
#'
#' @return A data frame of Kruskal-Wallis test results.
#' @export
as_alpha_tests <- function(x) {
  if (!inherits(x, "microeda_alpha_compare")) {
    stop("`x` must be a microeda_alpha_compare object.", call. = FALSE)
  }

  x$tests
}

#' Extract alpha diversity pairwise tests
#'
#' @param x A `microeda_alpha_compare` object.
#'
#' @return A data frame of pairwise Wilcoxon test results.
#' @export
as_alpha_pairwise <- function(x) {
  if (!inherits(x, "microeda_alpha_compare")) {
    stop("`x` must be a microeda_alpha_compare object.", call. = FALSE)
  }

  x$pairwise
}

alpha_numeric_metric_names <- function(alpha_table, group = NULL) {
  excluded <- c("sample_id", group)
  metric_names <- names(alpha_table)[vapply(alpha_table, is.numeric, logical(1))]
  setdiff(metric_names, excluded)
}

default_alpha_plot_metric <- function(metric_names) {
  preferred <- c(
    "observed",
    "chao1",
    "shannon",
    "simpson",
    "inverse_simpson",
    "hill_q0",
    "hill_q1",
    "hill_q2",
    "pielou_evenness",
    "goods_coverage"
  )
  metric <- intersect(preferred, metric_names)[1]

  if (is.na(metric)) {
    metric <- metric_names[1]
  }

  if (is.na(metric)) {
    stop("No numeric alpha metrics are available to plot.", call. = FALSE)
  }

  metric
}

calculate_alpha_indices <- function(counts, count_based_ok = TRUE) {
  metrics <- t(apply(counts, 1, alpha_one_sample, count_based_ok = count_based_ok))
  metrics <- as.data.frame(metrics, stringsAsFactors = FALSE)

  data.frame(
    sample_id = rownames(counts),
    metrics,
    row.names = rownames(counts),
    check.names = FALSE
  )
}

alpha_one_sample <- function(x, count_based_ok = TRUE) {
  n_reads <- sum(x)
  observed <- sum(x > 0)

  if (n_reads > 0) {
    p <- x[x > 0] / n_reads
    shannon <- -sum(p * log(p))
    simpson_d <- sum(p^2)
    simpson <- 1 - simpson_d
    inverse_simpson <- safe_ratio(1, simpson_d)
  } else {
    shannon <- NA_real_
    simpson <- NA_real_
    inverse_simpson <- NA_real_
  }

  hill_q0 <- observed
  hill_q1 <- if (!is.na(shannon)) exp(shannon) else NA_real_
  hill_q2 <- inverse_simpson
  pielou_evenness <- if (observed > 1) shannon / log(observed) else NA_real_

  if (count_based_ok) {
    singletons <- sum(x == 1)
    doubletons <- sum(x == 2)
    chao1 <- observed + singletons * (singletons - 1) / (2 * (doubletons + 1))
    goods_coverage <- if (n_reads > 0) 1 - singletons / n_reads else NA_real_
  } else {
    singletons <- NA_real_
    doubletons <- NA_real_
    chao1 <- NA_real_
    goods_coverage <- NA_real_
  }

  c(
    n_reads = n_reads,
    observed = observed,
    chao1 = chao1,
    shannon = shannon,
    simpson = simpson,
    inverse_simpson = inverse_simpson,
    hill_q0 = hill_q0,
    hill_q1 = hill_q1,
    hill_q2 = hill_q2,
    pielou_evenness = pielou_evenness,
    goods_coverage = goods_coverage,
    singletons = singletons,
    doubletons = doubletons
  )
}

insert_group_column <- function(alpha, group, values) {
  group_frame <- data.frame(values, stringsAsFactors = FALSE)
  names(group_frame) <- group

  out <- cbind(
    alpha["sample_id"],
    group_frame,
    alpha[setdiff(names(alpha), "sample_id")]
  )
  rownames(out) <- rownames(alpha)
  out
}

summarize_alpha_groups <- function(alpha, group = NULL) {
  if (is.null(group) || !group %in% names(alpha)) {
    return(data.frame())
  }

  numeric_columns <- names(alpha)[vapply(alpha, is.numeric, logical(1))]
  numeric_columns <- setdiff(numeric_columns, c("singletons", "doubletons"))
  group_values <- as.character(alpha[[group]])
  group_values[is.na(group_values)] <- "<NA>"

  summaries <- lapply(numeric_columns, function(index_name) {
    values <- alpha[[index_name]]
    by_group <- split(values, group_values, drop = TRUE)

    do.call(rbind, lapply(names(by_group), function(group_name) {
      group_index_values <- by_group[[group_name]]
      group_index_values <- group_index_values[!is.na(group_index_values)]

      if (length(group_index_values) == 0) {
        return(data.frame(
          group = group_name,
          index = index_name,
          n = 0,
          mean = NA_real_,
          sd = NA_real_,
          median = NA_real_,
          q1 = NA_real_,
          q3 = NA_real_,
          min = NA_real_,
          max = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      data.frame(
        group = group_name,
        index = index_name,
        n = length(group_index_values),
        mean = mean(group_index_values),
        sd = stats::sd(group_index_values),
        median = stats::median(group_index_values),
        q1 = stats::quantile(group_index_values, 0.25, names = FALSE),
        q3 = stats::quantile(group_index_values, 0.75, names = FALSE),
        min = min(group_index_values),
        max = max(group_index_values),
        stringsAsFactors = FALSE
      )
    }))
  })

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  summary
}

alpha_notes <- function(count_type, count_based_ok) {
  notes <- c(
    "Shannon and Simpson are also reported as Hill numbers: hill_q1 = exp(Shannon), hill_q2 = inverse Simpson."
  )

  if (!count_based_ok) {
    notes <- c(
      notes,
      "Chao1, singletons, doubletons, and Good's coverage are set to NA because the input does not look like raw integer counts."
    )
  }

  if (isTRUE(count_type$has_zero_sum_samples)) {
    notes <- c(
      notes,
      "Some samples have zero total reads; their proportion-based alpha metrics are NA."
    )
  }

  notes
}

coerce_to_alpha <- function(x, metadata = NULL, group = NULL, taxa_are_rows = TRUE) {
  if (inherits(x, "microeda_alpha")) {
    return(x)
  }

  microeda_alpha(
    x = x,
    metadata = metadata,
    group = group,
    taxa_are_rows = taxa_are_rows
  )
}

resolve_alpha_group <- function(alpha_table, group = NULL, stored_group = NULL) {
  if (is.null(group)) {
    group <- stored_group
  }

  if (is.null(group)) {
    candidate_groups <- names(alpha_table)[!vapply(alpha_table, is.numeric, logical(1))]
    candidate_groups <- setdiff(candidate_groups, "sample_id")

    if (length(candidate_groups) == 1) {
      group <- candidate_groups
    }
  }

  if (is.null(group) || !group %in% names(alpha_table)) {
    stop(
      "`group` is required and must name a group column in the alpha table.",
      call. = FALSE
    )
  }

  group
}

validate_alpha_indices <- function(alpha_table, indices) {
  missing_indices <- setdiff(indices, names(alpha_table))
  if (length(missing_indices) > 0) {
    stop(
      "Unknown alpha indices: ",
      paste(missing_indices, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric <- indices[!vapply(alpha_table[indices], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(
      "Alpha indices must be numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  indices
}

run_alpha_omnibus <- function(alpha_table, group, index) {
  test_data <- valid_alpha_test_data(alpha_table, group = group, index = index)
  group_counts <- table(test_data$group)

  base <- data.frame(
    index = index,
    method = "Kruskal-Wallis rank sum test",
    n = nrow(test_data),
    n_groups = length(group_counts),
    min_group_n = if (length(group_counts) > 0) unname(min(group_counts)) else NA_integer_,
    statistic = NA_real_,
    parameter = NA_real_,
    p_value = NA_real_,
    max_median_group = NA_character_,
    min_median_group = NA_character_,
    median_difference = NA_real_,
    stringsAsFactors = FALSE
  )

  if (nrow(test_data) == 0 || length(group_counts) < 2) {
    base$method <- "not tested"
    return(base)
  }

  medians <- tapply(test_data$value, test_data$group, stats::median)
  base$max_median_group <- names(medians)[which.max(medians)]
  base$min_median_group <- names(medians)[which.min(medians)]
  base$median_difference <- unname(max(medians) - min(medians))

  result <- tryCatch(
    suppressWarnings(stats::kruskal.test(value ~ group, data = test_data)),
    error = function(error) NULL
  )

  if (is.null(result)) {
    base$method <- "not tested"
    return(base)
  }

  base$statistic <- unname(result$statistic)
  base$parameter <- unname(result$parameter)
  base$p_value <- unname(result$p.value)
  base
}

run_alpha_pairwise <- function(alpha_table, group, index, p_adjust_method = "BH") {
  test_data <- valid_alpha_test_data(alpha_table, group = group, index = index)
  groups <- sort(unique(test_data$group))

  if (length(groups) < 2) {
    return(data.frame())
  }

  pairs <- utils::combn(groups, 2, simplify = FALSE)
  results <- do.call(rbind, lapply(pairs, function(pair) {
    group_1_values <- test_data$value[test_data$group == pair[1]]
    group_2_values <- test_data$value[test_data$group == pair[2]]

    p_value <- tryCatch(
      suppressWarnings(stats::wilcox.test(
        group_1_values,
        group_2_values,
        exact = FALSE
      )$p.value),
      error = function(error) NA_real_
    )

    data.frame(
      index = index,
      group_1 = pair[1],
      group_2 = pair[2],
      n_1 = length(group_1_values),
      n_2 = length(group_2_values),
      median_1 = stats::median(group_1_values),
      median_2 = stats::median(group_2_values),
      median_difference = stats::median(group_1_values) -
        stats::median(group_2_values),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }))

  results$p_value_adjusted <- stats::p.adjust(
    results$p_value,
    method = p_adjust_method
  )
  results$p_adjust_method <- p_adjust_method
  results$method <- "Wilcoxon rank sum test"
  results
}

valid_alpha_test_data <- function(alpha_table, group, index) {
  test_data <- data.frame(
    group = as.character(alpha_table[[group]]),
    value = alpha_table[[index]],
    stringsAsFactors = FALSE
  )

  keep <- !is.na(test_data$group) &
    trimws(test_data$group) != "" &
    !is.na(test_data$value) &
    is.finite(test_data$value)

  test_data <- test_data[keep, , drop = FALSE]
  test_data$group <- factor(test_data$group)
  test_data
}

alpha_compare_diagnostics <- function(alpha_table, group, p_adjust_method = "BH") {
  depth_test <- run_alpha_omnibus(alpha_table, group = group, index = "n_reads")

  coverage_test <- data.frame()
  if ("goods_coverage" %in% names(alpha_table)) {
    coverage_test <- run_alpha_omnibus(
      alpha_table,
      group = group,
      index = "goods_coverage"
    )
  }

  diagnostics <- rbind(depth_test, coverage_test)
  diagnostics$p_value_adjusted <- stats::p.adjust(
    diagnostics$p_value,
    method = p_adjust_method
  )
  diagnostics$p_adjust_method <- p_adjust_method
  rownames(diagnostics) <- NULL
  diagnostics
}

alpha_compare_notes <- function(alpha, tests, diagnostics) {
  notes <- c(
    "Alpha group tests are exploratory rank-based comparisons; inspect effect sizes, coverage, and sequencing depth before biological interpretation."
  )

  depth_row <- diagnostics[diagnostics$index == "n_reads", , drop = FALSE]
  if (nrow(depth_row) > 0 &&
      !is.na(depth_row$p_value_adjusted[1]) &&
      depth_row$p_value_adjusted[1] < 0.05) {
    notes <- c(
      notes,
      "Sequencing depth differs across groups after p-value correction; richness-like indices may partly reflect sampling effort."
    )
  }

  coverage_row <- diagnostics[diagnostics$index == "goods_coverage", , drop = FALSE]
  if (nrow(coverage_row) > 0 &&
      !is.na(coverage_row$p_value_adjusted[1]) &&
      coverage_row$p_value_adjusted[1] < 0.05) {
    notes <- c(
      notes,
      "Good's coverage differs across groups; consider rarefaction or coverage-based sensitivity analysis before strong alpha conclusions."
    )
  }

  if (any(tests$min_group_n < 5, na.rm = TRUE)) {
    notes <- c(
      notes,
      "At least one group has fewer than five non-missing observations for an alpha index; treat p-values as unstable."
    )
  }

  if (any(is.na(tests$p_value))) {
    notes <- c(
      notes,
      "Some alpha indices could not be tested, often because all values were missing or only one group had valid values."
    )
  }

  unique(c(notes, alpha$notes))
}
