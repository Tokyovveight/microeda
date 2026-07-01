#' Compare beta diversity distance methods
#'
#' `microeda_beta_compare()` computes multiple `microeda_beta` results from the
#' same input so distance summaries can be compared consistently.
#'
#' @inheritParams microeda_beta
#' @param methods Distance methods to compute. Supported values are `"bray"`,
#'   `"jaccard"`, and `"hellinger"`.
#'
#' @return A `microeda_beta_compare` object containing named `microeda_beta`
#'   results, method names, sample IDs, optional group information, count-type
#'   diagnostics, and matched call.
#' @examples
#' counts <- matrix(
#'   c(
#'     1, 2, 0,
#'     2, 1, 0,
#'     0, 0, 3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("S1", "S2", "S3")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
#' as_beta_compare_summary(beta_cmp)
#' @export
microeda_beta_compare <- function(x,
                                  metadata = NULL,
                                  group = NULL,
                                  taxa_are_rows = TRUE,
                                  methods = c("bray", "jaccard", "hellinger")) {
  methods <- validate_beta_compare_methods(methods)
  results <- lapply(methods, function(method) {
    microeda_beta(
      x = x,
      metadata = metadata,
      group = group,
      taxa_are_rows = taxa_are_rows,
      method = method
    )
  })
  names(results) <- methods

  first_result <- results[[1]]
  structure(
    list(
      results = results,
      methods = methods,
      sample_ids = first_result$sample_ids,
      group = first_result$group,
      group_values = first_result$group_values,
      count_type = first_result$count_type,
      call = match.call()
    ),
    class = "microeda_beta_compare"
  )
}

#' Summarize beta diversity method comparisons
#'
#' @param x A `microeda_beta_compare` object.
#'
#' @return A data frame with one row per method and distance summary columns.
#' @export
as_beta_compare_summary <- function(x) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  summaries <- lapply(x$results, beta_compare_summary_row)
  do.call(rbind, summaries)
}

#' Extract beta diversity method comparison distances
#'
#' @param x A `microeda_beta_compare` object.
#'
#' @return A data frame with one row per method and sample pair. Group columns
#'   are included when the comparison object has group metadata.
#' @export
as_beta_compare_distances <- function(x) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  rows <- lapply(x$methods, function(method) {
    beta_compare_distance_rows(
      x$results[[method]],
      group = x$group,
      group_values = x$group_values
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

#' Correlate beta diversity distance methods
#'
#' `as_beta_compare_distance_correlations()` computes descriptive correlations
#' between pairwise distance vectors from a `microeda_beta_compare` object.
#'
#' Distance-vector correlations are descriptive only. They do not identify the
#' correct method and do not replace PERMANOVA, dispersion checks, or
#' compositional diagnostics.
#'
#' @param x A `microeda_beta_compare` object.
#' @param correlation_method Correlation method. Supported values are
#'   `"pearson"`, `"spearman"`, and `"kendall"`.
#'
#' @return A data frame with one row per method pair.
#' @examples
#' counts <- matrix(
#'   c(
#'     1, 2, 0,
#'     2, 1, 0,
#'     0, 0, 3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("S1", "S2", "S3")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
#' as_beta_compare_distance_correlations(beta_cmp)
#' @export
as_beta_compare_distance_correlations <- function(x,
                                                  correlation_method = "spearman") {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  correlation_method <- validate_beta_correlation_method(correlation_method)
  if (length(x$methods) < 2) {
    return(beta_empty_distance_correlations())
  }

  distances <- as_beta_compare_distances(x)
  method_pairs <- utils::combn(x$methods, 2, simplify = FALSE)
  rows <- lapply(method_pairs, function(method_pair) {
    beta_distance_correlation_row(
      method_1 = method_pair[1],
      method_2 = method_pair[2],
      distances = distances,
      correlation_method = correlation_method
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

#' Summarize beta diversity method comparisons by group
#'
#' @param x A `microeda_beta_compare` object with group metadata.
#'
#' @return A data frame with one row per method and comparison level.
#' @export
as_beta_compare_group_summary <- function(x) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  if (is.null(x$group)) {
    stop("`x` must include group metadata.", call. = FALSE)
  }

  distances <- as_beta_compare_distances(x)
  distances <- distances[!is.na(distances$comparison), , drop = FALSE]
  if (nrow(distances) == 0) {
    return(beta_empty_group_summary())
  }

  rows <- list()
  comparisons <- c("within", "between")
  for (method in x$methods) {
    for (comparison in comparisons) {
      subset <- distances[
        distances$method == method & distances$comparison == comparison,
        ,
        drop = FALSE
      ]
      if (nrow(subset) > 0) {
        rows[[length(rows) + 1L]] <- beta_group_summary_row(
          method = method,
          comparison = comparison,
          distances = subset$distance
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

#' Build a compact beta diversity method comparison report
#'
#' `microeda_beta_compare_report()` combines beta method summaries into a
#' compact text report. Distance-method correlations are descriptive only. This
#' report does not run PERMANOVA or make formal method recommendations.
#'
#' @param x A `microeda_beta_compare` object.
#' @param correlation_method Correlation method passed to
#'   [as_beta_compare_distance_correlations()]. Supported values are
#'   `"pearson"`, `"spearman"`, and `"kendall"`.
#' @param digits Number of decimal places to use for numeric report values.
#'
#' @return A single character string.
#' @examples
#' counts <- matrix(
#'   c(
#'     1, 2, 0,
#'     2, 1, 0,
#'     0, 0, 3
#'   ),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#' rownames(counts) <- c("S1", "S2", "S3")
#' colnames(counts) <- paste0("ASV", 1:3)
#'
#' beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
#' microeda_beta_compare_report(beta_cmp)
#' @export
microeda_beta_compare_report <- function(x,
                                         correlation_method = "spearman",
                                         digits = 3) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  digits <- validate_beta_report_digits(digits)

  correlations <- as_beta_compare_distance_correlations(
    x,
    correlation_method = correlation_method
  )

  separator <- "========================================="
  lines <- c(
    separator,
    "Beta diversity method comparison",
    separator,
    "",
    paste0("Methods: ", paste(x$methods, collapse = ", ")),
    paste0("Samples: ", length(x$sample_ids)),
    paste0("Group: ", if (is.null(x$group)) "<none>" else x$group),
    "",
    "Method-level distance summary",
    beta_compare_report_table(as_beta_compare_summary(x), digits = digits),
    "",
    "Distance-method correlations"
  )

  if (nrow(correlations) == 0) {
    lines <- c(
      lines,
      "Distance-method correlations unavailable: fewer than two beta methods supplied."
    )
  } else {
    lines <- c(lines, beta_compare_report_table(correlations, digits = digits))
  }

  lines <- c(
    lines,
    "",
    "Group-level distance summary"
  )

  if (is.null(x$group)) {
    lines <- c(
      lines,
      "Group-level distance summary unavailable: no group metadata supplied."
    )
  } else {
    lines <- c(
      lines,
      beta_compare_report_table(
        as_beta_compare_group_summary(x),
        digits = digits
      )
    )
  }

  lines <- c(
    lines,
    "",
    "Notes",
    "- Bray-Curtis: abundance-sensitive distance.",
    "- Jaccard: binary presence/absence distance.",
    "- Hellinger: square-root relative abundance transform followed by Euclidean distance.",
    "- Distance-method correlations are descriptive only and do not identify the correct method.",
    "- PERMANOVA is not implemented in this report.",
    "- Formal method recommendation is not implemented yet."
  )

  paste(lines, collapse = "\n")
}

validate_beta_report_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1 || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("`digits` must be a single non-negative whole number.", call. = FALSE)
  }

  as.integer(digits)
}

validate_beta_compare_methods <- function(methods) {
  supported_methods <- supported_beta_methods()
  if (!is.character(methods) || length(methods) < 1 ||
      any(is.na(methods)) || any(!nzchar(methods)) ||
      any(!methods %in% supported_methods)) {
    stop(
      "`methods` must contain one or more supported methods: ",
      paste0("\"", supported_methods, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(methods)) {
    stop("`methods` cannot contain duplicate values.", call. = FALSE)
  }

  methods
}

validate_beta_correlation_method <- function(correlation_method) {
  supported_methods <- c("pearson", "spearman", "kendall")
  if (!is.character(correlation_method) ||
      length(correlation_method) != 1 ||
      is.na(correlation_method) ||
      !correlation_method %in% supported_methods) {
    stop(
      "`correlation_method` must be one of: ",
      paste0("\"", supported_methods, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  correlation_method
}

beta_compare_summary_row <- function(x) {
  distances <- as.numeric(as_beta_dist(x))
  if (length(distances) == 0) {
    min_distance <- NA_real_
    median_distance <- NA_real_
    max_distance <- NA_real_
  } else {
    min_distance <- min(distances)
    median_distance <- stats::median(distances)
    max_distance <- max(distances)
  }

  data.frame(
    method = x$method,
    n_samples = length(x$sample_ids),
    n_pairs = length(distances),
    min_distance = min_distance,
    median_distance = median_distance,
    max_distance = max_distance,
    stringsAsFactors = FALSE
  )
}

beta_compare_distance_rows <- function(x, group = NULL, group_values = NULL) {
  distance <- as_beta_dist(x)
  distances <- as.numeric(distance)
  sample_pairs <- beta_dist_sample_pairs(distance)
  out <- data.frame(
    method = rep(x$method, nrow(sample_pairs)),
    sample_1 = sample_pairs[, 1],
    sample_2 = sample_pairs[, 2],
    stringsAsFactors = FALSE
  )

  if (!is.null(group)) {
    group_1 <- unname(as.character(group_values[out$sample_1]))
    group_2 <- unname(as.character(group_values[out$sample_2]))
    out$group_1 <- group_1
    out$group_2 <- group_2
    out$comparison <- beta_group_comparison(group_1, group_2)
  }

  out$distance <- distances
  out
}

beta_distance_correlation_row <- function(method_1,
                                          method_2,
                                          distances,
                                          correlation_method) {
  distance_1 <- distances$distance[distances$method == method_1]
  distance_2 <- distances$distance[distances$method == method_2]
  correlation <- beta_distance_correlation(
    distance_1,
    distance_2,
    correlation_method
  )

  data.frame(
    method_1 = method_1,
    method_2 = method_2,
    n_pairs = length(distance_1),
    correlation = correlation,
    correlation_method = correlation_method,
    stringsAsFactors = FALSE
  )
}

beta_distance_correlation <- function(distance_1,
                                      distance_2,
                                      correlation_method) {
  if (length(distance_1) < 2 || length(distance_2) < 2) {
    return(NA_real_)
  }

  unname(suppressWarnings(stats::cor(
    distance_1,
    distance_2,
    method = correlation_method
  )))
}

beta_empty_distance_correlations <- function() {
  data.frame(
    method_1 = character(),
    method_2 = character(),
    n_pairs = integer(),
    correlation = numeric(),
    correlation_method = character(),
    stringsAsFactors = FALSE
  )
}

beta_dist_sample_pairs <- function(distance) {
  labels <- attr(distance, "Labels")
  if (length(distance) == 0) {
    return(matrix(
      character(),
      ncol = 2,
      dimnames = list(NULL, c("sample_1", "sample_2"))
    ))
  }

  pairs <- t(utils::combn(labels, 2))
  colnames(pairs) <- c("sample_1", "sample_2")
  pairs
}

beta_group_comparison <- function(group_1, group_2) {
  comparison <- ifelse(group_1 == group_2, "within", "between")
  comparison[is.na(group_1) | is.na(group_2)] <- NA_character_
  unname(comparison)
}

beta_group_summary_row <- function(method, comparison, distances) {
  data.frame(
    method = method,
    comparison = comparison,
    n_pairs = length(distances),
    min_distance = min(distances),
    median_distance = stats::median(distances),
    max_distance = max(distances),
    stringsAsFactors = FALSE
  )
}

beta_empty_group_summary <- function() {
  data.frame(
    method = character(),
    comparison = character(),
    n_pairs = integer(),
    min_distance = numeric(),
    median_distance = numeric(),
    max_distance = numeric(),
    stringsAsFactors = FALSE
  )
}

beta_compare_report_table <- function(x, digits) {
  values <- beta_compare_format_report_table(x, digits = digits)
  headers <- names(values)
  widths <- vapply(seq_along(values), function(i) {
    max(nchar(c(headers[i], values[[i]]), type = "width"), na.rm = TRUE)
  }, integer(1))

  format_row <- function(row_values) {
    paste(
      vapply(seq_along(row_values), function(i) {
        format(row_values[[i]], width = widths[i], justify = "left")
      }, character(1)),
      collapse = " "
    )
  }

  lines <- format_row(headers)
  if (nrow(values) == 0) {
    return(c(lines, "(no rows)"))
  }

  c(
    lines,
    vapply(seq_len(nrow(values)), function(i) {
      format_row(unname(values[i, , drop = TRUE]))
    }, character(1))
  )
}

beta_compare_format_report_table <- function(x, digits) {
  values <- as.data.frame(x, stringsAsFactors = FALSE)
  for (column in names(values)) {
    if (is.numeric(values[[column]])) {
      values[[column]] <- beta_compare_report_format_number(
        values[[column]],
        digits = digits
      )
    }
  }

  values[] <- lapply(values, as.character)
  values
}

beta_compare_report_format_number <- function(x, digits) {
  out <- rep(NA_character_, length(x))
  missing <- is.na(x)
  finite <- !missing & is.finite(x)
  whole <- finite & x == floor(x)
  decimal <- finite & !whole

  out[missing] <- "NA"
  out[whole] <- format(x[whole], trim = TRUE, scientific = FALSE)
  out[decimal] <- format(
    round(x[decimal], digits = digits),
    trim = TRUE,
    scientific = FALSE
  )
  out[!missing & !finite] <- as.character(x[!missing & !finite])
  out
}
