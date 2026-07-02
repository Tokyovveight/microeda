#' Test beta diversity groups with paired dispersion diagnostics
#'
#' `microeda_beta_test()` runs a PERMANOVA-style grouped distance test together
#' with betadisper-style dispersion diagnostics using the stored distance object
#' from a `microeda_beta` result. The test is exploratory: PERMANOVA can be
#' confounded by group dispersion differences and should not be interpreted
#' without the dispersion diagnostics.
#'
#' @param x A grouped `microeda_beta` object.
#' @param permutations Number of permutations passed to `vegan`.
#' @param seed Optional random seed used while running the permutation tests.
#'
#' @return A `microeda_beta_test` object with PERMANOVA results, dispersion
#'   diagnostics, caveats, parameters, and matched call.
#' @examples
#' counts <- matrix(
#'   c(
#'     10, 0, 0,
#'     8, 1, 0,
#'     0, 9, 1,
#'     0, 7, 2
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", 1:4)
#' colnames(counts) <- paste0("ASV", 1:3)
#' metadata <- data.frame(group = c("A", "A", "B", "B"), row.names = rownames(counts))
#'
#' beta <- microeda_beta(
#'   counts,
#'   metadata = metadata,
#'   group = "group",
#'   taxa_are_rows = FALSE
#' )
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)
#'   as_beta_test_summary(beta_test)
#'   cat(microeda_beta_test_report(beta_test))
#' }
#' @export
microeda_beta_test <- function(x, permutations = 999, seed = NULL) {
  if (!inherits(x, "microeda_beta")) {
    stop("`x` must be a microeda_beta object.", call. = FALSE)
  }

  if (is.null(x$group)) {
    stop("`x` must include group metadata.", call. = FALSE)
  }

  permutations <- validate_beta_test_permutations(permutations)
  seed <- validate_beta_test_seed(seed)
  group_values <- validate_beta_test_groups(x)

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop(
      "`microeda_beta_test()` requires the optional package `vegan`. ",
      "Install it with install.packages(\"vegan\").",
      call. = FALSE
    )
  }

  old_seed <- NULL
  old_seed_exists <- FALSE
  if (!is.null(seed)) {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (old_seed_exists) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit(beta_restore_random_seed(old_seed_exists, old_seed), add = TRUE)
    set.seed(seed)
  }

  distance <- as_beta_dist(x)
  group_factor <- factor(group_values)
  metadata <- data.frame(group = group_factor)

  permanova_result <- suppressMessages(beta_capture_warnings(
    vegan::adonis2(distance ~ group, data = metadata, permutations = permutations)
  ))
  dispersion_model <- suppressMessages(beta_capture_warnings(
    vegan::betadisper(distance, group_factor)
  ))
  dispersion_test <- suppressMessages(beta_capture_warnings(
    vegan::permutest(dispersion_model$value, permutations = permutations)
  ))

  permanova <- beta_permanova_table(permanova_result$value)
  dispersion <- list(
    test = beta_dispersion_test_table(dispersion_test$value),
    groups = beta_dispersion_group_table(dispersion_model$value, group_factor)
  )
  caveats <- beta_test_caveats(
    permanova = permanova,
    dispersion = dispersion$test,
    min_group_n = min(table(group_factor)),
    vegan_warnings = unique(c(
      permanova_result$warnings,
      dispersion_model$warnings,
      dispersion_test$warnings
    ))
  )

  structure(
    list(
      method = x$method,
      group = x$group,
      n_samples = length(x$sample_ids),
      n_groups = length(levels(group_factor)),
      min_group_n = min(table(group_factor)),
      permanova = permanova,
      dispersion = dispersion,
      caveats = caveats,
      params = list(
        permutations = permutations,
        seed = seed
      ),
      call = match.call()
    ),
    class = "microeda_beta_test"
  )
}

#' Summarize beta group test results
#'
#' @param x A `microeda_beta_test` object.
#'
#' @return A one-row data frame with compact PERMANOVA and dispersion columns.
#' @export
as_beta_test_summary <- function(x) {
  if (!inherits(x, "microeda_beta_test")) {
    stop("`x` must be a microeda_beta_test object.", call. = FALSE)
  }

  permanova <- beta_test_primary_row(x$permanova)
  dispersion <- beta_test_primary_row(x$dispersion$test)

  data.frame(
    method = x$method,
    group = x$group,
    n_samples = x$n_samples,
    n_groups = x$n_groups,
    min_group_n = x$min_group_n,
    permanova_r2 = permanova$r2,
    permanova_f = permanova$statistic,
    permanova_p = permanova$p_value,
    dispersion_f = dispersion$statistic,
    dispersion_p = dispersion$p_value,
    permutations = x$params$permutations,
    stringsAsFactors = FALSE
  )
}

#' Build a compact beta group test report
#'
#' `microeda_beta_test_report()` formats the paired PERMANOVA and dispersion
#' diagnostics from [microeda_beta_test()] as plain text. The report intentionally
#' presents PERMANOVA with dispersion diagnostics and caveats instead of a single
#' significance claim.
#'
#' @param x A `microeda_beta_test` object.
#' @param digits Number of decimal places to use for numeric report values.
#'
#' @return A single character string suitable for [cat()].
#' @export
microeda_beta_test_report <- function(x, digits = 3) {
  if (!inherits(x, "microeda_beta_test")) {
    stop("`x` must be a microeda_beta_test object.", call. = FALSE)
  }

  digits <- validate_beta_report_digits(digits)
  caveats <- x$caveats
  caveat_lines <- paste0("- [", caveats$severity, "] ", caveats$message)

  lines <- c(
    "=========================================",
    "Beta group test",
    "=========================================",
    "",
    paste0("Method: ", x$method),
    paste0("Group: ", x$group),
    paste0("Samples: ", x$n_samples),
    paste0("Groups: ", x$n_groups, " (min n = ", x$min_group_n, ")"),
    paste0("Permutations: ", x$params$permutations),
    "",
    "PERMANOVA",
    beta_compare_report_table(x$permanova, digits = digits),
    "",
    "Dispersion diagnostics",
    beta_compare_report_table(x$dispersion$test, digits = digits),
    "",
    "Mean distance to group centroid",
    beta_compare_report_table(x$dispersion$groups, digits = digits),
    "",
    "Caveats",
    caveat_lines
  )

  paste(lines, collapse = "\n")
}

#' Test beta diversity groups across compared distance methods
#'
#' `microeda_beta_compare_test()` runs [microeda_beta_test()] for each
#' `microeda_beta` object stored in a grouped `microeda_beta_compare` object.
#' The results are side-by-side exploratory diagnostics, not a formal ranking
#' of distance methods.
#'
#' @param x A grouped `microeda_beta_compare` object.
#' @param permutations Number of permutations passed to `vegan` for each
#'   method.
#' @param seed Optional random seed used while running each method's
#'   permutation tests.
#'
#' @return A `microeda_beta_compare_test` object containing one
#'   `microeda_beta_test` result per method, aggregated caveats, parameters,
#'   and matched call.
#' @examples
#' counts <- matrix(
#'   c(
#'     10, 0, 0,
#'     8, 1, 0,
#'     0, 9, 1,
#'     0, 7, 2
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#' rownames(counts) <- paste0("S", 1:4)
#' colnames(counts) <- paste0("ASV", 1:3)
#' metadata <- data.frame(group = c("A", "A", "B", "B"), row.names = rownames(counts))
#'
#' beta_cmp <- microeda_beta_compare(
#'   counts,
#'   metadata = metadata,
#'   group = "group",
#'   taxa_are_rows = FALSE,
#'   methods = c("bray", "jaccard")
#' )
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   beta_cmp_test <- microeda_beta_compare_test(beta_cmp, permutations = 99, seed = 1)
#'   as_beta_compare_test_summary(beta_cmp_test)
#'   cat(microeda_beta_compare_test_report(beta_cmp_test))
#' }
#' @export
microeda_beta_compare_test <- function(x, permutations = 999, seed = NULL) {
  if (!inherits(x, "microeda_beta_compare")) {
    stop("`x` must be a microeda_beta_compare object.", call. = FALSE)
  }

  if (is.null(x$group)) {
    stop("`x` must include group metadata.", call. = FALSE)
  }

  results <- lapply(x$methods, function(method) {
    microeda_beta_test(
      x$results[[method]],
      permutations = permutations,
      seed = seed
    )
  })
  names(results) <- x$methods
  first_result <- results[[1L]]

  structure(
    list(
      results = results,
      methods = x$methods,
      group = first_result$group,
      n_samples = first_result$n_samples,
      n_groups = first_result$n_groups,
      min_group_n = first_result$min_group_n,
      caveats = beta_compare_test_caveats(results, x$methods),
      params = list(
        permutations = first_result$params$permutations,
        seed = first_result$params$seed
      ),
      call = match.call()
    ),
    class = "microeda_beta_compare_test"
  )
}

#' Summarize beta group tests across compared methods
#'
#' @param x A `microeda_beta_compare_test` object.
#'
#' @return A data frame with one row per distance method and compact PERMANOVA
#'   and dispersion columns.
#' @export
as_beta_compare_test_summary <- function(x) {
  if (!inherits(x, "microeda_beta_compare_test")) {
    stop("`x` must be a microeda_beta_compare_test object.", call. = FALSE)
  }

  rows <- lapply(x$methods, function(method) {
    as_beta_test_summary(x$results[[method]])
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

#' Build a compact report for beta group tests across compared methods
#'
#' `microeda_beta_compare_test_report()` formats side-by-side PERMANOVA and
#' dispersion diagnostics from [microeda_beta_compare_test()] as plain text.
#' The report intentionally avoids method ranking.
#'
#' @param x A `microeda_beta_compare_test` object.
#' @param digits Number of decimal places to use for numeric report values.
#'
#' @return A single character string suitable for [cat()].
#' @export
microeda_beta_compare_test_report <- function(x, digits = 3) {
  if (!inherits(x, "microeda_beta_compare_test")) {
    stop("`x` must be a microeda_beta_compare_test object.", call. = FALSE)
  }

  digits <- validate_beta_report_digits(digits)
  separator <- "========================================="
  lines <- c(
    separator,
    "Beta comparison group tests",
    separator,
    "",
    paste0("Methods: ", paste(x$methods, collapse = ", ")),
    paste0("Group: ", x$group),
    paste0("Samples: ", x$n_samples),
    paste0("Groups: ", x$n_groups, " (min n = ", x$min_group_n, ")"),
    paste0("Permutations: ", x$params$permutations),
    "",
    "These side-by-side diagnostics are not a formal method ranking."
  )

  for (method in x$methods) {
    result <- x$results[[method]]
    lines <- c(
      lines,
      "",
      separator,
      paste0("Method: ", method),
      separator,
      "",
      "PERMANOVA",
      beta_compare_report_table(result$permanova, digits = digits),
      "",
      "Dispersion diagnostics",
      beta_compare_report_table(result$dispersion$test, digits = digits),
      "",
      "Mean distance to group centroid",
      beta_compare_report_table(result$dispersion$groups, digits = digits)
    )
  }

  caveat_lines <- paste0(
    "- [",
    x$caveats$method,
    " ",
    x$caveats$severity,
    "] ",
    x$caveats$message
  )
  lines <- c(
    lines,
    "",
    "Aggregated caveats",
    caveat_lines
  )

  paste(lines, collapse = "\n")
}

validate_beta_test_permutations <- function(permutations) {
  if (!is.numeric(permutations) || length(permutations) != 1 ||
      is.na(permutations) || !is.finite(permutations) ||
      permutations != floor(permutations) || permutations < 1) {
    stop("`permutations` must be a single positive whole number.", call. = FALSE)
  }

  as.integer(permutations)
}

validate_beta_test_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("`seed` must be NULL or a single finite whole number.", call. = FALSE)
  }

  as.integer(seed)
}

validate_beta_test_groups <- function(x) {
  group_values <- x$group_values
  group_labels <- as.character(group_values)

  if (length(group_labels) != length(x$sample_ids) ||
      any(is.na(group_labels)) ||
      any(!nzchar(group_labels))) {
    stop("Group labels must be present for all samples.", call. = FALSE)
  }

  n_groups <- length(unique(group_labels))
  if (n_groups < 2) {
    stop("`x` must contain at least two groups.", call. = FALSE)
  }

  group_values
}

beta_restore_random_seed <- function(old_seed_exists, old_seed) {
  if (isTRUE(old_seed_exists)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

beta_capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = unique(warnings))
}

beta_permanova_table <- function(x) {
  table <- as.data.frame(x)
  data.frame(
    term = row.names(table),
    df = table$Df,
    sum_of_squares = table$SumOfSqs,
    r2 = table$R2,
    statistic = table$F,
    p_value = table[["Pr(>F)"]],
    stringsAsFactors = FALSE
  )
}

beta_dispersion_test_table <- function(x) {
  table <- as.data.frame(x$tab)
  data.frame(
    term = row.names(table),
    df = table$Df,
    sum_of_squares = table[["Sum Sq"]],
    mean_square = table[["Mean Sq"]],
    statistic = table$F,
    p_value = table[["Pr(>F)"]],
    stringsAsFactors = FALSE
  )
}

beta_dispersion_group_table <- function(x, group_factor) {
  distances <- x$distances
  groups <- levels(group_factor)
  rows <- lapply(groups, function(group_name) {
    values <- distances[group_factor == group_name]
    data.frame(
      group = group_name,
      n = length(values),
      mean_distance_to_centroid = mean(values),
      median_distance_to_centroid = stats::median(values),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

beta_test_primary_row <- function(x) {
  x[1L, , drop = FALSE]
}

beta_test_caveats <- function(permanova,
                              dispersion,
                              min_group_n,
                              vegan_warnings = character()) {
  rows <- list(beta_test_caveat(
    caveat_id = "permanova_dispersion_confounding",
    severity = "warning",
    message = paste(
      "PERMANOVA can be confounded by group dispersion differences;",
      "inspect dispersion diagnostics and do not interpret it alone."
    )
  ))

  dispersion_p <- beta_test_primary_row(dispersion)$p_value
  if (!is.na(dispersion_p) && dispersion_p <= 0.05) {
    rows[[length(rows) + 1L]] <- beta_test_caveat(
      caveat_id = "dispersion_difference_detected",
      severity = "warning",
      message = paste(
        "Dispersion diagnostics suggest group dispersion differences;",
        "PERMANOVA location effects may be difficult to interpret."
      )
    )
  }

  if (min_group_n < 5) {
    rows[[length(rows) + 1L]] <- beta_test_caveat(
      caveat_id = "small_group_size",
      severity = "warning",
      message = paste(
        "At least one group has fewer than five samples;",
        "permutation p-values and dispersion diagnostics may be unstable."
      )
    )
  }

  if (length(vegan_warnings) > 0) {
    for (warning_message in vegan_warnings) {
      rows[[length(rows) + 1L]] <- beta_test_caveat(
        caveat_id = "vegan_warning",
        severity = "info",
        message = paste("vegan reported:", warning_message)
      )
    }
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

beta_test_caveat <- function(caveat_id, severity, message) {
  data.frame(
    caveat_id = caveat_id,
    severity = severity,
    message = message,
    stringsAsFactors = FALSE
  )
}

beta_compare_test_caveats <- function(results, methods) {
  rows <- lapply(methods, function(method) {
    caveats <- results[[method]]$caveats
    data.frame(
      method = method,
      caveat_id = caveats$caveat_id,
      severity = caveats$severity,
      message = caveats$message,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}
