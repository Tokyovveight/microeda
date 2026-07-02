#' @export
print.microeda_report <- function(x, ...) {
  diagnostics <- x$diagnostics

  cat("microeda report\n")
  cat("---------------\n")
  cat("Samples:  ", diagnostics$n_samples, "\n", sep = "")
  cat("Features: ", diagnostics$n_features, "\n", sep = "")
  cat("Input:    ", diagnostics$count_type$input_scale, "\n", sep = "")
  cat(
    "Zeros:    ",
    format_percent(diagnostics$sparsity$zero_fraction),
    "\n",
    sep = ""
  )
  cat(
    "Reads total: ",
    format_number(diagnostics$reads$total_reads),
    "\n",
    sep = ""
  )
  cat(
    "First ",
    diagnostics$reads$feature_read_n,
    " features by table order: ",
    format_number(diagnostics$reads$first_n_reads),
    " (",
    format_percent(diagnostics$reads$first_n_fraction),
    ")\n",
    sep = ""
  )
  cat(
    "Top ",
    diagnostics$reads$feature_read_n,
    " features by reads:       ",
    format_number(diagnostics$reads$top_n_reads),
    " (",
    format_percent(diagnostics$reads$top_n_fraction),
    ")\n",
    sep = ""
  )
  cat(
    "Library size median: ",
    format_number(diagnostics$library_size$median),
    "\n",
    sep = ""
  )
  cat(
    "Library size range:  ",
    format_number(diagnostics$library_size$min),
    " - ",
    format_number(diagnostics$library_size$max),
    "\n",
    sep = ""
  )
  cat(
    "Min reads sample: ",
    diagnostics$reads$min_sample$sample_id,
    " (#",
    diagnostics$reads$min_sample$sample_index,
    ", ",
    format_number(diagnostics$reads$min_sample$reads),
    " reads)\n",
    sep = ""
  )
  cat(
    "Max reads sample: ",
    diagnostics$reads$max_sample$sample_id,
    " (#",
    diagnostics$reads$max_sample$sample_index,
    ", ",
    format_number(diagnostics$reads$max_sample$reads),
    " reads)\n",
    sep = ""
  )

  if (!is.null(diagnostics$metadata$group)) {
    group_info <- diagnostics$metadata$group
    cat(
      "Groups:   ",
      group_info$n_groups,
      " (min n = ",
      group_info$min_n,
      ")\n",
      sep = ""
    )
  }

  if (isTRUE(diagnostics$taxonomy$provided)) {
    cat(
      "Taxa by rank: ",
      collapse_named_numbers(diagnostics$taxonomy$unique_taxa_by_rank),
      "\n",
      sep = ""
    )
  }

  recommendations <- x$recommendations
  cat("\nScreening notes: ", nrow(recommendations), "\n", sep = "")
  if (nrow(recommendations) > 0) {
    severity_counts <- table(recommendations$severity)
    cat(
      "By severity: ",
      paste0(names(severity_counts), "=", as.integer(severity_counts), collapse = "; "),
      "\n",
      sep = ""
    )
    cat(
      "Topics:      ",
      paste(sort(unique(recommendations$topic)), collapse = ", "),
      "\n",
      sep = ""
    )
  }
  cat("Use as_recommendations(x) to inspect broad screening notes.\n")

  invisible(x)
}

format_percent <- function(x) {
  paste0(format_number(100 * x), "%")
}

format_number <- function(x) {
  formatted <- format(round(x, 2), big.mark = ",", trim = TRUE, scientific = FALSE)
  formatted[is.na(x)] <- "NA"
  formatted
}

collapse_named_numbers <- function(x) {
  paste0(names(x), "=", format_number(x), collapse = "; ")
}

#' @export
print.microeda_alpha <- function(x, ...) {
  indices <- x$indices

  cat("microeda alpha diversity\n")
  cat("------------------------\n")
  cat("Samples: ", nrow(indices), "\n", sep = "")
  cat(
    "Indices: observed, chao1, shannon, simpson, inverse_simpson, ",
    "hill_q0, hill_q1, hill_q2, pielou_evenness, goods_coverage\n",
    sep = ""
  )

  if (nrow(x$group_summary) > 0) {
    groups <- unique(x$group_summary$group)
    cat("Groups:  ", paste(groups, collapse = ", "), "\n", sep = "")
  }

  cat("\nNotes:\n")
  for (note in x$notes) {
    cat("- ", note, "\n", sep = "")
  }

  cat("\nUse microeda_alpha_report(x) for a readable alpha report.\n")
  cat("Use as_alpha_table(x) for per-sample indices")
  if (nrow(x$group_summary) > 0) {
    cat(" and as_alpha_summary(x) for group summaries")
  }
  cat(".\n")

  invisible(x)
}

#' @export
print.microeda_alpha_compare <- function(x, ...) {
  cat("microeda alpha comparison\n")
  cat("-------------------------\n")
  cat("Group:   ", x$group, "\n", sep = "")
  cat("Indices: ", paste(x$tests$index, collapse = ", "), "\n", sep = "")

  if (nrow(x$tests) > 0) {
    cat("\nOmnibus tests:\n")
    for (i in seq_len(nrow(x$tests))) {
      row <- x$tests[i, , drop = FALSE]
      cat(
        "- ",
        row$index,
        ": p=",
        format_p(row$p_value),
        ", p_adj=",
        format_p(row$p_value_adjusted),
        ", max median group=",
        row$max_median_group,
        "\n",
        sep = ""
      )
    }
  }

  if (nrow(x$diagnostics) > 0) {
    cat("\nDepth / coverage diagnostics:\n")
    for (i in seq_len(nrow(x$diagnostics))) {
      row <- x$diagnostics[i, , drop = FALSE]
      cat(
        "- ",
        row$index,
        ": p_adj=",
        format_p(row$p_value_adjusted),
        "\n",
        sep = ""
      )
    }
  }

  cat("\nNotes:\n")
  for (note in x$notes) {
    cat("- ", note, "\n", sep = "")
  }

  cat(
    "\nUse microeda_alpha_report(alpha, alpha_compare = x) ",
    "for grouped test reporting.\n",
    sep = ""
  )
  cat("Use microeda_alpha_pairwise_report(x) for pairwise comparison reporting.\n")
  cat("Use as_alpha_tests(x) and as_alpha_pairwise(x) for result tables.\n")
  invisible(x)
}

#' @export
print.microeda_beta <- function(x, ...) {
  cat("<microeda_beta>\n")
  cat("Method:   ", x$method, "\n", sep = "")
  cat("Samples:  ", length(x$sample_ids), "\n", sep = "")
  cat("Group:    ", if (is.null(x$group)) "<none>" else x$group, "\n", sep = "")
  cat("Distance: ", class(x$distance)[1], " (", length(x$distance), " pairs)\n", sep = "")
  invisible(x)
}

#' @export
print.microeda_beta_compare <- function(x, ...) {
  cat("<microeda_beta_compare>\n")
  cat("Methods: ", paste(x$methods, collapse = ", "), "\n", sep = "")
  cat("Samples: ", length(x$sample_ids), "\n", sep = "")
  cat("Group:   ", if (is.null(x$group)) "<none>" else x$group, "\n", sep = "")
  cat("\nUse microeda_beta_compare_report() for a readable report.\n")
  cat("Use as_beta_compare_summary() for machine-readable distance summaries.\n")
  invisible(x)
}

#' @export
print.microeda_beta_test <- function(x, ...) {
  summary <- as_beta_test_summary(x)

  cat("<microeda_beta_test>\n")
  cat("Method:       ", x$method, "\n", sep = "")
  cat("Group:        ", x$group, "\n", sep = "")
  cat("Samples:      ", x$n_samples, "\n", sep = "")
  cat("Groups:       ", x$n_groups, " (min n = ", x$min_group_n, ")\n", sep = "")
  cat("Permutations: ", x$params$permutations, "\n", sep = "")
  cat(
    "PERMANOVA:    R2=",
    format_number(summary$permanova_r2),
    ", F=",
    format_number(summary$permanova_f),
    ", p=",
    format_p(summary$permanova_p),
    "\n",
    sep = ""
  )
  cat(
    "Dispersion:   F=",
    format_number(summary$dispersion_f),
    ", p=",
    format_p(summary$dispersion_p),
    "\n",
    sep = ""
  )

  caveats <- x$caveats
  cat("Caveats:      ", nrow(caveats), sep = "")
  if (nrow(caveats) > 0 && "severity" %in% names(caveats)) {
    severity_counts <- table(caveats$severity)
    cat(
      " (",
      paste0(
        names(severity_counts),
        "=",
        as.integer(severity_counts),
        collapse = "; "
      ),
      ")",
      sep = ""
    )
  }
  cat("\n")

  cat("\nUse microeda_beta_test_report(x) for the paired PERMANOVA/dispersion report.\n")
  cat("Use as_beta_test_summary(x) for machine-readable summary values.\n")
  invisible(x)
}

#' @export
print.microeda_beta_compare_ordination <- function(x, ...) {
  cat("<microeda_beta_compare_ordination>\n")
  cat("Methods:            ", paste(x$methods, collapse = ", "), "\n", sep = "")
  cat("Ordination method:  ", x$ordination_method, "\n", sep = "")
  cat("Dimensions:         ", x$dimensions, "\n", sep = "")
  cat("Samples:            ", length(x$sample_ids), "\n", sep = "")
  cat("Group:              ", if (is.null(x$group)) "<none>" else x$group, "\n", sep = "")
  cat("\nUse as_beta_compare_coordinates() to inspect coordinates.\n")
  invisible(x)
}

#' @export
print.microeda_beta_ordination <- function(x, ...) {
  cat("<microeda_beta_ordination>\n")
  cat("Method:          ", x$method, "\n", sep = "")
  cat("Distance method: ", x$distance_method, "\n", sep = "")
  cat("Samples:         ", length(x$sample_ids), "\n", sep = "")
  cat("Dimensions:      ", x$dimensions, "\n", sep = "")
  cat("Group:           ", if (is.null(x$group)) "<none>" else x$group, "\n", sep = "")
  invisible(x)
}

format_p <- function(x) {
  if (is.na(x)) {
    return("NA")
  }

  if (x < 0.001) {
    return("<0.001")
  }

  format(round(x, 4), trim = TRUE, scientific = FALSE)
}
#' @export
print.microeda_qc <- function(x, ...) {
  cat("microeda qc summary\n")
  cat("-------------------\n")
  cat("Per-sample:       ", nrow(x$per_sample), " samples\n", sep = "")
  cat("Per-feature:      ", nrow(x$per_feature), " features\n", sep = "")
  cat(
    "Reads total:      ",
    format_number(x$library_size_summary$total_reads),
    "\n",
    sep = ""
  )
  cat(
    "Library median:   ",
    format_number(x$library_size_summary$median),
    " (range ",
    format_number(x$library_size_summary$min),
    " - ",
    format_number(x$library_size_summary$max),
    ")\n",
    sep = ""
  )
  cat(
    "Zero libraries:   ",
    x$library_size_summary$zero_library_samples,
    " (",
    format_percent(x$sparsity_summary$zero_library_sample_fraction),
    ")\n",
    sep = ""
  )
  cat(
    "Overall zeros:    ",
    format_percent(x$sparsity_summary$overall_zero_fraction),
    "\n",
    sep = ""
  )
  cat(
    "Zero features:    ",
    x$sparsity_summary$zero_abundance_features,
    " (",
    format_percent(x$sparsity_summary$zero_abundance_feature_fraction),
    ")\n",
    sep = ""
  )
  cat(
    "Prevalence >= ",
    format_number(x$prevalence_summary$min_prevalence_threshold),
    ": ",
    x$prevalence_summary$n_features_above_threshold,
    " features (",
    format_percent(x$prevalence_summary$fraction_features_above_threshold),
    ")\n",
    sep = ""
  )
  cat(
    "One-sample features: ",
    x$prevalence_summary$n_features_detected_in_one_sample,
    " (",
    format_percent(
      x$prevalence_summary$fraction_features_detected_in_one_sample
    ),
    ")\n",
    sep = ""
  )
  cat(
    "Median prevalence: ",
    format_percent(x$prevalence_summary$median_prevalence),
    "\n",
    sep = ""
  )
  cat(
    "Top 10 feature reads: ",
    format_percent(x$feature_dominance$top_10_read_fraction),
    "\n",
    sep = ""
  )

  if (!is.null(x$per_rank)) {
    cat("Per-rank:         ", nrow(x$per_rank), " ranks\n", sep = "")
  } else {
    cat("Per-rank:         not provided\n")
  }

  if (!is.null(x$metadata_completeness)) {
    cat("Metadata columns: ", nrow(x$metadata_completeness), "\n", sep = "")
  } else {
    cat("Metadata:         not provided\n")
  }

  if (nrow(x$qc_flags) > 0) {
    cat("\nQC flags:\n")
    for (i in seq_len(nrow(x$qc_flags))) {
      flag <- x$qc_flags[i, , drop = FALSE]
      cat("- [", flag$severity, "] ", flag$message, "\n", sep = "")
    }
  } else {
    cat("\nQC flags: none\n")
  }

  cat(
    "Observations: ",
    nrow(x$qc_observations),
    " (",
    sum(x$qc_observations$severity == "warning"),
    " warning)\n",
    sep = ""
  )

  cat("\nUse microeda_qc_report(x) for a readable QC report.\n")
  cat(
    "Use x$per_sample, x$per_feature, x$library_size_summary, ",
    "x$sparsity_summary, x$prevalence_summary, x$feature_dominance, ",
    "x$per_rank, x$metadata_completeness, x$qc_observations.\n",
    sep = ""
  )
  invisible(x)
}
