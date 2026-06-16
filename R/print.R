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

  cat("\nRecommendations: ", nrow(x$recommendations), "\n", sep = "")
  for (i in seq_len(nrow(x$recommendations))) {
    row <- x$recommendations[i, , drop = FALSE]
    cat(
      "- [",
      row$severity,
      "] ",
      row$topic,
      ": ",
      row$recommendation,
      "\n",
      sep = ""
    )
  }

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

  cat("\nUse as_alpha_table(x) for per-sample indices")
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

  cat("\nUse as_alpha_tests(x) and as_alpha_pairwise(x) for result tables.\n")
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

  cat("\nUse x$per_sample, x$per_feature, x$per_rank, x$metadata_completeness.\n")
  invisible(x)
}
