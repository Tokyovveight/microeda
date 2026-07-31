.qc_validate_target_request <- function(target_kingdom, target_match) {
  target_match <- match.arg(target_match, c("normalized", "exact"))

  if (is.null(target_kingdom)) {
    return(list(
      target_kingdom = NULL,
      target_match = target_match
    ))
  }

  if (!is.character(target_kingdom) || length(target_kingdom) == 0 ||
      anyNA(target_kingdom) || any(!nzchar(trimws(target_kingdom)))) {
    stop(
      "`target_kingdom` must be a non-empty character vector without NA ",
      "or empty values.",
      call. = FALSE
    )
  }

  target_kingdom <- unique(trimws(target_kingdom))
  normalized <- .qc_normalize_kingdom(target_kingdom)
  if (any(.qc_is_unclassified_kingdom(target_kingdom, normalized))) {
    stop(
      "`target_kingdom` must contain classified kingdom labels.",
      call. = FALSE
    )
  }

  list(
    target_kingdom = target_kingdom,
    target_match = target_match
  )
}

.qc_validate_original_taxonomy_ids <- function(taxonomy) {
  if (is.null(taxonomy) || is.null(rownames(taxonomy))) {
    return(invisible(NULL))
  }

  if (anyDuplicated(rownames(taxonomy))) {
    stop(
      "`taxonomy` row names must be unique for target composition ",
      "diagnostics.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

.qc_resolve_kingdom_rank <- function(taxonomy, kingdom_rank = NULL) {
  if (is.null(taxonomy)) {
    stop(
      "Target composition diagnostics require a taxonomy table.",
      call. = FALSE
    )
  }

  ranks <- colnames(taxonomy)
  if (is.null(ranks) || length(ranks) == 0) {
    stop(
      "Target composition diagnostics require named taxonomy ranks.",
      call. = FALSE
    )
  }

  available <- paste(ranks, collapse = ", ")
  normalized_ranks <- tolower(trimws(ranks))

  if (!is.null(kingdom_rank)) {
    if (!is.character(kingdom_rank) || length(kingdom_rank) != 1 ||
        is.na(kingdom_rank) || !nzchar(trimws(kingdom_rank))) {
      stop(
        "`kingdom_rank` must be NULL or one non-empty character value.",
        call. = FALSE
      )
    }

    matched <- which(
      normalized_ranks == tolower(trimws(kingdom_rank))
    )
    if (length(matched) == 0) {
      stop(
        "`kingdom_rank` was not found in taxonomy. Available ranks: ",
        available,
        ".",
        call. = FALSE
      )
    }
    if (length(matched) > 1) {
      stop(
        "`kingdom_rank` matches multiple taxonomy columns ",
        "case-insensitively: ",
        paste(ranks[matched], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(ranks[[matched]])
  }

  for (candidate in c("kingdom", "domain", "superkingdom")) {
    matched <- which(normalized_ranks == candidate)
    if (length(matched) > 1) {
      stop(
        "Automatic kingdom-rank selection is ambiguous for `",
        candidate,
        "`: ",
        paste(ranks[matched], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if (length(matched) == 1) {
      return(ranks[[matched]])
    }
  }

  stop(
    "Could not find a Kingdom, Domain, or Superkingdom taxonomy rank. ",
    "Available ranks: ",
    available,
    ".",
    call. = FALSE
  )
}

.qc_normalize_kingdom <- function(x) {
  missing <- is.na(x)
  out <- trimws(as.character(x))
  out <- sub(
    "^(d__|k__|domain__|kingdom__|d_0__)",
    "",
    out,
    ignore.case = TRUE,
    perl = TRUE
  )
  out <- tolower(trimws(out))
  out[missing | !nzchar(out)] <- NA_character_
  out
}

.qc_is_unclassified_kingdom <- function(original, normalized) {
  unclassified_terms <- c(
    "na", "n/a", "n.a.", "none", "null", "not available",
    "unclassified", "unknown", "unidentified", "uncultured",
    "unassigned"
  )
  trimmed <- trimws(as.character(original))

  is.na(original) |
    !nzchar(trimmed) |
    is.na(normalized) |
    normalized %in% unclassified_terms
}

.qc_target_composition <- function(counts,
                                   taxonomy,
                                   target_kingdom,
                                   kingdom_rank = NULL,
                                   target_match = "normalized") {
  if (anyDuplicated(rownames(counts))) {
    stop(
      "Sample IDs must be unique for target composition diagnostics.",
      call. = FALSE
    )
  }
  if (anyDuplicated(colnames(counts))) {
    stop(
      "Feature IDs must be unique for target composition diagnostics.",
      call. = FALSE
    )
  }
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop(
      "Target composition diagnostics require finite, non-negative counts.",
      call. = FALSE
    )
  }

  rank <- .qc_resolve_kingdom_rank(taxonomy, kingdom_rank)
  if (anyDuplicated(rownames(taxonomy))) {
    stop(
      "`taxonomy` row names must be unique for target composition ",
      "diagnostics.",
      call. = FALSE
    )
  }
  if (!identical(rownames(taxonomy), colnames(counts))) {
    stop(
      "Taxonomy rows must align exactly with count feature IDs.",
      call. = FALSE
    )
  }

  original <- as.character(taxonomy[[rank]])
  trimmed <- trimws(original)
  normalized <- .qc_normalize_kingdom(original)
  unclassified <- .qc_is_unclassified_kingdom(original, normalized)
  target_key <- if (identical(target_match, "normalized")) {
    .qc_normalize_kingdom(target_kingdom)
  } else {
    trimws(target_kingdom)
  }
  feature_key <- if (identical(target_match, "normalized")) {
    normalized
  } else {
    trimmed
  }

  status <- rep("non_target", length(original))
  status[unclassified] <- "unclassified"
  status[!unclassified & feature_key %in% target_key] <- "target"

  feature_totals <- unname(colSums(counts))
  total_reads <- unname(sum(feature_totals))
  prevalence <- unname(colMeans(counts > 0))
  per_feature <- data.frame(
    feature_id = colnames(counts),
    original_kingdom = original,
    normalized_kingdom = normalized,
    target_status = status,
    total_reads = feature_totals,
    prevalence = prevalence,
    overall_read_proportion = .qc_safe_divide(
      feature_totals,
      rep(total_reads, length(feature_totals))
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  overall <- .qc_target_overall(per_feature, total_reads)
  by_kingdom <- .qc_target_by_kingdom(
    per_feature = per_feature,
    total_reads = total_reads,
    rank = rank
  )
  per_sample <- .qc_target_per_sample(
    counts = counts,
    per_feature = per_feature,
    by_kingdom = by_kingdom
  )
  observations <- .qc_target_observations(
    overall = overall,
    by_kingdom = by_kingdom,
    per_sample = per_sample
  )

  list(
    target_kingdom = target_kingdom,
    kingdom_rank = rank,
    target_match = target_match,
    overall = overall,
    by_kingdom = by_kingdom,
    per_sample = per_sample,
    per_feature = per_feature,
    observations = observations,
    caveats = c(
      paste(
        "This target-composition diagnostic is descriptive; reads assigned",
        "outside the target kingdom are not automatically evidence of",
        "contamination."
      ),
      paste(
        "Calculations use the supplied raw counts. No features or samples",
        "were filtered or removed."
      ),
      paste(
        "Unclassified kingdom assignments are reported separately from",
        "known non-target assignments."
      )
    )
  )
}

.qc_target_overall <- function(per_feature, total_reads) {
  statuses <- c("target", "non_target", "unclassified")
  feature_counts <- vapply(statuses, function(status) {
    sum(per_feature$target_status == status)
  }, integer(1))
  reads <- vapply(statuses, function(status) {
    sum(per_feature$total_reads[per_feature$target_status == status])
  }, numeric(1))
  total_features <- nrow(per_feature)

  data.frame(
    total_features = as.integer(total_features),
    total_reads = as.numeric(total_reads),
    target_features = unname(feature_counts[["target"]]),
    non_target_features = unname(feature_counts[["non_target"]]),
    unclassified_features = unname(feature_counts[["unclassified"]]),
    target_feature_proportion = .qc_safe_ratio(
      feature_counts[["target"]],
      total_features
    ),
    non_target_feature_proportion = .qc_safe_ratio(
      feature_counts[["non_target"]],
      total_features
    ),
    unclassified_feature_proportion = .qc_safe_ratio(
      feature_counts[["unclassified"]],
      total_features
    ),
    target_reads = unname(reads[["target"]]),
    non_target_reads = unname(reads[["non_target"]]),
    unclassified_reads = unname(reads[["unclassified"]]),
    target_read_proportion = .qc_safe_ratio(
      reads[["target"]],
      total_reads
    ),
    non_target_read_proportion = .qc_safe_ratio(
      reads[["non_target"]],
      total_reads
    ),
    unclassified_read_proportion = .qc_safe_ratio(
      reads[["unclassified"]],
      total_reads
    ),
    stringsAsFactors = FALSE
  )
}

.qc_target_by_kingdom <- function(per_feature, total_reads, rank) {
  kingdom_label <- trimws(per_feature$original_kingdom)
  kingdom_label[
    is.na(per_feature$original_kingdom) | !nzchar(kingdom_label)
  ] <- "<unclassified>"

  grouping <- data.frame(
    kingdom_label = kingdom_label,
    normalized_kingdom = per_feature$normalized_kingdom,
    target_status = per_feature$target_status,
    stringsAsFactors = FALSE
  )
  grouping_key <- paste(
    grouping$kingdom_label,
    ifelse(is.na(grouping$normalized_kingdom), "<NA>", grouping$normalized_kingdom),
    grouping$target_status,
    sep = "\r"
  )
  first <- !duplicated(grouping_key)
  groups <- grouping[first, , drop = FALSE]
  group_keys <- grouping_key[first]

  n_features <- vapply(group_keys, function(key) {
    sum(grouping_key == key)
  }, integer(1))
  reads <- vapply(group_keys, function(key) {
    sum(per_feature$total_reads[grouping_key == key])
  }, numeric(1))
  non_target_reads <- sum(
    per_feature$total_reads[per_feature$target_status == "non_target"]
  )
  within_non_target <- rep(NA_real_, length(group_keys))
  is_non_target <- groups$target_status == "non_target"
  if (non_target_reads > 0) {
    within_non_target[is_non_target] <- (
      reads[is_non_target] / non_target_reads
    )
  }

  data.frame(
    kingdom_label = groups$kingdom_label,
    normalized_kingdom = groups$normalized_kingdom,
    target_status = groups$target_status,
    n_features = unname(n_features),
    feature_proportion = .qc_safe_divide(
      n_features,
      rep(nrow(per_feature), length(n_features))
    ),
    reads = unname(reads),
    read_proportion = .qc_safe_divide(
      reads,
      rep(total_reads, length(reads))
    ),
    proportion_within_non_target_reads = within_non_target,
    rank = rep(rank, length(group_keys)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.qc_target_per_sample <- function(counts, per_feature, by_kingdom) {
  statuses <- c("target", "non_target", "unclassified")
  status_reads <- do.call(cbind, lapply(statuses, function(status) {
    rowSums(counts[, per_feature$target_status == status, drop = FALSE])
  }))
  colnames(status_reads) <- statuses
  total_reads <- unname(rowSums(counts))
  proportions <- .qc_safe_divide(
    status_reads,
    matrix(total_reads, nrow = nrow(counts), ncol = length(statuses))
  )
  colnames(proportions) <- statuses

  dominant_index <- max.col(status_reads, ties.method = "first")
  dominant_status <- statuses[dominant_index]
  dominant_status[total_reads == 0] <- NA_character_

  non_target <- per_feature$target_status == "non_target"
  non_target_keys <- unique(per_feature$normalized_kingdom[non_target])
  non_target_keys <- non_target_keys[!is.na(non_target_keys)]
  dominant_non_target <- rep(NA_character_, nrow(counts))
  dominant_non_target_reads <- rep(0, nrow(counts))

  if (length(non_target_keys) > 0) {
    non_target_matrix <- do.call(cbind, lapply(non_target_keys, function(key) {
      rowSums(counts[
        ,
        non_target & per_feature$normalized_kingdom == key,
        drop = FALSE
      ])
    }))
    colnames(non_target_matrix) <- non_target_keys
    dominant_non_target_index <- max.col(
      non_target_matrix,
      ties.method = "first"
    )
    dominant_non_target_reads <- unname(
      non_target_matrix[
        cbind(seq_len(nrow(non_target_matrix)), dominant_non_target_index)
      ]
    )
    display_labels <- vapply(non_target_keys, function(key) {
      by_kingdom$kingdom_label[
        which(
          by_kingdom$target_status == "non_target" &
            by_kingdom$normalized_kingdom == key
        )[[1]]
      ]
    }, character(1))
    dominant_non_target <- display_labels[dominant_non_target_index]
    dominant_non_target[dominant_non_target_reads == 0] <- NA_character_
  }

  data.frame(
    sample_id = rownames(counts),
    total_reads = total_reads,
    target_reads = unname(status_reads[, "target"]),
    non_target_reads = unname(status_reads[, "non_target"]),
    unclassified_reads = unname(status_reads[, "unclassified"]),
    target_read_proportion = unname(proportions[, "target"]),
    non_target_read_proportion = unname(proportions[, "non_target"]),
    unclassified_read_proportion = unname(proportions[, "unclassified"]),
    dominant_status = dominant_status,
    dominant_non_target_kingdom = dominant_non_target,
    dominant_non_target_reads = dominant_non_target_reads,
    dominant_non_target_read_proportion = .qc_safe_divide(
      dominant_non_target_reads,
      total_reads
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.qc_target_observations <- function(overall, by_kingdom, per_sample) {
  observations <- list()

  if (overall$target_features == 0) {
    observations[[length(observations) + 1]] <- .qc_observation(
      "target_kingdom_absent",
      "target_composition",
      "warning",
      "No features were assigned to the requested target kingdom(s)."
    )
  }
  if (overall$non_target_features == 0) {
    observations[[length(observations) + 1]] <- .qc_observation(
      "non_target_kingdom_absent",
      "target_composition",
      "info",
      "No features were assigned to a known non-target kingdom."
    )
  }
  if (overall$unclassified_features > 0) {
    observations[[length(observations) + 1]] <- .qc_observation(
      "unclassified_kingdom_signal",
      "target_composition",
      "warning",
      paste0(
        overall$unclassified_features,
        " feature(s) have an unclassified kingdom assignment."
      )
    )
  }

  non_target <- by_kingdom[
    by_kingdom$target_status == "non_target",
    ,
    drop = FALSE
  ]
  if (nrow(non_target) > 0 && sum(non_target$reads) > 0) {
    top <- non_target[order(-non_target$reads, seq_len(nrow(non_target))), , drop = FALSE][1, ]
    observations[[length(observations) + 1]] <- .qc_observation(
      "dominant_non_target_kingdom",
      "target_composition",
      "info",
      paste0(
        "The largest known non-target signal by reads is `",
        top$kingdom_label,
        "` (",
        .qc_format_target_percent(top$proportion_within_non_target_reads),
        " of non-target reads)."
      )
    )
  }

  observations[[length(observations) + 1]] <- .qc_target_max_sample_observation(
    per_sample,
    proportion = "non_target_read_proportion",
    observation_id = "max_sample_non_target_signal",
    label = "non-target read"
  )
  observations[[length(observations) + 1]] <- .qc_target_max_sample_observation(
    per_sample,
    proportion = "unclassified_read_proportion",
    observation_id = "max_sample_unclassified_signal",
    label = "unclassified kingdom read"
  )
  observations[[length(observations) + 1]] <- .qc_observation(
    "target_read_share",
    "target_composition",
    "info",
    paste0(
      "Target reads account for ",
      .qc_format_target_percent(overall$target_read_proportion),
      " of all reads."
    )
  )

  out <- do.call(rbind, observations)
  rownames(out) <- NULL
  out
}

.qc_target_max_sample_observation <- function(per_sample,
                                              proportion,
                                              observation_id,
                                              label) {
  values <- per_sample[[proportion]]
  if (all(is.na(values))) {
    return(.qc_observation(
      observation_id,
      "target_composition",
      "info",
      paste0(
        "The maximum sample-level ",
        label,
        " proportion is unavailable because sample read totals are zero."
      )
    ))
  }

  index <- which.max(replace(values, is.na(values), -Inf))
  .qc_observation(
    observation_id,
    "target_composition",
    "info",
    paste0(
      "Sample `",
      per_sample$sample_id[[index]],
      "` has the highest ",
      label,
      " proportion (",
      .qc_format_target_percent(values[[index]]),
      ")."
    )
  )
}

.qc_safe_divide <- function(numerator, denominator) {
  out <- numerator
  storage.mode(out) <- "double"
  out[] <- NA_real_
  valid <- !is.na(denominator) & denominator != 0
  out[valid] <- numerator[valid] / denominator[valid]
  out
}

#' Extract target-versus-non-target QC composition tables
#'
#' `as_qc_target_composition()` returns one target-composition table created
#' by [microeda_qc()]. The diagnostic uses raw counts and separates target,
#' known non-target, and unclassified kingdom assignments. It does not filter
#' features or establish contamination. In normalized matching mode, only
#' leading `d__`, `k__`, `domain__`, `kingdom__`, and `D_0__` prefixes are
#' removed before case-insensitive comparison; fuzzy and substring matching
#' are not used.
#'
#' The sample-level `dominant_non_target_read_proportion` is the dominant known
#' non-target kingdom's reads divided by all reads in that sample. A
#' zero-library sample therefore has `NA_real_` proportions. The kingdom-level
#' `proportion_within_non_target_reads` is defined only for known non-target
#' labels and is `NA_real_` when no non-target reads are present.
#'
#' @param x A `microeda_qc` object created with a non-`NULL`
#'   `target_kingdom`.
#' @param level Table to return: `"overall"`, `"kingdom"`, `"sample"`, or
#'   `"feature"`.
#'
#' @return A base data frame. `"overall"` has one row; the other levels have
#'   one row per kingdom label, sample, or feature.
#'
#' @seealso [microeda_qc()], [microeda_qc_report()], [microeda_qc_plot()]
#'
#' @export
as_qc_target_composition <- function(
    x,
    level = c("overall", "kingdom", "sample", "feature")) {
  if (!inherits(x, "microeda_qc")) {
    stop("`x` must be a microeda_qc object.", call. = FALSE)
  }
  level <- match.arg(level)
  if (is.null(x$target_composition)) {
    stop(
      "Target composition was not computed. Run ",
      "`microeda_qc(..., target_kingdom = \"Fungi\")` first.",
      call. = FALSE
    )
  }

  table <- switch(
    level,
    overall = x$target_composition$overall,
    kingdom = x$target_composition$by_kingdom,
    sample = x$target_composition$per_sample,
    feature = x$target_composition$per_feature
  )
  as.data.frame(table, stringsAsFactors = FALSE)
}

.qc_validate_target_top_n <- function(top_n, argument = "target_top_n") {
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) ||
      !is.finite(top_n) || top_n <= 0 || top_n != floor(top_n)) {
    stop(
      "`",
      argument,
      "` must be a single positive integer-like number.",
      call. = FALSE
    )
  }

  as.integer(top_n)
}

.qc_target_report_lines <- function(target_composition, top_n) {
  overall <- target_composition$overall
  by_kingdom <- target_composition$by_kingdom
  per_sample <- target_composition$per_sample
  target_labels <- paste(target_composition$target_kingdom, collapse = ", ")

  lines <- c(
    "",
    "Target composition",
    "------------------",
    paste0("Target kingdom(s): ", target_labels),
    paste0("Taxonomy rank: ", target_composition$kingdom_rank),
    paste0("Matching mode: ", target_composition$target_match),
    "",
    "Features:",
    .qc_target_category_lines(overall, "feature"),
    "",
    "Reads:",
    .qc_target_category_lines(overall, "read"),
    "",
    "Top known non-target kingdoms by reads:"
  )

  non_target <- by_kingdom[
    by_kingdom$target_status == "non_target",
    ,
    drop = FALSE
  ]
  if (nrow(non_target) == 0) {
    lines <- c(lines, "- No known non-target kingdoms were observed.")
  } else {
    order_index <- order(-non_target$reads, seq_len(nrow(non_target)))
    non_target <- non_target[order_index, , drop = FALSE]
    non_target <- utils::head(non_target, top_n)
    lines <- c(
      lines,
      vapply(seq_len(nrow(non_target)), function(i) {
        row <- non_target[i, , drop = FALSE]
        paste0(
          "- ",
          row$kingdom_label,
          ": ",
          .qc_format_target_count(row$reads),
          " reads (",
          .qc_format_target_percent(row$read_proportion),
          " of all reads; ",
          .qc_format_target_percent(row$proportion_within_non_target_reads),
          " of non-target reads)."
        )
      }, character(1))
    )
  }

  lines <- c(
    lines,
    "",
    "Samples with highest non-target read proportions:",
    .qc_target_sample_report_lines(
      per_sample,
      proportion = "non_target_read_proportion",
      reads = "non_target_reads",
      top_n = top_n,
      include_dominant = TRUE
    ),
    "",
    "Samples with highest unclassified read proportions:",
    .qc_target_sample_report_lines(
      per_sample,
      proportion = "unclassified_read_proportion",
      reads = "unclassified_reads",
      top_n = top_n,
      include_dominant = FALSE
    ),
    "",
    "Target-composition observations:",
    paste0(
      "- [",
      target_composition$observations$severity,
      "] ",
      target_composition$observations$message
    ),
    "",
    "Target-composition caveats:",
    paste0("- ", target_composition$caveats)
  )

  lines
}

.qc_target_category_lines <- function(overall, type) {
  statuses <- c("target", "non_target", "unclassified")
  labels <- c("Target", "Non-target", "Unclassified")
  if (identical(type, "feature")) {
    counts <- unlist(overall[paste0(statuses, "_features")], use.names = FALSE)
    proportions <- unlist(
      overall[paste0(statuses, "_feature_proportion")],
      use.names = FALSE
    )
    singular <- "feature"
    plural <- "features"
  } else {
    counts <- unlist(overall[paste0(statuses, "_reads")], use.names = FALSE)
    proportions <- unlist(
      overall[paste0(statuses, "_read_proportion")],
      use.names = FALSE
    )
    singular <- "read"
    plural <- "reads"
  }

  vapply(seq_along(statuses), function(i) {
    noun <- if (counts[[i]] == 1) singular else plural
    paste0(
      "- ",
      labels[[i]],
      ": ",
      .qc_format_target_count(counts[[i]]),
      " ",
      noun,
      " (",
      .qc_format_target_percent(proportions[[i]]),
      ")."
    )
  }, character(1))
}

.qc_target_sample_report_lines <- function(per_sample,
                                           proportion,
                                           reads,
                                           top_n,
                                           include_dominant) {
  values <- per_sample[[proportion]]
  order_index <- order(
    is.na(values),
    -replace(values, is.na(values), -Inf),
    seq_len(nrow(per_sample))
  )
  selected <- per_sample[utils::head(order_index, top_n), , drop = FALSE]

  vapply(seq_len(nrow(selected)), function(i) {
    row <- selected[i, , drop = FALSE]
    line <- paste0(
      "- ",
      row$sample_id,
      ": ",
      .qc_format_target_count(row[[reads]]),
      " / ",
      .qc_format_target_count(row$total_reads),
      " reads (",
      .qc_format_target_percent(row[[proportion]]),
      ")."
    )
    if (isTRUE(include_dominant) &&
        !is.na(row$dominant_non_target_kingdom)) {
      line <- paste0(
        line,
        " Dominant known non-target: ",
        row$dominant_non_target_kingdom,
        "."
      )
    }
    line
  }, character(1))
}

.qc_format_target_count <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  format(x, trim = TRUE, scientific = FALSE, big.mark = ",")
}

.qc_format_target_percent <- function(x) {
  if (length(x) != 1 || is.na(x)) {
    return("NA")
  }
  percent <- 100 * x
  if (percent == 0) {
    return("0%")
  }
  if (abs(percent) < 0.01) {
    return(paste0(
      sub("e\\+?", "e", formatC(percent, format = "e", digits = 2)),
      "%"
    ))
  }

  out <- formatC(percent, format = "f", digits = 2)
  out <- sub("\\.?0+$", "", out)
  paste0(out, "%")
}

.qc_wrap_report_lines <- function(lines, width = 105L) {
  wrapped <- lapply(lines, function(line) {
    if (!nzchar(line) || nchar(line, type = "width") <= width) {
      return(line)
    }

    prefix <- if (startsWith(line, "- ")) "- " else ""
    subsequent <- if (nzchar(prefix)) "  " else ""
    strwrap(
      line,
      width = width,
      prefix = "",
      initial = "",
      exdent = nchar(subsequent)
    )
  })
  unlist(wrapped, use.names = FALSE)
}

.qc_target_plot <- function(x, target_view, top_n, dots) {
  if (is.null(x$target_composition)) {
    stop(
      "Target composition was not computed. Run ",
      "`microeda_qc(..., target_kingdom = \"Fungi\")` first.",
      call. = FALSE
    )
  }

  target_view <- match.arg(target_view, c("samples", "kingdoms"))
  top_n <- .qc_validate_target_top_n(top_n, argument = "top_n")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (identical(target_view, "samples")) {
    plotting_data <- .qc_target_sample_plot_data(
      x$target_composition$per_sample,
      top_n
    )
    .qc_draw_target_sample_plot(plotting_data, dots)
  } else {
    plotting_data <- .qc_target_kingdom_plot_data(
      x$target_composition$by_kingdom,
      top_n
    )
    .qc_draw_target_kingdom_plot(plotting_data, dots)
  }

  invisible(plotting_data)
}

.qc_target_sample_plot_data <- function(per_sample, top_n) {
  values <- per_sample$non_target_read_proportion
  order_index <- order(
    is.na(values),
    -replace(values, is.na(values), -Inf),
    seq_len(nrow(per_sample))
  )
  per_sample[
    utils::head(order_index, min(top_n, nrow(per_sample))),
    ,
    drop = FALSE
  ]
}

.qc_draw_target_sample_plot <- function(plotting_data, dots) {
  proportions <- rbind(
    target = plotting_data$target_read_proportion,
    non_target = plotting_data$non_target_read_proportion,
    unclassified = plotting_data$unclassified_read_proportion
  )
  proportions[is.na(proportions)] <- 0
  colnames(proportions) <- plotting_data$sample_id
  colors <- c(
    target = "#3B7A57",
    non_target = "#C85C5C",
    unclassified = "#8A8A8A"
  )
  defaults <- list(
    height = proportions,
    names.arg = plotting_data$sample_id,
    col = colors,
    border = NA,
    ylim = c(0, 1),
    las = 2,
    xlab = "Sample",
    ylab = "Proportion of sample reads",
    main = "Target composition by sample"
  )
  do.call(graphics::barplot, .qc_plot_arguments(defaults, dots))
  graphics::legend(
    "topright",
    legend = c("Target", "Non-target", "Unclassified"),
    fill = colors,
    border = NA,
    bty = "n"
  )
}

.qc_target_kingdom_plot_data <- function(by_kingdom, top_n) {
  known_input <- by_kingdom[
    by_kingdom$target_status %in% c("target", "non_target"),
    ,
    drop = FALSE
  ]
  known <- data.frame(
    kingdom_label = character(),
    target_status = character(),
    reads = numeric(),
    stringsAsFactors = FALSE
  )
  if (nrow(known_input) > 0) {
    keys <- unique(paste(
      known_input$normalized_kingdom,
      known_input$target_status,
      sep = "\r"
    ))
    known <- do.call(rbind, lapply(keys, function(key) {
      selected <- paste(
        known_input$normalized_kingdom,
        known_input$target_status,
        sep = "\r"
      ) == key
      data.frame(
        kingdom_label = known_input$kingdom_label[which(selected)[[1]]],
        target_status = known_input$target_status[which(selected)[[1]]],
        reads = sum(known_input$reads[selected]),
        stringsAsFactors = FALSE
      )
    }))
    rownames(known) <- NULL
    known <- known[
      order(-known$reads, seq_len(nrow(known))),
      ,
      drop = FALSE
    ]
  }

  shown <- utils::head(known, min(top_n, nrow(known)))
  remaining <- if (nrow(known) > nrow(shown)) {
    known[-seq_len(nrow(shown)), , drop = FALSE]
  } else {
    known[0, , drop = FALSE]
  }
  other <- lapply(c("target", "non_target"), function(status) {
    selected <- remaining$target_status == status
    if (!any(selected)) {
      return(NULL)
    }
    data.frame(
      kingdom_label = paste0(
        "Other (",
        if (status == "target") "target" else "non-target",
        ")"
      ),
      target_status = status,
      reads = sum(remaining$reads[selected]),
      stringsAsFactors = FALSE
    )
  })
  other <- Filter(Negate(is.null), other)
  if (length(other) > 0) {
    shown <- rbind(shown, do.call(rbind, other))
  }

  unclassified_reads <- sum(
    by_kingdom$reads[by_kingdom$target_status == "unclassified"]
  )
  if (any(by_kingdom$target_status == "unclassified")) {
    shown <- rbind(
      shown,
      data.frame(
        kingdom_label = "Unclassified",
        target_status = "unclassified",
        reads = unclassified_reads,
        stringsAsFactors = FALSE
      )
    )
  }
  total_reads <- sum(by_kingdom$reads)
  shown$read_proportion <- .qc_safe_divide(
    shown$reads,
    rep(total_reads, nrow(shown))
  )
  rownames(shown) <- NULL
  shown
}

.qc_draw_target_kingdom_plot <- function(plotting_data, dots) {
  colors <- c(
    target = "#3B7A57",
    non_target = "#C85C5C",
    unclassified = "#8A8A8A"
  )
  defaults <- list(
    height = plotting_data$reads,
    names.arg = plotting_data$kingdom_label,
    col = unname(colors[plotting_data$target_status]),
    border = NA,
    las = 2,
    xlab = "Kingdom",
    ylab = "Total reads",
    main = "Target composition by kingdom"
  )
  do.call(graphics::barplot, .qc_plot_arguments(defaults, dots))
  present <- intersect(names(colors), unique(plotting_data$target_status))
  graphics::legend(
    "topright",
    legend = c(
      target = "Target",
      non_target = "Non-target",
      unclassified = "Unclassified"
    )[present],
    fill = colors[present],
    border = NA,
    bty = "n"
  )
}

.qc_plot_arguments <- function(defaults, dots) {
  if ("height" %in% names(dots)) {
    stop("`height` cannot be supplied through `...`.", call. = FALSE)
  }
  for (name in names(dots)) {
    defaults[[name]] <- dots[[name]]
  }
  defaults
}
