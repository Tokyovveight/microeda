test_that("microeda_beta_compare computes default beta methods", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)

  expect_s3_class(beta_cmp, "microeda_beta_compare")
  expect_equal(beta_cmp$methods, c("bray", "jaccard", "hellinger"))
  expect_named(beta_cmp$results, beta_cmp$methods)
  expect_true(all(vapply(
    beta_cmp$results,
    inherits,
    logical(1),
    "microeda_beta"
  )))
  expect_true(all(vapply(
    beta_cmp$results,
    function(result) identical(result$sample_ids, beta_cmp$sample_ids),
    logical(1)
  )))
})

test_that("microeda_beta_compare preserves requested method order", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = c("hellinger", "bray")
  )

  expect_equal(beta_cmp$methods, c("hellinger", "bray"))
  expect_named(beta_cmp$results, c("hellinger", "bray"))
  expect_equal(
    unname(vapply(beta_cmp$results, `[[`, character(1), "method")),
    beta_cmp$methods
  )
})

test_that("microeda_beta_compare carries group metadata", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "A", "B"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  expect_equal(beta_cmp$group, "group")
  expect_equal(unname(beta_cmp$group_values), metadata$group)
  expect_equal(names(beta_cmp$group_values), rownames(counts))
})

test_that("microeda_beta_compare validates methods", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )

  expect_error(
    microeda_beta_compare(
      counts,
      taxa_are_rows = FALSE,
      methods = c("bray", "bray")
    ),
    "duplicate"
  )
  expect_error(
    microeda_beta_compare(
      counts,
      taxa_are_rows = FALSE,
      methods = c("bray", "unknown")
    ),
    "bray.*jaccard.*hellinger"
  )
})

test_that("as_beta_compare_summary returns stable method summaries", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  summary <- as_beta_compare_summary(beta_cmp)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c(
      "method",
      "n_samples",
      "n_pairs",
      "min_distance",
      "median_distance",
      "max_distance"
    )
  )
  expect_equal(summary$method, beta_cmp$methods)
  expect_equal(nrow(summary), length(beta_cmp$methods))
  expect_equal(summary$n_samples, rep(3, length(beta_cmp$methods)))
  expect_equal(summary$n_pairs, rep(3, length(beta_cmp$methods)))
  expect_error(as_beta_compare_summary(data.frame()), "microeda_beta_compare")
})

test_that("as_beta_compare_distances returns stable ungrouped pairwise distances", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = c("hellinger", "bray")
  )
  distances <- as_beta_compare_distances(beta_cmp)

  expect_s3_class(distances, "data.frame")
  expect_named(distances, c("method", "sample_1", "sample_2", "distance"))
  expect_equal(nrow(distances), length(beta_cmp$methods) * 3)
  expect_equal(
    distances$method,
    rep(beta_cmp$methods, each = 3)
  )
  expect_equal(
    distances$sample_1[seq_len(3)],
    c("S1", "S1", "S2")
  )
  expect_equal(
    distances$sample_2[seq_len(3)],
    c("S2", "S3", "S3")
  )

  bray_rows <- distances$method == "bray"
  expect_equal(
    distances$distance[bray_rows],
    as.numeric(as_beta_dist(beta_cmp$results$bray))
  )
  expect_error(as_beta_compare_distances(data.frame()), "microeda_beta_compare")
})

test_that("as_beta_compare_distances includes grouped pairwise labels", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "B", "A"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    methods = "bray"
  )
  distances <- as_beta_compare_distances(beta_cmp)

  expect_named(
    distances,
    c(
      "method",
      "sample_1",
      "sample_2",
      "group_1",
      "group_2",
      "comparison",
      "distance"
    )
  )
  expect_equal(distances$group_1, c("A", "A", "B"))
  expect_equal(distances$group_2, c("B", "A", "A"))
  expect_equal(distances$comparison, c("between", "within", "between"))
})

test_that("as_beta_compare_distances marks missing group comparisons as missing", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", NA, "A"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    methods = "bray"
  )
  distances <- as_beta_compare_distances(beta_cmp)

  expect_equal(distances$comparison, c(NA_character_, "within", NA_character_))
})

test_that("as_beta_compare_distance_correlations returns stable correlations", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  correlations <- as_beta_compare_distance_correlations(beta_cmp)

  expect_s3_class(correlations, "data.frame")
  expect_named(
    correlations,
    c(
      "method_1",
      "method_2",
      "n_pairs",
      "correlation",
      "correlation_method"
    )
  )
  expect_equal(nrow(correlations), 3L)
  expect_equal(correlations$method_1, c("bray", "bray", "jaccard"))
  expect_equal(correlations$method_2, c("jaccard", "hellinger", "hellinger"))
  expect_equal(correlations$n_pairs, rep(3L, 3))
  expect_equal(correlations$correlation_method, rep("spearman", 3))

  manual <- stats::cor(
    as.numeric(as_beta_dist(beta_cmp$results$bray)),
    as.numeric(as_beta_dist(beta_cmp$results$jaccard)),
    method = "spearman"
  )
  expect_equal(correlations$correlation[1], unname(manual))
})

test_that("as_beta_compare_distance_correlations accepts supported methods", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  for (correlation_method in c("pearson", "spearman", "kendall")) {
    correlations <- as_beta_compare_distance_correlations(
      beta_cmp,
      correlation_method = correlation_method
    )
    expect_equal(
      correlations$correlation_method,
      rep(correlation_method, nrow(correlations))
    )
  }

  pearson <- as_beta_compare_distance_correlations(
    beta_cmp,
    correlation_method = "pearson"
  )
  manual <- stats::cor(
    as.numeric(as_beta_dist(beta_cmp$results$bray)),
    as.numeric(as_beta_dist(beta_cmp$results$jaccard)),
    method = "pearson"
  )
  expect_equal(pearson$correlation[1], unname(manual))
})

test_that("as_beta_compare_distance_correlations preserves subset method order", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = c("hellinger", "bray")
  )
  correlations <- as_beta_compare_distance_correlations(beta_cmp)

  expect_equal(nrow(correlations), 1L)
  expect_equal(correlations$method_1, "hellinger")
  expect_equal(correlations$method_2, "bray")
})

test_that("as_beta_compare_distance_correlations handles one method", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = "bray"
  )
  correlations <- as_beta_compare_distance_correlations(beta_cmp)

  expect_s3_class(correlations, "data.frame")
  expect_named(
    correlations,
    c(
      "method_1",
      "method_2",
      "n_pairs",
      "correlation",
      "correlation_method"
    )
  )
  expect_equal(nrow(correlations), 0L)
})

test_that("as_beta_compare_distance_correlations validates input", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)

  expect_error(
    as_beta_compare_distance_correlations(data.frame()),
    "microeda_beta_compare"
  )
  expect_error(
    as_beta_compare_distance_correlations(
      beta_cmp,
      correlation_method = "unknown"
    ),
    "pearson.*spearman.*kendall"
  )
})

test_that("as_beta_compare_distance_correlations handles undefined correlations", {
  counts <- matrix(
    0,
    nrow = 3,
    ncol = 2,
    dimnames = list(c("S1", "S2", "S3"), c("ASV1", "ASV2"))
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = c("bray", "jaccard")
  )
  expect_silent(
    correlations <- as_beta_compare_distance_correlations(beta_cmp)
  )

  expect_equal(correlations$n_pairs, 3L)
  expect_true(is.na(correlations$correlation))
})

test_that("as_beta_compare_group_summary returns stable grouped summaries", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "B", "A"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    methods = c("hellinger", "bray")
  )
  summary <- as_beta_compare_group_summary(beta_cmp)
  distances <- as_beta_compare_distances(beta_cmp)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c(
      "method",
      "comparison",
      "n_pairs",
      "min_distance",
      "median_distance",
      "max_distance"
    )
  )
  expect_equal(
    summary$method,
    rep(c("hellinger", "bray"), each = 2)
  )
  expect_equal(
    summary$comparison,
    rep(c("within", "between"), times = 2)
  )
  expect_equal(summary$n_pairs, c(1L, 2L, 1L, 2L))

  manual <- distances[
    distances$method == "bray" & distances$comparison == "between",
    ,
    drop = FALSE
  ]
  bray_between <- summary[
    summary$method == "bray" & summary$comparison == "between",
    ,
    drop = FALSE
  ]

  expect_equal(bray_between$n_pairs, nrow(manual))
  expect_equal(bray_between$min_distance, min(manual$distance))
  expect_equal(bray_between$median_distance, stats::median(manual$distance))
  expect_equal(bray_between$max_distance, max(manual$distance))
})

test_that("as_beta_compare_group_summary validates grouped input", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)

  expect_error(as_beta_compare_group_summary(data.frame()), "microeda_beta_compare")
  expect_error(as_beta_compare_group_summary(beta_cmp), "group metadata")
})

test_that("as_beta_compare_group_summary returns stable empty summaries", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c(NA_character_, NA_character_, NA_character_),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    methods = "bray"
  )
  summary <- as_beta_compare_group_summary(beta_cmp)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c(
      "method",
      "comparison",
      "n_pairs",
      "min_distance",
      "median_distance",
      "max_distance"
    )
  )
  expect_equal(nrow(summary), 0L)
})

test_that("microeda_beta_compare_report returns compact grouped reports", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "B", "A"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  output <- capture.output(report <- microeda_beta_compare_report(beta_cmp))

  expect_identical(output, character())
  expect_type(report, "character")
  expect_length(report, 1)
  expect_match(report, "Beta diversity method comparison", fixed = TRUE)
  expect_match(report, "Methods: bray, jaccard, hellinger", fixed = TRUE)
  expect_match(report, "Samples: 3", fixed = TRUE)
  expect_match(report, "Group: group", fixed = TRUE)
  expect_match(report, "Method-level distance summary", fixed = TRUE)
  expect_match(report, "Distance-method correlations", fixed = TRUE)
  expect_match(report, "spearman", fixed = TRUE)
  expect_match(report, "Group-level distance summary", fixed = TRUE)
  expect_match(report, "Bray-Curtis: abundance-sensitive distance.", fixed = TRUE)
  expect_match(report, "Jaccard: binary presence/absence distance.", fixed = TRUE)
  expect_match(report, "Hellinger: square-root relative abundance", fixed = TRUE)
  expect_match(
    report,
    "Distance-method correlations are descriptive only",
    fixed = TRUE
  )
  expect_match(report, "PERMANOVA is not implemented", fixed = TRUE)
  expect_match(report, "Formal method recommendation is not implemented yet.", fixed = TRUE)
})

test_that("microeda_beta_compare_report handles ungrouped comparisons", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  report <- microeda_beta_compare_report(beta_cmp)

  expect_match(report, "Group: <none>", fixed = TRUE)
  expect_match(report, "Distance-method correlations", fixed = TRUE)
  expect_match(
    report,
    "Group-level distance summary unavailable: no group metadata supplied.",
    fixed = TRUE
  )
})

test_that("microeda_beta_compare_report validates input", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)

  expect_error(
    microeda_beta_compare_report(data.frame()),
    "microeda_beta_compare"
  )
  expect_error(
    microeda_beta_compare_report(beta_cmp, correlation_method = "unknown"),
    "pearson.*spearman.*kendall"
  )
})

test_that("microeda_beta_compare_report handles one-method correlations", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(
    counts,
    taxa_are_rows = FALSE,
    methods = "bray"
  )
  report <- microeda_beta_compare_report(beta_cmp)

  expect_match(report, "Distance-method correlations", fixed = TRUE)
  expect_match(
    report,
    "Distance-method correlations unavailable: fewer than two beta methods supplied.",
    fixed = TRUE
  )
})

test_that("microeda_beta_compare prints compact summaries", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "A", "B"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  output <- capture.output(result <- withVisible(print(beta_cmp)))

  expect_false(result$visible)
  expect_identical(result$value, beta_cmp)
  expect_true(any(grepl("microeda_beta_compare", output)))
  expect_true(any(grepl("Methods: +bray, jaccard, hellinger", output)))
  expect_true(any(grepl("Samples: +3", output)))
  expect_true(any(grepl("Group: +group", output)))
  expect_true(any(grepl("as_beta_compare_summary", output)))
})

test_that("beta_compare_rule_context returns stable internal context", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  context <- microeda:::beta_compare_rule_context(beta_cmp)

  expect_type(context, "list")
  expect_named(context, c("summary", "methods", "caveats"))
  expect_s3_class(context$summary, "data.frame")
  expect_named(
    context$summary,
    c(
      "n_methods",
      "methods",
      "n_samples",
      "group",
      "has_group",
      "has_group_summary",
      "has_distance_correlations"
    )
  )
  expect_equal(context$summary$n_methods, length(beta_cmp$methods))
  expect_equal(context$summary$methods, paste(beta_cmp$methods, collapse = ", "))
  expect_equal(context$summary$n_samples, length(beta_cmp$sample_ids))
  expect_false(context$summary$has_group)
  expect_false(context$summary$has_group_summary)
  expect_true(context$summary$has_distance_correlations)
  expect_equal(context$methods$method, beta_cmp$methods)
  expect_equal(context$methods$n_samples, rep(3, length(beta_cmp$methods)))
})

test_that("beta_compare_rule_context documents shared caveat vocabulary", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  context <- microeda:::beta_compare_rule_context(beta_cmp)

  expect_s3_class(context$summary, "data.frame")
  expect_equal(nrow(context$summary), 1L)
  expect_s3_class(context$caveats, "data.frame")
  expect_true(all(c("context_id", "topic") %in% names(context$caveats)))
  expect_true(all(c("method", "severity", "message") %in% names(context$caveats)))
  expect_type(context$caveats$context_id, "character")
  expect_type(context$caveats$topic, "character")
  expect_type(context$caveats$method, "character")
  expect_type(context$caveats$severity, "character")
  expect_type(context$caveats$message, "character")
  expect_true(all(context$caveats$severity == "info"))
  expect_false(any(grepl("best method|rank", context$caveats$message, ignore.case = TRUE)))
})

test_that("beta_compare_rule_context detects grouped comparisons", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "B", "A"),
    row.names = rownames(counts)
  )

  beta_cmp <- microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  context <- microeda:::beta_compare_rule_context(beta_cmp)

  expect_true(context$summary$has_group)
  expect_equal(context$summary$group, "group")
  expect_true(context$summary$has_group_summary)
})

test_that("beta_compare_rule_context includes non-ranking caveats", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)
  caveats <- microeda:::beta_compare_rule_context(beta_cmp)$caveats

  expect_named(
    caveats,
    c("context_id", "topic", "method", "severity", "message")
  )
  expect_true("jaccard_incidence" %in% caveats$context_id)
  expect_true("hellinger_not_log_ratio" %in% caveats$context_id)
  expect_true("distance_correlations_descriptive" %in% caveats$context_id)
  expect_true("pcoa_axes_method_specific" %in% caveats$context_id)
  expect_true(any(grepl("incidence-based", caveats$message, fixed = TRUE)))
  expect_true(any(grepl("not a log-ratio", caveats$message, fixed = TRUE)))
  expect_true(any(grepl("do not validate a method", caveats$message, fixed = TRUE)))
  expect_true(any(grepl("method-specific", caveats$message, fixed = TRUE)))
  expect_false(any(grepl("best method", caveats$message, ignore.case = TRUE)))
})

test_that("beta_compare_rule_context validates input", {
  expect_error(
    microeda:::beta_compare_rule_context(data.frame()),
    "microeda_beta_compare"
  )
})
