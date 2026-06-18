test_that("microeda_check summarizes matrix input", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )

  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))

  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    batch = c("x", "y", "x", "y"),
    row.names = rownames(counts)
  )

  report <- microeda_check(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  expect_s3_class(report, "microeda_report")
  expect_equal(report$diagnostics$n_samples, 4)
  expect_equal(report$diagnostics$n_features, 4)
  expect_equal(report$diagnostics$reads$total_reads, 46)
  expect_equal(report$diagnostics$reads$min_sample$sample_id, "S3")
  expect_equal(report$diagnostics$reads$max_sample$sample_id, "S2")
  expect_true("compositional_caveat" %in% report$recommendations$rule_id)
  expect_true("low_group_n" %in% report$recommendations$rule_id)
})

test_that("relative input is detected", {
  counts <- matrix(
    c(
      0.5, 0.5,
      0.2, 0.8
    ),
    nrow = 2,
    byrow = TRUE
  )

  report <- microeda_check(counts, taxa_are_rows = FALSE)

  expect_true(report$diagnostics$count_type$looks_relative)
  expect_true("relative_input" %in% report$recommendations$rule_id)
})

test_that("taxonomy rank summaries count assigned and unique taxa", {
  counts <- matrix(
    c(
      10, 0, 0,
      0, 5, 1
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2")
  colnames(counts) <- c("ASV1", "ASV2", "ASV3")

  taxonomy <- data.frame(
    Phylum = c("Ascomycota", "Ascomycota", "Basidiomycota"),
    Genus = c("Fusarium", "unclassified", "Candida"),
    row.names = colnames(counts)
  )

  report <- microeda_check(
    counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE
  )

  expect_equal(report$diagnostics$taxonomy$unique_taxa_by_rank[["Phylum"]], 2)
  expect_equal(report$diagnostics$taxonomy$unique_taxa_by_rank[["Genus"]], 2)
  expect_equal(report$diagnostics$taxonomy$assigned_features_by_rank[["Genus"]], 2)
})

test_that("microeda_alpha computes classic and Hill alpha indices", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))

  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = rownames(counts)
  )

  alpha <- microeda_alpha(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  alpha_table <- as_alpha_table(alpha)

  expect_s3_class(alpha, "microeda_alpha")
  expect_equal(alpha_table$observed[1], 2)
  expect_equal(alpha_table$hill_q0[1], alpha_table$observed[1])
  expect_equal(alpha_table$hill_q1[1], exp(alpha_table$shannon[1]))
  expect_equal(alpha_table$hill_q2[1], alpha_table$inverse_simpson[1])
  expect_true(nrow(as_alpha_summary(alpha)) > 0)
})

test_that("microeda_alpha_plot draws alpha metric barplots invisibly", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))

  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = rownames(counts)
  )

  alpha <- microeda_alpha(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  result <- withVisible(microeda_alpha_plot(alpha))
  shannon <- withVisible(microeda_alpha_plot(alpha, metric = "shannon"))

  expect_false(result$visible)
  expect_true(is.numeric(result$value))
  expect_length(result$value, nrow(as_alpha_table(alpha)))
  expect_false(shannon$visible)
  expect_true(is.numeric(shannon$value))
  expect_length(shannon$value, nrow(as_alpha_table(alpha)))
})

test_that("microeda_alpha_plot validates input and metric", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))

  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = rownames(counts)
  )

  alpha <- microeda_alpha(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  expect_error(microeda_alpha_plot(data.frame()), "microeda_alpha")
  expect_error(microeda_alpha_plot(alpha, metric = NA_character_), "metric")
  expect_error(
    microeda_alpha_plot(alpha, metric = c("observed", "shannon")),
    "metric"
  )
  expect_error(microeda_alpha_plot(alpha, metric = "unknown"), "metric")
  expect_error(
    microeda_alpha_plot(alpha, metric = "group"),
    "numeric alpha metric"
  )
})

test_that("microeda_alpha_compare runs omnibus and pairwise tests", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1,
      10, 8, 7, 6,
      1, 0, 0, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(6))
  colnames(counts) <- paste0("ASV", seq_len(4))

  metadata <- data.frame(
    group = c("A", "A", "B", "B", "C", "C"),
    row.names = rownames(counts)
  )

  alpha <- microeda_alpha(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  comparison <- microeda_alpha_compare(
    alpha,
    indices = c("observed", "hill_q1")
  )

  expect_s3_class(comparison, "microeda_alpha_compare")
  expect_equal(as_alpha_tests(comparison)$index, c("observed", "hill_q1"))
  expect_true(all(c("p_value", "p_value_adjusted") %in% names(as_alpha_tests(comparison))))
  expect_true(nrow(as_alpha_pairwise(comparison)) > 0)
})
