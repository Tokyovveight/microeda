test_that("microeda_beta_compare_ordination computes default ordinations", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp)

  expect_s3_class(ord_cmp, "microeda_beta_compare_ordination")
  expect_equal(ord_cmp$methods, c("bray", "jaccard", "hellinger"))
  expect_named(ord_cmp$results, ord_cmp$methods)
  expect_true(all(vapply(
    ord_cmp$results,
    inherits,
    logical(1),
    "microeda_beta_ordination"
  )))
  expect_equal(ord_cmp$ordination_method, "pcoa")
  expect_equal(ord_cmp$dimensions, 2L)
  expect_equal(ord_cmp$sample_ids, beta_cmp$sample_ids)
})

test_that("microeda_beta_compare_ordination preserves method order", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp, dimensions = 1)

  expect_equal(ord_cmp$methods, c("hellinger", "bray"))
  expect_named(ord_cmp$results, c("hellinger", "bray"))
  expect_equal(ord_cmp$dimensions, 1L)
  expect_equal(
    unname(vapply(ord_cmp$results, `[[`, character(1), "distance_method")),
    ord_cmp$methods
  )
})

test_that("microeda_beta_compare_ordination carries group metadata", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp)

  expect_equal(ord_cmp$group, "group")
  expect_equal(unname(ord_cmp$group_values), metadata$group)
  expect_equal(names(ord_cmp$group_values), rownames(counts))
})

test_that("microeda_beta_compare_ordination validates input, method, and dimensions", {
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

  expect_error(
    microeda_beta_compare_ordination(data.frame()),
    "microeda_beta_compare"
  )
  expect_error(
    microeda_beta_compare_ordination(beta_cmp, method = "nmds"),
    "pcoa"
  )
  expect_error(
    microeda_beta_compare_ordination(beta_cmp, dimensions = 3),
    "dimensions"
  )
})

test_that("as_beta_compare_coordinates returns stable ungrouped coordinates", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp)
  coordinates <- as_beta_compare_coordinates(ord_cmp)

  expect_s3_class(coordinates, "data.frame")
  expect_named(coordinates, c("method", "sample_id", "Axis1", "Axis2"))
  expect_equal(nrow(coordinates), length(ord_cmp$methods) * length(ord_cmp$sample_ids))
  expect_equal(coordinates$method, rep(ord_cmp$methods, each = length(ord_cmp$sample_ids)))
  expect_equal(
    coordinates$sample_id,
    rep(ord_cmp$sample_ids, times = length(ord_cmp$methods))
  )
})

test_that("as_beta_compare_coordinates returns grouped coordinates", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp)
  coordinates <- as_beta_compare_coordinates(ord_cmp)

  expect_named(coordinates, c("method", "sample_id", "Axis1", "Axis2", "group"))
  expect_equal(
    coordinates$group,
    rep(metadata$group, times = length(ord_cmp$methods))
  )
})

test_that("as_beta_compare_coordinates supports additional axes", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 1, 3,
      3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3", "S4")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE, methods = "bray")
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp, dimensions = 3)
  coordinates <- as_beta_compare_coordinates(ord_cmp)

  expect_named(coordinates, c("method", "sample_id", "Axis1", "Axis2", "Axis3"))
  expect_equal(nrow(coordinates), length(ord_cmp$sample_ids))
})

test_that("as_beta_compare_coordinates validates input", {
  expect_error(
    as_beta_compare_coordinates(data.frame()),
    "microeda_beta_compare_ordination"
  )
})

test_that("microeda_beta_compare_ordination prints compact summaries", {
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
  ord_cmp <- microeda_beta_compare_ordination(beta_cmp)

  output <- capture.output(result <- withVisible(print(ord_cmp)))

  expect_false(result$visible)
  expect_identical(result$value, ord_cmp)
  expect_true(any(grepl("microeda_beta_compare_ordination", output)))
  expect_true(any(grepl("Methods: +bray, jaccard, hellinger", output)))
  expect_true(any(grepl("Dimensions: +2", output)))
  expect_true(any(grepl("Samples: +3", output)))
  expect_true(any(grepl("Group: +group", output)))
  expect_true(any(grepl("as_beta_compare_coordinates", output)))
})
