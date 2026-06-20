test_that("microeda_beta_ordination computes PCoA coordinates", {
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

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  ord <- microeda_beta_ordination(beta)
  coordinates <- as_beta_coordinates(ord)

  expect_s3_class(ord, "microeda_beta_ordination")
  expect_named(coordinates, c("sample_id", "Axis1", "Axis2"))
  expect_equal(nrow(coordinates), length(beta$sample_ids))
  expect_equal(coordinates$sample_id, beta$sample_ids)
  expect_equal(ord$method, "pcoa")
  expect_equal(ord$distance_method, beta$method)
  expect_equal(ord$dimensions, 2L)
  expect_equal(ord$sample_ids, beta$sample_ids)
  expect_type(ord$eigenvalues, "double")
  expect_type(ord$variance_explained, "double")
})

test_that("microeda_beta_ordination preserves Jaccard distance method", {
  counts <- matrix(
    c(
      1, 1, 0,
      1, 0, 1,
      0, 1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- c("A", "B", "C")

  beta <- microeda_beta(counts, taxa_are_rows = FALSE, method = "jaccard")
  ord <- microeda_beta_ordination(beta)

  expect_s3_class(ord, "microeda_beta_ordination")
  expect_equal(ord$distance_method, "jaccard")
})

test_that("microeda_beta_ordination preserves Hellinger distance method", {
  counts <- matrix(
    c(
      4, 0,
      0, 9,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- c("A", "B")

  beta <- microeda_beta(counts, taxa_are_rows = FALSE, method = "hellinger")
  ord <- microeda_beta_ordination(beta)

  expect_s3_class(ord, "microeda_beta_ordination")
  expect_equal(ord$distance_method, "hellinger")
})

test_that("microeda_beta_ordination carries group values into coordinates", {
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

  beta <- microeda_beta(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  ord <- microeda_beta_ordination(beta)
  coordinates <- as_beta_coordinates(ord)

  expect_named(coordinates, c("sample_id", "Axis1", "Axis2", "group"))
  expect_equal(coordinates$sample_id, beta$sample_ids)
  expect_equal(coordinates$group, metadata$group)
  expect_equal(ord$group, "group")
  expect_equal(unname(ord$group_values), metadata$group)
})

test_that("as_beta_coordinates extracts beta ordination coordinates", {
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

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  ord <- microeda_beta_ordination(beta)

  expect_identical(as_beta_coordinates(ord), ord$coordinates)
  expect_error(as_beta_coordinates(data.frame()), "microeda_beta_ordination")
})

test_that("microeda_beta_ordination validates input, method, and dimensions", {
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
  beta <- microeda_beta(counts, taxa_are_rows = FALSE)

  expect_error(microeda_beta_ordination(data.frame()), "microeda_beta")
  expect_error(microeda_beta_ordination(beta, method = "nmds"), "pcoa")
  expect_error(
    microeda_beta_ordination(beta, method = NA_character_),
    "pcoa"
  )
  expect_error(microeda_beta_ordination(beta, dimensions = NA_real_), "dimensions")
  expect_error(microeda_beta_ordination(beta, dimensions = 1.5), "dimensions")
  expect_error(microeda_beta_ordination(beta, dimensions = 0), "dimensions")
  expect_error(microeda_beta_ordination(beta, dimensions = 3), "dimensions")
})

test_that("microeda_beta_ordination default dimensions require enough samples", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  beta <- microeda_beta(counts, taxa_are_rows = FALSE)

  expect_error(microeda_beta_ordination(beta), "dimensions")
  expect_s3_class(microeda_beta_ordination(beta, dimensions = 1), "microeda_beta_ordination")
})

test_that("microeda_beta_ordination prints compact summaries", {
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

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  ord <- microeda_beta_ordination(beta)

  output <- capture.output(result <- withVisible(print(ord)))

  expect_false(result$visible)
  expect_identical(result$value, ord)
  expect_true(any(grepl("microeda_beta_ordination", output)))
  expect_true(any(grepl("Method: +pcoa", output)))
  expect_true(any(grepl("Distance method: +bray", output)))
  expect_true(any(grepl("Dimensions: +2", output)))
})
