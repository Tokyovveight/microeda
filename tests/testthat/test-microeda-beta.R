test_that("microeda_beta computes Bray-Curtis distances from a matrix", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  distance_matrix <- as.matrix(beta$distance)

  expect_s3_class(beta, "microeda_beta")
  expect_s3_class(beta$distance, "dist")
  expect_equal(beta$method, "bray")
  expect_equal(beta$sample_ids, rownames(counts))
  expect_equal(attr(beta$distance, "Labels"), rownames(counts))
  expect_equal(distance_matrix["S1", "S2"], 2 / 6)
  expect_equal(distance_matrix["S1", "S3"], 1)
})

test_that("microeda_beta accepts explicit Bray-Curtis method", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )

  beta <- microeda_beta(counts, taxa_are_rows = FALSE, method = "bray")

  expect_s3_class(beta, "microeda_beta")
  expect_equal(beta$method, "bray")
})

test_that("as_beta_dist extracts stored beta distance objects", {
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
  distance <- as_beta_dist(beta)

  expect_s3_class(distance, "dist")
  expect_equal(attr(distance, "Labels"), beta$sample_ids)
  expect_equal(distance, beta$distance)
  expect_error(as_beta_dist(data.frame()), "microeda_beta")
})

test_that("as_beta_matrix extracts square beta distance matrices", {
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
  distance_matrix <- as_beta_matrix(beta)

  expect_true(is.matrix(distance_matrix))
  expect_equal(dim(distance_matrix), c(3L, 3L))
  expect_equal(rownames(distance_matrix), beta$sample_ids)
  expect_equal(colnames(distance_matrix), beta$sample_ids)
  expect_equal(distance_matrix, as.matrix(as_beta_dist(beta)))
  expect_equal(unname(diag(distance_matrix)), rep(0, length(beta$sample_ids)))
  expect_error(as_beta_matrix(data.frame()), "microeda_beta")
})

test_that("microeda_beta prints compact summaries", {
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

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  grouped_beta <- microeda_beta(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  output <- capture.output(result <- withVisible(print(beta)))
  grouped_output <- capture.output(grouped_result <- withVisible(print(grouped_beta)))

  expect_false(result$visible)
  expect_identical(result$value, beta)
  expect_true(any(grepl("microeda_beta", output)))
  expect_true(any(grepl("Method: +bray", output)))
  expect_true(any(grepl("Samples: +3", output)))
  expect_true(any(grepl("Group: +<none>", output)))
  expect_false(grouped_result$visible)
  expect_identical(grouped_result$value, grouped_beta)
  expect_true(any(grepl("Group: +group", grouped_output)))
})

test_that("microeda_beta assigns zero distance to all-zero sample pairs", {
  counts <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )

  beta <- microeda_beta(counts, taxa_are_rows = FALSE)

  expect_equal(as.matrix(beta$distance)["S1", "S2"], 0)
})

test_that("microeda_beta stores optional group metadata", {
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

  expect_equal(beta$group, "group")
  expect_equal(unname(beta$group_values), metadata$group)
  expect_equal(names(beta$group_values), rownames(counts))
})

test_that("microeda_beta validates method, group, and counts", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  metadata <- data.frame(
    group = c("A", "B"),
    row.names = rownames(counts)
  )

  expect_error(
    microeda_beta(counts, taxa_are_rows = FALSE, method = "jaccard"),
    "method"
  )
  expect_error(
    microeda_beta(counts, taxa_are_rows = FALSE, method = NA_character_),
    "method"
  )
  expect_error(
    microeda_beta(counts, taxa_are_rows = FALSE, group = "group"),
    "metadata"
  )
  expect_error(
    microeda_beta(
      counts,
      metadata = metadata,
      group = "unknown",
      taxa_are_rows = FALSE
    ),
    "group"
  )
  expect_error(
    microeda_beta(
      counts,
      metadata = metadata,
      group = c("group", "batch"),
      taxa_are_rows = FALSE
    ),
    "group"
  )

  negative_counts <- counts
  negative_counts[1, 1] <- -1
  infinite_counts <- counts
  infinite_counts[1, 1] <- Inf

  expect_error(
    microeda_beta(negative_counts, taxa_are_rows = FALSE),
    "negative"
  )
  expect_error(
    microeda_beta(infinite_counts, taxa_are_rows = FALSE),
    "finite"
  )
})

test_that("microeda_beta handles phyloseq otu_table input", {
  skip_if_not_installed("phyloseq")

  counts <- matrix(
    c(1, 2, 0, 2, 1, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), paste0("ASV", seq_len(3)))
  )
  otu <- phyloseq::otu_table(counts, taxa_are_rows = FALSE)

  beta <- microeda_beta(otu)

  expect_s3_class(beta, "microeda_beta")
  expect_s3_class(beta$distance, "dist")
  expect_equal(beta$sample_ids, rownames(counts))
})
