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

test_that("microeda_beta computes binary Jaccard distances from a matrix", {
  counts <- matrix(
    c(
      1, 1, 0,
      1, 0, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("A", "B", "C"))
  )

  beta <- microeda_beta(counts, taxa_are_rows = FALSE, method = "jaccard")
  distance_matrix <- as_beta_matrix(beta)

  expect_s3_class(beta, "microeda_beta")
  expect_equal(beta$method, "jaccard")
  expect_s3_class(as_beta_dist(beta), "dist")
  expect_equal(beta$sample_ids, rownames(counts))
  expect_equal(attr(as_beta_dist(beta), "Labels"), rownames(counts))
  expect_equal(distance_matrix["S1", "S2"], 2 / 3)
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

test_that("as_beta_samples extracts sample IDs and optional groups", {
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

  samples <- as_beta_samples(beta)
  grouped_samples <- as_beta_samples(grouped_beta)

  expect_s3_class(samples, "data.frame")
  expect_named(samples, "sample_id")
  expect_equal(samples$sample_id, beta$sample_ids)
  expect_s3_class(grouped_samples, "data.frame")
  expect_named(grouped_samples, c("sample_id", "group"))
  expect_equal(grouped_samples$sample_id, grouped_beta$sample_ids)
  expect_equal(grouped_samples$group, metadata$group)
  expect_error(as_beta_samples(data.frame()), "microeda_beta")
})

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

test_that("microeda_beta_plot draws heatmaps invisibly", {
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
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  result <- withVisible(microeda_beta_plot(beta))
  explicit_type <- withVisible(microeda_beta_plot(beta, type = "heatmap"))
  custom_labels <- withVisible(microeda_beta_plot(
    beta,
    main = "Custom beta heatmap",
    xlab = "Custom x",
    ylab = "Custom y"
  ))

  expect_false(result$visible)
  expect_null(result$value)
  expect_false(explicit_type$visible)
  expect_null(explicit_type$value)
  expect_false(custom_labels$visible)
  expect_null(custom_labels$value)
})

test_that("microeda_beta_plot validates input and type", {
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), c("ASV1", "ASV2"))
  )
  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  supported_type_pattern <- "heatmap"

  expect_error(microeda_beta_plot(data.frame()), "microeda_beta")
  expect_error(
    microeda_beta_plot(beta, type = "unknown"),
    supported_type_pattern
  )
  expect_error(
    microeda_beta_plot(beta, type = NA_character_),
    supported_type_pattern
  )
  expect_error(
    microeda_beta_plot(beta, type = c("heatmap", "heatmap")),
    supported_type_pattern
  )
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
  jaccard_beta <- microeda_beta(counts, taxa_are_rows = FALSE, method = "jaccard")

  expect_equal(as.matrix(beta$distance)["S1", "S2"], 0)
  expect_equal(as_beta_matrix(jaccard_beta)["S1", "S2"], 0)
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
    microeda_beta(counts, taxa_are_rows = FALSE, method = "unknown"),
    "bray.*jaccard"
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
