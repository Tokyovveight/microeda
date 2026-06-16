test_that("microeda_qc builds per-sample and per-feature tables from a matrix", {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0,  4, 0, 0,
      2,  3, 0, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)

  expect_s3_class(qc, "microeda_qc")
  expect_named(
    qc$per_sample,
    c("sample_id", "library_size", "zero_fraction",
      "n_features_detected", "n_features_above_prevalence")
  )
  expect_named(
    qc$per_feature,
    c("feature_id", "total_reads", "prevalence",
      "n_samples_detected", "above_prevalence_threshold")
  )
  expect_equal(nrow(qc$per_sample), 4)
  expect_equal(nrow(qc$per_feature), 4)
  expect_equal(qc$per_sample$sample_id, paste0("S", 1:4))
  expect_equal(qc$per_sample$library_size, c(15, 21, 4, 6))
  expect_equal(qc$per_sample$n_features_detected, c(2, 2, 1, 3))
  expect_equal(qc$per_feature$total_reads, c(32, 7, 1, 6))
  expect_true(all(c(
    "library_size_summary", "sparsity_summary", "qc_flags"
  ) %in% names(qc)))
  expect_equal(qc$library_size_summary$n_samples, 4)
  expect_equal(qc$library_size_summary$total_reads, 46)
  expect_equal(qc$library_size_summary$min, 4)
  expect_equal(qc$library_size_summary$median, 10.5)
  expect_equal(qc$library_size_summary$nonzero_min, 4)
  expect_equal(qc$library_size_summary$max_to_median_ratio, 2)
  expect_equal(qc$library_size_summary$max_to_min_nonzero_ratio, 5.25)
  expect_equal(qc$sparsity_summary$overall_zero_fraction, 0.5)
  expect_equal(qc$sparsity_summary$zero_library_samples, 0)
  expect_equal(qc$sparsity_summary$zero_abundance_features, 0)
  expect_equal(qc$sparsity_summary$median_sample_zero_fraction, 0.5)
  expect_equal(qc$sparsity_summary$median_feature_zero_fraction, 0.5)
  expect_equal(nrow(qc$qc_flags), 0)
})

test_that("microeda_qc flags zero-library samples", {
  counts <- matrix(
    c(
      0, 0, 0,
      1, 0, 2,
      0, 3, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:3))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)

  expect_equal(qc$library_size_summary$zero_library_samples, 1)
  expect_equal(qc$sparsity_summary$zero_library_samples, 1)
  expect_equal(qc$sparsity_summary$zero_library_sample_fraction, 1 / 3)
  expect_true("zero_library_samples" %in% qc$qc_flags$flag_id)
})

test_that("microeda_qc flags zero-abundance features", {
  counts <- matrix(
    c(
      1, 0, 0, 2,
      3, 0, 4, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(paste0("s", 1:2), paste0("f", 1:4))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)

  expect_equal(qc$sparsity_summary$zero_abundance_features, 1)
  expect_equal(qc$sparsity_summary$zero_abundance_feature_fraction, 0.25)
  expect_true("zero_abundance_features" %in% qc$qc_flags$flag_id)
})

test_that("microeda_qc handles all-zero libraries without Inf or NaN ratios", {
  counts <- matrix(
    0,
    nrow = 3,
    ncol = 2,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:2))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)
  ratios <- c(
    qc$library_size_summary$max_to_median_ratio,
    qc$library_size_summary$max_to_min_nonzero_ratio
  )

  expect_equal(qc$library_size_summary$total_reads, 0)
  expect_equal(qc$library_size_summary$zero_library_samples, 3)
  expect_true(is.na(qc$library_size_summary$nonzero_min))
  expect_true(all(is.na(ratios)))
  expect_false(any(is.infinite(ratios)))
  expect_false(any(is.nan(ratios)))
  expect_equal(qc$sparsity_summary$overall_zero_fraction, 1)
  expect_equal(qc$sparsity_summary$zero_abundance_features, 2)
  expect_true("high_sparsity" %in% qc$qc_flags$flag_id)
})

test_that("microeda_qc returns NULL for missing taxonomy and metadata", {
  counts <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)

  expect_null(qc$per_rank)
  expect_null(qc$metadata_completeness)
})

test_that("microeda_qc summarizes taxonomy per rank with missing markers", {
  counts <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  tax <- data.frame(
    Phylum = c("Firmicutes", "Bacteroidota"),
    Genus = c("Faecalibacterium", "unclassified"),
    row.names = c("f1", "f2")
  )

  qc <- microeda_qc(counts, taxonomy = tax, taxa_are_rows = FALSE)

  expect_equal(qc$per_rank$rank, c("Phylum", "Genus"))
  expect_equal(qc$per_rank$n_assigned, c(2, 1))
  expect_equal(qc$per_rank$n_unique, c(2, 1))
  expect_equal(qc$per_rank$missing_fraction, c(0, 0.5))
})

test_that("microeda_qc summarizes metadata completeness with group flag", {
  counts <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  meta <- data.frame(
    group = c("A", "B"),
    batch = c("x", NA),
    constant = c("c", "c"),
    row.names = c("s1", "s2")
  )

  qc <- microeda_qc(counts, metadata = meta, group = "group", taxa_are_rows = FALSE)

  expect_equal(qc$metadata_completeness$column, c("group", "batch", "constant"))
  expect_equal(qc$metadata_completeness$missing_fraction, c(0, 0.5, 0))
  expect_equal(qc$metadata_completeness$n_unique, c(2, 1, 1))
  expect_equal(qc$metadata_completeness$is_constant, c(FALSE, TRUE, TRUE))
  expect_equal(qc$metadata_completeness$is_group, c(TRUE, FALSE, FALSE))
})

test_that("microeda_qc handles phyloseq otu_table input", {
  skip_if_not_installed("phyloseq")
  set.seed(1)
  counts <- matrix(
    sample(0:10, 12, replace = TRUE),
    nrow = 3,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:4))
  )
  otu <- phyloseq::otu_table(counts, taxa_are_rows = FALSE)

  qc <- microeda_qc(otu)

  expect_s3_class(qc, "microeda_qc")
  expect_equal(nrow(qc$per_sample), 3)
  expect_equal(nrow(qc$per_feature), 4)
})

test_that("microeda_qc handles full phyloseq input", {
  skip_if_not_installed("phyloseq")
  counts <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:4))
  )
  otu <- phyloseq::otu_table(counts, taxa_are_rows = FALSE)
  meta <- phyloseq::sample_data(data.frame(
    group = c("A", "A", "B"),
    row.names = paste0("s", 1:3)
  ))
  ps <- phyloseq::phyloseq(otu, meta)

  qc <- microeda_qc(ps, group = "group")

  expect_s3_class(qc, "microeda_qc")
  expect_equal(nrow(qc$per_sample), 3)
  expect_equal(nrow(qc$per_feature), 4)
  expect_equal(qc$metadata_completeness$column, "group")
  expect_true(qc$metadata_completeness$is_group)
})

test_that("microeda_qc applies the min_prevalence threshold", {
  counts <- matrix(
    c(
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 0, 0, 0,
      1, 1, 1, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("s", 1:4), paste0("f", 1:4))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE, min_prevalence = 0.5)

  expect_equal(qc$per_feature$above_prevalence_threshold, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(qc$per_sample$n_features_above_prevalence, c(1, 1, 1, 1))
})

test_that("microeda_qc validates min_prevalence", {
  counts <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_error(
    microeda_qc(counts, taxa_are_rows = FALSE, min_prevalence = -0.1),
    "min_prevalence"
  )
  expect_error(
    microeda_qc(counts, taxa_are_rows = FALSE, min_prevalence = c(0.1, 0.2)),
    "min_prevalence"
  )
  expect_error(
    microeda_qc(counts, taxa_are_rows = FALSE, min_prevalence = NA_real_),
    "min_prevalence"
  )
})
