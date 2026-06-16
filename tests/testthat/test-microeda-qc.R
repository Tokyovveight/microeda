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
    "library_size_summary", "sparsity_summary", "prevalence_summary",
    "qc_flags"
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
  expect_equal(qc$prevalence_summary$n_features, 4)
  expect_equal(qc$prevalence_summary$n_samples, 4)
  expect_equal(qc$prevalence_summary$min_prevalence_threshold, 0.05)
  expect_equal(qc$prevalence_summary$n_features_above_threshold, 4)
  expect_equal(qc$prevalence_summary$n_features_below_threshold, 0)
  expect_equal(qc$prevalence_summary$fraction_features_above_threshold, 1)
  expect_equal(qc$prevalence_summary$fraction_features_below_threshold, 0)
  expect_equal(qc$prevalence_summary$min_prevalence, 0.25)
  expect_equal(qc$prevalence_summary$median_prevalence, 0.5)
  expect_equal(qc$prevalence_summary$mean_prevalence, 0.5)
  expect_equal(qc$prevalence_summary$max_prevalence, 0.75)
  expect_equal(qc$prevalence_summary$n_features_detected_in_all_samples, 0)
  expect_equal(qc$prevalence_summary$n_features_detected_in_one_sample, 1)
  expect_equal(
    qc$prevalence_summary$fraction_features_detected_in_one_sample,
    0.25
  )
  expect_equal(nrow(qc$qc_flags), 0)
})

test_that("microeda_qc records parameters and call", {
  counts <- matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  metadata <- data.frame(
    group = c("A", "B"),
    row.names = c("s1", "s2")
  )

  qc <- microeda_qc(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    min_prevalence = 0.25
  )

  expect_true("params" %in% names(qc))
  expect_equal(qc$params$group, "group")
  expect_equal(qc$params$min_prevalence, 0.25)
  expect_true("call" %in% names(qc))
  expect_true(is.call(qc$call))
})

test_that("microeda_qc returns stable human-readable observations", {
  counts <- matrix(
    c(
      1, 0,
      0, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE)

  expect_named(
    qc$qc_observations,
    c("observation_id", "category", "severity", "message")
  )
  expect_equal(
    qc$qc_observations$observation_id[1:6],
    c(
      "input_dimensions",
      "total_reads",
      "overall_zero_fraction",
      "features_above_prevalence",
      "metadata_absent",
      "taxonomy_absent"
    )
  )
  expect_true(all(c("input", "library_size", "sparsity", "prevalence") %in%
    qc$qc_observations$category))
  expect_true(all(qc$qc_observations$severity %in% c("info", "warning")))
  expect_equal(sum(duplicated(qc$qc_observations$observation_id)), 0)
  expect_true(any(grepl(
    "2 sample\\(s\\) and 2 feature\\(s\\)",
    qc$qc_observations$message
  )))
  expect_true(any(grepl("Total reads across all samples: 3", qc$qc_observations$message)))
  expect_true(any(grepl("Overall zero fraction: 50%", qc$qc_observations$message)))
  expect_true(any(grepl("Metadata table is absent", qc$qc_observations$message)))
  expect_true(any(grepl("Taxonomy table is absent", qc$qc_observations$message)))
})

test_that("microeda_qc observations summarize metadata and taxonomy presence", {
  counts <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  metadata <- data.frame(
    group = c("A", "B"),
    batch = c("x", "y"),
    row.names = c("s1", "s2")
  )
  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Bacteroidota"),
    Genus = c("Faecalibacterium", "Bacteroides"),
    row.names = c("f1", "f2")
  )

  qc <- microeda_qc(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE
  )

  expect_true("metadata_columns" %in% qc$qc_observations$observation_id)
  expect_true("taxonomy_ranks" %in% qc$qc_observations$observation_id)
  expect_false("metadata_absent" %in% qc$qc_observations$observation_id)
  expect_false("taxonomy_absent" %in% qc$qc_observations$observation_id)
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
  expect_true("flag_zero_library_samples" %in% qc$qc_observations$observation_id)
  expect_equal(
    qc$qc_observations$severity[
      qc$qc_observations$observation_id == "flag_zero_library_samples"
    ],
    "warning"
  )
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
  expect_equal(qc$prevalence_summary$n_features_above_threshold, 3)
  expect_equal(qc$prevalence_summary$n_features_below_threshold, 1)
  expect_equal(qc$prevalence_summary$fraction_features_below_threshold, 0.25)
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
  expect_equal(qc$prevalence_summary$n_features_above_threshold, 0)
  expect_equal(qc$prevalence_summary$n_features_below_threshold, 2)
  expect_equal(qc$prevalence_summary$min_prevalence, 0)
  expect_equal(qc$prevalence_summary$median_prevalence, 0)
  expect_equal(qc$prevalence_summary$max_prevalence, 0)
  expect_equal(qc$prevalence_summary$n_features_detected_in_all_samples, 0)
  expect_equal(qc$prevalence_summary$n_features_detected_in_one_sample, 0)
  expect_equal(
    qc$prevalence_summary$fraction_features_detected_in_one_sample,
    0
  )
  expect_true("high_sparsity" %in% qc$qc_flags$flag_id)
  expect_true("many_features_below_prevalence" %in% qc$qc_flags$flag_id)
  expect_true("flag_high_sparsity" %in% qc$qc_observations$observation_id)
  expect_true("flag_many_features_below_prevalence" %in%
    qc$qc_observations$observation_id)
  expect_false(any(is.na(qc$qc_observations$message)))
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
  expect_equal(qc$prevalence_summary$min_prevalence_threshold, 0.5)
  expect_equal(qc$prevalence_summary$n_features_above_threshold, 1)
  expect_equal(qc$prevalence_summary$n_features_below_threshold, 3)
  expect_equal(qc$prevalence_summary$fraction_features_above_threshold, 0.25)
  expect_equal(qc$prevalence_summary$fraction_features_below_threshold, 0.75)
  expect_true("many_features_below_prevalence" %in% qc$qc_flags$flag_id)
  expect_true("many_single_sample_features" %in% qc$qc_flags$flag_id)
  expect_true("flag_many_features_below_prevalence" %in%
    qc$qc_observations$observation_id)
  expect_true("flag_many_single_sample_features" %in%
    qc$qc_observations$observation_id)
})

test_that("microeda_qc summarizes prevalence edge cases", {
  counts <- matrix(
    c(
      1, 0, 5,
      1, 0, 0,
      1, 0, 0,
      1, 2, 0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("s", 1:4), paste0("f", 1:3))
  )

  qc <- microeda_qc(counts, taxa_are_rows = FALSE, min_prevalence = 0.5)

  expect_equal(qc$per_feature$prevalence, c(1, 0.25, 0.25))
  expect_equal(qc$prevalence_summary$n_features_detected_in_all_samples, 1)
  expect_equal(qc$prevalence_summary$n_features_detected_in_one_sample, 2)
  expect_equal(
    qc$prevalence_summary$fraction_features_detected_in_one_sample,
    2 / 3
  )
  expect_equal(qc$prevalence_summary$q1_prevalence, 0.25)
  expect_equal(qc$prevalence_summary$median_prevalence, 0.25)
  expect_equal(qc$prevalence_summary$q3_prevalence, 0.625)
  expect_true("many_single_sample_features" %in% qc$qc_flags$flag_id)
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
