da_example_counts <- function() {
  counts <- matrix(
    c(
      10, 0, 1, 0,
      8, 2, 0, 1,
      0, 7, 2, 0,
      1, 6, 3, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(4))
  colnames(counts) <- paste0("ASV", seq_len(4))
  counts
}

da_example_metadata <- function(sample_ids = paste0("S", seq_len(4))) {
  data.frame(
    group = c("A", "A", "B", "B"),
    batch = c("x", "y", "x", "y"),
    row.names = sample_ids
  )
}

da_expected_result_columns <- function() {
  c(
    "feature_id",
    "taxon_label",
    "rank",
    "method",
    "contrast",
    "group1",
    "group2",
    "effect",
    "effect_type",
    "log_fold_change",
    "statistic",
    "standard_error",
    "ci_low",
    "ci_high",
    "p_value",
    "p_adjusted",
    "p_adjust_method",
    "p_adjust_scope",
    "significance",
    "direction",
    "method_note"
  )
}

test_that("da_prepare_context accepts matrix and data frame inputs", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Firmicutes", "Bacteroidota", "Bacteroidota"),
    Genus = c("A", "B", "C", "D"),
    row.names = colnames(counts)
  )

  context <- microeda:::da_prepare_context(
    as.data.frame(counts),
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = c("aldex2", "deseq2"),
    tax_rank = "Genus",
    prevalence_filter = 0.25,
    min_count = 2,
    p_adjust_method = "holm",
    taxa_are_rows = FALSE
  )

  expect_s3_class(context, "microeda_da_context")
  expect_equal(dim(context$counts), c(4L, 4L))
  expect_equal(context$feature_ids, colnames(counts))
  expect_equal(context$sample_ids, rownames(counts))
  expect_equal(context$group, "group")
  expect_equal(context$contrast, c("A", "B"))
  expect_equal(context$contrast_label, "A_vs_B")
  expect_equal(context$methods, c("aldex2", "deseq2"))
  expect_equal(context$tax_rank, "Genus")
  expect_equal(context$filters$prevalence_filter, 0.25)
  expect_equal(context$filters$min_count, 2L)
  expect_false(context$filters$applied)
  expect_equal(context$p_adjust_method, "holm")
  expect_equal(context$params$p_adjust_method, "holm")
})

test_that("da_prepare_context validates group metadata clearly", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  metadata_with_missing <- metadata
  metadata_with_missing$group[2] <- NA_character_

  expect_error(
    microeda:::da_prepare_context(
      counts,
      group = "group",
      contrast = c("A", "B"),
      taxa_are_rows = FALSE
    ),
    "metadata"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "missing",
      contrast = c("A", "B"),
      taxa_are_rows = FALSE
    ),
    "column"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata_with_missing,
      group = "group",
      contrast = c("A", "B"),
      taxa_are_rows = FALSE
    ),
    "Group labels"
  )
})

test_that("da_prepare_context validates contrasts", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))

  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "A",
      taxa_are_rows = FALSE
    ),
    "length-2"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "A"),
      taxa_are_rows = FALSE
    ),
    "different"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "C"),
      taxa_are_rows = FALSE
    ),
    "not found"
  )
})

test_that("da_prepare_context validates methods and p adjustment", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))

  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = c("aldex2", "aldex2"),
      taxa_are_rows = FALSE
    ),
    "duplicate"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "edgeR",
      taxa_are_rows = FALSE
    ),
    "Unknown DA method"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      p_adjust_method = "not_a_method",
      taxa_are_rows = FALSE
    ),
    "p_adjust_method"
  )

  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    p_adjust_method = "holm",
    taxa_are_rows = FALSE
  )
  expect_equal(context$p_adjust_method, "holm")
  expect_equal(context$params$p_adjust_method, "holm")
  expect_true("method_native_p_adjustment" %in% context$caveats$caveat_id)
})

test_that("da_prepare_context returns expected internal fields", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    taxa_are_rows = FALSE
  )

  expect_named(
    context,
    c(
      "counts",
      "metadata",
      "taxonomy",
      "group",
      "contrast",
      "contrast_label",
      "methods",
      "tax_rank",
      "filters",
      "p_adjust_method",
      "feature_ids",
      "sample_ids",
      "group_values",
      "caveats",
      "params",
      "call"
    )
  )
  expect_null(context$p_adjust_method)
  expect_null(context$params$p_adjust_method)
  expect_true("method_native_p_adjustment" %in% context$caveats$caveat_id)
  adjustment_note <- context$caveats[
    context$caveats$caveat_id == "method_native_p_adjustment",
    ,
    drop = FALSE
  ]
  expect_match(adjustment_note$message, "method-native p-value adjustment")
  expect_match(adjustment_note$message, "not globally re-adjust")
  expect_named(
    context$caveats,
    c("method", "caveat_id", "topic", "severity", "message")
  )
})

test_that("da_prepare_context records input caveats without requiring backends", {
  counts <- matrix(
    c(
      0.5, 0, 0, 0,
      0, 1.5, 0, 0,
      0, 0, 2.5, 0
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(3))
  colnames(counts) <- paste0("ASV", seq_len(4))
  metadata <- data.frame(
    group = c("A", "A", "B"),
    row.names = rownames(counts)
  )

  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )

  expect_true("small_group_size" %in% context$caveats$caveat_id)
  expect_true("high_sparsity" %in% context$caveats$caveat_id)
  expect_true("non_integer_counts" %in% context$caveats$caveat_id)
  expect_true("taxonomy_unavailable" %in% context$caveats$caveat_id)
})

test_that("standardized DA empty result has stable columns", {
  empty <- microeda:::da_empty_result()

  expect_s3_class(empty, "data.frame")
  expect_equal(nrow(empty), 0L)
  expect_named(empty, da_expected_result_columns())

  result <- microeda:::da_standard_result(
    feature_id = c("ASV1", "ASV2"),
    method = "aldex2",
    contrast = "A_vs_B"
  )
  expect_named(result, da_expected_result_columns())
  expect_equal(nrow(result), 2L)
  expect_equal(result$method, c("aldex2", "aldex2"))
  expect_equal(result$contrast, c("A_vs_B", "A_vs_B"))
})

test_that("DA method notes include DESeq2 sensitivity caveat", {
  notes <- microeda:::da_method_notes(c("aldex2", "ancombc2", "deseq2"))

  expect_s3_class(notes, "data.frame")
  expect_named(notes, c("method", "caveat_id", "topic", "severity", "message"))
  expect_equal(notes$method, c("aldex2", "ancombc2", "deseq2"))
  deseq2_note <- notes[notes$method == "deseq2", , drop = FALSE]
  expect_match(deseq2_note$message, "sensitivity", ignore.case = TRUE)
  expect_match(deseq2_note$message, "compositionality", ignore.case = TRUE)
})

test_that("DA skeleton adds no public exports or backend dependencies", {
  exports <- getNamespaceExports("microeda")
  expect_false("microeda_da" %in% exports)
  expect_false(any(c(
    "da_prepare_context",
    "da_standard_result",
    "da_empty_result",
    "da_method_notes"
  ) %in% exports))

  description <- utils::packageDescription("microeda")
  dependency_text <- paste(description$Imports, description$Suggests, collapse = ",")
  expect_false(grepl("\\bALDEx2\\b", dependency_text))
  expect_false(grepl("\\bANCOMBC\\b", dependency_text))
  expect_false(grepl("\\bDESeq2\\b", dependency_text))
})
