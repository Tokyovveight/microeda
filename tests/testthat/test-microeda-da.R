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

da_fake_backend_result <- function(method,
                                   features = c("ASV1", "ASV2"),
                                   raw_output = list(fixture = method),
                                   notes = NULL) {
  result <- microeda:::da_standard_result(
    feature_id = features,
    method = method,
    contrast = "A_vs_B",
    group1 = "A",
    group2 = "B",
    effect = seq_along(features),
    effect_type = "fixture",
    p_value = seq_along(features) / 10,
    p_adjusted = seq_along(features) / 5,
    p_adjust_method = "native",
    p_adjust_scope = "method"
  )

  microeda:::da_backend_result(
    method = method,
    results = result,
    raw_output = raw_output,
    notes = notes,
    params = list(fixture = TRUE)
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
  expect_named(
    context$contrast_plan,
    c("contrast", "group1", "group2", "contrast_type")
  )
  expect_equal(nrow(context$contrast_plan), 1L)
  expect_equal(context$contrast_plan$contrast, "A_vs_B")
  expect_equal(context$contrast_plan$group1, "A")
  expect_equal(context$contrast_plan$group2, "B")
  expect_equal(context$contrast_plan$contrast_type, "explicit")
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
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "all",
      taxa_are_rows = FALSE
    ),
    "pairwise"
  )
  expect_error(
    microeda:::da_prepare_context(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B", "C"),
      taxa_are_rows = FALSE
    ),
    "length-2"
  )
})

test_that("da_prepare_context creates pairwise contrast plans", {
  counts <- matrix(
    c(
      10, 0, 1,
      8, 2, 0,
      0, 7, 2,
      1, 6, 3,
      2, 1, 8,
      1, 2, 7
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(6))
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "A", "B", "B", "C", "C"),
    row.names = rownames(counts)
  )

  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    taxa_are_rows = FALSE
  )

  expect_equal(context$contrast, "pairwise")
  expect_equal(context$contrast_label, "pairwise")
  expect_named(
    context$contrast_plan,
    c("contrast", "group1", "group2", "contrast_type")
  )
  expect_equal(nrow(context$contrast_plan), 3L)
  expect_equal(context$contrast_plan$contrast, c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_equal(context$contrast_plan$group1, c("A", "A", "B"))
  expect_equal(context$contrast_plan$group2, c("B", "C", "C"))
  expect_equal(
    context$contrast_plan$contrast_type,
    rep("pairwise", 3)
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
      "contrast_plan",
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

test_that("da_backend_result creates valid backend result objects", {
  backend <- da_fake_backend_result("aldex2")

  expect_s3_class(backend, "microeda_da_backend_result")
  expect_named(
    backend,
    c("method", "results", "raw_output", "notes", "params")
  )
  expect_equal(backend$method, "aldex2")
  expect_named(backend$results, da_expected_result_columns())
  expect_type(backend$raw_output, "list")
  expect_s3_class(backend$notes, "data.frame")
  expect_type(backend$params, "list")
})

test_that("da_validate_backend_result rejects invalid method and result schema", {
  unknown_method <- structure(
    list(
      method = "edgeR",
      results = microeda:::da_empty_result(),
      raw_output = NULL,
      notes = microeda:::da_empty_caveats(),
      params = list()
    ),
    class = "microeda_da_backend_result"
  )
  missing_columns <- structure(
    list(
      method = "aldex2",
      results = microeda:::da_empty_result()[, -1, drop = FALSE],
      raw_output = NULL,
      notes = microeda:::da_empty_caveats(),
      params = list()
    ),
    class = "microeda_da_backend_result"
  )

  expect_error(
    microeda:::da_validate_backend_result(unknown_method),
    "Unknown DA method"
  )
  expect_error(
    microeda:::da_validate_backend_result(missing_columns),
    "standardized DA result columns"
  )
})

test_that("da_combine_method_results preserves method order and raw outputs", {
  backend_results <- list(
    da_fake_backend_result("deseq2"),
    da_fake_backend_result("aldex2")
  )

  combined <- microeda:::da_combine_method_results(backend_results)

  expect_equal(combined$methods, c("deseq2", "aldex2"))
  expect_named(combined$method_results, c("deseq2", "aldex2"))
  expect_named(combined$raw_outputs, c("deseq2", "aldex2"))
  expect_named(combined$results, da_expected_result_columns())
  expect_equal(unique(combined$results$method), c("deseq2", "aldex2"))
  expect_equal(combined$raw_outputs$deseq2$fixture, "deseq2")
  expect_equal(combined$raw_outputs$aldex2$fixture, "aldex2")
})

test_that("da_combine_method_results aggregates notes with method labels", {
  unlabeled_note <- microeda:::da_caveat(
    method = NA_character_,
    caveat_id = "fixture_note",
    topic = "differential_abundance",
    severity = "info",
    message = "Fixture backend note."
  )
  backend_results <- list(
    da_fake_backend_result("aldex2", notes = unlabeled_note),
    da_fake_backend_result("ancombc2")
  )

  combined <- microeda:::da_combine_method_results(backend_results)

  expect_named(
    combined$caveats,
    c("method", "caveat_id", "topic", "severity", "message")
  )
  expect_true("fixture_note" %in% combined$caveats$caveat_id)
  fixture_note <- combined$caveats[
    combined$caveats$caveat_id == "fixture_note",
    ,
    drop = FALSE
  ]
  expect_equal(fixture_note$method, "aldex2")
  expect_true(all(combined$caveats$method %in% c("aldex2", "ancombc2")))
})

test_that("da_build_result_object creates internal DA result skeleton", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Firmicutes", "Bacteroidota", "Bacteroidota"),
    Genus = c("A", "B", "C", "D"),
    row.names = colnames(counts)
  )
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = c("aldex2", "deseq2"),
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
  backend_results <- list(
    da_fake_backend_result("aldex2"),
    da_fake_backend_result("deseq2")
  )

  da_result <- microeda:::da_build_result_object(context, backend_results)

  expect_s3_class(da_result, "microeda_da")
  expect_named(
    da_result,
    c(
      "results",
      "method_results",
      "raw_outputs",
      "methods",
      "group",
      "contrast",
      "contrast_plan",
      "contrast_label",
      "tax_rank",
      "feature_metadata",
      "filters",
      "caveats",
      "params",
      "call"
    )
  )
  expect_named(da_result$results, da_expected_result_columns())
  expect_named(da_result$method_results, c("aldex2", "deseq2"))
  expect_named(da_result$raw_outputs, c("aldex2", "deseq2"))
  expect_equal(da_result$methods, c("aldex2", "deseq2"))
  expect_equal(da_result$contrast_plan$contrast, "A_vs_B")
  expect_equal(da_result$contrast_label, "A_vs_B")
  expect_named(da_result$feature_metadata, c("feature_id", "Phylum", "Genus"))
  expect_true("method_native_p_adjustment" %in% da_result$caveats$caveat_id)
  expect_true("deseq2_sensitivity_note" %in% da_result$caveats$caveat_id)
})

test_that("DA skeleton adds no public exports or backend dependencies", {
  exports <- getNamespaceExports("microeda")
  expect_false("microeda_da" %in% exports)
  expect_false(any(c(
    "da_prepare_context",
    "da_standard_result",
    "da_empty_result",
    "da_method_notes",
    "da_backend_result",
    "da_validate_backend_result",
    "da_standardize_backend_result",
    "da_combine_method_results",
    "da_build_result_object"
  ) %in% exports))

  description <- utils::packageDescription("microeda")
  dependency_text <- paste(description$Imports, description$Suggests, collapse = ",")
  expect_false(grepl("\\bALDEx2\\b", dependency_text))
  expect_false(grepl("\\bANCOMBC\\b", dependency_text))
  expect_false(grepl("\\bDESeq2\\b", dependency_text))

  da_function_names <- c(
    "da_prepare_context",
    "da_backend_result",
    "da_standardize_backend_result",
    "da_combine_method_results",
    "da_build_result_object",
    "da_validate_p_adjust_method"
  )
  da_code <- unlist(lapply(da_function_names, function(function_name) {
    deparse(get(function_name, envir = asNamespace("microeda")))
  }))
  expect_false(any(grepl("p\\.adjust\\s*\\(", da_code)))
})
