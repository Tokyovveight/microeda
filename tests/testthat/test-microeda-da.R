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

da_report_fixture <- function(pairwise = FALSE) {
  if (pairwise) {
    contrast <- c("A_vs_B", "A_vs_B", "A_vs_C", "A_vs_C")
    group1 <- c("A", "A", "A", "A")
    group2 <- c("B", "B", "C", "C")
    contrast_plan <- data.frame(
      contrast = c("A_vs_B", "A_vs_C"),
      group1 = c("A", "A"),
      group2 = c("B", "C"),
      contrast_type = "pairwise",
      stringsAsFactors = FALSE
    )
    contrast_label <- "pairwise"
    contrast_input <- "pairwise"
  } else {
    contrast <- rep("A_vs_B", 4)
    group1 <- rep("A", 4)
    group2 <- rep("B", 4)
    contrast_plan <- data.frame(
      contrast = "A_vs_B",
      group1 = "A",
      group2 = "B",
      contrast_type = "explicit",
      stringsAsFactors = FALSE
    )
    contrast_label <- "A_vs_B"
    contrast_input <- c("A", "B")
  }

  results <- microeda:::da_standard_result(
    feature_id = paste0("ASV", seq_len(4)),
    taxon_label = paste0("Genus", seq_len(4)),
    rank = "Genus",
    method = "aldex2",
    contrast = contrast,
    group1 = group1,
    group2 = group2,
    effect = c(0.4, -0.2, 3.1, 1.7),
    effect_type = "aldex2_effect",
    p_value = c(0.18, NA, 0.009, 0.031),
    p_adjusted = c(0.2, NA, 0.01, 0.04),
    p_adjust_method = "aldex2_native_BH",
    p_adjust_scope = "method_contrast",
    significance = c("ns", NA, "**", "*"),
    method_note = "ALDEx2 note"
  )

  structure(
    list(
      results = results,
      method_results = list(aldex2 = list(results = results)),
      raw_outputs = list(aldex2 = list(fixture = TRUE)),
      methods = "aldex2",
      group = "group",
      contrast = contrast_input,
      contrast_plan = contrast_plan,
      contrast_label = contrast_label,
      tax_rank = "Genus",
      feature_metadata = data.frame(
        feature_id = paste0("ASV", seq_len(4)),
        Genus = paste0("Genus", seq_len(4)),
        stringsAsFactors = FALSE
      ),
      filters = list(),
      caveats = microeda:::da_caveat(
        method = NA_character_,
        caveat_id = "method_native_p_adjustment",
        topic = "differential_abundance",
        severity = "info",
        message = "Use method-native adjustment."
      ),
      params = list(),
      call = quote(microeda_da())
    ),
    class = "microeda_da"
  )
}

da_aldex2_example_counts <- function() {
  counts <- matrix(
    c(
      40, 8, 20, 2,
      38, 9, 22, 3,
      42, 7, 19, 2,
      12, 35, 18, 4,
      10, 37, 17, 5,
      11, 33, 20, 4
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(6))
  colnames(counts) <- paste0("ASV", seq_len(4))
  counts
}

da_aldex2_example_metadata <- function(sample_ids = paste0("S", seq_len(6))) {
  data.frame(
    group = c("A", "A", "A", "B", "B", "B"),
    row.names = sample_ids
  )
}

da_aldex2_pairwise_counts <- function() {
  counts <- matrix(
    c(
      40, 8, 20, 2,
      38, 9, 22, 3,
      42, 7, 19, 2,
      12, 35, 18, 4,
      10, 37, 17, 5,
      11, 33, 20, 4,
      8, 15, 40, 12,
      7, 14, 38, 13,
      9, 16, 42, 11
    ),
    nrow = 9,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(9))
  colnames(counts) <- paste0("ASV", seq_len(4))
  counts
}

da_aldex2_pairwise_metadata <- function(sample_ids = paste0("S", seq_len(9))) {
  data.frame(
    group = c("A", "A", "A", "B", "B", "B", "C", "C", "C"),
    row.names = sample_ids
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

test_that("da_deduplicate_caveats collapses exact duplicates stably", {
  input_note <- microeda:::da_caveat(
    method = NA_character_,
    caveat_id = "input_note",
    topic = "input",
    severity = "warning",
    message = "Input-level fixture note."
  )
  method_note <- microeda:::da_caveat(
    method = "aldex2",
    caveat_id = "method_note",
    topic = "differential_abundance",
    severity = "info",
    message = "Method-level fixture note."
  )
  later_note <- microeda:::da_caveat(
    method = NA_character_,
    caveat_id = "later_note",
    topic = "taxonomy",
    severity = "info",
    message = "Later fixture note."
  )

  caveats <- rbind(input_note, method_note, method_note, later_note)
  deduplicated <- microeda:::da_deduplicate_caveats(caveats)

  expect_equal(nrow(deduplicated), 3L)
  expect_equal(
    deduplicated$caveat_id,
    c("input_note", "method_note", "later_note")
  )
  expect_true(is.na(deduplicated$method[1]))
  expect_equal(deduplicated$method[2], "aldex2")
  expect_true(is.na(deduplicated$method[3]))
})

test_that("da_run_aldex2 validates context and method before execution", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "deseq2",
    taxa_are_rows = FALSE
  )

  expect_error(
    microeda:::da_run_aldex2(list()),
    "microeda_da_context"
  )
  expect_error(
    microeda:::da_run_aldex2(context),
    "aldex2"
  )
})

test_that("da_run_aldex2 reports missing optional ALDEx2 clearly", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE
  )

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) FALSE,
    .package = "microeda"
  )

  expect_error(
    microeda:::da_run_aldex2(context),
    "requires the optional package `ALDEx2`",
    fixed = TRUE
  )
})

test_that("da_run_aldex2 executes pairwise contrast plans", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_pairwise_counts()
  metadata <- da_aldex2_pairwise_metadata(rownames(counts))
  taxonomy <- data.frame(
    Genus = paste0("Genus", seq_len(ncol(counts))),
    row.names = colnames(counts)
  )
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
  backend <- microeda:::da_run_aldex2(context, mc.samples = 16)

  expected_contrasts <- c("A_vs_B", "A_vs_C", "B_vs_C")
  expect_s3_class(backend, "microeda_da_backend_result")
  expect_equal(backend$method, "aldex2")
  expect_named(backend$results, da_expected_result_columns())
  expect_equal(nrow(backend$results), ncol(counts) * length(expected_contrasts))
  expect_equal(unique(backend$results$contrast), expected_contrasts)
  expect_equal(
    unique(backend$results[backend$results$contrast == "A_vs_B", "group1"]),
    "A"
  )
  expect_equal(
    unique(backend$results[backend$results$contrast == "A_vs_B", "group2"]),
    "B"
  )
  expect_equal(
    unique(backend$results[backend$results$contrast == "A_vs_C", "group1"]),
    "A"
  )
  expect_equal(
    unique(backend$results[backend$results$contrast == "A_vs_C", "group2"]),
    "C"
  )
  expect_equal(
    unique(backend$results[backend$results$contrast == "B_vs_C", "group1"]),
    "B"
  )
  expect_equal(
    unique(backend$results[backend$results$contrast == "B_vs_C", "group2"]),
    "C"
  )
  expect_equal(backend$results$feature_id, rep(colnames(counts), length(expected_contrasts)))
  expect_equal(backend$results$taxon_label, rep(taxonomy$Genus, length(expected_contrasts)))
  expect_named(
    backend$raw_output,
    c(
      "contrasts",
      "contrast_plan",
      "params",
      "input_orientation",
      "transposed_from_context"
    )
  )
  expect_named(backend$raw_output$contrasts, expected_contrasts)
  expect_equal(backend$raw_output$contrast_plan, context$contrast_plan)
  expect_equal(backend$raw_output$params$mc.samples, 16L)
  expect_equal(backend$raw_output$input_orientation, "feature_by_sample")
  expect_true(backend$raw_output$transposed_from_context)
  for (contrast in expected_contrasts) {
    raw_contrast <- backend$raw_output$contrasts[[contrast]]
    expect_named(
      raw_contrast,
      c(
        "clr",
        "ttest",
        "effect",
        "combined",
        "conditions",
        "contrast_row",
        "input_orientation",
        "transposed_from_context",
        "params",
        "warnings",
        "messages"
      )
    )
    expect_true(methods::is(raw_contrast$clr, "aldex.clr"))
    expect_s3_class(raw_contrast$ttest, "data.frame")
    expect_s3_class(raw_contrast$effect, "data.frame")
    expect_s3_class(raw_contrast$combined, "data.frame")
    expect_equal(raw_contrast$contrast_row$contrast, contrast)
    expect_equal(raw_contrast$params, backend$params)
  }

  da_result <- microeda:::da_build_result_object(context, list(aldex2 = backend))
  expect_s3_class(da_result, "microeda_da")
  expect_equal(nrow(da_result$results), nrow(backend$results))
  expect_named(da_result$raw_outputs, "aldex2")
  expect_named(da_result$method_results, "aldex2")
  expect_false(any(duplicated(
    da_result$caveats[c("method", "caveat_id", "topic", "severity", "message")]
  )))
})

test_that("da_run_aldex2 validates unsupported contrast plan types", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE
  )
  context$contrast_plan$contrast_type <- "unsupported"
  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    .package = "microeda"
  )

  expect_error(
    microeda:::da_run_aldex2(context, mc.samples = 16),
    "explicit or pairwise"
  )
})

test_that("da_run_aldex2 returns standardized backend results", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_example_counts()
  metadata <- da_aldex2_example_metadata(rownames(counts))
  taxonomy <- data.frame(
    Genus = paste0("Genus", seq_len(ncol(counts))),
    row.names = colnames(counts)
  )
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )

  backend <- microeda:::da_run_aldex2(context, mc.samples = 16)

  expect_s3_class(backend, "microeda_da_backend_result")
  expect_equal(backend$method, "aldex2")
  expect_named(backend$results, da_expected_result_columns())
  expect_equal(nrow(backend$results), ncol(counts))
  expect_equal(backend$results$feature_id, colnames(counts))
  expect_equal(backend$results$taxon_label, taxonomy$Genus)
  expect_equal(backend$results$rank, rep("Genus", ncol(counts)))
  expect_true(all(backend$results$method == "aldex2"))
  expect_true(all(backend$results$contrast == "A_vs_B"))
  expect_true(all(backend$results$group1 == "A"))
  expect_true(all(backend$results$group2 == "B"))
  expect_true(all(backend$results$effect_type == "aldex2_effect"))
  expect_true(all(is.na(backend$results$log_fold_change)))
  expect_true(all(is.na(backend$results$statistic)))
  expect_true(all(is.na(backend$results$standard_error)))
  expect_true(all(is.na(backend$results$ci_low)))
  expect_true(all(is.na(backend$results$ci_high)))
  expect_equal(backend$results$p_adjust_method, rep("aldex2_native_BH", ncol(counts)))
  expect_equal(backend$results$p_adjust_scope, rep("method_contrast", ncol(counts)))
  expect_equal(backend$results$direction, rep(NA_character_, ncol(counts)))
  expect_match(backend$results$method_note[1], "ALDEx2")
  expect_equal(
    backend$results$p_value,
    as.numeric(backend$raw_output$combined$we.ep),
    tolerance = 1e-12
  )
  expect_equal(
    backend$results$p_adjusted,
    as.numeric(backend$raw_output$combined$we.eBH),
    tolerance = 1e-12
  )
  expect_true(all(backend$results$significance %in% c("***", "**", "*", "ns", NA)))
})

test_that("da_run_aldex2 preserves raw output and params", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_example_counts()
  metadata <- da_aldex2_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE
  )

  backend <- microeda:::da_run_aldex2(
    context,
    mc.samples = 16,
    denom = "all",
    paired.test = FALSE
  )

  expect_named(
    backend$raw_output,
    c(
      "clr",
      "ttest",
      "effect",
      "combined",
      "conditions",
      "contrast_row",
      "input_orientation",
      "transposed_from_context",
      "params",
      "warnings",
      "messages"
    )
  )
  expect_true(methods::is(backend$raw_output$clr, "aldex.clr"))
  expect_s3_class(backend$raw_output$ttest, "data.frame")
  expect_s3_class(backend$raw_output$effect, "data.frame")
  expect_s3_class(backend$raw_output$combined, "data.frame")
  expect_equal(rownames(backend$raw_output$combined), colnames(counts))
  expect_equal(backend$raw_output$conditions, metadata$group)
  expect_equal(backend$raw_output$contrast_row$contrast, "A_vs_B")
  expect_equal(backend$raw_output$input_orientation, "feature_by_sample")
  expect_true(backend$raw_output$transposed_from_context)
  expect_equal(backend$raw_output$params$mc.samples, 16L)
  expect_equal(backend$raw_output$params$denom, "all")
  expect_false(backend$raw_output$params$paired.test)
  expect_equal(backend$params, backend$raw_output$params)
})

test_that("microeda_da returns public DA object for explicit ALDEx2 contrast", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_example_counts()
  metadata <- da_aldex2_example_metadata(rownames(counts))
  taxonomy <- data.frame(
    Genus = paste0("Genus", seq_len(ncol(counts))),
    row.names = colnames(counts)
  )

  da <- microeda_da(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )

  expect_s3_class(da, "microeda_da")
  expect_equal(da$methods, "aldex2")
  expect_equal(da$contrast_label, "A_vs_B")
  expect_named(da$raw_outputs, "aldex2")
  expect_named(da$method_results, "aldex2")
  expect_named(as_da_results(da), da_expected_result_columns())
  expect_equal(as_da_results(da), da$results)
  expect_equal(nrow(as_da_results(da)), ncol(counts))
  expect_true(all(as_da_results(da)$method == "aldex2"))
  expect_true(all(as_da_results(da)$contrast == "A_vs_B"))
})

test_that("microeda_da returns public DA object for pairwise ALDEx2 contrast", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_pairwise_counts()
  metadata <- da_aldex2_pairwise_metadata(rownames(counts))
  expected_contrasts <- c("A_vs_B", "A_vs_C", "B_vs_C")

  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )

  expect_s3_class(da, "microeda_da")
  expect_equal(da$contrast_label, "pairwise")
  expect_equal(da$contrast_plan$contrast, expected_contrasts)
  expect_equal(nrow(as_da_results(da)), ncol(counts) * length(expected_contrasts))
  expect_equal(unique(as_da_results(da)$contrast), expected_contrasts)
  expect_named(da$raw_outputs$aldex2$contrasts, expected_contrasts)
})

test_that("as_da_results validates input", {
  expect_error(
    as_da_results(list()),
    "microeda_da"
  )
})

test_that("print.microeda_da is compact and points to DA outputs", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_pairwise_counts()
  metadata <- da_aldex2_pairwise_metadata(rownames(counts))
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )

  output <- capture.output(returned <- print(da))

  expect_identical(returned, da)
  expect_true(any(grepl("<microeda_da>", output, fixed = TRUE)))
  expect_true(any(grepl("Methods:", output, fixed = TRUE)))
  expect_true(any(grepl("Group:", output, fixed = TRUE)))
  expect_true(any(grepl("Contrasts:", output, fixed = TRUE)))
  expect_true(any(grepl("Result rows:", output, fixed = TRUE)))
  expect_true(any(grepl("Caveats:", output, fixed = TRUE)))
  expect_true(any(grepl("as_da_results", output, fixed = TRUE)))
  expect_true(any(grepl("raw_outputs", output, fixed = TRUE)))
  expect_false(any(grepl("ASV1", output, fixed = TRUE)))
})

test_that("microeda_da_report returns compact character output", {
  da <- da_report_fixture()
  report <- microeda_da_report(da, top_n = 3, alpha = 0.05, digits = 3)

  expect_type(report, "character")
  expect_length(report, 1L)
  expect_match(report, "microeda differential representation report", fixed = TRUE)
  expect_match(report, "Methods: aldex2", fixed = TRUE)
  expect_match(report, "Group: group", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_B", fixed = TRUE)
  expect_match(report, "Result rows: 4", fixed = TRUE)
  expect_match(report, "Unique features: 4", fixed = TRUE)
  expect_match(report, "Contrasts: 1", fixed = TRUE)
  expect_match(report, "p_adjust_method: aldex2_native_BH", fixed = TRUE)
  expect_match(report, "Backend-native/method-specific adjustment is used", fixed = TRUE)
  expect_match(report, "Caveats:", fixed = TRUE)
  expect_match(
    report,
    "These rows are exploratory method outputs, not confirmed biological discoveries.",
    fixed = TRUE
  )
  expect_match(report, "x$raw_outputs", fixed = TRUE)
  expect_match(report, "as_da_results(x)", fixed = TRUE)
})

test_that("microeda_da_report sorts and limits top rows", {
  da <- da_report_fixture()
  report <- microeda_da_report(da, top_n = 2, alpha = 0.05, digits = 3)

  asv3_position <- regexpr("ASV3", report, fixed = TRUE)[[1]]
  asv4_position <- regexpr("ASV4", report, fixed = TRUE)[[1]]
  asv1_position <- regexpr("ASV1", report, fixed = TRUE)[[1]]
  asv2_position <- regexpr("ASV2", report, fixed = TRUE)[[1]]

  expect_gt(asv3_position, 0L)
  expect_gt(asv4_position, 0L)
  expect_lt(asv3_position, asv4_position)
  expect_equal(asv1_position, -1L)
  expect_equal(asv2_position, -1L)
})

test_that("microeda_da_report alpha threshold affects per-contrast counts", {
  da <- da_report_fixture()
  report_005 <- microeda_da_report(da, top_n = 0, alpha = 0.05, digits = 3)
  report_002 <- microeda_da_report(da, top_n = 0, alpha = 0.02, digits = 3)

  expect_match(report_005, "A_vs_B\\s+4\\s+2\\s+0\\.01")
  expect_match(report_002, "A_vs_B\\s+4\\s+1\\s+0\\.01")
  expect_match(report_005, "No top rows requested.", fixed = TRUE)
})

test_that("microeda_da_report works for explicit and pairwise DA objects", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_example_counts()
  metadata <- da_aldex2_example_metadata(rownames(counts))
  da_explicit <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  before <- as_da_results(da_explicit)
  explicit_report <- microeda_da_report(da_explicit, top_n = 2)

  expect_match(explicit_report, "Contrast: A_vs_B", fixed = TRUE)
  expect_equal(as_da_results(da_explicit), before)

  pairwise_counts <- da_aldex2_pairwise_counts()
  pairwise_metadata <- da_aldex2_pairwise_metadata(rownames(pairwise_counts))
  da_pairwise <- microeda_da(
    pairwise_counts,
    metadata = pairwise_metadata,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  pairwise_report <- microeda_da_report(da_pairwise, top_n = 2)

  expect_match(pairwise_report, "Pairwise contrasts: 3", fixed = TRUE)
  expect_match(pairwise_report, "A_vs_B", fixed = TRUE)
  expect_match(pairwise_report, "A_vs_C", fixed = TRUE)
  expect_match(pairwise_report, "B_vs_C", fixed = TRUE)
})

test_that("microeda_da_report validates inputs clearly", {
  da <- da_report_fixture()

  expect_error(microeda_da_report(list()), "microeda_da")
  expect_error(microeda_da_report(da, top_n = -1), "top_n")
  expect_error(microeda_da_report(da, top_n = 1.5), "top_n")
  expect_error(microeda_da_report(da, top_n = c(1, 2)), "top_n")
  expect_error(microeda_da_report(da, alpha = -0.1), "alpha")
  expect_error(microeda_da_report(da, alpha = 1.1), "alpha")
  expect_error(microeda_da_report(da, alpha = c(0.01, 0.05)), "alpha")
  expect_error(microeda_da_report(da, digits = -1), "digits")
  expect_error(microeda_da_report(da, digits = 1.2), "digits")
  expect_error(microeda_da_report(da, digits = c(2, 3)), "digits")
})

test_that("write_da_results writes standardized DA CSV only", {
  da <- da_report_fixture()
  before <- as_da_results(da)
  file <- tempfile(fileext = ".csv")

  returned <- write_da_results(da, file)
  after <- as_da_results(da)
  csv <- utils::read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = ""
  )
  raw_lines <- readLines(file, warn = FALSE)

  expect_identical(returned, file)
  expect_true(file.exists(file))
  expect_named(csv, da_expected_result_columns())
  expect_false(any(names(csv) %in% c("", "X", "row.names")))
  expect_equal(nrow(csv), nrow(before))
  expect_equal(csv$feature_id, before$feature_id)
  expect_equal(csv$method, before$method)
  expect_equal(csv$contrast, before$contrast)
  expect_equal(csv$taxon_label, before$taxon_label)
  expect_equal(csv$p_adjusted, before$p_adjusted)
  expect_false(any(grepl("raw_outputs", raw_lines, fixed = TRUE)))
  expect_false(any(grepl("fixture", raw_lines, fixed = TRUE)))
  expect_equal(after, before)

  invisible_file <- tempfile(fileext = ".csv")
  output <- capture.output(
    visible_result <- withVisible(write_da_results(da, invisible_file))
  )
  expect_equal(output, character())
  expect_false(visible_result$visible)
  expect_identical(visible_result$value, invisible_file)
})

test_that("write_da_results validates inputs clearly", {
  da <- da_report_fixture()
  file <- tempfile(fileext = ".csv")

  expect_error(write_da_results(list(), file), "microeda_da")
  expect_error(write_da_results(da, ""), "file")
  expect_error(write_da_results(da, NA_character_), "file")
  expect_error(write_da_results(da, c(file, file)), "file")
  expect_error(write_da_results(da, file, na = c("", "NA")), "na")
  expect_error(write_da_results(da, file, na = NA_character_), "na")
  expect_error(write_da_results(da, file, quote = c(TRUE, FALSE)), "quote")
  expect_error(write_da_results(da, file, quote = NA), "quote")
  expect_error(write_da_results(da, file, quote = "TRUE"), "quote")
})

test_that("write_da_results works for explicit and pairwise DA objects", {
  skip_if_not_installed("ALDEx2")

  counts <- da_aldex2_example_counts()
  metadata <- da_aldex2_example_metadata(rownames(counts))
  da_explicit <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  explicit_file <- tempfile(fileext = ".csv")
  explicit_before <- as_da_results(da_explicit)
  explicit_returned <- write_da_results(da_explicit, explicit_file)
  explicit_csv <- utils::read.csv(
    explicit_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = ""
  )

  expect_identical(explicit_returned, explicit_file)
  expect_named(explicit_csv, da_expected_result_columns())
  expect_equal(nrow(explicit_csv), nrow(explicit_before))
  expect_equal(explicit_csv$feature_id, explicit_before$feature_id)
  expect_equal(explicit_csv$contrast, explicit_before$contrast)
  expect_equal(as_da_results(da_explicit), explicit_before)

  pairwise_counts <- da_aldex2_pairwise_counts()
  pairwise_metadata <- da_aldex2_pairwise_metadata(rownames(pairwise_counts))
  da_pairwise <- microeda_da(
    pairwise_counts,
    metadata = pairwise_metadata,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  pairwise_file <- tempfile(fileext = ".csv")
  pairwise_before <- as_da_results(da_pairwise)
  pairwise_returned <- write_da_results(da_pairwise, pairwise_file)
  pairwise_csv <- utils::read.csv(
    pairwise_file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = ""
  )

  expect_identical(pairwise_returned, pairwise_file)
  expect_named(pairwise_csv, da_expected_result_columns())
  expect_equal(nrow(pairwise_csv), nrow(pairwise_before))
  expect_equal(unique(pairwise_csv$contrast), c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_false(any(grepl("raw_outputs", readLines(pairwise_file, warn = FALSE), fixed = TRUE)))
  expect_equal(as_da_results(da_pairwise), pairwise_before)
})

test_that("microeda_da rejects planned and unknown methods clearly", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE
    ),
    "Only methods = \"aldex2\" is implemented",
    fixed = TRUE
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "deseq2",
      taxa_are_rows = FALSE
    ),
    "Only methods = \"aldex2\" is implemented",
    fixed = TRUE
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "unknown",
      taxa_are_rows = FALSE
    ),
    "Unknown DA method"
  )
})

test_that("microeda_da reports missing optional ALDEx2 clearly", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) FALSE,
    .package = "microeda"
  )

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "aldex2",
      taxa_are_rows = FALSE
    ),
    "requires the optional package `ALDEx2`",
    fixed = TRUE
  )
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
  expect_true("small_group_size" %in% da_result$caveats$caveat_id)
  expect_true("aldex2_compositional_note" %in% da_result$caveats$caveat_id)
  expect_true("deseq2_sensitivity_note" %in% da_result$caveats$caveat_id)
  expect_equal(
    sum(da_result$caveats$caveat_id == "aldex2_compositional_note"),
    1L
  )
  expect_equal(
    sum(da_result$caveats$caveat_id == "deseq2_sensitivity_note"),
    1L
  )
  expect_equal(
    da_result$caveats$caveat_id[seq_len(nrow(context$caveats))],
    context$caveats$caveat_id
  )
})

test_that("da_build_result_object accepts named backend result lists", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE
  )
  backend <- da_fake_backend_result("aldex2")

  named_result <- microeda:::da_build_result_object(
    context,
    list(aldex2 = backend)
  )
  unnamed_result <- microeda:::da_build_result_object(
    context,
    list(backend)
  )

  expect_s3_class(named_result, "microeda_da")
  expect_s3_class(unnamed_result, "microeda_da")
  expect_equal(named_result$methods, "aldex2")
  expect_equal(unnamed_result$methods, "aldex2")
  expect_named(named_result$raw_outputs, "aldex2")
  expect_named(named_result$method_results, "aldex2")
  expect_named(unnamed_result$raw_outputs, "aldex2")
  expect_named(unnamed_result$method_results, "aldex2")
})

test_that("da_build_result_object still rejects method mismatches", {
  counts <- da_example_counts()
  metadata <- da_example_metadata(rownames(counts))
  context_one <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = "aldex2",
    taxa_are_rows = FALSE
  )
  context_two <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = c("aldex2", "deseq2"),
    taxa_are_rows = FALSE
  )

  expect_error(
    microeda:::da_build_result_object(
      context_one,
      list(aldex2 = da_fake_backend_result("deseq2"))
    ),
    "methods must match"
  )
  expect_error(
    microeda:::da_build_result_object(
      context_two,
      list(
        deseq2 = da_fake_backend_result("deseq2"),
        aldex2 = da_fake_backend_result("aldex2")
      )
    ),
    "methods must match"
  )
  expect_error(
    microeda:::da_build_result_object(
      context_two,
      list(aldex2 = da_fake_backend_result("aldex2"))
    ),
    "methods must match"
  )
})

test_that("DA public exports are limited and backend dependencies stay optional", {
  exports <- getNamespaceExports("microeda")
  expect_true("microeda_da" %in% exports)
  expect_true("as_da_results" %in% exports)
  expect_true("microeda_da_report" %in% exports)
  expect_true("write_da_results" %in% exports)
  expect_false(any(c(
    "da_prepare_context",
    "da_standard_result",
    "da_empty_result",
    "da_method_notes",
    "da_backend_result",
    "da_validate_backend_result",
    "da_standardize_backend_result",
    "da_combine_method_results",
    "da_build_result_object",
    "da_deduplicate_caveats",
    "da_run_aldex2",
    "da_run_aldex2_contrast",
    "da_standardize_aldex2_result",
    "da_report_contrast_summary",
    "da_report_top_rows",
    "da_report_table_lines"
  ) %in% exports))

  description <- utils::packageDescription("microeda")
  imports <- if (is.null(description$Imports)) "" else description$Imports
  suggests <- if (is.null(description$Suggests)) "" else description$Suggests
  dependency_text <- paste(imports, suggests, collapse = ",")
  expect_false(grepl("\\bALDEx2\\b", imports))
  expect_true(grepl("\\bALDEx2\\b", suggests))
  expect_false(grepl("\\bANCOMBC\\b", dependency_text))
  expect_false(grepl("\\bDESeq2\\b", dependency_text))

  da_function_names <- c(
    "microeda_da",
    "as_da_results",
    "microeda_da_report",
    "write_da_results",
    "da_prepare_context",
    "da_run_aldex2",
    "da_run_aldex2_contrast",
    "da_standardize_aldex2_result",
    "da_backend_result",
    "da_standardize_backend_result",
    "da_combine_method_results",
    "da_build_result_object",
    "da_deduplicate_caveats",
    "da_p_significance",
    "da_validate_p_adjust_method"
  )
  da_code <- unlist(lapply(da_function_names, function(function_name) {
    deparse(get(function_name, envir = asNamespace("microeda")))
  }))
  expect_false(any(grepl("p\\.adjust\\s*\\(", da_code)))
})
