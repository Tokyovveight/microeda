da_ancombc2_counts <- function() {
  counts <- t(rbind(
    F_A = c(310, 300, 320, 305, 315, 295, 20, 22, 18, 21, 19, 23),
    F_B = c(18, 20, 22, 19, 21, 17, 300, 310, 295, 305, 315, 290),
    F_stable = c(100, 98, 102, 101, 99, 103, 100, 102, 98, 101, 99, 103),
    F_other = c(50, 48, 52, 49, 51, 47, 50, 49, 51, 48, 52, 47)
  ))
  rownames(counts) <- c(paste0("A", seq_len(6)), paste0("B", seq_len(6)))
  counts
}

da_ancombc2_metadata <- function(sample_ids = NULL) {
  metadata <- data.frame(
    group = c(rep("A", 6), rep("B", 6)),
    row.names = c(paste0("A", seq_len(6)), paste0("B", seq_len(6)))
  )
  if (is.null(sample_ids)) {
    return(metadata)
  }
  metadata[sample_ids, , drop = FALSE]
}

da_ancombc2_taxonomy <- function() {
  data.frame(
    Genus = c("Genus A", "Genus B", "Stable genus", "Other genus"),
    row.names = c("F_A", "F_B", "F_stable", "F_other")
  )
}

da_ancombc2_expected_result_columns <- function() {
  names(microeda:::da_empty_result())
}

da_ancombc2_context <- function(counts = da_ancombc2_counts(),
                                metadata = da_ancombc2_metadata(rownames(counts)),
                                taxonomy = da_ancombc2_taxonomy()) {
  microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "ancombc2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
}

da_ancombc2_fake_native <- function(data,
                                    meta_data,
                                    params,
                                    exclude_last = FALSE) {
  model_columns <- colnames(stats::model.matrix(
    stats::reformulate(params$fix_formula),
    data = meta_data
  ))
  coefficient <- setdiff(model_columns, "(Intercept)")
  stopifnot(length(coefficient) == 1)

  feature_ids <- rownames(data)
  if (exclude_last) {
    feature_ids <- head(feature_ids, -1)
  }
  n <- length(feature_ids)
  lfc <- rep(c(-2.5, 2.25, 0.1, -0.1), length.out = n)
  se <- rep(c(0.2, 0.25, 0.3, 0.35), length.out = n)
  statistic <- lfc / se
  p_value <- rep(c(0.001, 0.002, 0.8, 0.9), length.out = n)
  p_adjusted <- rep(c(0.004, 0.006, 1, 1), length.out = n)

  res <- data.frame(taxon = feature_ids, stringsAsFactors = FALSE)
  res[[paste0("lfc_", coefficient)]] <- lfc
  res[[paste0("se_", coefficient)]] <- se
  res[[paste0("W_", coefficient)]] <- statistic
  res[[paste0("p_", coefficient)]] <- p_value
  res[[paste0("q_", coefficient)]] <- p_adjusted
  res[[paste0("diff_", coefficient)]] <- p_adjusted <= 0.05

  list(
    feature_table = data,
    bias_correct_log_table = NULL,
    ss_tab = data.frame(taxon = feature_ids),
    zero_ind = NULL,
    samp_frac = NULL,
    delta_em = NULL,
    delta_wls = NULL,
    res = res,
    res_global = NULL,
    res_pair = NULL,
    res_dunn = NULL,
    res_trend = NULL
  )
}

da_mock_ancombc2_backend <- function(exclude_last = FALSE,
                                     emit_conditions = FALSE) {
  force(exclude_last)
  force(emit_conditions)
  function(data, meta_data, params) {
    if (emit_conditions) {
      warning("fixture backend warning", call. = FALSE)
      message("fixture backend message")
    }
    da_ancombc2_fake_native(
      data = data,
      meta_data = meta_data,
      params = params,
      exclude_last = exclude_last
    )
  }
}

test_that("microeda_da validates ANCOM-BC2 dispatch boundaries", {
  counts <- da_ancombc2_counts()
  metadata <- da_ancombc2_metadata(rownames(counts))

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "pairwise",
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "exactly one explicit contrast"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE
    ),
    "one method per call"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE,
      paired.test = TRUE
    ),
    "only supported by the ALDEx2 backend"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE,
      pair_id = "pair"
    ),
    "only supported by the paired ALDEx2 backend"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE,
      mc.samples = 16
    ),
    "ALDEx2-specific"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE,
      denom = "all"
    ),
    "ALDEx2-specific"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "aldex2",
      taxa_are_rows = FALSE,
      ancombc2_p_adj_method = "holm"
    ),
    "can only be supplied"
  )
})

test_that("ANCOM-BC2 optional package and adjustment errors are clear", {
  context <- da_ancombc2_context()
  contrast_row <- context$contrast_plan[1, , drop = FALSE]

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) FALSE,
    .package = "microeda"
  )
  expect_error(
    microeda:::da_run_ancombc2(
      context,
      contrast_row,
      params = list(p_adj_method = "holm")
    ),
    'BiocManager::install("ANCOMBC")',
    fixed = TRUE
  )

  expect_error(
    microeda:::da_validate_ancombc2_params(
      list(p_adj_method = "not-a-method")
    ),
    "stats::p.adjust.methods"
  )
})

test_that("ANCOM-BC2 validates counts without modifying them", {
  counts <- da_ancombc2_counts()
  metadata <- da_ancombc2_metadata(rownames(counts))

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    .package = "microeda"
  )

  non_integer <- counts
  non_integer[1, 1] <- non_integer[1, 1] + 0.5
  expect_error(
    microeda_da(
      non_integer,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "integer-like"
  )

  negative <- counts
  negative[1, 1] <- -1
  expect_error(
    microeda_da(
      negative,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "negative"
  )

  non_finite <- counts
  non_finite[1, 1] <- Inf
  expect_error(
    microeda_da(
      non_finite,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "finite"
  )

  zero_library <- counts
  zero_library[1, ] <- 0
  expect_error(
    microeda_da(
      zero_library,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "zero-library sample"
  )

  all_zero_feature <- counts
  all_zero_feature[, "F_other"] <- 0
  expect_error(
    microeda_da(
      all_zero_feature,
      metadata = metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      taxa_are_rows = FALSE
    ),
    "all-zero feature"
  )
})

test_that("ANCOM-BC2 input is aligned and uses safe two-column metadata", {
  counts <- da_ancombc2_counts()
  metadata <- da_ancombc2_metadata(rownames(counts))
  metadata$microeda_ancombc2_group <- "collision"
  metadata$microeda_ancombc2_metadata <- "collision"
  context <- da_ancombc2_context(counts, metadata)
  input <- microeda:::da_prepare_ancombc2_input(
    context,
    context$contrast_plan[1, , drop = FALSE]
  )

  expect_equal(dim(input$counts), c(ncol(counts), nrow(counts)))
  expect_equal(colnames(input$counts), rownames(counts))
  expect_equal(rownames(input$counts), colnames(counts))
  expect_identical(input$sample_order, rownames(counts))
  expect_identical(input$feature_order, colnames(counts))
  expect_equal(ncol(input$metadata), 2L)
  expect_false(input$backend_group_column %in% colnames(metadata))
  expect_false(input$inert_metadata_column %in% colnames(metadata))
  expect_false(input$backend_group_column == input$inert_metadata_column)
  expect_identical(
    levels(input$metadata[[input$backend_group_column]]),
    c("A", "B")
  )
  expect_identical(
    input$coefficient,
    paste0(input$backend_group_column, "B")
  )
  expect_false(input$inert_metadata_column %in%
                 microeda:::da_ancombc2_call_params(
                   input$backend_group_column,
                   "holm"
                 )$fix_formula)
  expect_equal(metadata$microeda_ancombc2_group, rep("collision", nrow(metadata)))

  misaligned <- context
  rownames(misaligned$metadata) <- rev(rownames(misaligned$metadata))
  expect_error(
    microeda:::da_prepare_ancombc2_input(
      misaligned,
      misaligned$contrast_plan[1, , drop = FALSE]
    ),
    "aligned"
  )

  duplicated <- context
  duplicated$sample_ids[2] <- duplicated$sample_ids[1]
  expect_error(
    microeda:::da_prepare_ancombc2_input(
      duplicated,
      duplicated$contrast_plan[1, , drop = FALSE]
    ),
    "unique"
  )
})

test_that("ANCOM-BC2 standardizer maps exact native coefficient columns", {
  context <- da_ancombc2_context()
  contrast_row <- context$contrast_plan[1, , drop = FALSE]
  input <- microeda:::da_prepare_ancombc2_input(context, contrast_row)
  native <- da_ancombc2_fake_native(
    data = input$counts,
    meta_data = input$metadata,
    params = microeda:::da_ancombc2_call_params(
      input$backend_group_column,
      "holm"
    )
  )
  results <- microeda:::da_standardize_ancombc2_result(
    res = native$res,
    context = context,
    contrast_row = contrast_row,
    coefficient = input$coefficient,
    p_adj_method = "holm"
  )

  expect_named(results, da_ancombc2_expected_result_columns())
  expect_equal(ncol(results), 21L)
  expect_equal(results$feature_id, rownames(input$counts))
  expect_equal(results$taxon_label, da_ancombc2_taxonomy()$Genus)
  expect_equal(results$rank, rep("Genus", nrow(results)))
  expect_equal(results$method, rep("ancombc2", nrow(results)))
  expect_equal(results$contrast, rep("A_vs_B", nrow(results)))
  expect_equal(results$group1, rep("A", nrow(results)))
  expect_equal(results$group2, rep("B", nrow(results)))
  expect_equal(results$effect, native$res[[paste0("lfc_", input$coefficient)]])
  expect_equal(results$log_fold_change, results$effect)
  expect_equal(
    results$effect_type,
    rep("ancombc2_log_fold_change_group2_vs_group1", nrow(results))
  )
  expect_equal(
    results$statistic,
    native$res[[paste0("W_", input$coefficient)]]
  )
  expect_equal(
    results$standard_error,
    native$res[[paste0("se_", input$coefficient)]]
  )
  expect_equal(results$p_value, native$res[[paste0("p_", input$coefficient)]])
  expect_equal(
    results$p_adjusted,
    native$res[[paste0("q_", input$coefficient)]]
  )
  expect_equal(results$p_adjust_method, rep("holm", nrow(results)))
  expect_equal(results$p_adjust_scope, rep("method_contrast", nrow(results)))
  expect_true(all(is.na(results$ci_low)))
  expect_true(all(is.na(results$ci_high)))
  expect_true(all(is.na(results$direction)))
  expect_type(results$ci_low, "double")
  expect_type(results$ci_high, "double")
  expect_type(results$direction, "character")
  expect_match(results$method_note[1], "natural-log")
  expect_lt(results$effect[results$feature_id == "F_A"], 0)
  expect_gt(results$effect[results$feature_id == "F_B"], 0)

  broken <- native$res
  broken[[paste0("q_", input$coefficient)]] <- NULL
  expect_error(
    microeda:::da_standardize_ancombc2_result(
      broken,
      context,
      contrast_row,
      input$coefficient,
      "holm"
    ),
    "missing required coefficient"
  )
})

test_that("ANCOM-BC2 runner preserves native output and captures diagnostics", {
  context <- da_ancombc2_context()
  contrast_row <- context$contrast_plan[1, , drop = FALSE]

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_call_ancombc2 = da_mock_ancombc2_backend(
      exclude_last = TRUE,
      emit_conditions = TRUE
    ),
    da_ancombc2_package_version = function() "2.14.0",
    .package = "microeda"
  )
  backend <- microeda:::da_run_ancombc2(
    context,
    contrast_row,
    params = list(p_adj_method = "holm")
  )

  expect_s3_class(backend, "microeda_da_backend_result")
  expect_equal(backend$method, "ancombc2")
  expect_named(backend$results, da_ancombc2_expected_result_columns())
  expect_equal(nrow(backend$results), ncol(context$counts) - 1L)
  expect_named(
    backend$raw_output,
    c(
      "native_result",
      "res",
      "res_global",
      "res_pair",
      "res_dunn",
      "res_trend",
      "zero_ind",
      "ss_tab",
      "coefficient",
      "reference_level",
      "factor_levels",
      "original_conditions",
      "backend_conditions",
      "sample_order",
      "feature_order",
      "backend_excluded_features",
      "backend_excluded_feature_count",
      "contrast_row",
      "input_orientation",
      "transposed_from_context",
      "backend_group_column",
      "inert_metadata_column",
      "backend_metadata",
      "params",
      "warnings",
      "messages",
      "package_version"
    )
  )
  expect_identical(backend$raw_output$res, backend$raw_output$native_result$res)
  expect_null(backend$raw_output$res_pair)
  expect_null(backend$raw_output$zero_ind)
  expect_equal(backend$raw_output$reference_level, "A")
  expect_equal(backend$raw_output$factor_levels, c("A", "B"))
  expect_equal(backend$raw_output$input_orientation, "feature_by_sample")
  expect_true(backend$raw_output$transposed_from_context)
  expect_equal(backend$raw_output$package_version, "2.14.0")
  expect_equal(backend$raw_output$params$p_adj_method, "holm")
  expect_false(backend$raw_output$params$pairwise)
  expect_false(backend$raw_output$params$global)
  expect_false(backend$raw_output$params$struc_zero)
  expect_equal(backend$raw_output$backend_excluded_features, "F_other")
  expect_equal(backend$raw_output$backend_excluded_feature_count, 1L)
  expect_match(backend$raw_output$warnings, "fixture backend warning")
  expect_true(any(grepl(
    "fixture backend message",
    backend$raw_output$messages,
    fixed = TRUE
  )))
  expect_true("ancombc2_backend_warning_1" %in% backend$notes$caveat_id)
  expect_true(
    "ancombc2_backend_feature_exclusion" %in% backend$notes$caveat_id
  )

  da <- microeda:::da_build_result_object(
    context,
    list(ancombc2 = backend)
  )
  expect_s3_class(da, "microeda_da")
  expect_named(da$raw_outputs, "ancombc2")
  expect_named(da$method_results, "ancombc2")
  expect_equal(
    sum(da$caveats$caveat_id == "ancombc2_compositional_note"),
    1L
  )
})

test_that("ANCOM-BC2 works with DA summary, report, writer, and print", {
  context <- da_ancombc2_context()
  contrast_row <- context$contrast_plan[1, , drop = FALSE]

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_call_ancombc2 = da_mock_ancombc2_backend(),
    da_ancombc2_package_version = function() "2.14.0",
    .package = "microeda"
  )
  backend <- microeda:::da_run_ancombc2(
    context,
    contrast_row,
    params = list(p_adj_method = "holm")
  )
  da <- microeda:::da_build_result_object(
    context,
    list(ancombc2 = backend)
  )
  before <- as_da_results(da)

  summary <- as_da_summary(da)
  expect_equal(summary$method, "ancombc2")
  expect_equal(summary$contrast, "A_vs_B")
  expect_equal(summary$tested_features, ncol(context$counts))

  report <- microeda_da_report(da, top_n = 2)
  expect_match(report, "Methods: ancombc2", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_B", fixed = TRUE)
  expect_match(report, "natural-log estimates", fixed = TRUE)
  expect_match(report, "without additional adjustment", fixed = TRUE)
  expect_match(report, "exploratory/sensitivity", fixed = TRUE)

  output <- capture.output(visible <- withVisible(print(da)))
  expect_identical(visible$value, da)
  expect_false(visible$visible)
  expect_true(any(grepl("Methods:     ancombc2", output, fixed = TRUE)))

  file <- tempfile(fileext = ".csv")
  write_da_results(da, file)
  csv <- utils::read.csv(file, check.names = FALSE)
  expect_named(csv, da_ancombc2_expected_result_columns())
  expect_equal(nrow(csv), nrow(before))
  expect_false(any(grepl(
    "native_result",
    readLines(file, warn = FALSE),
    fixed = TRUE
  )))
  expect_equal(as_da_results(da), before)
})

test_that("real ANCOM-BC2 explicit backend preserves sign and schema", {
  skip_if_not_installed("ANCOMBC")

  counts <- da_ancombc2_counts()
  metadata <- da_ancombc2_metadata(rownames(counts))
  taxonomy <- da_ancombc2_taxonomy()
  set.seed(1)
  expect_silent(
    da <- microeda_da(
      counts,
      metadata = metadata,
      taxonomy = taxonomy,
      group = "group",
      contrast = c("A", "B"),
      methods = "ancombc2",
      tax_rank = "Genus",
      taxa_are_rows = FALSE,
      ancombc2_p_adj_method = "holm"
    )
  )

  results <- as_da_results(da)
  expect_s3_class(da, "microeda_da")
  expect_named(results, da_ancombc2_expected_result_columns())
  expect_equal(ncol(results), 21L)
  expect_equal(nrow(results), ncol(counts))
  expect_true(all(results$method == "ancombc2"))
  expect_true(all(results$contrast == "A_vs_B"))
  expect_lt(results$effect[results$feature_id == "F_A"], 0)
  expect_gt(results$effect[results$feature_id == "F_B"], 0)
  expect_equal(results$effect, results$log_fold_change)
  expect_equal(results$p_adjust_method, rep("holm", nrow(results)))
  expect_equal(results$p_adjust_scope, rep("method_contrast", nrow(results)))
  expect_identical(
    results$p_adjusted,
    da$raw_outputs$ancombc2$res[[
      paste0("q_", da$raw_outputs$ancombc2$coefficient)
    ]]
  )
  expect_identical(
    results$p_value,
    da$raw_outputs$ancombc2$res[[
      paste0("p_", da$raw_outputs$ancombc2$coefficient)
    ]]
  )
  expect_identical(
    results$standard_error,
    da$raw_outputs$ancombc2$res[[
      paste0("se_", da$raw_outputs$ancombc2$coefficient)
    ]]
  )
  expect_identical(
    results$statistic,
    da$raw_outputs$ancombc2$res[[
      paste0("W_", da$raw_outputs$ancombc2$coefficient)
    ]]
  )
  expect_identical(
    levels(da$raw_outputs$ancombc2$backend_conditions),
    c("A", "B")
  )
  expect_false(da$raw_outputs$ancombc2$params$pairwise)
  expect_null(da$raw_outputs$ancombc2$res_pair)

  set.seed(1)
  expect_silent(
    reversed <- microeda_da(
      counts,
      metadata = metadata,
      taxonomy = taxonomy,
      group = "group",
      contrast = c("B", "A"),
      methods = "ancombc2",
      tax_rank = "Genus",
      taxa_are_rows = FALSE,
      ancombc2_p_adj_method = "holm"
    )
  )
  reversed_results <- as_da_results(reversed)
  expect_gt(reversed_results$effect[reversed_results$feature_id == "F_A"], 0)
  expect_lt(reversed_results$effect[reversed_results$feature_id == "F_B"], 0)
  aligned_reverse <- reversed_results$effect[
    match(results$feature_id, reversed_results$feature_id)
  ]
  expect_equal(aligned_reverse, -results$effect, tolerance = 1e-08)
  expect_identical(
    levels(reversed$raw_outputs$ancombc2$backend_conditions),
    c("B", "A")
  )
})
