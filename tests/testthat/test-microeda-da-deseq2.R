da_deseq2_fixture <- function(three_groups = TRUE) {
  set.seed(20260731)
  group_levels <- if (three_groups) c("A", "B", "C") else c("A", "B")
  groups <- rep(group_levels, each = 6)
  n_samples <- length(groups)
  n_features <- 220L
  sample_ids <- paste0("S", seq_len(n_samples))
  feature_ids <- paste0("F", seq_len(n_features))
  library_scale <- rep(c(0.75, 0.9, 1.05, 1.2, 1.35, 0.85),
                       length(group_levels))
  feature_means <- c(runif(150, 15, 90), runif(70, 0.2, 4))
  means <- outer(feature_means, library_scale)
  means[seq_len(8), groups == "A"] <- means[seq_len(8), groups == "A"] * 7
  means[9:16, groups == "B"] <- means[9:16, groups == "B"] * 7
  if (three_groups) {
    means[17:24, groups == "C"] <- means[17:24, groups == "C"] * 7
  }

  feature_by_sample <- matrix(
    stats::rnbinom(
      n_features * n_samples,
      mu = as.vector(means),
      size = 15
    ),
    nrow = n_features,
    dimnames = list(feature_ids, sample_ids)
  )
  feature_by_sample <- rbind(
    F_A_high = c(
      rep(450L, 6),
      rep(18L, 6),
      if (three_groups) rep(20L, 6) else integer()
    ),
    F_B_high = c(
      rep(18L, 6),
      rep(450L, 6),
      if (three_groups) rep(20L, 6) else integer()
    ),
    F_C_high = c(
      rep(20L, 12),
      if (three_groups) rep(450L, 6) else integer()
    ),
    F_stable = rep(c(95L, 100L, 105L, 98L, 102L, 97L),
                   length(group_levels)),
    F_all_zero = rep(0L, n_samples),
    feature_by_sample
  )
  storage.mode(feature_by_sample) <- "integer"
  counts <- t(feature_by_sample)
  metadata <- data.frame(
    group = groups,
    microeda_deseq2_group = rep("user_column", n_samples),
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  taxonomy <- data.frame(
    Genus = paste("Genus", rownames(feature_by_sample)),
    row.names = rownames(feature_by_sample),
    stringsAsFactors = FALSE
  )

  list(counts = counts, metadata = metadata, taxonomy = taxonomy)
}

da_deseq2_context <- function(contrast = c("A", "B"),
                              methods = "deseq2",
                              three_groups = identical(contrast, "pairwise")) {
  fixture <- da_deseq2_fixture(three_groups = three_groups)
  microeda:::da_prepare_context(
    fixture$counts,
    metadata = fixture$metadata,
    taxonomy = fixture$taxonomy,
    group = "group",
    contrast = contrast,
    methods = methods,
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
}

da_deseq2_params <- function() {
  list(
    sf_type = "poscounts",
    fit_type = "parametric",
    independent_filtering = TRUE,
    alpha = 0.05
  )
}

da_deseq2_fake_backend <- function(context, method, params) {
  rows <- lapply(seq_len(nrow(context$contrast_plan)), function(i) {
    contrast_row <- context$contrast_plan[i, , drop = FALSE]
    feature_ids <- head(context$feature_ids, 2)
    adjusted <- c(0.01, NA_real_)
    microeda:::da_standard_result(
      feature_id = feature_ids,
      method = method,
      contrast = contrast_row$contrast,
      group1 = contrast_row$group1,
      group2 = contrast_row$group2,
      effect = c(-1, 1),
      effect_type = paste0(method, "_fixture"),
      p_value = c(0.005, NA_real_),
      p_adjusted = adjusted,
      p_adjust_method = if (identical(method, "aldex2")) {
        "aldex2_native_BH"
      } else if (identical(method, "ancombc2")) {
        params$p_adj_method
      } else {
        "BH"
      },
      p_adjust_scope = "method_contrast",
      significance = microeda:::da_p_significance(adjusted),
      method_note = microeda:::da_method_note(method)$message
    )
  })
  results <- do.call(rbind, rows)
  row.names(results) <- NULL
  raw <- lapply(context$contrast_plan$contrast, function(contrast) {
    list(method = method, contrast = contrast, fixture = TRUE)
  })
  names(raw) <- context$contrast_plan$contrast
  if (nrow(context$contrast_plan) == 1) {
    raw <- raw[[1]]
  }

  microeda:::da_backend_result(
    method = method,
    results = results,
    raw_output = raw,
    notes = microeda:::da_method_notes(method),
    params = params
  )
}

da_deseq2_mock_runner <- function(context, contrast_row = NULL, params) {
  da_deseq2_fake_backend(context, "deseq2", params)
}

da_deseq2_mock_aldex2 <- function(context,
                                  mc.samples,
                                  denom,
                                  paired.test) {
  da_deseq2_fake_backend(
    context,
    "aldex2",
    list(
      mc.samples = mc.samples,
      denom = denom,
      paired.test = paired.test
    )
  )
}

da_deseq2_mock_ancombc2 <- function(context,
                                    contrast_row = NULL,
                                    params) {
  da_deseq2_fake_backend(context, "ancombc2", params)
}

test_that("DESeq2 parameters and public boundaries validate clearly", {
  expect_equal(
    microeda:::da_validate_deseq2_params(da_deseq2_params()),
    da_deseq2_params()
  )
  bad <- da_deseq2_params()
  bad$sf_type <- "bad"
  expect_error(
    microeda:::da_validate_deseq2_params(bad),
    "deseq2_sf_type"
  )
  bad <- da_deseq2_params()
  bad$fit_type <- "bad"
  expect_error(
    microeda:::da_validate_deseq2_params(bad),
    "deseq2_fit_type"
  )
  bad <- da_deseq2_params()
  bad$independent_filtering <- NA
  expect_error(
    microeda:::da_validate_deseq2_params(bad),
    "deseq2_independent_filtering"
  )
  bad <- da_deseq2_params()
  bad$alpha <- 1
  expect_error(
    microeda:::da_validate_deseq2_params(bad),
    "strictly between"
  )

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) FALSE,
    .package = "microeda"
  )
  bad <- da_deseq2_params()
  bad$fit_type <- "glmGamPoi"
  expect_error(
    microeda:::da_validate_deseq2_params(bad),
    "requires the optional package"
  )
})

test_that("DESeq2 count validation rejects unsafe input without rounding", {
  counts <- matrix(
    c(10, 5, 8, 7, 6, 9),
    nrow = 3,
    dimnames = list(paste0("S", 1:3), c("F1", "F2"))
  )
  expect_invisible(microeda:::da_validate_deseq2_counts(counts))

  fractional <- counts
  fractional[1, 1] <- 1.5
  expect_error(
    microeda:::da_validate_deseq2_counts(fractional),
    "does not round"
  )
  negative <- counts
  negative[1, 1] <- -1
  expect_error(
    microeda:::da_validate_deseq2_counts(negative),
    "negative"
  )
  non_finite <- counts
  non_finite[1, 1] <- Inf
  expect_error(
    microeda:::da_validate_deseq2_counts(non_finite),
    "finite"
  )
  zero_library <- counts
  zero_library[1, ] <- 0
  expect_error(
    microeda:::da_validate_deseq2_counts(zero_library),
    "zero-library sample"
  )
  all_zero_feature <- cbind(counts, F_zero = 0)
  expect_invisible(
    microeda:::da_validate_deseq2_counts(all_zero_feature)
  )
})

test_that("DESeq2 input preparation aligns metadata and uses safe names", {
  context <- da_deseq2_context()
  input <- microeda:::da_prepare_deseq2_input(
    context,
    context$contrast_plan[1, , drop = FALSE]
  )

  expect_identical(colnames(input$counts), input$sample_order)
  expect_identical(rownames(input$metadata), input$sample_order)
  expect_identical(
    levels(input$metadata[[input$backend_group_column]]),
    c("A", "B")
  )
  expect_identical(
    as.character(input$design),
    c("~", "microeda_deseq2_group_1")
  )
  expect_equal(dim(input$counts), c(length(context$feature_ids), 12L))
  expect_identical(rownames(input$counts), context$feature_ids)
  expect_false(input$backend_group_column %in% colnames(context$metadata))

  duplicate_samples <- context
  duplicate_samples$sample_ids[2] <- duplicate_samples$sample_ids[1]
  expect_error(
    microeda:::da_prepare_deseq2_input(
      duplicate_samples,
      duplicate_samples$contrast_plan[1, , drop = FALSE]
    ),
    "unique"
  )
  duplicate_features <- context
  duplicate_features$feature_ids[2] <- duplicate_features$feature_ids[1]
  expect_error(
    microeda:::da_prepare_deseq2_input(
      duplicate_features,
      duplicate_features$contrast_plan[1, , drop = FALSE]
    ),
    "unique"
  )
  misaligned <- context
  misaligned$metadata <- misaligned$metadata[rev(rownames(misaligned$metadata)), ,
                                             drop = FALSE]
  expect_error(
    microeda:::da_prepare_deseq2_input(
      misaligned,
      misaligned$contrast_plan[1, , drop = FALSE]
    ),
    "aligned"
  )
})

test_that("DESeq2 standardizer preserves the 21-column mapping and native NA", {
  counts <- matrix(
    c(40, 20, 10, 42, 18, 11, 12, 35, 9, 10, 38, 8),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("S", 1:4), c("F1", "F2", "F3"))
  )
  metadata <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = rownames(counts)
  )
  taxonomy <- data.frame(
    Genus = c("G1", "G2", "G3"),
    row.names = colnames(counts)
  )
  context <- microeda:::da_prepare_context(
    counts,
    metadata = metadata,
    taxonomy = taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "deseq2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )
  native <- data.frame(
    baseMean = c(20, 25, 0),
    log2FoldChange = c(-2, 3, NA_real_),
    lfcSE = c(0.2, 0.3, NA_real_),
    stat = c(-10, 10, NA_real_),
    pvalue = c(0.001, 0.002, NA_real_),
    padj = c(0.003, 0.004, NA_real_),
    row.names = c("F1", "F2", "F3")
  )

  result <- microeda:::da_standardize_deseq2_result(
    results_table = native,
    context = context,
    contrast_row = context$contrast_plan[1, , drop = FALSE],
    feature_order = context$feature_ids
  )

  expect_named(result, names(microeda:::da_empty_result()))
  expect_equal(ncol(result), 21L)
  expect_equal(result$feature_id, c("F1", "F2", "F3"))
  expect_equal(result$taxon_label, c("G1", "G2", "G3"))
  expect_equal(result$method, rep("deseq2", 3))
  expect_equal(result$effect, native$log2FoldChange)
  expect_equal(result$log_fold_change, native$log2FoldChange)
  expect_equal(result$statistic, native$stat)
  expect_equal(result$standard_error, native$lfcSE)
  expect_equal(result$p_value, native$pvalue)
  expect_equal(result$p_adjusted, native$padj)
  expect_true(all(is.na(result$ci_low)))
  expect_true(all(is.na(result$ci_high)))
  expect_true(all(is.na(result$direction)))
  expect_identical(result$p_adjust_method, rep("BH", 3))
  expect_identical(
    result$effect_type,
    rep("deseq2_log2_fold_change_group2_vs_group1", 3)
  )
  expect_true(is.na(result$significance[[3]]))
})

test_that("DESeq2 feature diagnostics distinguish native missing-value classes", {
  results <- data.frame(
    pvalue = c(NA, NA, 0.2, NA),
    padj = c(NA, NA, NA, NA),
    row.names = c("all_zero", "cook", "filtered", "other")
  )
  dds_mcols <- data.frame(
    allZero = c(TRUE, FALSE, FALSE, FALSE),
    maxCooks = c(NA, 100, 1, 1)
  )
  diagnostics <- microeda:::da_deseq2_feature_diagnostics(
    results_table = results,
    dds_mcols = dds_mcols,
    cooks_cutoff = 10,
    independent_filtering = TRUE
  )

  expect_equal(which(diagnostics$all_zero), 1L)
  expect_equal(which(diagnostics$cooks_outlier), 2L)
  expect_equal(which(diagnostics$independently_filtered), 3L)
  expect_equal(which(diagnostics$other_na), 4L)
  expect_named(
    diagnostics,
    c(
      "feature_id",
      "all_zero",
      "cooks_outlier",
      "independently_filtered",
      "other_na",
      "p_value_missing",
      "p_adjusted_missing",
      "max_cooks"
    )
  )
})

test_that("DESeq2 pairwise runner preserves plan order and aborts on failure", {
  context <- da_deseq2_context(contrast = "pairwise")
  calls <- character()
  fake_contrast <- function(context, contrast_row, params) {
    calls <<- c(calls, contrast_row$contrast)
    if (identical(contrast_row$contrast, "A_vs_C") &&
        isTRUE(getOption("microeda.deseq2.fail", FALSE))) {
      stop("fixture failure", call. = FALSE)
    }
    result <- microeda:::da_standard_result(
      feature_id = head(context$feature_ids, 2),
      method = "deseq2",
      contrast = contrast_row$contrast,
      group1 = contrast_row$group1,
      group2 = contrast_row$group2,
      effect = c(-1, 1),
      p_value = c(0.1, 0.2),
      p_adjusted = c(0.2, 0.3),
      p_adjust_method = "BH",
      p_adjust_scope = "method_contrast"
    )
    microeda:::da_backend_result(
      method = "deseq2",
      results = result,
      raw_output = list(contrast = contrast_row$contrast),
      params = params
    )
  }

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_deseq2_contrast = fake_contrast,
    .package = "microeda"
  )
  backend <- microeda:::da_run_deseq2(context, params = da_deseq2_params())
  expected <- c("A_vs_B", "A_vs_C", "B_vs_C")
  expect_s3_class(backend, "microeda_da_backend_result")
  expect_equal(unique(backend$results$contrast), expected)
  expect_named(backend$raw_output, expected)
  expect_equal(calls, expected)

  calls <- character()
  old_options <- options(microeda.deseq2.fail = TRUE)
  on.exit(options(old_options), add = TRUE)
  expect_error(
    microeda:::da_run_deseq2(context, params = da_deseq2_params()),
    "DESeq2 pairwise execution failed for contrast `A_vs_C`"
  )
  expect_equal(calls, c("A_vs_B", "A_vs_C"))
})

test_that("public DESeq2 dispatch validates availability and paired boundaries", {
  fixture <- da_deseq2_fixture(three_groups = FALSE)
  calls <- 0L
  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) {
      !identical(package, "DESeq2")
    },
    da_run_aldex2 = function(...) {
      calls <<- calls + 1L
      stop("should not run")
    },
    .package = "microeda"
  )
  expect_error(
    microeda_da(
      fixture$counts,
      metadata = fixture$metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = c("aldex2", "deseq2"),
      taxa_are_rows = FALSE
    ),
    "BiocManager::install\\(\"DESeq2\"\\)"
  )
  expect_equal(calls, 0L)

  expect_error(
    microeda_da(
      fixture$counts,
      metadata = fixture$metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "deseq2",
      taxa_are_rows = FALSE,
      paired.test = TRUE
    ),
    "Paired/repeated DESeq2"
  )
  expect_error(
    microeda_da(
      fixture$counts,
      metadata = fixture$metadata,
      group = "group",
      contrast = c("A", "B"),
      methods = "deseq2",
      taxa_are_rows = FALSE,
      pair_id = "pair"
    ),
    "Paired/repeated DESeq2"
  )
})

test_that("ordered dispatch supports DESeq2 method combinations", {
  fixture <- da_deseq2_fixture()
  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_deseq2_mock_aldex2,
    da_run_ancombc2 = da_deseq2_mock_ancombc2,
    da_run_deseq2 = da_deseq2_mock_runner,
    .package = "microeda"
  )
  combinations <- list(
    "deseq2",
    c("aldex2", "deseq2"),
    c("ancombc2", "deseq2"),
    c("aldex2", "ancombc2", "deseq2"),
    c("deseq2", "ancombc2", "aldex2")
  )

  for (methods in combinations) {
    da <- microeda_da(
      fixture$counts,
      metadata = fixture$metadata,
      group = "group",
      contrast = "pairwise",
      methods = methods,
      taxa_are_rows = FALSE
    )
    expect_equal(da$methods, methods)
    expect_named(da$method_results, methods)
    expect_named(da$raw_outputs, methods)
    expect_equal(unique(as_da_results(da)$method), methods)
    expect_equal(
      unique(as_da_results(da)$contrast),
      c("A_vs_B", "A_vs_C", "B_vs_C")
    )
  }

  da <- microeda_da(
    fixture$counts,
    metadata = fixture$metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = c("aldex2", "ancombc2", "deseq2"),
    taxa_are_rows = FALSE
  )
  expect_named(da$params$backend, c("aldex2", "ancombc2", "deseq2"))
  expect_false("sf_type" %in% names(da$params$backend$aldex2))
  expect_false("sf_type" %in% names(da$params$backend$ancombc2))
  expect_equal(da$params$backend$deseq2$sf_type, "poscounts")
})

test_that("real explicit DESeq2 backend preserves native contracts", {
  skip_if_not_installed("DESeq2")
  fixture <- da_deseq2_fixture(three_groups = FALSE)
  before_counts <- fixture$counts
  set.seed(1)
  da <- microeda_da(
    fixture$counts,
    metadata = fixture$metadata,
    taxonomy = fixture$taxonomy,
    group = "group",
    contrast = c("A", "B"),
    methods = "deseq2",
    tax_rank = "Genus",
    taxa_are_rows = FALSE
  )

  results <- as_da_results(da)
  raw <- da$raw_outputs$deseq2
  expect_s3_class(da, "microeda_da")
  expect_named(results, names(microeda:::da_empty_result()))
  expect_equal(ncol(results), 21L)
  expect_equal(nrow(results), ncol(fixture$counts))
  expect_lt(results$effect[results$feature_id == "F_A_high"], 0)
  expect_gt(results$effect[results$feature_id == "F_B_high"], 0)
  expect_equal(results$taxon_label, fixture$taxonomy[results$feature_id, "Genus"])
  expect_equal(results$p_value, raw$results_table$pvalue)
  expect_equal(results$p_adjusted, raw$results_table$padj)
  expect_equal(results$standard_error, raw$results_table$lfcSE)
  expect_equal(results$statistic, raw$results_table$stat)
  expect_true(is.na(results$p_value[results$feature_id == "F_all_zero"]))
  expect_equal(fixture$counts, before_counts)

  required_raw <- c(
    "native_dds",
    "native_results",
    "results_table",
    "results_metadata",
    "results_mcols",
    "dds_mcols",
    "results_names",
    "requested_contrast",
    "reference_level",
    "factor_levels",
    "sample_order",
    "feature_order",
    "size_factors",
    "normalization_factors",
    "dispersions",
    "design",
    "requested_fit_type",
    "actual_fit_type",
    "deseq_params",
    "results_params",
    "cooks",
    "cooks_cutoff",
    "replacement_diagnostics",
    "feature_diagnostics",
    "warnings",
    "messages",
    "package_version"
  )
  expect_true(all(required_raw %in% names(raw)))
  expect_s4_class(raw$native_dds, "DESeqDataSet")
  expect_s4_class(raw$native_results, "DESeqResults")
  expect_equal(raw$requested_contrast, c(
    raw$backend_group_column,
    "B",
    "A"
  ))
  expect_equal(raw$reference_level, "A")
  expect_equal(raw$factor_levels, c("A", "B"))
  expect_false(raw$replacement_diagnostics$replace_counts_present)
  expect_true(raw$feature_diagnostics$all_zero[
    raw$feature_diagnostics$feature_id == "F_all_zero"
  ])
  expect_identical(raw$deseq_params$minReplicatesForReplace, Inf)
  expect_identical(raw$results_params$pAdjustMethod, "BH")
  expect_identical(raw$results_params$independentFiltering, TRUE)
  expect_equal(raw$package_version, as.character(packageVersion("DESeq2")))

  reversed <- DESeq2::results(
    raw$native_dds,
    contrast = c(raw$backend_group_column, "A", "B"),
    lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    independentFiltering = TRUE,
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  expect_equal(
    results$effect,
    -as.numeric(reversed$log2FoldChange),
    tolerance = 1e-8
  )
  expect_equal(results$p_value, as.numeric(reversed$pvalue))
  expect_equal(results$p_adjusted, as.numeric(reversed$padj))

  summary <- as_da_summary(da)
  expect_equal(summary$tested_features, nrow(results))
  expect_equal(summary$n_p_value, sum(!is.na(results$p_value)))
  expect_equal(summary$n_p_adjusted, sum(!is.na(results$p_adjusted)))
  report <- microeda_da_report(da, top_n = 2)
  expect_match(report, "DESeq2 coefficients are unshrunk log2")
  expect_match(report, "independent filtering")
  expect_match(report, "sensitivity/comparison")
  output <- capture.output(visible <- withVisible(print(da)))
  expect_false(visible$visible)
  expect_identical(visible$value, da)
  expect_true(any(grepl("deseq2", output, fixed = TRUE)))

  file <- tempfile(fileext = ".csv")
  write_da_results(da, file)
  csv <- utils::read.csv(file, check.names = FALSE, na.strings = "")
  expect_named(csv, names(results))
  expect_equal(nrow(csv), nrow(results))
  expect_false(any(grepl(
    "native_dds",
    readLines(file, warn = FALSE),
    fixed = TRUE
  )))
})

test_that("real pairwise DESeq2 refits every contrast separately", {
  skip_if_not_installed("DESeq2")
  fixture <- da_deseq2_fixture(three_groups = TRUE)
  set.seed(2)
  da <- microeda_da(
    fixture$counts,
    metadata = fixture$metadata,
    group = "group",
    contrast = "pairwise",
    methods = "deseq2",
    taxa_are_rows = FALSE
  )

  expected <- c("A_vs_B", "A_vs_C", "B_vs_C")
  results <- as_da_results(da)
  expect_equal(unique(results$contrast), expected)
  expect_equal(
    as.integer(table(factor(results$contrast, levels = expected))),
    rep(ncol(fixture$counts), 3)
  )
  expect_named(da$raw_outputs$deseq2, expected)
  expect_equal(
    unname(vapply(
      da$raw_outputs$deseq2,
      function(raw) ncol(raw$native_dds),
      integer(1)
    )),
    rep(12L, 3)
  )
  expect_equal(
    unname(lapply(
      da$raw_outputs$deseq2,
      function(raw) raw$factor_levels
    )),
    list(c("A", "B"), c("A", "C"), c("B", "C"))
  )
  expect_false(identical(
    da$raw_outputs$deseq2$A_vs_B$size_factors,
    da$raw_outputs$deseq2$A_vs_C$size_factors
  ))
  expect_false(identical(
    da$raw_outputs$deseq2$A_vs_B$dispersions,
    da$raw_outputs$deseq2$A_vs_C$dispersions
  ))
  expect_lt(
    results$effect[
      results$contrast == "A_vs_B" & results$feature_id == "F_A_high"
    ],
    0
  )
  expect_gt(
    results$effect[
      results$contrast == "A_vs_B" & results$feature_id == "F_B_high"
    ],
    0
  )
  expect_equal(nrow(as_da_summary(da)), 3L)
})

test_that("DESeq2 backend remains internal and adds no p.adjust call", {
  exports <- getNamespaceExports("microeda")
  internal <- c(
    "da_run_deseq2",
    "da_run_deseq2_contrast",
    "da_prepare_deseq2_input",
    "da_standardize_deseq2_result",
    "da_deseq2_feature_diagnostics"
  )
  expect_false(any(internal %in% exports))

  code <- unlist(lapply(internal, function(name) {
    deparse(get(name, envir = asNamespace("microeda")))
  }))
  expect_false(any(grepl("p\\.adjust\\s*\\(", code)))

  description <- utils::packageDescription("microeda")
  imports <- if (is.null(description$Imports)) "" else description$Imports
  suggests <- if (is.null(description$Suggests)) "" else description$Suggests
  expect_false(grepl("\\bDESeq2\\b", imports))
  expect_true(grepl("\\bDESeq2\\b", suggests))
})
