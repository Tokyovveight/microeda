da_performance_counts <- function() {
  counts <- matrix(
    c(
      20, 4, 8, 2,
      22, 5, 7, 3,
      5, 21, 9, 2,
      4, 23, 8, 3,
      8, 6, 24, 2,
      7, 5, 26, 3
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(6))
  colnames(counts) <- paste0("F", seq_len(4))
  counts
}

da_performance_metadata <- function() {
  data.frame(
    group = rep(c("A", "B", "C"), each = 2),
    row.names = paste0("S", seq_len(6))
  )
}

da_performance_standardized <- function(context, method, contrast_row) {
  feature_ids <- context$feature_ids
  adjusted <- c(1e-300, NA_real_, 0.04, 0.4)
  adjustment <- switch(
    method,
    aldex2 = "aldex2_native_BH",
    ancombc2 = "holm",
    deseq2 = "BH"
  )
  effect_type <- switch(
    method,
    aldex2 = "aldex2_effect",
    ancombc2 = "ancombc2_log_fold_change_group2_vs_group1",
    deseq2 = "deseq2_log2_fold_change_group2_vs_group1"
  )

  microeda:::da_standard_result(
    feature_id = feature_ids,
    taxon_label = c(
      paste(rep("Very long taxonomy label for report wrapping", 4), collapse = " "),
      "Genus_2",
      "Genus_3",
      "Genus_4"
    ),
    rank = "Genus",
    method = method,
    contrast = contrast_row$contrast,
    group1 = contrast_row$group1,
    group2 = contrast_row$group2,
    effect = c(2, -1, 0, 0.5),
    effect_type = effect_type,
    log_fold_change = if (identical(method, "aldex2")) {
      NA_real_
    } else {
      c(2, -1, 0, 0.5)
    },
    statistic = c(4, 3, 0, 1),
    standard_error = rep(0.2, 4),
    p_value = c(1e-300, NA_real_, 0.02, 0.3),
    p_adjusted = adjusted,
    p_adjust_method = adjustment,
    p_adjust_scope = "method_contrast",
    significance = microeda:::da_p_significance(adjusted),
    direction = NA_character_,
    method_note = microeda:::da_method_note(method)$message
  )
}

da_performance_raw <- function(context, method, contrast_row, params) {
  feature_ids <- context$feature_ids
  sample_ids <- context$sample_ids
  native_table <- data.frame(
    feature_id = feature_ids,
    effect = c(2, -1, 0, 0.5),
    p_value = c(1e-300, NA_real_, 0.02, 0.3),
    p_adjusted = c(1e-300, NA_real_, 0.04, 0.4),
    row.names = feature_ids,
    stringsAsFactors = FALSE
  )
  heavy <- matrix(seq_len(160000), nrow = 400)

  if (identical(method, "aldex2")) {
    return(list(
      clr = heavy,
      ttest = native_table[c("p_value", "p_adjusted")],
      effect = native_table["effect"],
      combined = native_table,
      conditions = rep(c("A", "B"), each = 2),
      backend_conditions = rep(c("microeda_group1", "microeda_group2"), each = 2),
      condition_mapping = data.frame(group = c("A", "B")),
      sample_order = head(sample_ids, 4),
      pair_id = NULL,
      pairing = NULL,
      contrast_row = contrast_row,
      input_orientation = "feature_by_sample",
      transposed_from_context = TRUE,
      params = params,
      warnings = "fixture warning",
      messages = "fixture message",
      package_version = "fixture"
    ))
  }

  if (identical(method, "ancombc2")) {
    return(list(
      native_result = list(res = native_table, heavy = heavy),
      res = native_table,
      res_global = heavy,
      res_pair = heavy,
      res_dunn = NULL,
      res_trend = NULL,
      zero_ind = NULL,
      ss_tab = data.frame(feature_id = feature_ids),
      coefficient = "microeda_groupmicroeda_group2",
      reference_level = contrast_row$group1,
      factor_levels = c(contrast_row$group1, contrast_row$group2),
      original_conditions = rep(c("A", "B"), each = 2),
      backend_conditions = factor(rep(c("A", "B"), each = 2)),
      sample_order = head(sample_ids, 4),
      feature_order = feature_ids,
      backend_excluded_features = character(),
      backend_excluded_feature_count = 0L,
      contrast_row = contrast_row,
      input_orientation = "feature_by_sample",
      transposed_from_context = TRUE,
      backend_group_column = "microeda_group",
      inert_metadata_column = "microeda_inert",
      params = params,
      warnings = "fixture warning",
      messages = "fixture message",
      package_version = "fixture"
    ))
  }

  list(
    native_dds = heavy,
    native_results = heavy,
    results_table = native_table,
    results_metadata = list(filterThreshold = 0.1),
    results_mcols = data.frame(type = "results"),
    dds_mcols = heavy,
    results_names = "microeda_group_B_vs_A",
    requested_contrast = c("microeda_group", "B", "A"),
    reference_level = contrast_row$group1,
    factor_levels = c(contrast_row$group1, contrast_row$group2),
    original_conditions = rep(c("A", "B"), each = 2),
    backend_conditions = factor(rep(c("A", "B"), each = 2)),
    sample_order = head(sample_ids, 4),
    feature_order = feature_ids,
    size_factors = rep(1, 4),
    normalization_factors = heavy,
    dispersions = rep(0.1, 4),
    design = ~microeda_group,
    backend_group_column = "microeda_group",
    backend_metadata = data.frame(microeda_group = rep(c("A", "B"), each = 2)),
    contrast_row = contrast_row,
    input_orientation = "feature_by_sample",
    transposed_from_context = TRUE,
    requested_fit_type = "parametric",
    actual_fit_type = "parametric",
    deseq_params = list(test = "Wald"),
    results_params = list(pAdjustMethod = "BH"),
    cooks = heavy,
    cooks_cutoff = 10,
    replacement_diagnostics = list(replaced = FALSE),
    feature_diagnostics = data.frame(
      feature_id = feature_ids,
      all_zero = FALSE,
      cooks_outlier = c(FALSE, TRUE, FALSE, FALSE),
      independently_filtered = c(FALSE, FALSE, TRUE, FALSE),
      other_na = FALSE
    ),
    warnings = "fixture warning",
    messages = "fixture message",
    params = params,
    package_version = "fixture"
  )
}

da_performance_backend <- function(context, method, params) {
  contrast_results <- lapply(seq_len(nrow(context$contrast_plan)), function(i) {
    contrast_row <- context$contrast_plan[i, , drop = FALSE]
    microeda:::da_execute_contrast(
      context = context,
      method = method,
      contrast_row = contrast_row,
      runner = function() {
        results <- da_performance_standardized(context, method, contrast_row)
        raw <- da_performance_raw(context, method, contrast_row, params)
        notes <- microeda:::da_method_notes(method)
        if (identical(method, "ancombc2")) {
          notes <- rbind(
            notes,
            microeda:::da_caveat(
              method,
              "fixture_warning",
              "backend",
              "warning",
              paste(rep("Long native warning retained for report wrapping", 5),
                    collapse = " ")
            )
          )
        }
        microeda:::da_backend_result(
          method = method,
          results = results,
          raw_output = raw,
          notes = notes,
          params = params
        )
      }
    )
  })

  results <- do.call(rbind, lapply(contrast_results, `[[`, "results"))
  row.names(results) <- NULL
  raw <- lapply(contrast_results, `[[`, "raw_output")
  names(raw) <- context$contrast_plan$contrast
  if (nrow(context$contrast_plan) == 1 &&
      identical(context$contrast_plan$contrast_type, "explicit")) {
    raw <- raw[[1]]
  } else if (identical(method, "aldex2")) {
    raw <- list(
      contrasts = raw,
      contrast_plan = context$contrast_plan,
      params = params,
      input_orientation = "feature_by_sample",
      transposed_from_context = TRUE
    )
  }

  notes <- do.call(rbind, lapply(contrast_results, `[[`, "notes"))
  notes <- microeda:::da_deduplicate_caveats(notes)
  microeda:::da_backend_result(
    method = method,
    results = results,
    raw_output = raw,
    notes = notes,
    params = params
  )
}

da_performance_mock_aldex2 <- function(context,
                                       mc.samples,
                                       denom,
                                       paired.test) {
  da_performance_backend(
    context,
    "aldex2",
    list(
      mc.samples = mc.samples,
      denom = denom,
      paired.test = paired.test,
      pair_id = context$pair_id
    )
  )
}

da_performance_mock_ancombc2 <- function(context,
                                         contrast_row = NULL,
                                         params) {
  da_performance_backend(context, "ancombc2", params)
}

da_performance_mock_deseq2 <- function(context,
                                       contrast_row = NULL,
                                       params) {
  da_performance_backend(context, "deseq2", params)
}

da_performance_run <- function(methods = "aldex2",
                               contrast = c("A", "B"),
                               raw_storage = "full",
                               progress = FALSE) {
  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_performance_mock_aldex2,
    da_run_ancombc2 = da_performance_mock_ancombc2,
    da_run_deseq2 = da_performance_mock_deseq2,
    .package = "microeda"
  )
  microeda_da(
    da_performance_counts(),
    metadata = da_performance_metadata(),
    group = "group",
    contrast = contrast,
    methods = methods,
    taxa_are_rows = FALSE,
    mc.samples = 16,
    raw_storage = raw_storage,
    progress = progress
  )
}

da_capture_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

test_that("raw storage modes preserve standardized results and reduce size", {
  full <- da_performance_run()
  compact <- da_performance_run(raw_storage = "compact")
  none <- da_performance_run(raw_storage = "none")

  expect_equal(as_da_results(full), as_da_results(compact))
  expect_equal(as_da_results(full), as_da_results(none))
  expect_equal(full$params$raw_storage, "full")
  expect_equal(compact$params$raw_storage, "compact")
  expect_equal(none$params$raw_storage, "none")

  expect_true("clr" %in% names(as_da_raw_output(full)))
  compact_raw <- as_da_raw_output(compact)
  expect_equal(compact_raw$raw_storage, "compact")
  expect_true(all(c("ttest", "effect", "combined", "timing") %in%
                  names(compact_raw)))
  expect_false("clr" %in% names(compact_raw))
  expect_null(none$raw_outputs$aldex2)
  expect_null(none$method_results$aldex2$raw_output)
  expect_error(as_da_raw_output(none), "raw_storage = \"none\"")

  sizes <- vapply(
    list(none = none, compact = compact, full = full),
    object.size,
    numeric(1)
  )
  expect_lt(sizes[["none"]], sizes[["compact"]])
  expect_lt(sizes[["compact"]], sizes[["full"]])

  files <- vapply(
    list(full, compact, none),
    function(x) {
      file <- tempfile(fileext = ".csv")
      write_da_results(x, file)
      file
    },
    character(1)
  )
  csv <- lapply(files, utils::read.csv, check.names = FALSE)
  expect_equal(csv[[1]], csv[[2]])
  expect_equal(csv[[1]], csv[[3]])

  expect_error(
    da_performance_run(raw_storage = "invalid"),
    "should be one of"
  )
  expect_error(da_performance_run(progress = NA), "progress")
})

test_that("compact storage retains diagnostics but omits heavy backend objects", {
  da <- da_performance_run(
    methods = c("aldex2", "ancombc2", "deseq2"),
    contrast = "pairwise",
    raw_storage = "compact"
  )

  for (method in da$methods) {
    for (contrast in da$contrast_plan$contrast) {
      raw <- as_da_raw_output(da, method = method, contrast = contrast)
      expect_equal(raw$raw_storage, "compact")
      expect_s3_class(raw$timing, "data.frame")
      expect_equal(raw$timing$contrast, contrast)
    }
  }

  aldex <- as_da_raw_output(da, "aldex2", "A_vs_B")
  expect_true(all(c("ttest", "effect", "combined", "conditions",
                    "sample_order", "params", "warnings", "messages") %in%
                  names(aldex)))
  expect_false("clr" %in% names(aldex))

  ancombc <- as_da_raw_output(da, "ancombc2", "A_vs_B")
  expect_true(all(c("res", "zero_ind", "ss_tab", "coefficient",
                    "backend_excluded_features", "params") %in%
                  names(ancombc)))
  expect_false(any(c("native_result", "res_global", "res_pair") %in%
                   names(ancombc)))

  deseq <- as_da_raw_output(da, "deseq2", "A_vs_B")
  expect_true(all(c("results_table", "results_metadata", "results_names",
                    "size_factors", "dispersions", "feature_diagnostics") %in%
                  names(deseq)))
  expect_false(any(c("native_dds", "native_results", "cooks",
                     "normalization_factors") %in% names(deseq)))
})

test_that("progress and timings follow method and contrast order", {
  quiet <- da_capture_messages(da_performance_run(progress = FALSE))
  expect_length(quiet$messages, 0L)

  captured <- da_capture_messages(da_performance_run(
    methods = c("deseq2", "aldex2", "ancombc2"),
    contrast = "pairwise",
    raw_storage = "none",
    progress = TRUE
  ))
  messages <- captured$messages
  expect_match(messages[[1]], "Method 1/3: deseq2", fixed = TRUE)
  expect_true(any(grepl("Contrast 1/3: A_vs_B - started", messages, fixed = TRUE)))
  expect_true(any(grepl("Contrast 3/3: B_vs_C - finished in", messages, fixed = TRUE)))
  expect_true(any(grepl("Method 2/3: aldex2", messages, fixed = TRUE)))
  expect_true(any(grepl("Method 3/3: ancombc2", messages, fixed = TRUE)))

  timings <- captured$value$params$timings
  expect_named(
    timings,
    c("method", "contrast", "elapsed_seconds", "status")
  )
  expect_equal(
    timings$method,
    rep(c("deseq2", "aldex2", "ancombc2"), each = 3)
  )
  expect_equal(
    timings$contrast,
    rep(c("A_vs_B", "A_vs_C", "B_vs_C"), 3)
  )
  expect_true(all(is.finite(timings$elapsed_seconds)))
  expect_true(all(timings$elapsed_seconds >= 0))
  expect_true(all(timings$status == "success"))
  expect_true(is.finite(captured$value$params$total_elapsed_seconds))
  expect_gte(captured$value$params$total_elapsed_seconds, 0)
})

test_that("DA reports use compact lines and preserve tiny probabilities", {
  da <- da_performance_run(
    methods = c("deseq2", "aldex2", "ancombc2"),
    raw_storage = "compact"
  )
  comparison <- microeda_da_comparison_report(
    da,
    alpha = 0.05,
    max_features = 2
  )
  report <- microeda_da_report(da, top_n = 2)

  comparison_lines <- strsplit(comparison, "\n", fixed = TRUE)[[1]]
  report_lines <- strsplit(report, "\n", fixed = TRUE)[[1]]
  expect_lte(max(nchar(comparison_lines)), 105L)
  expect_lte(max(nchar(report_lines)), 105L)
  expect_match(comparison, "Feature: F1", fixed = TRUE)
  expect_match(comparison, "Method    Effect", fixed = TRUE)
  expect_false(grepl("deseq2_effect\\s+aldex2_effect", comparison))
  expect_match(comparison, "display filter, not a consensus call", fixed = TRUE)
  expect_match(comparison, "FALSE = returned", fixed = TRUE)
  expect_match(comparison, "NA = feature\\s+absent")
  expect_false(grepl("winner|best method|votes", comparison, ignore.case = TRUE))
  expect_match(report, "1.000e-300", fixed = TRUE)

  formatted <- microeda:::da_report_probability(
    c(NA, 0, 1e-300, 1e-20, 0.0009, 0.001, 0.05, 1),
    digits = 3
  )
  expect_equal(
    formatted,
    c(
      "NA", "<2.2e-308", "1.000e-300", "1.000e-20",
      "9.000e-04", "0.001", "0.05", "1"
    )
  )
  expect_false(any(formatted[c(3, 4, 5)] == "0"))

  wide <- as_da_comparison(da)
  expect_equal(ncol(wide), 30L)
  expect_true(all(c("deseq2_effect", "aldex2_effect",
                    "ancombc2_effect") %in% names(wide)))
})

test_that("print shows storage and total timing without expanding details", {
  da <- da_performance_run(raw_storage = "compact")
  before <- as_da_results(da)
  output <- capture.output(visible <- withVisible(print(da)))

  expect_false(visible$visible)
  expect_identical(visible$value, da)
  expect_true(any(grepl("Raw storage: compact", output, fixed = TRUE)))
  expect_true(any(grepl("Elapsed:", output, fixed = TRUE)))
  expect_false(any(grepl("elapsed_seconds|A_vs_B.*success", output)))
  expect_equal(as_da_results(da), before)
})
