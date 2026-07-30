da_multimethod_counts <- function() {
  counts <- matrix(
    c(
      80, 10, 30,
      82, 11, 29,
      78, 9, 31,
      12, 75, 28,
      11, 78, 30,
      13, 72, 29,
      10, 12, 85,
      9, 11, 88,
      12, 10, 82
    ),
    nrow = 9,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(9))
  colnames(counts) <- paste0("F", seq_len(3))
  counts
}

da_multimethod_metadata <- function() {
  data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    pair = rep(paste0("P", seq_len(3)), 3),
    row.names = paste0("S", seq_len(9))
  )
}

da_multimethod_backend <- function(context, method, params) {
  contrast_results <- lapply(seq_len(nrow(context$contrast_plan)), function(i) {
    contrast_row <- context$contrast_plan[i, , drop = FALSE]
    n <- length(context$feature_ids)
    adjusted <- seq_len(n) / (n + 1)
    p_adjust_method <- if (identical(method, "aldex2")) {
      "aldex2_native_BH"
    } else {
      params$p_adj_method
    }

    microeda:::da_standard_result(
      feature_id = context$feature_ids,
      method = method,
      contrast = contrast_row$contrast,
      group1 = contrast_row$group1,
      group2 = contrast_row$group2,
      effect = seq_len(n),
      effect_type = paste0(method, "_fixture_effect"),
      p_value = adjusted / 2,
      p_adjusted = adjusted,
      p_adjust_method = p_adjust_method,
      p_adjust_scope = "method_contrast",
      significance = microeda:::da_p_significance(adjusted),
      method_note = microeda:::da_method_note(method)$message
    )
  })
  results <- do.call(rbind, contrast_results)
  row.names(results) <- NULL

  raw_contrasts <- lapply(seq_len(nrow(context$contrast_plan)), function(i) {
    list(
      method = method,
      contrast = context$contrast_plan$contrast[[i]],
      params = params,
      fixture_native_output = TRUE
    )
  })
  names(raw_contrasts) <- context$contrast_plan$contrast
  if (nrow(context$contrast_plan) == 1) {
    raw_output <- raw_contrasts[[1]]
  } else if (identical(method, "aldex2")) {
    raw_output <- list(
      contrasts = raw_contrasts,
      contrast_plan = context$contrast_plan,
      params = params,
      input_orientation = "feature_by_sample",
      transposed_from_context = TRUE
    )
  } else {
    raw_output <- raw_contrasts
  }

  microeda:::da_backend_result(
    method = method,
    results = results,
    raw_output = raw_output,
    notes = microeda:::da_method_notes(method),
    params = params
  )
}

da_multimethod_mock_aldex2 <- function(context,
                                       mc.samples,
                                       denom,
                                       paired.test) {
  da_multimethod_backend(
    context,
    method = "aldex2",
    params = list(
      mc.samples = as.integer(mc.samples),
      denom = denom,
      paired.test = paired.test,
      pair_id = context$pair_id
    )
  )
}

da_multimethod_mock_ancombc2 <- function(context,
                                         contrast_row = NULL,
                                         params) {
  da_multimethod_backend(
    context,
    method = "ancombc2",
    params = list(p_adj_method = params$p_adj_method)
  )
}

test_that("multi-method explicit dispatch preserves method order and contracts", {
  counts <- da_multimethod_counts()[seq_len(6), , drop = FALSE]
  metadata <- da_multimethod_metadata()[seq_len(6), , drop = FALSE]

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_multimethod_mock_aldex2,
    da_run_ancombc2 = da_multimethod_mock_ancombc2,
    .package = "microeda"
  )
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = c("A", "B"),
    methods = c("aldex2", "ancombc2"),
    taxa_are_rows = FALSE,
    mc.samples = 16,
    denom = "all",
    ancombc2_p_adj_method = "holm"
  )

  expect_s3_class(da, "microeda_da")
  expect_equal(da$methods, c("aldex2", "ancombc2"))
  expect_named(da$method_results, c("aldex2", "ancombc2"))
  expect_named(da$raw_outputs, c("aldex2", "ancombc2"))
  expect_named(da$params$backend, c("aldex2", "ancombc2"))
  expect_equal(
    as_da_results(da)$method,
    c(rep("aldex2", 3), rep("ancombc2", 3))
  )
  expect_true(all(as_da_results(da)$contrast == "A_vs_B"))
  expect_equal(da$params$backend$aldex2$mc.samples, 16L)
  expect_equal(da$params$backend$aldex2$denom, "all")
  expect_false("p_adj_method" %in% names(da$params$backend$aldex2))
  expect_equal(da$params$backend$ancombc2$p_adj_method, "holm")
  expect_false("mc.samples" %in% names(da$params$backend$ancombc2))
  expect_equal(
    sum(da$caveats$caveat_id == "aldex2_compositional_note"),
    1L
  )
  expect_equal(
    sum(da$caveats$caveat_id == "ancombc2_compositional_note"),
    1L
  )
})

test_that("multi-method pairwise rows follow methods then contrast plan", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_multimethod_mock_aldex2,
    da_run_ancombc2 = da_multimethod_mock_ancombc2,
    .package = "microeda"
  )
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = c("aldex2", "ancombc2"),
    taxa_are_rows = FALSE,
    mc.samples = 32,
    ancombc2_p_adj_method = "BH"
  )

  expected_contrasts <- c("A_vs_B", "A_vs_C", "B_vs_C")
  expected_methods <- rep(c("aldex2", "ancombc2"), each = 9)
  expected_rows <- rep(expected_contrasts, each = 3)
  expected_rows <- rep(expected_rows, 2)
  results <- as_da_results(da)

  expect_named(results, names(microeda:::da_empty_result()))
  expect_equal(ncol(results), 21L)
  expect_equal(results$method, expected_methods)
  expect_equal(results$contrast, expected_rows)
  expect_named(da$raw_outputs, c("aldex2", "ancombc2"))
  expect_named(da$raw_outputs$aldex2$contrasts, expected_contrasts)
  expect_named(da$raw_outputs$ancombc2, expected_contrasts)
  expect_named(da$method_results, c("aldex2", "ancombc2"))
  expect_equal(unique(results$p_adjust_scope), "method_contrast")
  expect_equal(
    unique(results$p_adjust_method[results$method == "aldex2"]),
    "aldex2_native_BH"
  )
  expect_equal(
    unique(results$p_adjust_method[results$method == "ancombc2"]),
    "BH"
  )
})

test_that("reverse multi-method order is preserved deterministically", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_multimethod_mock_aldex2,
    da_run_ancombc2 = da_multimethod_mock_ancombc2,
    .package = "microeda"
  )
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = c("ancombc2", "aldex2"),
    taxa_are_rows = FALSE,
    mc.samples = 16
  )

  expect_equal(da$methods, c("ancombc2", "aldex2"))
  expect_named(da$method_results, c("ancombc2", "aldex2"))
  expect_named(da$raw_outputs, c("ancombc2", "aldex2"))
  expect_equal(
    unique(as_da_results(da)$method),
    c("ancombc2", "aldex2")
  )
  expect_equal(
    unique(as_da_results(da)$contrast[
      as_da_results(da)$method == "ancombc2"
    ]),
    c("A_vs_B", "A_vs_C", "B_vs_C")
  )
})

test_that("paired restrictions are checked before multi-method backends run", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()
  calls <- character()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = function(...) {
      calls <<- c(calls, "aldex2")
      stop("should not run")
    },
    da_run_ancombc2 = function(...) {
      calls <<- c(calls, "ancombc2")
      stop("should not run")
    },
    .package = "microeda"
  )

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "pairwise",
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE,
      paired.test = TRUE,
      pair_id = "pair"
    ),
    "Paired ANCOM-BC2 execution is not implemented"
  )
  expect_length(calls, 0L)

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "pairwise",
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE,
      pair_id = "pair"
    ),
    "Paired ANCOM-BC2 execution is not implemented"
  )
  expect_length(calls, 0L)
})

test_that("paired ALDEx2-only dispatch remains available", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()
  captured <- NULL

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = function(context, mc.samples, denom, paired.test) {
      captured <<- list(
        paired.test = paired.test,
        pair_id = context$pair_id
      )
      da_multimethod_mock_aldex2(
        context,
        mc.samples,
        denom,
        paired.test
      )
    },
    .package = "microeda"
  )
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = "aldex2",
    taxa_are_rows = FALSE,
    paired.test = TRUE,
    pair_id = "pair",
    mc.samples = 16
  )

  expect_s3_class(da, "microeda_da")
  expect_true(captured$paired.test)
  expect_equal(captured$pair_id, "pair")
})

test_that("multi-method availability is checked before backend execution", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()
  calls <- character()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) {
      identical(package, "ALDEx2")
    },
    da_run_aldex2 = function(...) {
      calls <<- c(calls, "aldex2")
      stop("should not run")
    },
    da_run_ancombc2 = function(...) {
      calls <<- c(calls, "ancombc2")
      stop("should not run")
    },
    .package = "microeda"
  )

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "pairwise",
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE
    ),
    "requires the optional package `ANCOMBC`"
  )
  expect_length(calls, 0L)
})

test_that("backend failure does not return a partial multi-method object", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()
  calls <- character()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = function(context, mc.samples, denom, paired.test) {
      calls <<- c(calls, "aldex2")
      da_multimethod_mock_aldex2(
        context,
        mc.samples,
        denom,
        paired.test
      )
    },
    da_run_ancombc2 = function(...) {
      calls <<- c(calls, "ancombc2")
      stop("fixture ANCOM-BC2 failure", call. = FALSE)
    },
    .package = "microeda"
  )

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      contrast = "pairwise",
      methods = c("aldex2", "ancombc2"),
      taxa_are_rows = FALSE
    ),
    "fixture ANCOM-BC2 failure"
  )
  expect_equal(calls, c("aldex2", "ancombc2"))
})

test_that("multi-method helpers remain descriptive and raw-output safe", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()

  testthat::local_mocked_bindings(
    da_optional_package_available = function(package) TRUE,
    da_run_aldex2 = da_multimethod_mock_aldex2,
    da_run_ancombc2 = da_multimethod_mock_ancombc2,
    .package = "microeda"
  )
  da <- microeda_da(
    counts,
    metadata = metadata,
    group = "group",
    contrast = "pairwise",
    methods = c("aldex2", "ancombc2"),
    taxa_are_rows = FALSE,
    mc.samples = 16
  )
  before <- as_da_results(da)

  summary <- as_da_summary(da)
  expect_equal(nrow(summary), 6L)
  expect_equal(
    summary$method,
    rep(c("aldex2", "ancombc2"), each = 3)
  )
  expect_equal(
    summary$contrast,
    rep(c("A_vs_B", "A_vs_C", "B_vs_C"), 2)
  )
  expect_equal(summary$tested_features, rep(3L, 6))

  report <- microeda_da_report(da, top_n = 1)
  expect_match(report, "Method: aldex2", fixed = TRUE)
  expect_match(report, "Method: ancombc2", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_B", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_C", fixed = TRUE)
  expect_match(report, "Contrast: B_vs_C", fixed = TRUE)
  expect_match(report, "side-by-side sensitivity views", fixed = TRUE)
  expect_match(report, "does not create a consensus or rank", fixed = TRUE)
  expect_match(report, "aldex2_native_BH", fixed = TRUE)
  expect_match(report, "holm", fixed = TRUE)
  expect_match(report, "natural-log", fixed = TRUE)

  output <- capture.output(visible <- withVisible(print(da)))
  expect_identical(visible$value, da)
  expect_false(visible$visible)
  expect_true(any(grepl(
    "Methods:     aldex2, ancombc2",
    output,
    fixed = TRUE
  )))
  expect_true(any(grepl("Contrasts:   3 pairwise", output, fixed = TRUE)))

  file <- tempfile(fileext = ".csv")
  write_da_results(da, file)
  csv <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(csv$method, before$method)
  expect_equal(csv$contrast, before$contrast)
  expect_false(any(grepl(
    "fixture_native_output",
    readLines(file, warn = FALSE),
    fixed = TRUE
  )))
  expect_equal(as_da_results(da), before)
})

test_that("duplicate and unknown public methods still error clearly", {
  counts <- da_multimethod_counts()
  metadata <- da_multimethod_metadata()

  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      methods = c("aldex2", "aldex2"),
      taxa_are_rows = FALSE
    ),
    "duplicate"
  )
  expect_error(
    microeda_da(
      counts,
      metadata = metadata,
      group = "group",
      methods = c("aldex2", "unknown"),
      taxa_are_rows = FALSE
    ),
    "Unknown DA method"
  )
})
