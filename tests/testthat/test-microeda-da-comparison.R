da_comparison_result_rows <- function(method,
                                      contrast,
                                      group2,
                                      feature_id,
                                      effect,
                                      p_value,
                                      p_adjusted) {
  taxon <- paste0("Genus_", substring(feature_id, 2))
  effect_type <- switch(
    method,
    aldex2 = "aldex2_effect",
    ancombc2 = "ancombc2_log_fold_change_group2_vs_group1",
    deseq2 = "deseq2_log2_fold_change_group2_vs_group1"
  )
  adjustment <- switch(
    method,
    aldex2 = "aldex2_native_BH",
    ancombc2 = "holm",
    deseq2 = "BH"
  )

  microeda:::da_standard_result(
    feature_id = feature_id,
    taxon_label = taxon,
    rank = "Genus",
    method = method,
    contrast = contrast,
    group1 = "A",
    group2 = group2,
    effect = effect,
    effect_type = effect_type,
    log_fold_change = if (identical(method, "aldex2")) {
      NA_real_
    } else {
      effect
    },
    statistic = seq_along(feature_id),
    standard_error = rep(0.2, length(feature_id)),
    p_value = p_value,
    p_adjusted = p_adjusted,
    p_adjust_method = adjustment,
    p_adjust_scope = "method_contrast",
    significance = microeda:::da_p_significance(p_adjusted),
    direction = NA_character_,
    method_note = microeda:::da_method_note(method)$message
  )
}

da_comparison_fixture <- function(pairwise = TRUE,
                                  methods = c(
                                    "deseq2",
                                    "aldex2",
                                    "ancombc2"
                                  )) {
  plan <- if (pairwise) {
    data.frame(
      contrast = c("A_vs_B", "A_vs_C"),
      group1 = c("A", "A"),
      group2 = c("B", "C"),
      contrast_type = rep("pairwise", 2),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      contrast = "A_vs_B",
      group1 = "A",
      group2 = "B",
      contrast_type = "explicit",
      stringsAsFactors = FALSE
    )
  }

  specs <- list(
    A_vs_B = list(
      deseq2 = list(
        feature = c("F1", "F2"),
        effect = c(2, 0),
        p = c(0.001, 0.2),
        padj = c(0.003, NA_real_)
      ),
      aldex2 = list(
        feature = c("F1", "F2"),
        effect = c(1.5, -0.5),
        p = c(0.01, NA_real_),
        padj = c(0.02, NA_real_)
      ),
      ancombc2 = list(
        feature = c("F1", "F3"),
        effect = c(0.7, -1.1),
        p = c(0.03, 0.02),
        padj = c(0.06, 0.04)
      )
    ),
    A_vs_C = list(
      deseq2 = list(
        feature = c("F1", "F3"),
        effect = c(-2.2, 0.4),
        p = c(0.004, 0.5),
        padj = c(0.008, 0.6)
      ),
      aldex2 = list(
        feature = c("F1", "F2"),
        effect = c(-1.4, 0.2),
        p = c(0.02, 0.7),
        padj = c(0.04, 0.8)
      ),
      ancombc2 = list(
        feature = c("F2", "F3"),
        effect = c(0, 1.2),
        p = c(0.9, 0.03),
        padj = c(0.9, 0.07)
      )
    )
  )

  result_rows <- list()
  for (method in methods) {
    for (i in seq_len(nrow(plan))) {
      contrast <- plan$contrast[[i]]
      spec <- specs[[contrast]][[method]]
      result_rows[[length(result_rows) + 1L]] <-
        da_comparison_result_rows(
          method = method,
          contrast = contrast,
          group2 = plan$group2[[i]],
          feature_id = spec$feature,
          effect = spec$effect,
          p_value = spec$p,
          p_adjusted = spec$padj
        )
    }
  }
  results <- do.call(rbind, result_rows)
  row.names(results) <- NULL

  raw_by_method <- lapply(methods, function(method) {
    contrast_raw <- lapply(seq_len(nrow(plan)), function(i) {
      contrast <- plan$contrast[[i]]
      raw <- list(
        marker = paste(method, contrast, sep = ":"),
        contrast_row = plan[i, , drop = FALSE]
      )
      if (identical(method, "deseq2")) {
        raw$feature_diagnostics <- data.frame(
          feature_id = c("F1", "F2", "F3"),
          all_zero = c(FALSE, TRUE, FALSE),
          cooks_outlier = c(FALSE, FALSE, TRUE),
          independently_filtered = c(FALSE, TRUE, FALSE),
          other_na = c(FALSE, FALSE, FALSE),
          stringsAsFactors = FALSE
        )
      }
      if (identical(method, "ancombc2")) {
        raw$backend_excluded_feature_count <- 1L
        raw$backend_excluded_features <- "F4"
      }
      raw
    })
    names(contrast_raw) <- plan$contrast

    if (!pairwise) {
      return(contrast_raw[[1]])
    }
    if (identical(method, "aldex2")) {
      return(list(
        contrasts = contrast_raw,
        contrast_plan = plan
      ))
    }
    contrast_raw
  })
  names(raw_by_method) <- methods

  caveats <- do.call(
    rbind,
    lapply(methods, microeda:::da_method_note)
  )
  row.names(caveats) <- NULL

  structure(
    list(
      results = results,
      method_results = setNames(vector("list", length(methods)), methods),
      raw_outputs = raw_by_method,
      methods = methods,
      group = "group",
      contrast = if (pairwise) "pairwise" else c("A", "B"),
      contrast_plan = plan,
      contrast_label = if (pairwise) "pairwise" else "A_vs_B",
      tax_rank = "Genus",
      feature_metadata = data.frame(
        feature_id = c("F1", "F2", "F3", "F4"),
        Genus = paste0("Genus_", seq_len(4)),
        stringsAsFactors = FALSE
      ),
      filters = list(applied = FALSE),
      caveats = caveats,
      params = list(),
      call = quote(microeda_da())
    ),
    class = "microeda_da"
  )
}

test_that("as_da_raw_output infers only unambiguous selections", {
  for (method in c("aldex2", "ancombc2", "deseq2")) {
    explicit <- da_comparison_fixture(
      pairwise = FALSE,
      methods = method
    )
    expect_identical(
      as_da_raw_output(explicit),
      explicit$raw_outputs[[method]]
    )
  }

  multi <- da_comparison_fixture(pairwise = FALSE)
  expect_error(
    as_da_raw_output(multi),
    "method.*Available methods"
  )
  expect_identical(
    as_da_raw_output(multi, method = "deseq2"),
    multi$raw_outputs$deseq2
  )

  pairwise <- da_comparison_fixture(pairwise = TRUE, methods = "deseq2")
  expect_error(
    as_da_raw_output(pairwise),
    "contrast.*Available contrasts"
  )
})

test_that("as_da_raw_output adapts all pairwise backend nestings", {
  da <- da_comparison_fixture()

  expect_identical(
    as_da_raw_output(da, "aldex2", "A_vs_C"),
    da$raw_outputs$aldex2$contrasts$A_vs_C
  )
  expect_identical(
    as_da_raw_output(da, "ancombc2", "A_vs_C"),
    da$raw_outputs$ancombc2$A_vs_C
  )
  expect_identical(
    as_da_raw_output(da, "deseq2", "A_vs_C"),
    da$raw_outputs$deseq2$A_vs_C
  )

  expect_error(as_da_raw_output(list()), "microeda_da")
  expect_error(as_da_raw_output(da, "unknown", "A_vs_B"), "Available methods")
  expect_error(as_da_raw_output(da, "aldex2", "unknown"), "Available contrasts")
  expect_error(as_da_raw_output(da, c("aldex2", "deseq2")), "single")
})

test_that("as_da_raw_output rejects malformed nesting clearly", {
  aldex <- da_comparison_fixture(methods = "aldex2")
  aldex$raw_outputs$aldex2$contrasts$A_vs_B <- NULL
  expect_error(
    as_da_raw_output(aldex, contrast = "A_vs_B"),
    "Malformed ALDEx2"
  )

  deseq <- da_comparison_fixture(methods = "deseq2")
  deseq$raw_outputs$deseq2$A_vs_B <- NULL
  expect_error(
    as_da_raw_output(deseq, contrast = "A_vs_B"),
    "Malformed deseq2"
  )

  missing <- da_comparison_fixture(methods = "ancombc2")
  missing$raw_outputs <- list()
  expect_error(as_da_raw_output(missing, contrast = "A_vs_B"), "is missing")
})

test_that("as_da_comparison aligns method-specific rows without consensus", {
  da <- da_comparison_fixture()
  comparison <- as_da_comparison(da)

  expect_s3_class(comparison, "data.frame")
  expect_equal(nrow(comparison), 6L)
  expect_equal(
    comparison$contrast,
    rep(c("A_vs_B", "A_vs_C"), each = 3)
  )
  expect_equal(comparison$feature_id, rep(c("F1", "F2", "F3"), 2))
  expect_equal(
    names(comparison)[7:14],
    paste0(
      "deseq2_",
      c(
        "effect",
        "effect_type",
        "p_value",
        "p_adjusted",
        "significance",
        "effect_sign",
        "tested",
        "p_adjusted_le_alpha"
      )
    )
  )
  expect_equal(comparison$taxon_label, paste0(
    "Genus_",
    rep(seq_len(3), 2)
  ))
  expect_true(all(comparison$rank == "Genus"))
  expect_equal(attr(comparison, "alpha"), 0.05)
  expect_equal(attr(comparison, "methods"), da$methods)
})

test_that("as_da_comparison distinguishes missing and untested rows", {
  da <- da_comparison_fixture()
  comparison <- as_da_comparison(da, contrast = "A_vs_B")
  f2 <- comparison[comparison$feature_id == "F2", , drop = FALSE]

  expect_identical(f2$aldex2_tested, FALSE)
  expect_true(is.na(f2$ancombc2_tested))
  expect_true(is.na(f2$ancombc2_effect))
  expect_true(is.na(f2$ancombc2_p_value))
  expect_identical(f2$deseq2_tested, TRUE)
  expect_true(is.na(f2$deseq2_p_adjusted))
  expect_equal(f2$deseq2_effect_sign, "zero")

  f3 <- comparison[comparison$feature_id == "F3", , drop = FALSE]
  expect_equal(f3$ancombc2_effect, -1.1)
  expect_equal(f3$ancombc2_p_value, 0.02)
  expect_equal(f3$ancombc2_p_adjusted, 0.04)
  expect_equal(f3$ancombc2_effect_sign, "negative")
  expect_identical(f3$ancombc2_p_adjusted_le_alpha, TRUE)

  lower_alpha <- as_da_comparison(
    da,
    contrast = "A_vs_B",
    alpha = 0.01
  )
  expect_identical(
    lower_alpha$ancombc2_p_adjusted_le_alpha[
      lower_alpha$feature_id == "F3"
    ],
    FALSE
  )
})

test_that("as_da_comparison preserves object order under filtering", {
  da <- da_comparison_fixture()
  comparison <- as_da_comparison(
    da,
    contrast = c("A_vs_C", "A_vs_B"),
    methods = c("ancombc2", "deseq2"),
    features = c("F3", "F1")
  )

  expect_equal(unique(comparison$contrast), c("A_vs_B", "A_vs_C"))
  expect_equal(attr(comparison, "methods"), c("deseq2", "ancombc2"))
  expect_equal(unique(comparison$feature_id), c("F1", "F3"))
  expect_false(any(grepl("^aldex2_", names(comparison))))
  expect_true(all(c("deseq2_effect", "ancombc2_effect") %in%
                    names(comparison)))

  expect_error(as_da_comparison(da, methods = "unknown"), "Available methods")
  expect_error(as_da_comparison(da, contrast = "unknown"), "Available contrasts")
  expect_error(as_da_comparison(da, features = "unknown"), "Unknown `features`")
  expect_error(
    as_da_comparison(da, methods = c("aldex2", "aldex2")),
    "duplicates"
  )
  expect_error(
    as_da_comparison(da, features = c("F1", "F1")),
    "duplicates"
  )
  expect_error(as_da_comparison(da, alpha = 0), "strictly between")
})

test_that("comparison schema has stable types and no consensus fields", {
  comparison <- as_da_comparison(da_comparison_fixture())
  for (method in attr(comparison, "methods")) {
    expect_type(comparison[[paste0(method, "_effect")]], "double")
    expect_type(comparison[[paste0(method, "_effect_type")]], "character")
    expect_type(comparison[[paste0(method, "_p_value")]], "double")
    expect_type(comparison[[paste0(method, "_p_adjusted")]], "double")
    expect_type(comparison[[paste0(method, "_significance")]], "character")
    expect_type(comparison[[paste0(method, "_effect_sign")]], "character")
    expect_type(comparison[[paste0(method, "_tested")]], "logical")
    expect_type(
      comparison[[paste0(method, "_p_adjusted_le_alpha")]],
      "logical"
    )
  }
  expect_false(any(grepl(
    "consensus|vote|winner|best|agreement|n_methods",
    names(comparison),
    ignore.case = TRUE
  )))
})

test_that("comparison report is descriptive and contrast-specific", {
  da <- da_comparison_fixture()
  report <- microeda_da_comparison_report(
    da,
    alpha = 0.05,
    max_features = 1
  )

  expect_type(report, "character")
  expect_length(report, 1L)
  expect_match(report, "comparison report", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_B", fixed = TRUE)
  expect_match(report, "Contrast: A_vs_C", fixed = TRUE)
  expect_match(report, "Method: aldex2", fixed = TRUE)
  expect_match(report, "Method: ancombc2", fixed = TRUE)
  expect_match(report, "Method: deseq2", fixed = TRUE)
  expect_match(report, "Display criterion", fixed = TRUE)
  expect_match(report, "not a consensus call", fixed = TRUE)
  expect_match(report, "Effect scales differ", fixed = TRUE)
  expect_match(report, "not directly comparable", fixed = TRUE)
  expect_match(report, "Missing values are not evidence", fixed = TRUE)
  expect_match(report, "not proof", fixed = TRUE)
  expect_match(report, "not a compositional correction", fixed = TRUE)
  expect_match(report, "DESeq2 diagnostics", fixed = TRUE)
  expect_match(report, "ANCOM-BC2 backend-excluded features", fixed = TRUE)
  expect_match(report, "hidden by `max_features`", fixed = TRUE)

  ab_report <- microeda_da_comparison_report(
    da,
    contrast = "A_vs_B",
    alpha = 0.05,
    max_features = Inf
  )
  expect_match(ab_report, "F1", fixed = TRUE)
  expect_match(ab_report, "F3", fixed = TRUE)
  expect_false(grepl("F2", ab_report, fixed = TRUE))

  no_features <- microeda_da_comparison_report(
    da,
    alpha = 0.0001,
    max_features = Inf
  )
  expect_match(no_features, "No features meet the display criterion", fixed = TRUE)

  requested <- microeda_da_comparison_report(
    da,
    contrast = "A_vs_B",
    features = "F2"
  )
  expect_match(requested, "explicitly requested feature IDs", fixed = TRUE)
  expect_match(requested, "F2", fixed = TRUE)
  expect_error(
    microeda_da_comparison_report(da, max_features = 0),
    "positive"
  )
})

test_that("write_da_comparison exports only the wide table", {
  da <- da_comparison_fixture()
  file <- tempfile(fileext = ".csv")
  expected <- as_da_comparison(da, contrast = "A_vs_B")
  visible <- withVisible(write_da_comparison(
    da,
    file,
    contrast = "A_vs_B",
    na = "NA"
  ))

  expect_false(visible$visible)
  expect_identical(visible$value, file)
  expect_true(file.exists(file))
  csv <- utils::read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
  expect_named(csv, names(expected))
  expect_equal(nrow(csv), nrow(expected))
  expect_equal(csv$feature_id, expected$feature_id)
  expect_equal(csv$deseq2_effect, expected$deseq2_effect)
  expect_false(any(grepl(
    "raw_outputs|method_results|native_dds|native_results",
    readLines(file, warn = FALSE)
  )))
  expect_false(names(csv)[[1]] %in% c("X", "row.names"))

  expect_error(write_da_comparison(da, ""), "non-empty")
  expect_error(
    write_da_comparison(da, tempfile(), row.names = TRUE),
    "Unsupported CSV argument"
  )
  expect_error(
    write_da_comparison(
      x = da,
      file = tempfile(),
      contrast = NULL,
      methods = NULL,
      alpha = 0.05,
      features = NULL,
      1
    ),
    "uniquely named"
  )
})

test_that("comparison helpers are exported and do not adjust p-values", {
  exports <- getNamespaceExports("microeda")
  expected <- c(
    "as_da_raw_output",
    "as_da_comparison",
    "microeda_da_comparison_report",
    "write_da_comparison"
  )
  expect_true(all(expected %in% exports))

  functions <- c(
    expected,
    "da_comparison_row",
    "da_comparison_method_block",
    "da_comparison_display_rows"
  )
  code <- unlist(lapply(functions, function(name) {
    deparse(get(name, envir = asNamespace("microeda")))
  }))
  expect_false(any(grepl("p\\.adjust\\s*\\(", code)))
})
