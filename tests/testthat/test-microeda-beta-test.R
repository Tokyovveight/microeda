beta_test_example <- function() {
  counts <- matrix(
    c(
      10, 0, 0, 1,
      8, 2, 0, 1,
      6, 1, 1, 0,
      0, 9, 1, 0,
      1, 7, 2, 0,
      0, 6, 3, 1,
      1, 0, 8, 3,
      0, 2, 6, 2,
      2, 1, 5, 4
    ),
    nrow = 9,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(9))
  colnames(counts) <- paste0("ASV", seq_len(4))
  metadata <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    row.names = rownames(counts)
  )

  microeda_beta(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    method = "bray"
  )
}

beta_compare_test_example <- function(methods = c("bray", "jaccard")) {
  counts <- matrix(
    c(
      10, 0, 0, 1,
      8, 2, 0, 1,
      6, 1, 1, 0,
      0, 9, 1, 0,
      1, 7, 2, 0,
      0, 6, 3, 1,
      1, 0, 8, 3,
      0, 2, 6, 2,
      2, 1, 5, 4
    ),
    nrow = 9,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(9))
  colnames(counts) <- paste0("ASV", seq_len(4))
  metadata <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    row.names = rownames(counts)
  )

  microeda_beta_compare(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE,
    methods = methods
  )
}

test_that("microeda_beta_test validates input before running vegan", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", "A", "B"),
    row.names = rownames(counts)
  )
  beta <- microeda_beta(counts, taxa_are_rows = FALSE)
  grouped_beta <- microeda_beta(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  expect_error(microeda_beta_test(data.frame()), "microeda_beta")
  expect_error(microeda_beta_test(beta), "group metadata")
  expect_error(microeda_beta_test(grouped_beta, permutations = 0), "permutations")
  expect_error(microeda_beta_test(grouped_beta, seed = NA_real_), "seed")
})

test_that("microeda_beta_test requires complete group labels", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  metadata <- data.frame(
    group = c("A", NA_character_, "B"),
    row.names = rownames(counts)
  )
  beta <- microeda_beta(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )

  expect_error(microeda_beta_test(beta), "Group labels")
})

test_that("microeda_beta_test returns paired beta group diagnostics", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()

  beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)

  expect_s3_class(beta_test, "microeda_beta_test")
  expect_named(
    beta_test,
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova",
      "dispersion",
      "caveats",
      "params",
      "call"
    )
  )
  expect_equal(beta_test$method, "bray")
  expect_equal(beta_test$group, "group")
  expect_equal(beta_test$n_samples, 9L)
  expect_equal(beta_test$n_groups, 3L)
  expect_equal(beta_test$min_group_n, 3L)
  expect_s3_class(beta_test$permanova, "data.frame")
  expect_named(
    beta_test$permanova,
    c("term", "df", "sum_of_squares", "r2", "statistic", "p_value")
  )
  expect_type(beta_test$dispersion, "list")
  expect_named(beta_test$dispersion, c("test", "groups"))
  expect_s3_class(beta_test$dispersion$test, "data.frame")
  expect_s3_class(beta_test$dispersion$groups, "data.frame")
  expect_s3_class(beta_test$caveats, "data.frame")
  expect_equal(beta_test$params$permutations, 99L)
  expect_equal(beta_test$params$seed, 1L)
})

test_that("as_beta_test_summary returns stable compact columns", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()
  beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)

  summary <- as_beta_test_summary(beta_test)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova_r2",
      "permanova_f",
      "permanova_p",
      "dispersion_f",
      "dispersion_p",
      "permutations"
    )
  )
  expect_equal(nrow(summary), 1L)
  expect_equal(summary$method, "bray")
  expect_equal(summary$permutations, 99L)
  expect_error(as_beta_test_summary(data.frame()), "microeda_beta_test")
})

test_that("print.microeda_beta_test gives compact report and helper hints", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()
  beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)

  visible <- NULL
  output <- capture.output(visible <- withVisible(print(beta_test)))
  output <- paste(output, collapse = "\n")

  expect_false(visible$visible)
  expect_identical(visible$value, beta_test)
  expect_match(output, "<microeda_beta_test>", fixed = TRUE)
  expect_match(output, "Method:", fixed = TRUE)
  expect_match(output, "Group:", fixed = TRUE)
  expect_match(output, "Samples:", fixed = TRUE)
  expect_match(output, "Groups:", fixed = TRUE)
  expect_match(output, "Permutations:", fixed = TRUE)
  expect_match(output, "PERMANOVA:", fixed = TRUE)
  expect_match(output, "Dispersion:", fixed = TRUE)
  expect_match(output, "microeda_beta_test_report", fixed = TRUE)
  expect_match(output, "as_beta_test_summary", fixed = TRUE)
  expect_named(
    as_beta_test_summary(beta_test),
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova_r2",
      "permanova_f",
      "permanova_p",
      "dispersion_f",
      "dispersion_p",
      "permutations"
    )
  )
})

test_that("microeda_beta_test_report includes paired diagnostics and caveats", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()
  beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)

  report <- microeda_beta_test_report(beta_test)

  expect_type(report, "character")
  expect_length(report, 1)
  expect_match(report, "Beta group test", fixed = TRUE)
  expect_match(report, "PERMANOVA", fixed = TRUE)
  expect_match(report, "Dispersion diagnostics", fixed = TRUE)
  expect_match(report, "Mean distance to group centroid", fixed = TRUE)
  expect_match(report, "Caveats", fixed = TRUE)
  expect_match(report, "PERMANOVA can be confounded", fixed = TRUE)
  expect_match(report, "do not interpret it alone", fixed = TRUE)
  expect_error(microeda_beta_test_report(data.frame()), "microeda_beta_test")
})

test_that("microeda_beta_compare_test validates compare input", {
  counts <- matrix(
    c(
      1, 2, 0,
      2, 1, 0,
      0, 0, 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- c("S1", "S2", "S3")
  colnames(counts) <- paste0("ASV", seq_len(3))
  beta_cmp <- microeda_beta_compare(counts, taxa_are_rows = FALSE)

  expect_error(microeda_beta_compare_test(data.frame()), "microeda_beta_compare")
  expect_error(microeda_beta_compare_test(beta_cmp), "group metadata")
})

test_that("microeda_beta_compare_test returns one beta test per method", {
  skip_if_not_installed("vegan")
  beta_cmp <- beta_compare_test_example(methods = c("jaccard", "bray"))

  beta_cmp_test <- microeda_beta_compare_test(
    beta_cmp,
    permutations = 99,
    seed = 1
  )

  expect_s3_class(beta_cmp_test, "microeda_beta_compare_test")
  expect_named(
    beta_cmp_test,
    c(
      "results",
      "methods",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "caveats",
      "params",
      "call"
    )
  )
  expect_equal(beta_cmp_test$methods, c("jaccard", "bray"))
  expect_named(beta_cmp_test$results, c("jaccard", "bray"))
  expect_true(all(vapply(
    beta_cmp_test$results,
    inherits,
    logical(1),
    what = "microeda_beta_test"
  )))
  expect_equal(beta_cmp_test$group, "group")
  expect_equal(beta_cmp_test$n_samples, 9L)
  expect_equal(beta_cmp_test$n_groups, 3L)
  expect_equal(beta_cmp_test$min_group_n, 3L)
  expect_equal(beta_cmp_test$params$permutations, 99L)
  expect_equal(beta_cmp_test$params$seed, 1L)
  expect_named(
    beta_cmp_test$caveats,
    c("method", "caveat_id", "severity", "message")
  )
  expect_true(all(beta_cmp_test$caveats$method %in% c("jaccard", "bray")))
})

test_that("as_beta_compare_test_summary returns stable rows by method", {
  skip_if_not_installed("vegan")
  beta_cmp <- beta_compare_test_example()
  beta_cmp_test <- microeda_beta_compare_test(
    beta_cmp,
    permutations = 99,
    seed = 1
  )

  summary <- as_beta_compare_test_summary(beta_cmp_test)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova_r2",
      "permanova_f",
      "permanova_p",
      "dispersion_f",
      "dispersion_p",
      "permutations"
    )
  )
  expect_equal(summary$method, beta_cmp$methods)
  expect_equal(nrow(summary), length(beta_cmp$methods))
  expect_equal(summary$permutations, rep(99L, length(beta_cmp$methods)))
  expect_error(
    as_beta_compare_test_summary(data.frame()),
    "microeda_beta_compare_test"
  )
})

test_that("microeda_beta_compare_test_report formats method blocks and caveats", {
  skip_if_not_installed("vegan")
  beta_cmp <- beta_compare_test_example()
  beta_cmp_test <- microeda_beta_compare_test(
    beta_cmp,
    permutations = 99,
    seed = 1
  )

  report <- microeda_beta_compare_test_report(beta_cmp_test)

  expect_type(report, "character")
  expect_length(report, 1)
  expect_match(report, "Beta comparison group tests", fixed = TRUE)
  expect_match(report, "Method: bray", fixed = TRUE)
  expect_match(report, "Method: jaccard", fixed = TRUE)
  expect_match(report, "PERMANOVA", fixed = TRUE)
  expect_match(report, "Dispersion diagnostics", fixed = TRUE)
  expect_match(report, "Aggregated caveats", fixed = TRUE)
  expect_match(
    report,
    "These side-by-side diagnostics are not a formal method ranking.",
    fixed = TRUE
  )
  expect_match(report, "[bray warning]", fixed = TRUE)
  expect_match(report, "[jaccard warning]", fixed = TRUE)
  expect_false(grepl("best method", report, fixed = TRUE))
  expect_false(grepl("preferred method", report, fixed = TRUE))
  expect_false(grepl("most significant", report, fixed = TRUE))
  expect_false(grepl("winner", report, fixed = TRUE))
  expect_error(
    microeda_beta_compare_test_report(data.frame()),
    "microeda_beta_compare_test"
  )
})

test_that("print.microeda_beta_compare_test gives compact summary and helper hints", {
  skip_if_not_installed("vegan")
  beta_cmp <- beta_compare_test_example()
  beta_cmp_test <- microeda_beta_compare_test(
    beta_cmp,
    permutations = 99,
    seed = 1
  )

  visible <- NULL
  output <- capture.output(visible <- withVisible(print(beta_cmp_test)))
  output <- paste(output, collapse = "\n")

  expect_false(visible$visible)
  expect_identical(visible$value, beta_cmp_test)
  expect_match(output, "<microeda_beta_compare_test>", fixed = TRUE)
  expect_match(output, "Methods:", fixed = TRUE)
  expect_match(output, "Group:", fixed = TRUE)
  expect_match(output, "Samples:", fixed = TRUE)
  expect_match(output, "Groups:", fixed = TRUE)
  expect_match(output, "Permutations:", fixed = TRUE)
  expect_match(output, "bray", fixed = TRUE)
  expect_match(output, "jaccard", fixed = TRUE)
  expect_match(output, "PERMANOVA", fixed = TRUE)
  expect_match(output, "Dispersion", fixed = TRUE)
  expect_match(output, "microeda_beta_compare_test_report", fixed = TRUE)
  expect_match(output, "as_beta_compare_test_summary", fixed = TRUE)
  expect_named(
    as_beta_compare_test_summary(beta_cmp_test),
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova_r2",
      "permanova_f",
      "permanova_p",
      "dispersion_f",
      "dispersion_p",
      "permutations"
    )
  )
})

test_that("compare testing keeps single-method beta test behavior unchanged", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()
  beta_cmp <- beta_compare_test_example()
  beta_test <- microeda_beta_test(beta, permutations = 99, seed = 1)

  invisible(microeda_beta_compare_test(beta_cmp, permutations = 99, seed = 1))

  expect_s3_class(beta_test, "microeda_beta_test")
  expect_named(
    beta_test,
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova",
      "dispersion",
      "caveats",
      "params",
      "call"
    )
  )
  expect_named(
    as_beta_test_summary(beta_test),
    c(
      "method",
      "group",
      "n_samples",
      "n_groups",
      "min_group_n",
      "permanova_r2",
      "permanova_f",
      "permanova_p",
      "dispersion_f",
      "dispersion_p",
      "permutations"
    )
  )
})

test_that("beta testing does not change existing beta extractor columns", {
  skip_if_not_installed("vegan")
  beta <- beta_test_example()
  beta_cmp <- microeda_beta_compare(
    matrix(
      c(
        10, 0, 0, 1,
        8, 2, 0, 1,
        6, 1, 1, 0,
        0, 9, 1, 0,
        1, 7, 2, 0,
        0, 6, 3, 1,
        1, 0, 8, 3,
        0, 2, 6, 2,
        2, 1, 5, 4
      ),
      nrow = 9,
      byrow = TRUE,
      dimnames = list(paste0("S", seq_len(9)), paste0("ASV", seq_len(4)))
    ),
    metadata = data.frame(
      group = rep(c("A", "B", "C"), each = 3),
      row.names = paste0("S", seq_len(9))
    ),
    group = "group",
    taxa_are_rows = FALSE
  )

  invisible(microeda_beta_test(beta, permutations = 99, seed = 1))

  expect_named(as_beta_samples(beta), c("sample_id", "group"))
  expect_s3_class(as_beta_dist(beta), "dist")
  expect_true(is.matrix(as_beta_matrix(beta)))
  expect_named(
    as_beta_compare_summary(beta_cmp),
    c(
      "method",
      "n_samples",
      "n_pairs",
      "min_distance",
      "median_distance",
      "max_distance"
    )
  )
  expect_named(
    as_beta_compare_distances(beta_cmp),
    c(
      "method",
      "sample_1",
      "sample_2",
      "group_1",
      "group_2",
      "comparison",
      "distance"
    )
  )
  expect_named(
    as_beta_compare_group_summary(beta_cmp),
    c(
      "method",
      "comparison",
      "n_pairs",
      "min_distance",
      "median_distance",
      "max_distance"
    )
  )
})
