alpha_report_example <- function() {
  counts <- matrix(
    c(
      10, 0, 0, 5,
      20, 0, 1, 0,
      0, 4, 0, 0,
      2, 3, 0, 1,
      10, 8, 7, 6,
      1, 0, 0, 0
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(counts) <- paste0("S", seq_len(6))
  colnames(counts) <- paste0("ASV", seq_len(4))
  metadata <- data.frame(
    group = c("A", "A", "B", "B", "C", "C"),
    row.names = rownames(counts)
  )

  alpha <- microeda_alpha(
    counts,
    metadata = metadata,
    group = "group",
    taxa_are_rows = FALSE
  )
  alpha_cmp <- microeda_alpha_compare(
    alpha,
    indices = c("observed", "shannon", "simpson")
  )

  list(alpha = alpha, alpha_cmp = alpha_cmp)
}

test_that("microeda_alpha_report returns compact alpha summaries", {
  example <- alpha_report_example()
  report <- microeda_alpha_report(example$alpha)

  expect_type(report, "character")
  expect_length(report, 1)
  expect_match(report, "Alpha diversity summary", fixed = TRUE)
  expect_match(report, "Index: observed", fixed = TRUE)
  expect_match(report, "Index: shannon", fixed = TRUE)
  expect_match(report, "Index: simpson", fixed = TRUE)
  expect_match(report, "group", fixed = TRUE)
  expect_match(report, "n", fixed = TRUE)
  expect_match(report, "mean", fixed = TRUE)
  expect_match(report, "median", fixed = TRUE)
})

test_that("microeda_alpha_report includes compact alpha group tests", {
  example <- alpha_report_example()
  report <- microeda_alpha_report(
    example$alpha,
    alpha_compare = example$alpha_cmp
  )

  expect_match(report, "Alpha group tests", fixed = TRUE)
  expect_match(report, "p", fixed = TRUE)
  expect_match(report, "p.adj", fixed = TRUE)
  expect_match(report, "p.adj.signif", fixed = TRUE)
  expect_match(report, "ns", fixed = TRUE)

  lines <- strsplit(report, "\n", fixed = TRUE)[[1]]
  expect_true(any(grepl(
    "^index\\s+method\\s+n\\s+n_groups\\s+statistic\\s+p\\s+p\\.adj\\s+p\\.adj\\.signif$",
    lines
  )))
})

test_that("alpha report leaves extractor tables unchanged", {
  example <- alpha_report_example()
  summary <- as_alpha_summary(example$alpha)
  tests <- as_alpha_tests(example$alpha_cmp)

  expect_s3_class(summary, "data.frame")
  expect_named(
    summary,
    c("group", "index", "n", "mean", "sd", "median", "q1", "q3", "min", "max")
  )
  expect_s3_class(tests, "data.frame")
  expect_named(
    tests,
    c(
      "index",
      "method",
      "n",
      "n_groups",
      "min_group_n",
      "statistic",
      "parameter",
      "p_value",
      "max_median_group",
      "min_median_group",
      "median_difference",
      "p_value_adjusted",
      "p_adjust_method"
    )
  )
})

test_that("microeda_alpha_report validates inputs", {
  example <- alpha_report_example()

  expect_error(microeda_alpha_report(data.frame()), "microeda_alpha")
  expect_error(
    microeda_alpha_report(example$alpha, alpha_compare = data.frame()),
    "microeda_alpha_compare"
  )
  expect_error(microeda_alpha_report(example$alpha, digits = NA_real_), "digits")
  expect_error(microeda_alpha_report(example$alpha, digits = 1.5), "digits")
})
