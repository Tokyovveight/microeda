target_qc_fixture <- function() {
  counts <- matrix(
    c(
      10, 0, 5, 0, 5, 0,
      0, 20, 5, 5, 0, 0,
      5, 0, 0, 5, 10, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("s", 1:3), paste0("f", 1:6))
  )
  taxonomy <- data.frame(
    Kingdom = c(
      "Fungi",
      "k__fungi",
      "Metazoa",
      NA,
      "Plantae",
      "Fungi_unclassified"
    ),
    row.names = paste0("f", 1:6),
    stringsAsFactors = FALSE
  )
  list(counts = counts, taxonomy = taxonomy)
}

target_qc <- function(target_match = "normalized") {
  fixture <- target_qc_fixture()
  microeda_qc(
    fixture$counts,
    taxonomy = fixture$taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi",
    target_match = target_match
  )
}

test_that("target_kingdom NULL preserves the existing QC path", {
  fixture <- target_qc_fixture()
  qc <- microeda_qc(
    fixture$counts,
    taxonomy = fixture$taxonomy,
    taxa_are_rows = FALSE
  )

  expect_null(qc$target_composition)
  expect_null(qc$params$target_kingdom)
  expect_null(qc$params$kingdom_rank)
  expect_equal(qc$params$target_match, "normalized")
  expect_named(
    qc$per_sample,
    c(
      "sample_id", "library_size", "zero_fraction",
      "n_features_detected", "n_features_above_prevalence"
    )
  )
  expect_named(
    qc$per_feature,
    c(
      "feature_id", "total_reads", "prevalence",
      "n_samples_detected", "above_prevalence_threshold"
    )
  )
})

test_that("microeda_qc validates target taxonomy requests", {
  fixture <- target_qc_fixture()

  expect_error(
    microeda_qc(
      fixture$counts,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "require a taxonomy"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = character()
    ),
    "target_kingdom"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = NA_character_
    ),
    "target_kingdom"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = " "
    ),
    "target_kingdom"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "unclassified"
    ),
    "classified kingdom"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi",
      target_match = "fuzzy"
    ),
    "arg"
  )
})

test_that("kingdom rank selection is case-insensitive and explicit", {
  fixture <- target_qc_fixture()
  explicit <- microeda_qc(
    fixture$counts,
    taxonomy = fixture$taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi",
    kingdom_rank = "kINGdom"
  )
  expect_equal(explicit$target_composition$kingdom_rank, "Kingdom")

  domain_taxonomy <- fixture$taxonomy
  names(domain_taxonomy) <- "Domain"
  automatic <- microeda_qc(
    fixture$counts,
    taxonomy = domain_taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  expect_equal(automatic$target_composition$kingdom_rank, "Domain")

  unknown_rank <- fixture$taxonomy
  names(unknown_rank) <- "Phylum"
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = unknown_rank,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "Available ranks: Phylum"
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi",
      kingdom_rank = "Phylum"
    ),
    "Available ranks: Kingdom"
  )

  ambiguous <- data.frame(
    Kingdom = fixture$taxonomy$Kingdom,
    kingdom = fixture$taxonomy$Kingdom,
    row.names = rownames(fixture$taxonomy),
    check.names = FALSE
  )
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = ambiguous,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "ambiguous"
  )
})

test_that("normalized matching handles only documented prefixes", {
  labels <- c(
    "Fungi", "fungi", "k__Fungi", "Kingdom__Fungi", "D_0__Fungi",
    "Fungi_unclassified", "Fungus", NA, "", "Unclassified"
  )
  counts <- matrix(
    1,
    nrow = 2,
    ncol = length(labels),
    dimnames = list(c("s1", "s2"), paste0("f", seq_along(labels)))
  )
  taxonomy <- data.frame(
    Kingdom = labels,
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )

  qc <- microeda_qc(
    counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  feature <- as_qc_target_composition(qc, "feature")

  expect_equal(
    feature$target_status,
    c(
      rep("target", 5),
      "non_target", "non_target",
      "unclassified", "unclassified", "unclassified"
    )
  )
  expect_equal(feature$original_kingdom, labels)
  expect_equal(feature$normalized_kingdom[1:5], rep("fungi", 5))
  expect_equal(feature$normalized_kingdom[[6]], "fungi_unclassified")
})

test_that("exact matching trims but preserves case and prefixes", {
  qc <- target_qc(target_match = "exact")
  feature <- as_qc_target_composition(qc, "feature")

  expect_equal(feature$target_status[1:2], c("target", "non_target"))
  expect_equal(qc$target_composition$target_match, "exact")
  expect_equal(qc$target_composition$overall$target_features, 1L)
})

test_that("multiple target kingdoms are supported", {
  fixture <- target_qc_fixture()
  qc <- microeda_qc(
    fixture$counts,
    taxonomy = fixture$taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = c("Fungi", "Metazoa")
  )

  expect_equal(
    qc$target_composition$target_kingdom,
    c("Fungi", "Metazoa")
  )
  expect_equal(qc$target_composition$overall$target_features, 3L)
  expect_equal(qc$target_composition$overall$target_reads, 45)
})

test_that("target composition calculates exact overall and kingdom totals", {
  qc <- target_qc()
  overall <- as_qc_target_composition(qc, "overall")
  kingdom <- as_qc_target_composition(qc, "kingdom")

  expect_named(
    overall,
    c(
      "total_features", "total_reads", "target_features",
      "non_target_features", "unclassified_features",
      "target_feature_proportion", "non_target_feature_proportion",
      "unclassified_feature_proportion", "target_reads",
      "non_target_reads", "unclassified_reads", "target_read_proportion",
      "non_target_read_proportion", "unclassified_read_proportion"
    )
  )
  expect_equal(
    unname(unlist(overall[
      c("target_features", "non_target_features", "unclassified_features")
    ])),
    c(2, 3, 1)
  )
  expect_equal(
    unname(unlist(overall[
      c("target_reads", "non_target_reads", "unclassified_reads")
    ])),
    c(35, 30, 10)
  )
  expect_equal(
    sum(unlist(overall[
      c(
        "target_feature_proportion",
        "non_target_feature_proportion",
        "unclassified_feature_proportion"
      )
    ])),
    1
  )
  expect_equal(
    sum(unlist(overall[
      c(
        "target_read_proportion",
        "non_target_read_proportion",
        "unclassified_read_proportion"
      )
    ])),
    1
  )
  expect_named(
    kingdom,
    c(
      "kingdom_label", "normalized_kingdom", "target_status",
      "n_features", "feature_proportion", "reads", "read_proportion",
      "proportion_within_non_target_reads", "rank"
    )
  )
  expect_equal(sum(kingdom$n_features), 6)
  expect_equal(sum(kingdom$reads), 75)
  expect_equal(
    sum(
      kingdom$proportion_within_non_target_reads[
        kingdom$target_status == "non_target"
      ]
    ),
    1
  )
  expect_true(all(is.na(
    kingdom$proportion_within_non_target_reads[
      kingdom$target_status != "non_target"
    ]
  )))
})

test_that("target composition calculates per-sample and per-feature values", {
  qc <- target_qc()
  sample <- as_qc_target_composition(qc, "sample")
  feature <- as_qc_target_composition(qc, "feature")

  expect_named(
    sample,
    c(
      "sample_id", "total_reads", "target_reads", "non_target_reads",
      "unclassified_reads", "target_read_proportion",
      "non_target_read_proportion", "unclassified_read_proportion",
      "dominant_status", "dominant_non_target_kingdom",
      "dominant_non_target_reads", "dominant_non_target_read_proportion"
    )
  )
  expect_equal(sample$total_reads, c(20, 30, 25))
  expect_equal(sample$target_reads, c(10, 20, 5))
  expect_equal(sample$non_target_reads, c(10, 5, 15))
  expect_equal(sample$unclassified_reads, c(0, 5, 5))
  expect_equal(sample$dominant_status, c("target", "target", "non_target"))
  expect_equal(
    sample$dominant_non_target_kingdom,
    c("Metazoa", "Metazoa", "Plantae")
  )
  expect_equal(sample$dominant_non_target_reads, c(5, 5, 10))
  expect_equal(
    sample$dominant_non_target_read_proportion,
    c(0.25, 1 / 6, 0.4)
  )

  expect_named(
    feature,
    c(
      "feature_id", "original_kingdom", "normalized_kingdom",
      "target_status", "total_reads", "prevalence",
      "overall_read_proportion"
    )
  )
  expect_equal(feature$total_reads, c(15, 20, 10, 10, 15, 5))
  expect_equal(feature$prevalence, c(2 / 3, 1 / 3, 2 / 3, 2 / 3, 2 / 3, 1 / 3))
  expect_equal(sum(feature$overall_read_proportion), 1)
  expect_false(any(names(feature) %in% c("Phylum", "Genus")))
})

test_that("zero read denominators produce typed NA values", {
  counts <- matrix(
    0,
    nrow = 2,
    ncol = 3,
    dimnames = list(c("s1", "s2"), c("f1", "f2", "f3"))
  )
  taxonomy <- data.frame(
    Kingdom = c("Fungi", "Metazoa", NA),
    row.names = colnames(counts)
  )
  qc <- microeda_qc(
    counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  overall <- qc$target_composition$overall
  sample <- qc$target_composition$per_sample
  kingdom <- qc$target_composition$by_kingdom

  expect_true(all(is.na(unlist(overall[
    c(
      "target_read_proportion",
      "non_target_read_proportion",
      "unclassified_read_proportion"
    )
  ]))))
  expect_true(all(is.na(sample$target_read_proportion)))
  expect_true(all(is.na(sample$non_target_read_proportion)))
  expect_true(all(is.na(sample$unclassified_read_proportion)))
  expect_true(all(is.na(sample$dominant_status)))
  expect_true(all(is.na(sample$dominant_non_target_read_proportion)))
  expect_true(all(is.na(kingdom$read_proportion)))
  expect_true(all(is.na(kingdom$proportion_within_non_target_reads)))
  expect_false(any(is.nan(unlist(overall))))
  expect_false(any(is.infinite(unlist(overall))))
})

test_that("target composition handles absent category edge cases", {
  counts <- matrix(
    c(5, 0, 0, 5),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("f1", "f2"))
  )
  all_target_tax <- data.frame(
    Kingdom = c("Fungi", "Fungi"),
    row.names = colnames(counts)
  )
  all_target <- microeda_qc(
    counts,
    taxonomy = all_target_tax,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  expect_equal(all_target$target_composition$overall$target_read_proportion, 1)
  expect_equal(all_target$target_composition$overall$non_target_reads, 0)
  expect_true("non_target_kingdom_absent" %in%
    all_target$target_composition$observations$observation_id)

  no_target_tax <- data.frame(
    Kingdom = c("Metazoa", "Plantae"),
    row.names = colnames(counts)
  )
  no_target <- microeda_qc(
    counts,
    taxonomy = no_target_tax,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  expect_equal(no_target$target_composition$overall$target_features, 0)
  expect_true("target_kingdom_absent" %in%
    no_target$target_composition$observations$observation_id)

  unknown_tax <- data.frame(
    Kingdom = c("unknown", NA),
    row.names = colnames(counts)
  )
  unknown <- microeda_qc(
    counts,
    taxonomy = unknown_tax,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  expect_equal(unknown$target_composition$overall$unclassified_features, 2L)
  expect_equal(unknown$target_composition$overall$non_target_features, 0L)
})

test_that("samples containing one status retain explicit missing semantics", {
  counts <- rbind(
    target_only = c(10, 0, 0),
    non_target_only = c(0, 10, 0),
    unclassified_only = c(0, 0, 10),
    zero_library = c(0, 0, 0)
  )
  colnames(counts) <- c("target", "non_target", "unknown")
  taxonomy <- data.frame(
    Kingdom = c("Fungi", "Metazoa", "unknown"),
    row.names = colnames(counts)
  )
  qc <- microeda_qc(
    counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  sample <- qc$target_composition$per_sample

  expect_equal(
    sample$dominant_status,
    c("target", "non_target", "unclassified", NA)
  )
  expect_equal(sample$target_read_proportion[1:3], c(1, 0, 0))
  expect_equal(sample$non_target_read_proportion[1:3], c(0, 1, 0))
  expect_equal(sample$unclassified_read_proportion[1:3], c(0, 0, 1))
  expect_true(all(is.na(sample[4, c(
    "target_read_proportion",
    "non_target_read_proportion",
    "unclassified_read_proportion"
  )])))
})

test_that("taxonomy rows are aligned and extra rows do not affect results", {
  fixture <- target_qc_fixture()
  shuffled <- fixture$taxonomy[
    c("f6", "f4", "f2", "f1", "f5", "f3"),
    ,
    drop = FALSE
  ]
  extra <- rbind(
    shuffled,
    data.frame(Kingdom = "Metazoa", row.names = "extra")
  )
  qc <- microeda_qc(
    fixture$counts,
    taxonomy = extra,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )

  expect_equal(qc$target_composition$per_feature$feature_id, paste0("f", 1:6))
  expect_equal(
    qc$target_composition$per_feature$original_kingdom,
    fixture$taxonomy$Kingdom
  )
})

test_that("target composition rejects missing or duplicated identifiers", {
  fixture <- target_qc_fixture()
  missing <- fixture$taxonomy[-1, , drop = FALSE]
  expect_error(
    microeda_qc(
      fixture$counts,
      taxonomy = missing,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "missing rows for features"
  )

  duplicate_counts <- fixture$counts
  colnames(duplicate_counts)[[2]] <- colnames(duplicate_counts)[[1]]
  expect_error(
    microeda_qc(
      duplicate_counts,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "Feature IDs must be unique"
  )

  duplicate_samples <- fixture$counts
  rownames(duplicate_samples)[[2]] <- rownames(duplicate_samples)[[1]]
  expect_error(
    microeda_qc(
      duplicate_samples,
      taxonomy = fixture$taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "Sample IDs must be unique"
  )

  duplicate_taxonomy <- matrix(
    c("Fungi", "Metazoa"),
    ncol = 1,
    dimnames = list(c("f1", "f1"), "Kingdom")
  )
  small_counts <- fixture$counts[, 1:2, drop = FALSE]
  expect_error(
    microeda_qc(
      small_counts,
      taxonomy = duplicate_taxonomy,
      taxa_are_rows = FALSE,
      target_kingdom = "Fungi"
    ),
    "row names must be unique"
  )
})

test_that("matrix orientations and data-frame counts produce the same diagnostic", {
  fixture <- target_qc_fixture()
  sample_rows <- microeda_qc(
    as.data.frame(fixture$counts),
    taxonomy = fixture$taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  feature_rows <- microeda_qc(
    t(fixture$counts),
    taxonomy = fixture$taxonomy,
    taxa_are_rows = TRUE,
    target_kingdom = "Fungi"
  )

  expect_equal(
    sample_rows$target_composition,
    feature_rows$target_composition
  )
})

test_that("phyloseq target composition respects object orientation", {
  skip_if_not_installed("phyloseq")
  fixture <- target_qc_fixture()
  otu <- phyloseq::otu_table(t(fixture$counts), taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(as.matrix(fixture$taxonomy))
  ps <- phyloseq::phyloseq(otu, tax)

  qc <- microeda_qc(ps, target_kingdom = "Fungi")

  expect_equal(qc$target_composition$overall$total_features, 6L)
  expect_equal(qc$target_composition$overall$total_reads, 75)
  expect_equal(qc$target_composition$per_sample$sample_id, paste0("s", 1:3))
})

test_that("as_qc_target_composition returns stable data frames", {
  qc <- target_qc()
  levels <- c("overall", "kingdom", "sample", "feature")

  for (level in levels) {
    expect_s3_class(as_qc_target_composition(qc, level), "data.frame")
  }
  expect_equal(nrow(as_qc_target_composition(qc, "overall")), 1)
  expect_type(as_qc_target_composition(qc, "overall")$total_features, "integer")
  expect_type(as_qc_target_composition(qc, "overall")$total_reads, "double")
  expect_type(
    as_qc_target_composition(qc, "sample")$dominant_status,
    "character"
  )
  expect_type(
    as_qc_target_composition(qc, "feature")$prevalence,
    "double"
  )

  expect_error(as_qc_target_composition(list()), "microeda_qc")
  fixture <- target_qc_fixture()
  no_target <- microeda_qc(fixture$counts, taxa_are_rows = FALSE)
  expect_error(
    as_qc_target_composition(no_target),
    "target_kingdom = \"Fungi\""
  )
  expect_error(as_qc_target_composition(qc, "unknown"), "arg")
})

test_that("target observations add only meaningful issue rows", {
  qc <- target_qc()
  issues <- as_qc_issues(qc)

  expect_true("unclassified_kingdom_signal" %in% issues$issue_id)
  expect_false("dominant_non_target_kingdom" %in% issues$issue_id)
  expect_false("max_sample_non_target_signal" %in% issues$issue_id)
  expect_false(any(grepl("contamination confirmed", issues$message, ignore.case = TRUE)))

  fixture <- target_qc_fixture()
  taxonomy <- fixture$taxonomy
  taxonomy$Kingdom <- rep("Metazoa", nrow(taxonomy))
  no_target <- microeda_qc(
    fixture$counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  expect_true("target_kingdom_absent" %in% as_qc_issues(no_target)$issue_id)
  expect_false(any(grepl("remove", as_qc_issues(no_target)$message, ignore.case = TRUE)))
})

test_that("microeda_qc_report includes compact target composition", {
  qc <- target_qc()
  report <- microeda_qc_report(qc, target_top_n = 1)
  lines <- strsplit(report, "\n", fixed = TRUE)[[1]]

  expect_type(report, "character")
  expect_length(report, 1)
  expect_match(report, "Target composition", fixed = TRUE)
  expect_match(report, "Target kingdom(s): Fungi", fixed = TRUE)
  expect_match(report, "Taxonomy rank: Kingdom", fixed = TRUE)
  expect_match(report, "Matching mode: normalized", fixed = TRUE)
  expect_match(report, "Target: 2 features", fixed = TRUE)
  expect_match(report, "Non-target: 30 reads", fixed = TRUE)
  expect_match(report, "Top known non-target kingdoms", fixed = TRUE)
  expect_match(report, "Samples with highest non-target", fixed = TRUE)
  expect_match(report, "automatically evidence of contamination", fixed = TRUE)
  expect_lte(max(nchar(lines, type = "width")), 110)
  kingdom_heading <- match("Top known non-target kingdoms by reads:", lines)
  kingdom_end <- which(
    seq_along(lines) > kingdom_heading & !nzchar(lines)
  )[[1]]
  kingdom_lines <- lines[
    seq.int(kingdom_heading + 1L, kingdom_end - 1L)
  ]
  expect_length(kingdom_lines, 1)

  no_section <- microeda_qc_report(
    qc,
    include_target_composition = FALSE
  )
  expect_false(grepl("Target composition", no_section, fixed = TRUE))

  fixture <- target_qc_fixture()
  ordinary <- microeda_qc(fixture$counts, taxa_are_rows = FALSE)
  expect_false(grepl(
    "Target composition",
    microeda_qc_report(ordinary),
    fixed = TRUE
  ))
})

test_that("target report validates controls and formats tiny proportions", {
  qc <- target_qc()
  expect_error(
    microeda_qc_report(qc, include_target_composition = NA),
    "include_target_composition"
  )
  expect_error(microeda_qc_report(qc, target_top_n = 0), "target_top_n")
  expect_error(microeda_qc_report(qc, target_top_n = 1.5), "target_top_n")

  counts <- matrix(
    c(1, 1e8),
    nrow = 1,
    dimnames = list("s1", c("target", "other"))
  )
  taxonomy <- data.frame(
    Kingdom = c("Fungi", "Metazoa"),
    row.names = colnames(counts)
  )
  tiny <- microeda_qc(
    counts,
    taxonomy = taxonomy,
    taxa_are_rows = FALSE,
    target_kingdom = "Fungi"
  )
  report <- microeda_qc_report(tiny)

  expect_match(report, "e-06%", fixed = TRUE)
  expect_false(grepl("Target: 1 read (0%).", report, fixed = TRUE))
})

test_that("microeda_qc_write_report includes the same target section", {
  qc <- target_qc()
  path <- tempfile(fileext = ".txt")
  expected <- microeda_qc_report(qc, target_top_n = 2)

  returned <- withVisible(microeda_qc_write_report(
    qc,
    path,
    target_top_n = 2
  ))
  written <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_false(returned$visible)
  expect_equal(returned$value, path)
  expect_identical(written, expected)
})

test_that("microeda_qc_plot draws target sample composition and restores par", {
  qc <- target_qc()
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)

  expect_silent(
    result <- withVisible(microeda_qc_plot(
      qc,
      type = "target_composition",
      target_view = "samples",
      top_n = 2
    ))
  )

  expect_false(result$visible)
  expect_s3_class(result$value, "data.frame")
  expect_equal(nrow(result$value), 2)
  expect_equal(result$value$sample_id, c("s3", "s1"))
  expect_equal(graphics::par(no.readonly = TRUE), old_par)
})

test_that("microeda_qc_plot draws kingdom composition without changing tables", {
  qc <- target_qc()
  before <- qc$target_composition$by_kingdom
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(
    result <- withVisible(microeda_qc_plot(
      qc,
      type = "target_composition",
      target_view = "kingdoms",
      top_n = 1
    ))
  )

  expect_false(result$visible)
  expect_s3_class(result$value, "data.frame")
  expect_true("Unclassified" %in% result$value$kingdom_label)
  expect_true(any(grepl("^Other", result$value$kingdom_label)))
  expect_identical(qc$target_composition$by_kingdom, before)
})

test_that("target composition plot validates diagnostic and controls", {
  fixture <- target_qc_fixture()
  ordinary <- microeda_qc(fixture$counts, taxa_are_rows = FALSE)
  qc <- target_qc()

  expect_error(
    microeda_qc_plot(ordinary, type = "target_composition"),
    "target_kingdom = \"Fungi\""
  )
  expect_error(
    microeda_qc_plot(
      qc,
      type = "target_composition",
      target_view = "unknown"
    ),
    "arg"
  )
  expect_error(
    microeda_qc_plot(qc, type = "target_composition", top_n = 0),
    "top_n"
  )
})

test_that("print.microeda_qc adds compact target hints", {
  qc <- target_qc()
  before <- qc$target_composition
  output <- capture.output(result <- withVisible(print(qc)))

  expect_false(result$visible)
  expect_identical(result$value, qc)
  expect_true(any(grepl("Target kingdom:", output, fixed = TRUE)))
  expect_true(any(grepl("Target reads:", output, fixed = TRUE)))
  expect_true(any(grepl("Non-target reads:", output, fixed = TRUE)))
  expect_true(any(grepl("as_qc_target_composition", output, fixed = TRUE)))
  expect_identical(qc$target_composition, before)
})
