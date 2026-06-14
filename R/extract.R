microeda_extract <- function(x, metadata = NULL, taxonomy = NULL,
                             taxa_are_rows = TRUE) {
  if (inherits(x, "phyloseq")) {
    return(extract_phyloseq(x))
  }

  counts <- orient_count_matrix(x, taxa_are_rows = taxa_are_rows)
  metadata <- normalize_metadata(metadata, sample_names = rownames(counts))
  taxonomy <- normalize_taxonomy(taxonomy, feature_names = colnames(counts))

  list(
    counts = counts,
    metadata = metadata,
    taxonomy = taxonomy
  )
}

extract_phyloseq <- function(x) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop(
      "`phyloseq` input requires the phyloseq package to be installed.",
      call. = FALSE
    )
  }

  otu <- phyloseq::otu_table(x)
  counts <- as(otu, "matrix")

  if (phyloseq::taxa_are_rows(otu)) {
    counts <- t(counts)
  }

  metadata <- NULL
  if (!is.null(phyloseq::sample_data(x, errorIfNULL = FALSE))) {
    metadata <- as.data.frame(phyloseq::sample_data(x))
  }
  metadata <- normalize_metadata(metadata, sample_names = rownames(counts))

  taxonomy <- NULL
  if (!is.null(phyloseq::tax_table(x, errorIfNULL = FALSE))) {
    taxonomy <- as.data.frame(as(phyloseq::tax_table(x), "matrix"))
  }
  taxonomy <- normalize_taxonomy(taxonomy, feature_names = colnames(counts))

  list(
    counts = counts,
    metadata = metadata,
    taxonomy = taxonomy
  )
}

orient_count_matrix <- function(x, taxa_are_rows = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop("`x` must be a phyloseq object, matrix, or data frame.", call. = FALSE)
  }

  suppressWarnings(storage.mode(x) <- "double")

  if (any(is.na(x))) {
    stop("`x` must be a numeric count matrix without missing values.", call. = FALSE)
  }

  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("`x` must contain at least one sample and one feature.", call. = FALSE)
  }

  if (taxa_are_rows) {
    x <- t(x)
  }

  if (is.null(rownames(x))) {
    rownames(x) <- paste0("sample_", seq_len(nrow(x)))
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("feature_", seq_len(ncol(x)))
  }

  x
}

normalize_metadata <- function(metadata, sample_names) {
  if (is.null(metadata)) {
    return(NULL)
  }

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)

  if (has_automatic_rownames(metadata)) {
    if (nrow(metadata) != length(sample_names)) {
      stop(
        "`metadata` needs row names matching samples, or the same number ",
        "of rows as samples.",
        call. = FALSE
      )
    }
    rownames(metadata) <- sample_names
  }

  missing_samples <- setdiff(sample_names, rownames(metadata))
  if (length(missing_samples) > 0) {
    stop(
      "`metadata` is missing rows for samples: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  metadata[sample_names, , drop = FALSE]
}

normalize_taxonomy <- function(taxonomy, feature_names) {
  if (is.null(taxonomy)) {
    return(NULL)
  }

  taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)

  if (has_automatic_rownames(taxonomy)) {
    if (nrow(taxonomy) != length(feature_names)) {
      stop(
        "`taxonomy` needs row names matching features, or the same number ",
        "of rows as features.",
        call. = FALSE
      )
    }
    rownames(taxonomy) <- feature_names
  }

  missing_features <- setdiff(feature_names, rownames(taxonomy))
  if (length(missing_features) > 0) {
    stop(
      "`taxonomy` is missing rows for features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }

  taxonomy[feature_names, , drop = FALSE]
}

has_automatic_rownames <- function(x) {
  row_info <- .row_names_info(x, type = 0L)
  length(row_info) == 2 && is.na(row_info[1]) && row_info[2] < 0
}
