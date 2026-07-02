#' Inspect microbiome count data and return EDA diagnostics
#'
#' `microeda_check()` accepts either a `phyloseq` object or a plain count
#' matrix/data frame. Internally, counts are oriented as samples by features.
#'
#' @param x A `phyloseq` object, matrix, or data frame containing counts.
#' @param metadata Optional sample metadata. Required for non-phyloseq inputs
#'   when group-level diagnostics are desired.
#' @param taxonomy Optional taxonomy table with features as rows.
#' @param group Optional metadata column used for group-size diagnostics and
#'   broad group-aware screening notes.
#' @param taxa_are_rows For matrix/data frame inputs, whether rows are taxa.
#'   Ignored for `phyloseq` input because the orientation is read from the
#'   object.
#' @param min_prevalence Features below this sample prevalence are treated as
#'   low-prevalence features.
#' @param feature_read_n Number of features used for the first-N and top-N read
#'   summaries.
#'
#' @return A `microeda_report` object with diagnostics and broad screening notes.
#' @examples
#' counts <- matrix(c(10, 0, 0, 5, 20, 0, 1, 0), nrow = 2, byrow = TRUE)
#' rownames(counts) <- c("S1", "S2")
#' colnames(counts) <- paste0("ASV", 1:4)
#' microeda_check(counts, taxa_are_rows = FALSE)
#' @export
microeda_check <- function(x,
                           metadata = NULL,
                           taxonomy = NULL,
                           group = NULL,
                           taxa_are_rows = TRUE,
                           min_prevalence = 0.05,
                           feature_read_n = 50) {
  extracted <- microeda_extract(x, metadata, taxonomy, taxa_are_rows)

  diagnostics <- diagnose_microbiome_data(
    counts = extracted$counts,
    metadata = extracted$metadata,
    taxonomy = extracted$taxonomy,
    group = group,
    min_prevalence = min_prevalence,
    feature_read_n = feature_read_n
  )

  recommendations <- build_recommendations(diagnostics, group = group)

  structure(
    list(
      diagnostics = diagnostics,
      recommendations = recommendations,
      call = match.call()
    ),
    class = "microeda_report"
  )
}

#' Extract broad screening notes from a microeda report
#'
#' Extract the broad screening notes table from the result of
#' [microeda_check()]. These notes are general caveats, not contextual method
#' rankings or final workflow recommendations.
#'
#' @param x A `microeda_report` object.
#'
#' @return A data frame of recommendations.
#' @export
as_recommendations <- function(x) {
  if (!inherits(x, "microeda_report")) {
    stop("`x` must be a microeda_report.", call. = FALSE)
  }

  x$recommendations
}
