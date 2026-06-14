\name{microeda_check}
\alias{microeda_check}
\title{Inspect microbiome count data and return evidence-guided EDA advice}
\usage{
microeda_check(
  x,
  metadata = NULL,
  taxonomy = NULL,
  group = NULL,
  taxa_are_rows = TRUE,
  min_prevalence = 0.05,
  feature_read_n = 50
)
}
\arguments{
\item{x}{A \code{phyloseq} object, matrix, or data frame containing counts.}
\item{metadata}{Optional sample metadata. Required for non-phyloseq inputs when group-level diagnostics are desired.}
\item{taxonomy}{Optional taxonomy table with features as rows.}
\item{group}{Optional metadata column used for group-size diagnostics and group-aware recommendations.}
\item{taxa_are_rows}{For matrix/data frame inputs, whether rows are taxa. Ignored for \code{phyloseq} input because the orientation is read from the object.}
\item{min_prevalence}{Features below this sample prevalence are treated as low-prevalence features.}
\item{feature_read_n}{Number of features used for the first-N and top-N read summaries.}
}
\value{
A \code{microeda_report} object with diagnostics and recommendations.
}
\description{
\code{microeda_check()} accepts either a \code{phyloseq} object or a plain count matrix/data frame. Internally, counts are oriented as samples by features.
}
\examples{
counts <- matrix(c(10, 0, 0, 5, 20, 0, 1, 0), nrow = 2, byrow = TRUE)
rownames(counts) <- c("S1", "S2")
colnames(counts) <- paste0("ASV", 1:4)
microeda_check(counts, taxa_are_rows = FALSE)
}
