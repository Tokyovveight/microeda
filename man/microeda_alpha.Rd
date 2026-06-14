\name{microeda_alpha}
\alias{microeda_alpha}
\title{Compute alpha diversity indices}
\usage{
microeda_alpha(x, metadata = NULL, group = NULL, taxa_are_rows = TRUE)
}
\arguments{
\item{x}{A \code{phyloseq} object, matrix, or data frame containing counts.}
\item{metadata}{Optional sample metadata. Required for non-phyloseq inputs when group-level summaries are desired.}
\item{group}{Optional metadata column used for group summaries.}
\item{taxa_are_rows}{For matrix/data frame inputs, whether rows are taxa. Ignored for \code{phyloseq} input because the orientation is read from the object.}
}
\value{
A \code{microeda_alpha} object with per-sample indices and optional group summaries.
}
\description{
\code{microeda_alpha()} calculates per-sample alpha diversity metrics from raw
count-like microbiome tables. It reports classic indices and their Hill number
interpretation so Shannon and Simpson can be read as effective taxa.
}
\examples{
counts <- matrix(c(10, 0, 0, 5, 20, 0, 1, 0), nrow = 2, byrow = TRUE)
rownames(counts) <- c("S1", "S2")
colnames(counts) <- paste0("ASV", 1:4)
alpha <- microeda_alpha(counts, taxa_are_rows = FALSE)
as_alpha_table(alpha)
}
