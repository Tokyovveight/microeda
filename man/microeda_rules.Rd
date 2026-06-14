\name{microeda_rules}
\alias{microeda_rules}
\title{Return the built-in evidence rules}
\usage{
microeda_rules()
}
\value{
A data frame containing rule metadata.
}
\description{
Return the built-in evidence rule table used by \code{microeda_check()}.
The same rule map is mirrored as \code{inst/extdata/evidence_rules.yml} so it
can be inspected or edited outside R.
}
