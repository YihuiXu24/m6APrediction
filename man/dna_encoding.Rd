\name{dna_encoding}
\alias{dna_encoding}
\title{Encode DNA Contexts as Nucleotide Factors}
\description{
Converts equal-length DNA strings into one categorical predictor per nucleotide
position for the course random forest model. This is an internal helper.
}
\usage{
dna_encoding(dna_strings)
}
\arguments{
  \item{dna_strings}{A non-empty character vector of equal-length, uppercase DNA
  strings containing only \code{A}, \code{T}, \code{C}, and \code{G}, without
  missing values. For the bundled prediction model, each string must have five
  bases and a central \code{A}.}
}
\value{
A data frame with one row per input string and columns named \code{nt_pos1}
through \code{nt_posN}, where \code{N} is the string length. Each column is a
factor with levels \code{c("A", "T", "C", "G")} in that order.
}
\details{
This encoding uses nucleotide factors, not one-hot indicator columns. Input
order is retained. The helper does not repair invalid or unequal-length
sequences; callers must supply valid inputs. It is documented for completeness
but is not exported from the package namespace.
}
\seealso{\code{\link{prediction_multiple}}, \code{\link{prediction_single}}}
\examples{
m6APrediction:::dna_encoding(c("GGACA", "GAACC"))
}
\keyword{internal}
