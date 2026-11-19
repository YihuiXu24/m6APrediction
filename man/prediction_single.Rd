\name{prediction_single}
\alias{prediction_single}
\title{Predict an m6A Score for a Single Candidate Site}
\description{
Builds a one-row feature table from individual feature values and uses the
fitted random forest model to predict an m6A score and status.
}
\usage{
prediction_single(ml_fit, gc_content, RNA_type, RNA_region,
  exon_length, distance_to_junction, evolutionary_conservation,
  DNA_5mer, positive_threshold = 0.5)
}
\arguments{
  \item{ml_fit}{A fitted \code{randomForest} classification model with a
  \code{Positive} class and the same feature encoding as the bundled model.
  The example below loads the supplied \code{rf_fit.rds} resource.}
  \item{gc_content}{One numeric GC-content proportion between 0 and 1.}
  \item{RNA_type}{One RNA type: \code{mRNA}, \code{lincRNA}, \code{lncRNA}, or
  \code{pseudogene}.}
  \item{RNA_region}{One RNA region: \code{CDS}, \code{intron}, \code{3'UTR}, or
  \code{5'UTR}.}
  \item{exon_length}{One numeric exon-length feature already on the log2 scale
  used during model training.}
  \item{distance_to_junction}{One numeric distance-to-junction feature already
  on the log2 scale used during model training.}
  \item{evolutionary_conservation}{One numeric conservation score on the same
  scale as the course training data.}
  \item{DNA_5mer}{One character string of exactly five uppercase DNA bases
  from \code{A}, \code{T}, \code{C}, and \code{G}, with a central \code{A}; for
  example, \code{"GGACA"}.}
  \item{positive_threshold}{One numeric cutoff between 0 and 1, default
  \code{0.5}. Only a score strictly greater than this cutoff is
  \code{Positive}; equality is \code{Negative}.}
}
\value{
A named character vector of length two. The element
\code{predicted_m6A_prob} is the positive-class vote proportion converted to
text; \code{predicted_m6A_status} is \code{Positive} or \code{Negative}. Use
\code{as.numeric(result["predicted_m6A_prob"])} if a numeric score is needed.
}
\details{
All feature arguments must be non-missing scalar values. The function calls
\code{\link{prediction_multiple}} and applies the same categorical encoding
and strict threshold rule. It does not transform the numeric inputs. In
particular, \code{exon_length = 12} and \code{distance_to_junction = 5} in the
example are already transformed feature values, not raw base-pair lengths.

Random forest vote proportions are model scores and need not be calibrated
biological probabilities. The function is intended for educational
exploration; a predicted status is not experimental confirmation of m6A.
}
\seealso{\code{\link{prediction_multiple}}}
\examples{
library(m6APrediction)
ml_fit <- readRDS(system.file(
  "extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE
))
result <- prediction_single(
  ml_fit,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "GGACA",
  positive_threshold = 0.5
)
result
as.numeric(result["predicted_m6A_prob"])
}
