\name{prediction_multiple}
\alias{prediction_multiple}
\title{Predict m6A Scores for Multiple Candidate Sites}
\description{
Applies a fitted random forest classifier to a feature table and returns a
positive-class vote proportion and threshold-based status for each candidate
m6A site.
}
\usage{
prediction_multiple(ml_fit, feature_df, positive_threshold = 0.5)
}
\arguments{
  \item{ml_fit}{A fitted classification model of class \code{randomForest} with
  a class named \code{Positive}, trained using the same features and factor
  levels as the bundled course model. Load the supplied model with
  \code{readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction",
  mustWork = TRUE))}.}
  \item{feature_df}{A non-empty data frame with one candidate site per row and
  these seven required columns:
  \describe{
    \item{\code{gc_content}}{Numeric GC-content proportion between 0 and 1.}
    \item{\code{RNA_type}}{Character or factor values from \code{mRNA},
    \code{lincRNA}, \code{lncRNA}, and \code{pseudogene}.}
    \item{\code{RNA_region}}{Character or factor values from \code{CDS},
    \code{intron}, \code{3'UTR}, and \code{5'UTR}.}
    \item{\code{exon_length}}{Numeric exon-length feature already transformed
    to the log2 scale used during model training.}
    \item{\code{distance_to_junction}}{Numeric distance-to-junction feature
    already transformed to the log2 scale used during model training.}
    \item{\code{evolutionary_conservation}}{Numeric conservation score on the
    same scale as the course training data.}
    \item{\code{DNA_5mer}}{Character strings of exactly five uppercase DNA
    bases from \code{A}, \code{T}, \code{C}, and \code{G}, with \code{A} at the
    central position.}
  }
  Required features must not contain missing values. Additional columns are
  retained in the returned table.}
  \item{positive_threshold}{A single numeric cutoff between 0 and 1; the
  default is \code{0.5}. A score strictly greater than the cutoff is labelled
  \code{Positive}; a score equal to or below it is labelled \code{Negative}.}
}
\value{
A data frame preserving input row order and columns, with \code{RNA_type} and
\code{RNA_region} converted to factors with the documented levels, and two
added columns: \code{predicted_m6A_prob}, a numeric positive-class vote
proportion, and \code{predicted_m6A_status}, a character vector containing
\code{Positive} or \code{Negative}. Existing columns with either prediction
name are replaced.
}
\details{
The DNA context is encoded as five categorical nucleotide predictors. The
function does not apply a log transformation to the numeric features: supply
the same preprocessing used during training. Do not pass raw lengths or
distances in place of the transformed course features.

The function checks required column names; the caller is responsible for the
remaining input constraints. The returned score is obtained from
\code{predict(ml_fit, type = "prob")} and is not necessarily a calibrated
biological probability. Predictions are for educational exploration and do
not establish that a site is experimentally methylated.

The bundled 100-row feature table demonstrates software use. It has no outcome
labels and must not be treated as an independent accuracy benchmark.
}
\seealso{\code{\link{prediction_single}},
\code{\link[randomForest:predict.randomForest]{predict.randomForest}}}
\examples{
library(m6APrediction)
ml_fit <- readRDS(system.file(
  "extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE
))
feature_df <- read.csv(system.file(
  "extdata", "m6A_input_example.csv",
  package = "m6APrediction", mustWork = TRUE
), stringsAsFactors = FALSE)
result <- prediction_multiple(
  ml_fit, feature_df[1:3, ], positive_threshold = 0.6
)
result[, c("predicted_m6A_prob", "predicted_m6A_status")]
}
