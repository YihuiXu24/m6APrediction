library(m6APrediction)

ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE))
feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction", mustWork = TRUE), stringsAsFactors = FALSE)

batch <- prediction_multiple(ml_fit, feature_df, positive_threshold = 0.6)
stopifnot(
  nrow(batch) == 100L,
  identical(names(batch), c(names(feature_df), "predicted_m6A_prob", "predicted_m6A_status")),
  isTRUE(all.equal(batch$predicted_m6A_prob[1:3], c(0.876, 0, 0.498))),
  identical(batch$predicted_m6A_status[1:3], c("Positive", "Negative", "Negative")),
  all(is.finite(batch$predicted_m6A_prob)),
  all(batch$predicted_m6A_prob >= 0 & batch$predicted_m6A_prob <= 1),
  identical(batch$predicted_m6A_status, ifelse(batch$predicted_m6A_prob > 0.6, "Positive", "Negative"))
)

single <- prediction_single(
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
stopifnot(identical(single, c(predicted_m6A_prob = "0.51", predicted_m6A_status = "Positive")))

equal_threshold <- prediction_multiple(ml_fit, feature_df[1, , drop = FALSE], positive_threshold = 0.876)
stopifnot(identical(equal_threshold$predicted_m6A_status, "Negative"))

for (i in seq_len(10L)) {
  row <- feature_df[i, , drop = FALSE]
  one <- do.call(prediction_single, c(list(ml_fit = ml_fit), as.list(row), list(positive_threshold = 0.6)))
  stopifnot(
    isTRUE(all.equal(as.numeric(one[["predicted_m6A_prob"]]), batch$predicted_m6A_prob[i])),
    identical(one[["predicted_m6A_status"]], batch$predicted_m6A_status[i])
  )
}

stopifnot(setequal(getNamespaceExports("m6APrediction"), c("prediction_multiple", "prediction_single")))
