dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", seq_len(nn))
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  seq_df
}

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(all(c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  ) %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(
    feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  feature_df$RNA_region <- factor(
    feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )
  inferred_prob <- predict(
    ml_fit,
    newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)),
    type = "prob"
  )
  feature_df$predicted_m6A_prob <- inferred_prob[, "Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[
    feature_df$predicted_m6A_prob > positive_threshold
  ] <- "Positive"
  feature_df
}

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = data.frame(
      gc_content = gc_content,
      RNA_type = factor(
        RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
      ),
      RNA_region = factor(
        RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")
      ),
      exon_length = exon_length,
      distance_to_junction = distance_to_junction,
      evolutionary_conservation = evolutionary_conservation,
      DNA_5mer = DNA_5mer,
      stringsAsFactors = FALSE
    ),
    positive_threshold = positive_threshold
  )
  returned_vector <- c(
    pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status
  )
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  returned_vector
}
