args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args) >= 1L) args[1L] else "."
output_dir <- if (length(args) >= 2L) args[2L] else file.path(package_dir, "man", "figures")
model_path <- file.path(package_dir, "inst", "extdata", "rf_fit.rds")
ml_fit <- readRDS(model_path)
scores <- ml_fit$votes[, "Positive"]
labels <- as.character(ml_fit$y) == "Positive"
stopifnot(length(scores) == length(labels), all(is.finite(scores)),
          all(scores >= 0 & scores <= 1), !anyNA(labels),
          all(ml_fit$oob.times > 0), any(labels), any(!labels))
ordering <- order(scores, decreasing = TRUE)
ordered_scores <- scores[ordering]
ordered_labels <- labels[ordering]
group_ends <- c(which(diff(ordered_scores) != 0), length(ordered_scores))
true_positive <- cumsum(ordered_labels)[group_ends]
false_positive <- cumsum(!ordered_labels)[group_ends]
n_positive <- sum(labels)
n_negative <- sum(!labels)
recall <- true_positive / n_positive
false_positive_rate <- false_positive / n_negative
precision <- true_positive / (true_positive + false_positive)
roc_x <- c(0, false_positive_rate)
roc_y <- c(0, recall)
auroc <- sum(diff(roc_x) * (head(roc_y, -1L) + tail(roc_y, -1L)) / 2)
rank_auroc <- (sum(rank(scores, ties.method = "average")[labels]) -
                n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
stopifnot(isTRUE(all.equal(auroc, rank_auroc, tolerance = 1e-12)))
average_precision <- sum(diff(c(0, recall)) * precision)
prevalence <- mean(labels)
metrics <- data.frame(
  metric = c("n", "positive", "negative", "prevalence", "auroc", "average_precision",
             "trees", "min_oob_trees", "max_oob_trees", "stored_oob_error"),
  value = c(length(labels), n_positive, n_negative, prevalence, auroc, average_precision,
            ml_fit$ntree, min(ml_fit$oob.times), max(ml_fit$oob.times),
            tail(ml_fit$err.rate[, "OOB"], 1L))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, file.path(output_dir, "oob_metrics.csv"), row.names = FALSE)
subtitle <- sprintf("%s sites  |  %s positive  |  %s negative",
                    format(length(labels), big.mark = ",", trim = TRUE),
                    format(n_positive, big.mark = ",", trim = TRUE),
                    format(n_negative, big.mark = ",", trim = TRUE))
footer <- "Stored out-of-bag votes | Internal assessment, not an independent test"
ink <- "#172B3A"
muted <- "#596A78"
grid_colour <- "#E7EDF1"
accent <- "#087F8C"
baseline <- "#8B99A5"
draw_frame <- function(x_label, y_label, heading) {
  par(mar = c(6.5, 5.6, 5.5, 1.6), mgp = c(3.1, 0.85, 0),
      family = "sans", col.axis = muted, col.lab = ink, col.main = ink,
      cex.axis = 1.05, cex.lab = 1.15, xaxs = "i", yaxs = "i", las = 1)
  plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(0, 1),
       xlab = x_label, ylab = y_label, axes = FALSE)
  abline(h = seq(0, 1, 0.2), v = seq(0, 1, 0.2), col = grid_colour, lwd = 1)
  axis(1, at = seq(0, 1, 0.2), col = baseline, col.ticks = baseline)
  axis(2, at = seq(0, 1, 0.2), col = baseline, col.ticks = baseline)
  box(col = baseline)
  title(main = heading, line = 3.0, adj = 0, cex.main = 1.45, font.main = 2)
  mtext(subtitle, side = 3, line = 1.35, adj = 0, col = muted, cex = 0.92)
  mtext(footer, side = 1, line = 5.0, col = muted, cex = 0.8)
}
png(file.path(output_dir, "oob-roc.png"), width = 1800, height = 1400, res = 180,
    bg = "white")
draw_frame("False positive rate (1 - specificity)", "True positive rate (sensitivity)",
           "OOB receiver operating characteristic")
abline(0, 1, col = baseline, lty = 2, lwd = 1.8)
lines(roc_x, roc_y, col = accent, lwd = 3)
legend("bottomright", inset = 0.035,
       legend = c(sprintf("AUROC = %.3f", auroc), "Chance baseline = 0.500"),
       col = c(accent, baseline), lty = c(1, 2), lwd = c(3, 1.8),
       bty = "n", text.col = ink, cex = 1.05, y.intersp = 1.3)
dev.off()
png(file.path(output_dir, "oob-prc.png"), width = 1800, height = 1400, res = 180,
    bg = "white")
draw_frame("Recall (sensitivity)", "Precision (positive predictive value)",
           "OOB precision-recall curve")
abline(h = prevalence, col = baseline, lty = 2, lwd = 1.8)
lines(c(0, recall), c(1, precision), col = accent, lwd = 2.7, type = "S")
legend("topright", inset = 0.035,
       legend = c(sprintf("Average precision = %.3f", average_precision),
                  sprintf("Prevalence baseline = %.3f", prevalence)),
       col = c(accent, baseline), lty = c(1, 2), lwd = c(2.7, 1.8),
       bty = "n", text.col = ink, cex = 1.05, y.intersp = 1.3)
dev.off()
print(metrics, row.names = FALSE, digits = 10)
