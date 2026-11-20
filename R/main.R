#' Batch prediction for multiple m6A candidate sites
#'
#' This function takes a pre-trained random forest model and a feature table,
#' returns prediction probability and classification result.
#'
#' @param ml_fit Fitted random forest model object
#' @param feature_df Input feature data frame
#' @param positive_threshold Numeric threshold for classification (0~1)
#' @return Data frame with original features + prediction probability and label
#' @export
#' @import randomForest
#'
#' @examples
#' \dontrun{
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#' pred_result <- prediction_multiple(rf_fit, input_df, positive_threshold = 0.5)
#' }
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  required_cols <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                     "distance_to_splice", "evolutionary_conservation", "DNA_5mer")
  if (!all(required_cols %in% colnames(feature_df))) {
    stop("feature_df is missing required columns!")
  }
  if (!is.numeric(positive_threshold) || positive_threshold <0 || positive_threshold >1) {
    stop("positive_threshold must be a numeric value between 0 and 1!")
  }

  pred_prob <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]
  pred_label <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")
  out_df <- feature_df
  out_df$pred_prob <- pred_prob
  out_df$pred_label <- pred_label
  return(out_df)
}

#' Single-site prediction for one m6A candidate
#'
#' Predict m6A modification for one single site
#'
#' @param ml_fit Fitted random forest model
#' @param gc_content Numeric gc content
#' @param RNA_type Character
#' @param RNA_region Character
#' @param exon_length Numeric
#' @param distance_to_splice Numeric
#' @param evolutionary_conservation Numeric
#' @param DNA_5mer Character
#' @param positive_threshold Numeric threshold (0~1)
#' @return List with probability and prediction label
#' @export
#' @import randomForest
#'
#' @examples
#' \dontrun{
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_splice = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#' }
prediction_single <- function(ml_fit,
                              gc_content,
                              RNA_type,
                              RNA_region,
                              exon_length,
                              distance_to_splice,
                              evolutionary_conservation,
                              DNA_5mer,
                              positive_threshold = 0.5) {
  if (!is.numeric(gc_content) || !is.numeric(exon_length) || !is.numeric(distance_to_splice) || !is.numeric(evolutionary_conservation)) {
    stop("Numeric inputs must be numeric values!")
  }
  if (!is.character(RNA_type) || !is.character(RNA_region) || !is.character(DNA_5mer)) {
    stop("RNA_type, RNA_region and DNA_5mer must be character strings!")
  }
  if (!is.numeric(positive_threshold) || positive_threshold <0 || positive_threshold >1) {
    stop("positive_threshold must be a numeric value between 0 and 1!")
  }

  new_row <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_splice = distance_to_splice,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  pred_prob <- predict(ml_fit, newdata = new_row, type = "prob")[, "Positive"]
  pred_label <- ifelse(pred_prob >= positive_threshold, "Positive", "Negative")
  res <- list(probability = pred_prob, label = pred_label)
  return(res)
}
