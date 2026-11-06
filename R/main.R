#' Convert DNA sequences to a data frame
#' This internal use function takes a string of DNA bases and separates them
#' in different columns to prepare the data for the model.
#' @param dna_strings A string of DNA bases (e.g. "GGACA").
#' @return A data frame where each column contains a base letter.
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A sites for a table of samples
#' This function takes a table of features, puts them into a trained ML (machine
#' learning) model, and calculates the probability to be an m6A site.
#' @param ml_fit The trained model.
#' @param feature_df A data frame with the information to predict.
#' @param positive_threshold A cutoff number (default is 0.5). If the probability is higher, it prints "Positive".
#' @return The original table containing two new columns with the prediction.
#' @examples
#' model_origin <- system.file("extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE)
#' ml_fit <- readRDS(model_origin)
#' file_origin <- system.file("extdata", "m6A_input_example.csv",
#'                             package = "m6APrediction", mustWork = TRUE)
#' example_df <- read.csv(file_origin)
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df)
}

#' Predict an m6A site for just one sample
#' This predicts a single case in which you may provide the
#' numbers and letters, and it gives an answer by calling the
#' table function behind the scenes.
#' @param ml_fit The trained model.
#' @param gc_content A number for the GC content.
#' @param RNA_type The type of RNA (e.g. "mRNA").
#' @param RNA_region The part of the RNA (e.g. "CDS").
#' @param exon_length A number regarding how long the exon is.
#' @param distance_to_junction A number regarding the distance.
#' @param evolutionary_conservation A number regarding how conserved it is.
#' @param DNA_5mer The 5-letter DNA sequence (e.g. "GGACA").
#' @param positive_threshold The cutoff number (default is 0.5).
#' @return An answer including the probability and if it's "Positive" or "Negative".
#' @examples
#' model_origin <- system.file("extdata", "rf_fit.rds", package = "m6APrediction", mustWork = TRUE)
#' ml_fit <- readRDS(model_origin)
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}
