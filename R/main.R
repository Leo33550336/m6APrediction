#' Encode DNA sequences by position
#'
#' Convert equal-length DNA strings into nucleotide factors for each position.
#'
#' @param dna_strings A character vector of equal-length DNA sequences
#'   containing A, T, C and G.
#' @return A data frame with one factor column per nucleotide position.
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A status for multiple sites
#'
#' Predict scores and statuses for multiple input rows using a fitted model.
#'
#' @param ml_fit A fitted random forest model.
#' @param feature_df A data frame containing columns gc_content, RNA_type,
#'   RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation, and DNA_5mer. Length and distance values
#'   should already be log2-transformed.
#' @param positive_threshold Numeric classification threshold; defaults to 0.5.
#'   Scores strictly greater than this value are classified as Positive.
#' @return The input data frame with predicted_m6A_prob and
#'   predicted_m6A_status columns added.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, head(example_df, 3))
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}
#' Predict m6A status for a single site
#'
#' Predict the score and status for one site using a fitted model.
#'
#' @param ml_fit A fitted random forest model.
#' @param gc_content Numeric GC content.
#' @param RNA_type RNA type: mRNA, lincRNA, lncRNA, or pseudogene.
#' @param RNA_region RNA region: CDS, intron, 3'UTR, or 5'UTR.
#' @param exon_length Numeric, log2-transformed exon length.
#' @param distance_to_junction Numeric, log2-transformed distance to junction.
#' @param evolutionary_conservation Numeric evolutionary conservation value.
#' @param DNA_5mer A five-base DNA sequence, such as GGACA.
#' @param positive_threshold Numeric classification threshold; defaults to 0.5.
#'   Scores strictly greater than this value are classified as Positive.
#' @return A named character vector containing predicted_m6A_prob and
#'   predicted_m6A_status.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA")
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
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
