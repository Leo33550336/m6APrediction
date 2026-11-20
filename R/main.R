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
#'   DNA_5mer values must contain five uppercase A/T/C/G bases with
#'   A at the central position. Required features must not contain missing values.
#' @param positive_threshold A single finite numeric value between 0 and 1;
#'   defaults to 0.5. Scores strictly greater than this value are Positive.
#' @return The input data frame with predicted_m6A_prob and
#'   predicted_m6A_status columns added.
#' @import randomForest
#' @importFrom stats predict
#' @details Input must be a non-empty data frame. RNA_type must be mRNA,
#'   lincRNA, lncRNA, or pseudogene; RNA_region must be CDS, intron,
#'   3'UTR, or 5'UTR. gc_content, exon_length, distance_to_junction,
#'   and evolutionary_conservation must contain finite numeric values.
#'   DNA_5mer may be character or factor. Additional input columns are
#'   preserved in the output but are not used as model predictors.
#' @details Invalid inputs cause an error. Prediction also stops if the
#'   model returns a missing or non-finite score.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, head(example_df, 3))
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  if (!is.numeric(positive_threshold) ||
      length(positive_threshold) != 1L ||
      !is.finite(positive_threshold) ||
      positive_threshold < 0 || positive_threshold > 1) {
    stop("positive_threshold must be one finite number between 0 and 1.")
  }

  required <- c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  )

  if (!is.data.frame(feature_df) || nrow(feature_df) == 0L) {
    stop("feature_df must be a non-empty data frame.")
  }

  if (!all(required %in% names(feature_df))) {
    stop("feature_df is missing one or more required columns.")
  }

  if (anyNA(feature_df[, required, drop = FALSE])) {
    stop("Required input features must not contain missing values.")
  }

  sequences <- as.character(feature_df$DNA_5mer)
  if (any(!grepl("^[ATCG]{5}$", sequences))) {
    stop("DNA_5mer must contain exactly five uppercase A/T/C/G bases.")
  }

  if (any(substr(sequences, 3, 3) != "A")) {
    stop("DNA_5mer must have A at the central position.")
  }

  if (any(!feature_df$RNA_type %in%
          c("mRNA", "lincRNA", "lncRNA", "pseudogene"))) {
    stop("RNA_type contains an unsupported category.")
  }

  if (any(!feature_df$RNA_region %in%
          c("CDS", "intron", "3'UTR", "5'UTR"))) {
    stop("RNA_region contains an unsupported category.")
  }

  numeric_features <- c(
    "gc_content", "exon_length", "distance_to_junction",
    "evolutionary_conservation"
  )

  valid_numeric <- vapply(
    feature_df[, numeric_features, drop = FALSE],
    function(x) is.numeric(x) && all(is.finite(x)),
    logical(1)
  )

  if (!all(valid_numeric)) {
    stop("Numeric input features must contain finite numeric values.")
  }
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  model_input <- cbind(
    feature_df[, setdiff(required, "DNA_5mer"), drop = FALSE],
    dna_encoding(sequences)
  )

  infered_prob <- predict(
    ml_fit, newdata = model_input, type = "prob"
  )
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  if (any(!is.finite(feature_df$predicted_m6A_prob))) {
    stop("Prediction produced a missing or non-finite score; no classification returned.")
  }
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
#' @param DNA_5mer A five-base uppercase DNA sequence containing only
#'   A/T/C/G, with A at the central position, such as GGACA.
#' @param positive_threshold A single finite numeric value between 0 and 1;
#'   defaults to 0.5. Scores strictly greater than this value are Positive.
#' @return A named character vector containing predicted_m6A_prob and
#'   predicted_m6A_status.
#' @details Each site input argument must contain exactly one value.
#'   Use prediction_multiple() for multiple sites.
#' @details Invalid inputs cause an error. Prediction also stops if the
#'   model returns a missing or non-finite score.
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
  site_inputs <- list(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )

  if (any(lengths(site_inputs) != 1L)) {
    stop("Each site input must contain exactly one value.")
  }
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
