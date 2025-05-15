#' Z-score normalization for omics data
#'
#' @description Apply Z-score normalization to omics data at the gene level
#' @param mat A matrix or data frame with genes in rows and samples in columns
#' @return A normalized matrix with Z-scores
#' @export
zscoreNormalize <- function(mat) {
  # Apply z-score normalization to each row (gene)
  normalized <- t(scale(t(mat)))
  # Handle any potential NaN values (e.g., if SD was 0)
  normalized[is.nan(normalized)] <- 0
  return(normalized)
}

# Commented out alternative implementation
# zscoreNormalizem <- function(mat) {
#   # Apply z-score normalization to each row (drug) independently
#   # and then transform to 0-1 scale centered at 0.5
#   # Initialize the result matrix with the same dimensions as the input
#   normalized <- as.matrix(mat)
#
#   # Process each row (drug) independently
#   for (i in 1:nrow(normalized)) {
#     # Get the current row
#     row_data <- normalized[i, ]
#
#     # Calculate mean and standard deviation for this drug
#     row_mean <- mean(row_data, na.rm = TRUE)
#     row_sd <- sd(row_data, na.rm = TRUE)
#
#     # Apply z-score normalization only if SD is not zero
#     if (row_sd > 0) {
#       # First apply z-score normalization
#       z_scores <- (row_data - row_mean) / row_sd
#
#       # Then transform z-scores to 0-1 scale centered at 0.5
#       # We'll use a sigmoid-like transformation
#       # This will map most values to the 0-1 range with 0.5 as center
#       normalized[i, ] <- 1 / (1 + exp(-z_scores))
#     } else {
#       # If SD is zero, all values are the same, so set to 0.5
#       normalized[i, ] <- 0.5
#     }
#     print(i)
#   }
#
#   return(normalized)
# }

#' Apply Z-score normalization to all continuous data
#'
#' @description Applies Z-score normalization to all loaded omics and drug datasets in the global environment
#' @return None, directly modifies objects in the global environment
#' @export
applyZscoreNormalization <- function() {
  # mRNA datasets
  if (exists("ccle_mRNA", envir = .GlobalEnv)) assign("ccle_mRNA", zscoreNormalize(base::get("ccle_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc_mRNA", envir = .GlobalEnv)) assign("gdsc_mRNA", zscoreNormalize(base::get("gdsc_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("NCI60_mRNA", envir = .GlobalEnv)) assign("NCI60_mRNA", zscoreNormalize(base::get("NCI60_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO1_mRNA", envir = .GlobalEnv)) assign("UMPDO1_mRNA", zscoreNormalize(base::get("UMPDO1_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO2_mRNA", envir = .GlobalEnv)) assign("UMPDO2_mRNA", zscoreNormalize(base::get("UMPDO2_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO3_mRNA", envir = .GlobalEnv)) assign("UMPDO3_mRNA", zscoreNormalize(base::get("UMPDO3_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("tavor_mRNA", envir = .GlobalEnv)) assign("tavor_mRNA", zscoreNormalize(base::get("tavor_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("Xeva_mRNA", envir = .GlobalEnv)) assign("Xeva_mRNA", zscoreNormalize(base::get("Xeva_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)

  # CNV datasets
  if (exists("ccle_cnv", envir = .GlobalEnv)) assign("ccle_cnv", zscoreNormalize(base::get("ccle_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc_cnv", envir = .GlobalEnv)) assign("gdsc_cnv", zscoreNormalize(base::get("gdsc_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gCSI_cnv", envir = .GlobalEnv)) assign("gCSI_cnv", zscoreNormalize(base::get("gCSI_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("Xeva_cnv", envir = .GlobalEnv)) assign("Xeva_cnv", zscoreNormalize(base::get("Xeva_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Methylation dataset
  if (exists("ccle_meth", envir = .GlobalEnv)) assign("ccle_meth", zscoreNormalize(base::get("ccle_meth", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Protein datasets
  if (exists("ccle_proteinms", envir = .GlobalEnv)) assign("ccle_proteinms", zscoreNormalize(base::get("ccle_proteinms", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ccle_proteinrppa", envir = .GlobalEnv)) assign("ccle_proteinrppa", zscoreNormalize(base::get("ccle_proteinrppa", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Drug datasets
  # PDC
  if (exists("tavor_drug", envir = .GlobalEnv)) assign("tavor_drug", zscoreNormalize(base::get("tavor_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("PDTXBreast_drug", envir = .GlobalEnv)) assign("PDTXBreast_drug", zscoreNormalize(base::get("PDTXBreast_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # CellLine
  if (exists("ccle_drug", envir = .GlobalEnv)) assign("ccle_drug", zscoreNormalize(base::get("ccle_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ctrp1_drug", envir = .GlobalEnv)) assign("ctrp1_drug", zscoreNormalize(base::get("ctrp1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ctrp2_drug", envir = .GlobalEnv)) assign("ctrp2_drug", zscoreNormalize(base::get("ctrp2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc1_drug", envir = .GlobalEnv)) assign("gdsc1_drug", zscoreNormalize(base::get("gdsc1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc2_drug", envir = .GlobalEnv)) assign("gdsc2_drug", zscoreNormalize(base::get("gdsc2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gCSI_drug", envir = .GlobalEnv)) assign("gCSI_drug", zscoreNormalize(base::get("gCSI_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("prism_drug", envir = .GlobalEnv)) assign("prism_drug", zscoreNormalize(base::get("prism_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("FIMM_drug", envir = .GlobalEnv)) assign("FIMM_drug", zscoreNormalize(base::get("FIMM_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UHNBreast_drug", envir = .GlobalEnv)) assign("UHNBreast_drug", zscoreNormalize(base::get("UHNBreast_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("GRAY_drug", envir = .GlobalEnv)) assign("GRAY_drug", zscoreNormalize(base::get("GRAY_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("NCI60_drug", envir = .GlobalEnv)) assign("NCI60_drug", zscoreNormalize(base::get("NCI60_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # PDO
  if (exists("UMPDO1_drug", envir = .GlobalEnv)) assign("UMPDO1_drug", zscoreNormalize(base::get("UMPDO1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO2_drug", envir = .GlobalEnv)) assign("UMPDO2_drug", zscoreNormalize(base::get("UMPDO2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO3_drug", envir = .GlobalEnv)) assign("UMPDO3_drug", zscoreNormalize(base::get("UMPDO3_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # PDX
  if (exists("Xeva_drug", envir = .GlobalEnv)) assign("Xeva_drug", zscoreNormalize(base::get("Xeva_drug", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Set normalization state
  assign("normalization_state", TRUE, envir = .GlobalEnv)
}

