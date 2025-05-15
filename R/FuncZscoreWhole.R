# Z-score normalization function for omics data (gene-level)
zscore_normalize <- function(mat) {
  # Apply z-score normalization to each row (gene)
  normalized <- t(scale(t(mat)))
  # Handle any potential NaN values (e.g., if SD was 0)
  normalized[is.nan(normalized)] <- 0
  return(normalized)
}

# Z-score normalization function for drug data (drug-level)
zscore_normalize_drug <- function(mat) {
  # Apply z-score normalization to each row (drug) independently
  # Initialize the result matrix with the same dimensions as the input
  normalized <- as.matrix(mat)
  
  # Process each row (drug) independently
  for (i in 1:nrow(normalized)) {
    # Get the current row
    row_data <- normalized[i, ]
    
    # Calculate mean and standard deviation for this drug
    row_mean <- mean(row_data, na.rm = TRUE)
    row_sd <- sd(row_data, na.rm = TRUE)
    
    # Apply z-score normalization only if SD is not zero
    if (row_sd > 0) {
      normalized[i, ] <- (row_data - row_mean) / row_sd
    } else {
      # If SD is zero, just center the data
      normalized[i, ] <- row_data - row_mean
    }
  }
  return(as.data.frame(normalized))
}

# zscore_normalize_drug <- function(mat) {
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

# Function to apply z-score normalization to all continuous data
apply_zscore_normalization <- function() {
  # mRNA datasets
  if (exists("ccle_mRNA", envir = .GlobalEnv)) assign("ccle_mRNA", zscore_normalize(base::get("ccle_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc_mRNA", envir = .GlobalEnv)) assign("gdsc_mRNA", zscore_normalize(base::get("gdsc_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("NCI60_mRNA", envir = .GlobalEnv)) assign("NCI60_mRNA", zscore_normalize(base::get("NCI60_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO1_mRNA", envir = .GlobalEnv)) assign("UMPDO1_mRNA", zscore_normalize(base::get("UMPDO1_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO2_mRNA", envir = .GlobalEnv)) assign("UMPDO2_mRNA", zscore_normalize(base::get("UMPDO2_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO3_mRNA", envir = .GlobalEnv)) assign("UMPDO3_mRNA", zscore_normalize(base::get("UMPDO3_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("tavor_mRNA", envir = .GlobalEnv)) assign("tavor_mRNA", zscore_normalize(base::get("tavor_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("Xeva_mRNA", envir = .GlobalEnv)) assign("Xeva_mRNA", zscore_normalize(base::get("Xeva_mRNA", envir = .GlobalEnv)), envir = .GlobalEnv)
  
  # CNV datasets
  if (exists("ccle_cnv", envir = .GlobalEnv)) assign("ccle_cnv", zscore_normalize(base::get("ccle_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc_cnv", envir = .GlobalEnv)) assign("gdsc_cnv", zscore_normalize(base::get("gdsc_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gCSI_cnv", envir = .GlobalEnv)) assign("gCSI_cnv", zscore_normalize(base::get("gCSI_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("Xeva_cnv", envir = .GlobalEnv)) assign("Xeva_cnv", zscore_normalize(base::get("Xeva_cnv", envir = .GlobalEnv)), envir = .GlobalEnv)
  
  # Methylation dataset
  if (exists("ccle_meth", envir = .GlobalEnv)) assign("ccle_meth", zscore_normalize(base::get("ccle_meth", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Protein datasets
  if (exists("ccle_proteinms", envir = .GlobalEnv)) assign("ccle_proteinms", zscore_normalize(base::get("ccle_proteinms", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ccle_proteinrppa", envir = .GlobalEnv)) assign("ccle_proteinrppa", zscore_normalize(base::get("ccle_proteinrppa", envir = .GlobalEnv)), envir = .GlobalEnv)

  # Drug datasets
  # PDC
  if (exists("tavor_drug", envir = .GlobalEnv)) assign("tavor_drug", zscore_normalize_drug(base::get("tavor_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("PDTXBreast_drug", envir = .GlobalEnv)) assign("PDTXBreast_drug", zscore_normalize_drug(base::get("PDTXBreast_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # CellLine
  if (exists("ccle_drug", envir = .GlobalEnv)) assign("ccle_drug", zscore_normalize_drug(base::get("ccle_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ctrp1_drug", envir = .GlobalEnv)) assign("ctrp1_drug", zscore_normalize_drug(base::get("ctrp1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("ctrp2_drug", envir = .GlobalEnv)) assign("ctrp2_drug", zscore_normalize_drug(base::get("ctrp2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc1_drug", envir = .GlobalEnv)) assign("gdsc1_drug", zscore_normalize_drug(base::get("gdsc1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gdsc2_drug", envir = .GlobalEnv)) assign("gdsc2_drug", zscore_normalize_drug(base::get("gdsc2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("gCSI_drug", envir = .GlobalEnv)) assign("gCSI_drug", zscore_normalize_drug(base::get("gCSI_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("prism_drug", envir = .GlobalEnv)) assign("prism_drug", zscore_normalize_drug(base::get("prism_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("FIMM_drug", envir = .GlobalEnv)) assign("FIMM_drug", zscore_normalize_drug(base::get("FIMM_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UHNBreast_drug", envir = .GlobalEnv)) assign("UHNBreast_drug", zscore_normalize_drug(base::get("UHNBreast_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("GRAY_drug", envir = .GlobalEnv)) assign("GRAY_drug", zscore_normalize_drug(base::get("GRAY_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("NCI60_drug", envir = .GlobalEnv)) assign("NCI60_drug", zscore_normalize_drug(base::get("NCI60_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # PDO
  if (exists("UMPDO1_drug", envir = .GlobalEnv)) assign("UMPDO1_drug", zscore_normalize_drug(base::get("UMPDO1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO2_drug", envir = .GlobalEnv)) assign("UMPDO2_drug", zscore_normalize_drug(base::get("UMPDO2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  if (exists("UMPDO3_drug", envir = .GlobalEnv)) assign("UMPDO3_drug", zscore_normalize_drug(base::get("UMPDO3_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  # PDX
  if (exists("Xeva_drug", envir = .GlobalEnv)) assign("Xeva_drug", zscore_normalize_drug(base::get("Xeva_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
  
  # Set normalization state
  assign("normalization_state", TRUE, envir = .GlobalEnv)
}

# Function to reset data to original values
reset_to_original <- function() {
  # Reload original data
  source("Modules/LoadData_2.R", local = FALSE)
  # Set normalization state
  assign("normalization_state", FALSE, envir = .GlobalEnv)
} 