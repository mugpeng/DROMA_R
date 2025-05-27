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
