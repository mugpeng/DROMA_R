# Data Loading Functions with Z-score Normalization ----

#' Load molecular profiles with optional z-score normalization
#'
#' @description Wrapper function for loadMolecularProfiles that applies z-score normalization by default
#' @param dromaset_object A DromaSet object
#' @param molecular_type Type of molecular data to load (e.g., "mRNA", "cnv", "meth", "proteinrppa", "proteinms")
#' @param features Optional vector of specific features to load. If NULL, loads all features
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @param return_data Logical, whether to return the data directly (default: TRUE)
#' @return Matrix or data frame with molecular profile data, optionally z-score normalized
#' @export
#' @examples
#' \dontrun{
#' # Load mRNA data with z-score normalization (default)
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' mrna_data <- loadMolecularProfilesNormalized(gCSI, "mRNA", features = "ABCB1")
#'
#' # Load without z-score normalization
#' mrna_raw <- loadMolecularProfilesNormalized(gCSI, "mRNA", features = "ABCB1", zscore = FALSE)
#'
#' # Load all mRNA features with normalization
#' all_mrna <- loadMolecularProfilesNormalized(gCSI, "mRNA")
#' }
loadMolecularProfilesNormalized <- function(dromaset_object,
                                          molecular_type,
                                          features = NULL,
                                          data_type = "all",
                                          tumor_type = "all",
                                          zscore = TRUE,
                                          return_data = TRUE) {

  # Validate input object
  if (!inherits(dromaset_object, "DromaSet")) {
    stop("Input must be a DromaSet object from DROMA.Set package")
  }

  # Load molecular profiles using DROMA.Set function
  molecular_data <- loadMolecularProfiles(
    dromaset_object,
    molecular_type = molecular_type,
    features = features,
    data_type = data_type,
    tumor_type = tumor_type,
    return_data = return_data
  )

  # Apply z-score normalization if requested and data is continuous
  if (zscore && !is.null(molecular_data)) {
    # Check if data is continuous (matrix format) and suitable for normalization
    continuous_types <- c("mRNA", "cnv", "meth", "proteinrppa", "proteinms")

    if (molecular_type %in% continuous_types && is.matrix(molecular_data)) {
      # Apply z-score normalization
      molecular_data <- zscoreNormalize(molecular_data)

      # Add attribute to indicate normalization was applied
      attr(molecular_data, "zscore_normalized") <- TRUE
    } else if (zscore && !molecular_type %in% continuous_types) {
      warning(paste("Z-score normalization not applicable for molecular type:", molecular_type))
    }
  }

  return(molecular_data)
}

#' Load treatment response data with optional z-score normalization
#'
#' @description Wrapper function for loadTreatmentResponse that applies z-score normalization by default
#' @param dromaset_object A DromaSet object
#' @param drugs Optional vector of specific drugs to load. If NULL, loads all drugs
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @param return_data Logical, whether to return the data directly (default: TRUE)
#' @return Matrix with treatment response data, optionally z-score normalized
#' @export
#' @examples
#' \dontrun{
#' # Load drug data with z-score normalization (default)
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' drug_data <- loadTreatmentResponseNormalized(gCSI, drugs = "Paclitaxel")
#'
#' # Load without z-score normalization
#' drug_raw <- loadTreatmentResponseNormalized(gCSI, drugs = "Paclitaxel", zscore = FALSE)
#'
#' # Load all drug data with normalization
#' all_drugs <- loadTreatmentResponseNormalized(gCSI)
#' }
loadTreatmentResponseNormalized <- function(dromaset_object,
                                          drugs = NULL,
                                          data_type = "all",
                                          tumor_type = "all",
                                          zscore = TRUE,
                                          return_data = TRUE) {

  # Validate input object
  if (!inherits(dromaset_object, "DromaSet")) {
    stop("Input must be a DromaSet object from DROMA.Set package")
  }

  # Load treatment response using DROMA.Set function
  treatment_data <- loadTreatmentResponse(
    dromaset_object,
    drugs = drugs,
    data_type = data_type,
    tumor_type = tumor_type,
    return_data = return_data
  )

  # Apply z-score normalization if requested
  if (zscore && !is.null(treatment_data) && is.matrix(treatment_data)) {
    # Apply z-score normalization
    treatment_data <- zscoreNormalize(treatment_data)

    # Add attribute to indicate normalization was applied
    attr(treatment_data, "zscore_normalized") <- TRUE
  }

  return(treatment_data)
}

#' Load multi-project molecular profiles with optional z-score normalization
#'
#' @description Wrapper function for loadMultiProjectMolecularProfiles that applies z-score normalization by default
#' @param multidromaset_object A MultiDromaSet object
#' @param molecular_type Type of molecular data to load (e.g., "mRNA", "cnv", "meth", "proteinrppa", "proteinms")
#' @param features Optional vector of specific features to load. If NULL, loads all features
#' @param overlap_only Logical, whether to use only overlapping samples (default: FALSE)
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @return List of matrices/data frames with molecular profile data from each project, optionally z-score normalized
#' @export
#' @examples
#' \dontrun{
#' # Load mRNA data across projects with z-score normalization (default)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")
#' mrna_data <- loadMultiProjectMolecularProfilesNormalized(multi_set, "mRNA", features = "ABCB1")
#'
#' # Load without z-score normalization
#' mrna_raw <- loadMultiProjectMolecularProfilesNormalized(multi_set, "mRNA",
#'                                                         features = "ABCB1", zscore = FALSE)
#' }
loadMultiProjectMolecularProfilesNormalized <- function(multidromaset_object,
                                                      molecular_type,
                                                      features = NULL,
                                                      overlap_only = FALSE,
                                                      data_type = "all",
                                                      tumor_type = "all",
                                                      zscore = TRUE) {

  # Validate input object
  if (!inherits(multidromaset_object, "MultiDromaSet")) {
    stop("Input must be a MultiDromaSet object from DROMA.Set package")
  }

  # Load molecular profiles using DROMA.Set function
  molecular_data_list <- loadMultiProjectMolecularProfiles(
    multidromaset_object,
    molecular_type = molecular_type,
    features = features,
    overlap_only = overlap_only,
    data_type = data_type,
    tumor_type = tumor_type
  )

  # Apply z-score normalization if requested and data is continuous
  if (zscore && !is.null(molecular_data_list)) {
    # Check if data is continuous (matrix format) and suitable for normalization
    continuous_types <- c("mRNA", "cnv", "meth", "proteinrppa", "proteinms")

    if (molecular_type %in% continuous_types) {
      # Apply z-score normalization to each project's data
      molecular_data_list <- lapply(molecular_data_list, function(project_data) {
        if (is.matrix(project_data)) {
          normalized_data <- zscoreNormalize(project_data)
          # Add attribute to indicate normalization was applied
          attr(normalized_data, "zscore_normalized") <- TRUE
          return(normalized_data)
        } else {
          return(project_data)
        }
      })
    } else if (zscore && !molecular_type %in% continuous_types) {
      warning(paste("Z-score normalization not applicable for molecular type:", molecular_type))
    }
  }

  return(molecular_data_list)
}

#' Load multi-project treatment response data with optional z-score normalization
#'
#' @description Wrapper function for loadMultiProjectTreatmentResponse that applies z-score normalization by default
#' @param multidromaset_object A MultiDromaSet object
#' @param drugs Optional vector of specific drugs to load. If NULL, loads all drugs
#' @param overlap_only Logical, whether to use only overlapping samples (default: FALSE)
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @return List of matrices with treatment response data from each project, optionally z-score normalized
#' @export
#' @examples
#' \dontrun{
#' # Load drug data across projects with z-score normalization (default)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")
#' drug_data <- loadMultiProjectTreatmentResponseNormalized(multi_set, drugs = "Paclitaxel")
#'
#' # Load without z-score normalization
#' drug_raw <- loadMultiProjectTreatmentResponseNormalized(multi_set, drugs = "Paclitaxel",
#'                                                        zscore = FALSE)
#' }
loadMultiProjectTreatmentResponseNormalized <- function(multidromaset_object,
                                                      drugs = NULL,
                                                      overlap_only = FALSE,
                                                      data_type = "all",
                                                      tumor_type = "all",
                                                      zscore = TRUE) {

  # Validate input object
  if (!inherits(multidromaset_object, "MultiDromaSet")) {
    stop("Input must be a MultiDromaSet object from DROMA.Set package")
  }

  # Load treatment response using DROMA.Set function
  treatment_data_list <- loadMultiProjectTreatmentResponse(
    multidromaset_object,
    drugs = drugs,
    overlap_only = overlap_only,
    data_type = data_type,
    tumor_type = tumor_type
  )

  # Apply z-score normalization if requested
  if (zscore && !is.null(treatment_data_list)) {
    # Apply z-score normalization to each project's data
    treatment_data_list <- lapply(treatment_data_list, function(project_data) {
      if (is.matrix(project_data)) {
        normalized_data <- zscoreNormalize(project_data)
        # Add attribute to indicate normalization was applied
        attr(normalized_data, "zscore_normalized") <- TRUE
        return(normalized_data)
      } else {
        return(project_data)
      }
    })
  }

  return(treatment_data_list)
}

#' Apply z-score normalization to existing data
#'
#' @description Convenience function to apply z-score normalization to already loaded data
#' @param data Matrix or data frame with features in rows and samples in columns
#' @param check_type Logical, whether to check if data appears to be continuous (default: TRUE)
#' @return Z-score normalized matrix
#' @export
#' @examples
#' \dontrun{
#' # Apply normalization to existing data
#' normalized_data <- applyZscoreNormalization(raw_data)
#'
#' # Skip type checking
#' normalized_data <- applyZscoreNormalization(raw_data, check_type = FALSE)
#' }
applyZscoreNormalization <- function(data, check_type = TRUE) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }

  # Convert to matrix if data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # Check if data appears to be continuous
  if (check_type) {
    # Simple heuristic: if more than 50% of values are unique, assume continuous
    unique_ratio <- length(unique(as.vector(data))) / length(as.vector(data))
    if (unique_ratio < 0.1) {
      warning("Data appears to be discrete. Z-score normalization may not be appropriate.")
    }
  }

  # Apply z-score normalization
  normalized_data <- zscoreNormalize(data)

  # Add attribute to indicate normalization was applied
  attr(normalized_data, "zscore_normalized") <- TRUE

  return(normalized_data)
}

#' Check if data has been z-score normalized
#'
#' @description Utility function to check if data has the z-score normalization attribute
#' @param data Matrix or data frame to check
#' @return Logical indicating whether data has been z-score normalized
#' @export
#' @examples
#' \dontrun{
#' # Check if data is normalized
#' is_normalized <- isZscoreNormalized(my_data)
#' }
isZscoreNormalized <- function(data) {
  return(!is.null(attr(data, "zscore_normalized")) && attr(data, "zscore_normalized"))
}
