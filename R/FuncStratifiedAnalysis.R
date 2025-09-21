# Stratified Analysis Functions ----

# Core principle: Analysis should be performed within biologically meaningful subgroups
# to uncover context-dependent associations that may be masked in pooled analysis

#' Create stratified DromaSet objects based on drug response
#'
#' @description Creates stratified datasets by dividing samples into subgroups based on
#' their response to a stratification drug. This enables context-dependent analysis
#' where biological associations may differ between resistant and sensitive populations.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Name of the drug used for stratification
#' @param strata_quantile Quantile threshold for creating strata (default: 0.33 for tertiles)
#' @param min_samples Minimum number of samples required in each stratum (default: 5)
#' @return List containing sensitive and resistant DromaSet objects
#' @export
createStratifiedDromaSets <- function(dromaset_object,
                                    stratification_drug,
                                    strata_quantile = 0.33,
                                    min_samples = 5) {

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object")
  }

  # Load stratification drug response data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    strat_data <- loadTreatmentResponseNormalized(dromaset_object,
                                                 drugs = stratification_drug,
                                                 data_type = "all",
                                                 tumor_type = "all",
                                                 return_data = TRUE)

    if (!is.matrix(strat_data) || !stratification_drug %in% rownames(strat_data)) {
      stop("Stratification drug not found in the dataset")
    }

    # Extract drug response vector
    drug_response <- as.numeric(strat_data[stratification_drug, ])
    names(drug_response) <- colnames(strat_data)
    drug_response <- drug_response[!is.na(drug_response)]

    if (length(drug_response) < (2 * min_samples)) {
      stop("Insufficient samples for stratification")
    }

    # Calculate quantiles
    lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
    upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

    # Create strata
    sensitive_samples <- names(drug_response[drug_response <= lower_quantile])
    resistant_samples <- names(drug_response[drug_response >= upper_quantile])

    # Validate sample sizes
    if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
      stop("Insufficient samples in one or more strata after stratification")
    }

    # Create stratified DromaSet objects
    sensitive_set <- createDromaSetFromDromaSet(
      dromaset_object,
      sample_subset = sensitive_samples,
      new_name = paste0(dromaset_object@name, "_Sensitive")
    )

    resistant_set <- createDromaSetFromDromaSet(
      dromaset_object,
      sample_subset = resistant_samples,
      new_name = paste0(dromaset_object@name, "_Resistant")
    )

    result <- list(
      sensitive = sensitive_set,
      resistant = resistant_set,
      stratification_info = list(
        drug = stratification_drug,
        quantile = strata_quantile,
        lower_threshold = lower_quantile,
        upper_threshold = upper_quantile,
        n_sensitive = length(sensitive_samples),
        n_resistant = length(resistant_samples)
      )
    )

  } else {
    # MultiDromaSet - handle each project separately
    strat_data_list <- loadMultiProjectTreatmentResponseNormalized(
      dromaset_object,
      drugs = stratification_drug,
      overlap_only = FALSE,
      data_type = "all",
      tumor_type = "all"
    )

    sensitive_sets <- list()
    resistant_sets <- list()
    stratification_info <- list()

    for (project_name in names(strat_data_list)) {
      strat_data <- strat_data_list[[project_name]]

      if (!is.matrix(strat_data) || !stratification_drug %in% rownames(strat_data)) {
        next
      }

      # Extract drug response vector
      drug_response <- as.numeric(strat_data[stratification_drug, ])
      names(drug_response) <- colnames(strat_data)
      drug_response <- drug_response[!is.na(drug_response)]

      if (length(drug_response) < (2 * min_samples)) {
        next
      }

      # Calculate quantiles
      lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
      upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

      # Create strata
      sensitive_samples <- names(drug_response[drug_response <= lower_quantile])
      resistant_samples <- names(drug_response[drug_response >= upper_quantile])

      # Skip if insufficient samples
      if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
        next
      }

      # Get original DromaSet
      original_set <- dromaset_object@DromaSets[[project_name]]

      # Create stratified DromaSet objects
      sensitive_sets[[project_name]] <- createDromaSetFromDromaSet(
        original_set,
        sample_subset = sensitive_samples,
        new_name = paste0(project_name, "_Sensitive")
      )

      resistant_sets[[project_name]] <- createDromaSetFromDromaSet(
        original_set,
        sample_subset = resistant_samples,
        new_name = paste0(project_name, "_Resistant")
      )

      stratification_info[[project_name]] <- list(
        drug = stratification_drug,
        quantile = strata_quantile,
        lower_threshold = lower_quantile,
        upper_threshold = upper_quantile,
        n_sensitive = length(sensitive_samples),
        n_resistant = length(resistant_samples)
      )
    }

    if (length(sensitive_sets) == 0) {
      stop("No projects had sufficient samples for stratification")
    }

    # Create MultiDromaSet objects
    sensitive_multi <- createMultiDromaSetFromList(sensitive_sets, "Sensitive_Group")
    resistant_multi <- createMultiDromaSetFromList(resistant_sets, "Resistant_Group")

    result <- list(
      sensitive = sensitive_multi,
      resistant = resistant_multi,
      stratification_info = stratification_info
    )
  }

  return(result)
}

#' Create DromaSet from existing DromaSet with sample subset
#'
#' @description Helper function to create a new DromaSet object from an existing one
#' with only a subset of samples. This maintains all metadata while filtering data.
#' @param original_set Original DromaSet object
#' @param sample_subset Vector of sample IDs to include
#' @param new_name Name for the new DromaSet
#' @return New DromaSet object with subset of samples
createDromaSetFromDromaSet <- function(original_set, sample_subset, new_name) {

  # Validate samples - check across all molecular profiles
  valid_samples <- character(0)
  for (mol_type in names(original_set@molecularProfiles)) {
    if (is.matrix(original_set@molecularProfiles[[mol_type]])) {
      valid_samples <- unique(c(valid_samples,
                              colnames(original_set@molecularProfiles[[mol_type]])))
    }
  }
  valid_samples <- intersect(sample_subset, valid_samples)

  if (length(valid_samples) == 0) {
    stop("No valid samples found for subset creation")
  }

  # Create new DromaSet object
  new_set <- new("DromaSet")

  # Copy basic attributes
  new_set@name <- new_name
  new_set@description <- paste0("Subset of ", original_set@name, " (", length(valid_samples), " samples)")

  # Copy treatment response data and filter samples
  new_set@treatmentResponse <- list()
  for (drug_type in names(original_set@treatmentResponse)) {
    if (is.matrix(original_set@treatmentResponse[[drug_type]]) &&
        ncol(original_set@treatmentResponse[[drug_type]]) > 0) {
      # Filter samples that exist in this matrix
      valid_samples_matrix <- intersect(valid_samples,
                                       colnames(original_set@treatmentResponse[[drug_type]]))
      if (length(valid_samples_matrix) > 0) {
        new_set@treatmentResponse[[drug_type]] <-
          original_set@treatmentResponse[[drug_type]][, valid_samples_matrix, drop = FALSE]
      }
    }
  }

  # Copy molecular profiles data and filter samples
  new_set@molecularProfiles <- list()
  for (mol_type in names(original_set@molecularProfiles)) {
    if (is.matrix(original_set@molecularProfiles[[mol_type]]) &&
        ncol(original_set@molecularProfiles[[mol_type]]) > 0) {
      # Filter samples that exist in this matrix
      valid_samples_matrix <- intersect(valid_samples,
                                       colnames(original_set@molecularProfiles[[mol_type]]))
      if (length(valid_samples_matrix) > 0) {
        new_set@molecularProfiles[[mol_type]] <-
          original_set@molecularProfiles[[mol_type]][, valid_samples_matrix, drop = FALSE]
      }
    }
  }

  # Copy and filter sample metadata
  if (nrow(original_set@sampleMetadata) > 0) {
    new_set@sampleMetadata <- original_set@sampleMetadata[
      original_set@sampleMetadata$sample_id %in% valid_samples,
    ]
  } else {
    new_set@sampleMetadata <- data.frame(sample_id = valid_samples)
  }

  # Copy metadata unchanged
  new_set@treatmentMetadata <- original_set@treatmentMetadata
  new_set@datasetType <- original_set@datasetType
  new_set@db_info <- original_set@db_info

  return(new_set)
}

#' Create MultiDromaSet from list of DromaSet objects
#'
#' @description Helper function to create a MultiDromaSet from a list of DromaSet objects
#' @param dromaset_list List of DromaSet objects
#' @param group_name Name for the group
#' @return MultiDromaSet object
createMultiDromaSetFromList <- function(dromaset_list, group_name) {

  # Create MultiDromaSet object
  multi_set <- new("MultiDromaSet")

  # Set basic attributes
  multi_set@name <- names(dromaset_list)
  multi_set@description <- paste0("Multi-project analysis for ", group_name)

  # Store DromaSet objects
  multi_set@DromaSets <- dromaset_list

  # Merge sample metadata
  all_sample_metadata <- do.call(rbind, lapply(dromaset_list, function(set) {
    if (nrow(set@sampleMetadata) > 0) {
      return(set@sampleMetadata)
    } else {
      return(data.frame(sample_id = colnames(set@molecularProfiles$mRNA)))
    }
  }))
  multi_set@sampleMetadata <- all_sample_metadata

  # Copy treatment metadata from first set
  multi_set@treatmentMetadata <- dromaset_list[[1]]@treatmentMetadata

  # Set dataset types
  multi_set@datasetType <- unique(sapply(dromaset_list, function(set) set@datasetType))

  # Copy db info from first set
  multi_set@db_info <- dromaset_list[[1]]@db_info

  return(multi_set)
}

#' Analyze stratified drug-omic associations
#'
#' @description Performs drug-omic association analysis within strata defined by
#' response to a stratification drug. This approach reveals context-dependent
#' associations that may differ between biological subgroups.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Drug used to define strata (e.g., "cisplatin")
#' @param strata_quantile Quantile threshold for stratification (default: 0.33)
#' @param select_omics_type Type of omics data to analyze
#' @param select_omics Name of the specific omics feature
#' @param select_drugs Name of the target drug for association analysis
#' @param min_samples Minimum samples per stratum (default: 5)
#' @param ... Additional parameters passed to analyzeDrugOmicPair
#' @return List containing stratified analysis results and between-stratum comparisons
#' @export
#' @examples
#' \dontrun{
#' # Stratified analysis of ERCC1 expression and bortezomib response
#' # based on cisplatin sensitivity
#' result <- analyzeStratifiedDrugOmic(
#'   dromaset_object = multi_set,
#'   stratification_drug = "cisplatin",
#'   strata_quantile = 0.33,
#'   select_omics_type = "mRNA",
#'   select_omics = "ERCC1",
#'   select_drugs = "bortezomib",
#'   data_type = "CellLine",
#'   tumor_type = "all"
#' )
#' }
analyzeStratifiedDrugOmic <- function(dromaset_object,
                                    stratification_drug,
                                    strata_quantile = 0.33,
                                    select_omics_type,
                                    select_omics,
                                    select_drugs,
                                    min_samples = 5,
                                    ...) {

  # Validate inputs
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object")
  }

  if (missing(select_omics_type) || missing(select_omics) || missing(select_drugs)) {
    stop("select_omics_type, select_omics, and select_drugs are required")
  }

  # Create stratified datasets
  stratified_sets <- createStratifiedDromaSets(
    dromaset_object = dromaset_object,
    stratification_drug = stratification_drug,
    strata_quantile = strata_quantile,
    min_samples = min_samples
  )

  # Extract additional parameters
  extra_params <- list(...)

  # Perform analysis in sensitive group
  sensitive_result <- do.call(analyzeDrugOmicPair, c(
    list(
      dromaset_object = stratified_sets$sensitive,
      select_omics_type = select_omics_type,
      select_omics = select_omics,
      select_drugs = select_drugs
    ),
    extra_params
  ))

  # Perform analysis in resistant group
  resistant_result <- do.call(analyzeDrugOmicPair, c(
    list(
      dromaset_object = stratified_sets$resistant,
      select_omics_type = select_omics_type,
      select_omics = select_omics,
      select_drugs = select_drugs
    ),
    extra_params
  ))

  # Compare results between strata
  comparison_result <- compareStratifiedResults(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    stratification_info = stratified_sets$stratification_info,
    select_omics_type = select_omics_type,
    select_omics = select_omics,
    select_drugs = select_drugs
  )

  # Compile final results
  result <- list(
    sensitive = sensitive_result,
    resistant = resistant_result,
    comparison = comparison_result,
    stratification_info = stratified_sets$stratification_info,
    analysis_parameters = list(
      stratification_drug = stratification_drug,
      strata_quantile = strata_quantile,
      select_omics_type = select_omics_type,
      select_omics = select_omics,
      select_drugs = select_drugs,
      min_samples = min_samples,
      extra_params = extra_params
    )
  )

  class(result) <- "StratifiedDrugOmicResult"

  return(result)
}

#' Compare results between strata
#'
#' @description Compares association results between sensitive and resistant groups
#' @param sensitive_result Analysis result from sensitive group
#' @param resistant_result Analysis result from resistant group
#' @param stratification_info Information about stratification
#' @param select_omics_type Type of omics data analyzed
#' @param select_omics Name of the omics feature
#' @param select_drugs Name of the drug analyzed
#' @return Comparison results with statistical tests
compareStratifiedResults <- function(sensitive_result, resistant_result, stratification_info,
                                   select_omics_type, select_omics, select_drugs) {

  comparison <- list()

  # Extract correlation coefficients if available
  if (!is.null(sensitive_result$meta) && !is.null(resistant_result$meta)) {
    if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
      # Continuous-continuous analysis
      sens_cor <- sensitive_result$meta$TE.random
      sens_se <- sensitive_result$meta$seTE.random
      res_cor <- resistant_result$meta$TE.random
      res_se <- resistant_result$meta$seTE.random

      # Test for difference in correlations
      z_diff <- (sens_cor - res_cor) / sqrt(sens_se^2 + res_se^2)
      p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)

      comparison$correlation_comparison <- data.frame(
        Group = c("Sensitive", "Resistant", "Difference"),
        Correlation = c(sens_cor, res_cor, sens_cor - res_cor),
        SE = c(sens_se, res_se, sqrt(sens_se^2 + res_se^2)),
        P_value = c(NA, NA, p_diff)
      )
    }
  }

  # Create comparison plot
  comparison$comparison_plot <- createStratifiedComparisonPlot(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    select_omics = select_omics,
    select_drugs = select_drugs
  )

  return(comparison)
}

#' Create stratified comparison plot
#'
#' @description Creates a visualization comparing results between strata
#' @param sensitive_result Results from sensitive group
#' @param resistant_result Results from resistant group
#' @param select_omics Name of omics feature
#' @param select_drugs Name of drug
#' @return ggplot object
createStratifiedComparisonPlot <- function(sensitive_result, resistant_result,
                                          select_omics, select_drugs) {

  # This function would create a comparison plot
  # Implementation would depend on specific visualization requirements

  return(NULL)
}