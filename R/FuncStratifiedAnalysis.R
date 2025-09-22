# Stratified Analysis Functions for ACC Data ----

# Core principle: Stratified analysis should follow the same pairing logic as
# the original analysis functions, filtering samples within each pairing operation
# rather than creating pre-filtered datasets
#
# NOTE: This implementation is designed for ACC (Activity Area Above the Curve) data
# where higher values indicate greater drug sensitivity

#' Load stratified drug response data
#'
#' @description Loads drug response data and stratifies samples into sensitive/resistant groups.
#' For ACC (Activity Area Above the Curve) data where higher values indicate greater sensitivity.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Name of the drug used for stratification
#' @param strata_quantile Quantile threshold for creating strata (default: 0.33)
#' @param min_samples Minimum samples required in each stratum (default: 5)
#' @return List containing stratification information for each project
#' @export
getStratificationInfo <- function(dromaset_object,
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
      stop("Stratification drug '", stratification_drug, "' not found in dataset")
    }

    # Extract drug response vector
    drug_response <- as.numeric(strat_data[stratification_drug, ])
    names(drug_response) <- colnames(strat_data)
    drug_response <- drug_response[!is.na(drug_response)]

    if (length(drug_response) < (2 * min_samples)) {
      stop("Insufficient samples (", length(drug_response), ") for stratification. Minimum required: ", 2 * min_samples)
    }

    # Calculate quantiles
    lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
    upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

    # Create strata - ACC data: higher values = more sensitive
    sensitive_samples <- names(drug_response[drug_response >= upper_quantile])
    resistant_samples <- names(drug_response[drug_response <= lower_quantile])

    # Validate sample sizes
    if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
      stop("Insufficient samples in strata. Sensitive: ", length(sensitive_samples),
           ", Resistant: ", length(resistant_samples), ". Minimum required: ", min_samples)
    }

    result <- list()
    result[[dromaset_object@name]] <- list(
      sensitive = sensitive_samples,
      resistant = resistant_samples,
      thresholds = c(lower = lower_quantile, upper = upper_quantile),
      drug_response = drug_response
    )

    attr(result, "type") <- "single"

  } else {
    # MultiDromaSet - handle each project separately
    strat_data_list <- loadMultiProjectTreatmentResponseNormalized(
      dromaset_object,
      drugs = stratification_drug,
      overlap_only = FALSE,
      data_type = "all",
      tumor_type = "all"
    )

    # Check if any projects have the stratification drug
    projects_with_drug <- names(strat_data_list)[
      sapply(strat_data_list, function(mat) is.matrix(mat) && stratification_drug %in% rownames(mat))
    ]

    if (length(projects_with_drug) == 0) {
      stop("Stratification drug '", stratification_drug, "' not found in any project")
    }

    # Warn about projects without the drug
    all_projects <- names(dromaset_object@DromaSets)
    if (length(all_projects) > length(projects_with_drug)) {
      projects_without <- setdiff(all_projects, projects_with_drug)
      warning("Stratification drug not found in projects: ",
              paste(projects_without, collapse = ", "))
    }

    stratification_info <- list()

    for (project_name in names(strat_data_list)) {
      strat_data <- strat_data_list[[project_name]]

      if (!is.matrix(strat_data) || !stratification_drug %in% rownames(strat_data)) {
        # Skip projects without the drug
        next
      }

      # Extract drug response vector
      drug_response <- as.numeric(strat_data[stratification_drug, ])
      names(drug_response) <- colnames(strat_data)
      drug_response <- drug_response[!is.na(drug_response)]

      if (length(drug_response) < (2 * min_samples)) {
        warning("Project ", project_name, ": Insufficient samples (", length(drug_response),
                ") for stratification. Minimum required: ", 2 * min_samples)
        next
      }

      # Calculate quantiles
      lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
      upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

      # Create strata - ACC data: higher values = more sensitive
      sensitive_samples <- names(drug_response[drug_response >= upper_quantile])
      resistant_samples <- names(drug_response[drug_response <= lower_quantile])

      # Skip if insufficient samples
      if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
        warning("Project ", project_name, ": Insufficient samples in strata. Sensitive: ",
                length(sensitive_samples), ", Resistant: ", length(resistant_samples),
                ". Minimum required: ", min_samples)
        next
      }

      stratification_info[[project_name]] <- list(
        sensitive = sensitive_samples,
        resistant = resistant_samples,
        thresholds = c(lower = lower_quantile, upper = upper_quantile),
        drug_response = drug_response
      )
    }

    if (length(stratification_info) == 0) {
      stop("No projects had sufficient samples for stratification")
    }

    result <- stratification_info
    attr(result, "type") <- "multi"
  }

  return(result)
}

#' Load data for stratified analysis project by project
#' @description This function is deprecated and no longer used in the stratified analysis workflow.
#' Data loading is now handled directly within analyzeStratifiedDrugOmic using the standard
#' loading functions from FuncPairDrugOmicPair.R
#' @export
loadStratifiedDataByProject <- function(...) {
  warning("loadStratifiedDataByProject is deprecated and no longer used.",
          " Data loading is handled within analyzeStratifiedDrugOmic.")
  return(NULL)
}

#' Stratified version of pairContinuousFeatures
#'
#' @description Creates paired datasets within specified strata (sensitive/resistant)
#' @param dataset1 First dataset list with named vectors
#' @param dataset2 Second dataset list with named vectors
#' @param stratification_info List containing sample stratification information
#' @param stratum Which stratum to analyze ("sensitive" or "resistant")
#' @param merged Logical, if TRUE, creates an additional merged dataset
#' @return List of paired data for the specified stratum
#' @export
pairContinuousFeaturesStratified <- function(dataset1, dataset2, stratification_info,
                                           stratum = c("sensitive", "resistant"),
                                           merged = FALSE) {

  stratum <- match.arg(stratum)

  pair_list2 <- lapply(seq_along(dataset1), function(x) {
    feat1_sel <- dataset1[[x]]
    study_name <- names(dataset1)[x]

    # Get stratum-specific samples for this study
    if (study_name %in% names(stratification_info)) {
      stratum_samples <- stratification_info[[study_name]][[stratum]]
    } else {
      # If study not in stratification info, use all samples
      stratum_samples <- names(feat1_sel)
    }

    # Filter feature to only include stratum samples
    feat1_sel_strat <- feat1_sel[names(feat1_sel) %in% stratum_samples]

    pair_list <- lapply(1:length(dataset2), function(y) {
      feat2_sel <- dataset2[[y]]
      study_name2 <- names(dataset2)[y]

      # Only pair within the same study
      if (study_name != study_name2) {
        return(NULL)
      }

      # Filter feature to only include stratum samples
      feat2_sel_strat <- feat2_sel[names(feat2_sel) %in% stratum_samples]

      feat1_sel_strat <- na.omit(feat1_sel_strat)
      feat2_sel_strat <- na.omit(feat2_sel_strat)

      if(length(feat1_sel_strat) == 0 | length(feat2_sel_strat) == 0) {
        return(NULL)
      }

      # Find intersection within the stratum
      intersected_cells <- intersect(names(feat1_sel_strat), names(feat2_sel_strat))
      if(length(intersected_cells) < 3) {
        return(NULL)
      }

      feat1_sel_strat <- feat1_sel_strat[match(intersected_cells, names(feat1_sel_strat))]
      feat2_sel_strat <- feat2_sel_strat[match(intersected_cells, names(feat2_sel_strat))]

      list(feature1 = feat1_sel_strat,
           feature2 = feat2_sel_strat)
    })

    names(pair_list) <- paste0(names(dataset1)[x], "_",
                               names(dataset2))
    pair_list
  })

  pair_list2 <- unlist(pair_list2, recursive = FALSE)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]

  if(length(pair_list2) < 1) {
    warning(paste("No valid pairs found for", stratum, "stratum"))
    return(NULL)
  }

  # If merged is TRUE, create a merged dataset
  if(merged & length(pair_list2) > 1) {
    combined_list <- list(
      # Combine all feature1 vectors
      unlist(lapply(pair_list2, function(x) x$feature1)),

      # Combine all feature2 vectors
      unlist(lapply(pair_list2, function(x) x$feature2))
    )

    pair_list2[["merged_dataset"]] <- list(
      "feature1" = combined_list[[1]],
      "feature2" = combined_list[[2]]
    )
  }

  pair_list2
}

#' Stratified version of pairDiscreteFeatures
#'
#' @description Creates paired datasets of discrete and continuous features within strata
#' @param discrete_dataset List with discrete feature data
#' @param continuous_dataset List with continuous feature data
#' @param stratification_info List containing sample stratification information
#' @param stratum Which stratum to analyze ("sensitive" or "resistant")
#' @param merged Logical, if TRUE, creates an additional merged dataset
#' @return List of paired data for the specified stratum
#' @export
pairDiscreteFeaturesStratified <- function(discrete_dataset, continuous_dataset,
                                          stratification_info,
                                          stratum = c("sensitive", "resistant"),
                                          merged = FALSE) {

  stratum <- match.arg(stratum)

  pair_list2 <- lapply(1:length(discrete_dataset), function(x) {
    discrete_sel <- discrete_dataset[[x]]
    study_name <- names(discrete_dataset)[x]

    # Get stratum-specific samples for this study
    if (study_name %in% names(stratification_info)) {
      stratum_samples <- stratification_info[[study_name]][[stratum]]
    } else {
      # If study not in stratification info, use all samples
      stratum_samples <- NULL
    }

    pair_list <- lapply(1:length(continuous_dataset), function(y) {
      continuous_sel <- continuous_dataset[[y]]
      study_name2 <- names(continuous_dataset)[y]

      # Only pair within the same study
      if (study_name != study_name2) {
        return(NULL)
      }

      # Filter continuous data to stratum if specified
      if (!is.null(stratum_samples)) {
        continuous_sel <- continuous_sel[names(continuous_sel) %in% stratum_samples]
      }

      yes_values <- na.omit(continuous_sel[names(continuous_sel) %in% discrete_sel])
      no_values <- na.omit(continuous_sel[!names(continuous_sel) %in% discrete_sel])

      if(length(yes_values) < 3 | length(no_values) < 3) {
        return(NULL)
      }

      list(yes = yes_values,
           no = no_values)
    })

    names(pair_list) <- paste0(names(discrete_dataset)[x], "_",
                               names(continuous_dataset))
    pair_list
  })

  pair_list2 <- unlist(pair_list2, recursive = FALSE)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]

  if(length(pair_list2) < 1) {
    warning(paste("No valid pairs found for", stratum, "stratum"))
    return(NULL)
  }

  # If merged is TRUE, create a merged dataset
  if(merged & length(pair_list2) > 1) {
    combined_list <- list(
      # Combine all yes vectors
      unlist(lapply(pair_list2, function(x) x$yes)),

      # Combine all no vectors
      unlist(lapply(pair_list2, function(x) x$no))
    )

    pair_list2[["merged_dataset"]] <- list(
      "yes" = combined_list[[1]],
      "no" = combined_list[[2]]
    )
  }

  return(pair_list2)
}

#' Analyze stratified drug-omic associations
#'
#' @description Performs drug-omic association analysis within strata defined by
#' response to a stratification drug. This approach follows the pairing principle
#' by filtering samples within each pairing operation rather than creating pre-filtered datasets.
#' Designed for ACC (Activity Area Above the Curve) data where higher values indicate greater sensitivity.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Drug used to define strata (e.g., "cisplatin")
#' @param strata_quantile Quantile threshold for stratification (default: 0.33)
#' @param select_omics_type Type of omics data to analyze
#' @param select_omics Name of the specific omics feature
#' @param select_drugs Name of the target drug for association analysis
#' @param min_samples Minimum samples per stratum (default: 5)
#' @param ... Additional parameters passed to the analysis functions
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

  # Extract additional parameters
  extra_params <- list(...)
  data_type <- extra_params$data_type %||% "all"
  tumor_type <- extra_params$tumor_type %||% "all"
  overlap_only <- extra_params$overlap_only %||% FALSE
  merged_enabled <- extra_params$merged_enabled %||% TRUE
  meta_enabled <- extra_params$meta_enabled %||% TRUE

  # Step 1: Get stratification information (sensitive/resistant sample names for each project)
  cat("Getting stratification info...\n")
  stratification_info <- getStratificationInfo(
    dromaset_object = dromaset_object,
    stratification_drug = stratification_drug,
    strata_quantile = strata_quantile,
    min_samples = min_samples
  )

  # Step 2: Load original data (same as analyzeDrugOmicPair)
  cat("Loading drug data...\n")
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    drug_data <- loadTreatmentResponseNormalized(dromaset_object,
                                              drugs = select_drugs,
                                              data_type = data_type,
                                              tumor_type = tumor_type,
                                              return_data = TRUE)

    # Convert matrix to list format
    if (is.matrix(drug_data) && select_drugs %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[select_drugs, ])
      names(drug_vector) <- colnames(drug_data)
      drug_data_list <- list()
      drug_data_list[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    } else {
      drug_data_list <- list()
    }
  } else {
    # MultiDromaSet
    drug_data_list <- loadMultiProjectTreatmentResponseNormalized(
      dromaset_object,
      drugs = select_drugs,
      overlap_only = overlap_only,
      data_type = data_type,
      tumor_type = tumor_type
    )

    # Extract specific drug from each project
    drug_data_list <- lapply(drug_data_list, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })
    drug_data_list <- drug_data_list[!sapply(drug_data_list, is.null)]
  }

  # Load omics data
  cat("Loading omics data...\n")
  if (select_omics_type %in% c("drug", "drug_raw")) {
    # Another drug
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadTreatmentResponseNormalized(dromaset_object,
                                                   drugs = select_omics,
                                                   data_type = data_type,
                                                   tumor_type = tumor_type,
                                                   return_data = TRUE)

      if (is.matrix(omics_data) && select_omics %in% rownames(omics_data)) {
        omics_vector <- as.numeric(omics_data[select_omics, ])
        names(omics_vector) <- colnames(omics_data)
        omics_data_list <- list()
        omics_data_list[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
      } else {
        omics_data_list <- list()
      }
    } else {
      omics_data_list <- loadMultiProjectTreatmentResponseNormalized(
        dromaset_object,
        drugs = select_omics,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type
      )

      omics_data_list <- lapply(omics_data_list, function(drug_matrix) {
        if (is.matrix(drug_matrix) && select_omics %in% rownames(drug_matrix)) {
          drug_vector <- as.numeric(drug_matrix[select_omics, ])
          names(drug_vector) <- colnames(drug_matrix)
          return(drug_vector[!is.na(drug_vector)])
        }
        return(NULL)
      })
      omics_data_list <- omics_data_list[!sapply(omics_data_list, is.null)]
    }
  } else {
    # Molecular profiles
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadMolecularProfilesNormalized(dromaset_object,
                                                  molecular_type = select_omics_type,
                                                  features = select_omics,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  return_data = TRUE)

      if (is.matrix(omics_data) && select_omics %in% rownames(omics_data)) {
        omics_data_list <- list()
        if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
          # Continuous
          omics_vector <- as.numeric(omics_data[select_omics, ])
          names(omics_vector) <- colnames(omics_data)
          omics_data_list[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
        } else {
          # Discrete
          gene_row <- omics_data[select_omics, ]
          present_samples <- names(gene_row)[gene_row != 0]
          omics_data_list[[dromaset_object@name]] <- present_samples
        }
      } else {
        omics_data_list <- list()
      }
    } else {
      omics_data_list <- loadMultiProjectMolecularProfilesNormalized(
        dromaset_object,
        molecular_type = select_omics_type,
        features = select_omics,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type
      )

      omics_data_list <- lapply(omics_data_list, function(omics_matrix) {
        if (is.matrix(omics_matrix) && select_omics %in% rownames(omics_matrix)) {
          if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
            # Continuous
            omics_vector <- as.numeric(omics_matrix[select_omics, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          } else {
            # Discrete
            gene_row <- omics_matrix[select_omics, ]
            present_samples <- names(gene_row)[gene_row != 0]
            return(present_samples)
          }
        }
        return(NULL)
      })
      omics_data_list <- omics_data_list[!sapply(omics_data_list, is.null)]
    }
  }

  # Check if we have data
  if (length(drug_data_list) == 0 || length(omics_data_list) == 0) {
    stop("No data found for the specified drug-omics pair")
  }

  # Step 3: Filter data for each stratum and analyze
  cat("Analyzing sensitive group...\n")
  sensitive_loaded_data <- list(drugs = drug_data_list, omics = omics_data_list)
  sensitive_data <- filterStratifiedData(sensitive_loaded_data, stratification_info, "sensitive")
  sensitive_result <- analyzeStratifiedData(
    omics_data = sensitive_data$omics,
    drug_data = sensitive_data$drugs,
    select_omics_type = select_omics_type,
    merged_enabled = merged_enabled,
    meta_enabled = meta_enabled
  )

  cat("Analyzing resistant group...\n")
  resistant_loaded_data <- list(drugs = drug_data_list, omics = omics_data_list)
  resistant_data <- filterStratifiedData(resistant_loaded_data, stratification_info, "resistant")
  resistant_result <- analyzeStratifiedData(
    omics_data = resistant_data$omics,
    drug_data = resistant_data$drugs,
    select_omics_type = select_omics_type,
    merged_enabled = merged_enabled,
    meta_enabled = meta_enabled
  )

  # Compare results between strata
  comparison_result <- compareStratifiedResults(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    stratification_info = stratification_info,
    select_omics_type = select_omics_type,
    select_omics = select_omics,
    select_drugs = select_drugs
  )

  # Compile final results
  result <- list(
    sensitive = sensitive_result,
    resistant = resistant_result,
    comparison = comparison_result,
    stratification_info = stratification_info,
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

  cat("Stratified analysis complete!\n")
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
    } else {
      # Discrete analysis - compare effect sizes
      sens_es <- sensitive_result$meta$TE.random
      sens_se <- sensitive_result$meta$seTE.random
      res_es <- resistant_result$meta$TE.random
      res_se <- resistant_result$meta$seTE.random

      # Test for difference in effect sizes
      z_diff <- (sens_es - res_es) / sqrt(sens_se^2 + res_se^2)
      p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)

      comparison$effect_comparison <- data.frame(
        Group = c("Sensitive", "Resistant", "Difference"),
        Effect_Size = c(sens_es, res_es, sens_es - res_es),
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
    select_drugs = select_drugs,
    select_omics_type = select_omics_type
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
#' @param select_omics_type Type of omics data
#' @return ggplot object
createStratifiedComparisonPlot <- function(sensitive_result, resistant_result,
                                          select_omics, select_drugs,
                                          select_omics_type) {

  # This function creates a comparison plot
  # For now, return NULL - can be implemented based on specific needs

  return(NULL)
}

#' Filter data for stratified analysis
#'
#' @description Filters loaded omics and drug data to keep only samples from specified stratum
#' @param dromaset_object DromaSet or MultiDromaSet object
#' @param stratification_info Stratification information from getStratificationInfo
#' @param stratum Which stratum ("sensitive" or "resistant")
#' @return List with filtered omics and drug data
#' @export
filterStratifiedData <- function(loaded_data, stratification_info, stratum) {

  # Validate stratum
  stratum <- match.arg(stratum, c("sensitive", "resistant"))

  # Filter data to keep only samples from specified stratum
  filtered_omics <- list()
  filtered_drugs <- list()

  for (project_name in names(loaded_data$omics)) {
    if (project_name %in% names(stratification_info)) {
      # Get stratum samples for this project
      stratum_samples <- stratification_info[[project_name]][[stratum]]

      # Filter omics data
      omics_data <- loaded_data$omics[[project_name]]
      if (is.numeric(omics_data)) {
        # Continuous omics data
        filtered_omics_data <- omics_data[names(omics_data) %in% stratum_samples]
        if (length(filtered_omics_data) >= 3) {
          filtered_omics[[project_name]] <- filtered_omics_data
        }
      } else if (is.character(omics_data)) {
        # Discrete omics data - keep samples that are both in omics and in stratum
        filtered_omics_data <- omics_data[omics_data %in% stratum_samples]
        if (length(filtered_omics_data) >= 3) {
          filtered_omics[[project_name]] <- filtered_omics_data
        }
      }

      # Filter drug data
      drug_data <- loaded_data$drugs[[project_name]]
      if (!is.null(drug_data)) {
        filtered_drug_data <- drug_data[names(drug_data) %in% stratum_samples]
        if (length(filtered_drug_data) >= 3) {
          filtered_drugs[[project_name]] <- filtered_drug_data
        }
      }
    }
  }

  return(list(omics = filtered_omics, drugs = filtered_drugs))
}

#' Backward compatibility wrapper for loadAndFilterStratifiedData
#' @description This function is deprecated. Use filterStratifiedData instead.
#' @export
loadAndFilterStratifiedData <- function(dromaset_object, select_omics_type, select_omics,
                                       select_drugs, stratification_info, stratum,
                                       data_type = "all", tumor_type = "all",
                                       overlap_only = FALSE) {

  # Load original data using standard functions
  loaded_data <- list()

  # Load drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    drug_data <- loadTreatmentResponseNormalized(dromaset_object,
                                              drugs = select_drugs,
                                              data_type = data_type,
                                              tumor_type = tumor_type,
                                              return_data = TRUE)

    if (is.matrix(drug_data) && select_drugs %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[select_drugs, ])
      names(drug_vector) <- colnames(drug_data)
      loaded_data$drugs <- list()
      loaded_data$drugs[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    }
  } else {
    # MultiDromaSet
    loaded_data$drugs <- loadMultiProjectTreatmentResponseNormalized(
      dromaset_object,
      drugs = select_drugs,
      overlap_only = overlap_only,
      data_type = data_type,
      tumor_type = tumor_type
    )

    # Extract specific drug
    loaded_data$drugs <- lapply(loaded_data$drugs, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })
    loaded_data$drugs <- loaded_data$drugs[!sapply(loaded_data$drugs, is.null)]
  }

  # Load omics data
  if (select_omics_type %in% c("drug", "drug_raw")) {
    # Another drug
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadTreatmentResponseNormalized(dromaset_object,
                                                   drugs = select_omics,
                                                   data_type = data_type,
                                                   tumor_type = tumor_type,
                                                   return_data = TRUE)

      if (is.matrix(omics_data) && select_omics %in% rownames(omics_data)) {
        omics_vector <- as.numeric(omics_data[select_omics, ])
        names(omics_vector) <- colnames(omics_data)
        loaded_data$omics <- list()
        loaded_data$omics[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
      }
    } else {
      loaded_data$omics <- loadMultiProjectTreatmentResponseNormalized(
        dromaset_object,
        drugs = select_omics,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type
      )

      loaded_data$omics <- lapply(loaded_data$omics, function(omics_matrix) {
        if (is.matrix(omics_matrix) && select_omics %in% rownames(omics_matrix)) {
          omics_vector <- as.numeric(omics_matrix[select_omics, ])
          names(omics_vector) <- colnames(omics_matrix)
          return(omics_vector[!is.na(omics_vector)])
        }
        return(NULL)
      })
      loaded_data$omics <- loaded_data$omics[!sapply(loaded_data$omics, is.null)]
    }
  } else {
    # Molecular profiles
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadMolecularProfilesNormalized(dromaset_object,
                                                  molecular_type = select_omics_type,
                                                  features = select_omics,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  return_data = TRUE)

      if (is.matrix(omics_data) && select_omics %in% rownames(omics_data)) {
        if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
          # Continuous
          omics_vector <- as.numeric(omics_data[select_omics, ])
          names(omics_vector) <- colnames(omics_data)
          loaded_data$omics <- list()
          loaded_data$omics[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
        } else {
          # Discrete
          gene_row <- omics_data[select_omics, ]
          present_samples <- names(gene_row)[gene_row != 0]
          loaded_data$omics <- list()
          loaded_data$omics[[dromaset_object@name]] <- present_samples
        }
      }
    } else {
      loaded_data$omics <- loadMultiProjectMolecularProfilesNormalized(
        dromaset_object,
        molecular_type = select_omics_type,
        features = select_omics,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type
      )

      loaded_data$omics <- lapply(loaded_data$omics, function(omics_matrix) {
        if (is.matrix(omics_matrix) && select_omics %in% rownames(omics_matrix)) {
          if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
            # Continuous
            omics_vector <- as.numeric(omics_matrix[select_omics, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          } else {
            # Discrete
            gene_row <- omics_matrix[select_omics, ]
            present_samples <- names(gene_row)[gene_row != 0]
            return(present_samples)
          }
        }
        return(NULL)
      })
      loaded_data$omics <- loaded_data$omics[!sapply(loaded_data$omics, is.null)]
    }
  }

  # Now filter the loaded data
  return(filterStratifiedData(loaded_data, stratification_info, stratum))
}

#' Analyze stratified data
#'
#' @description Performs drug-omics analysis on pre-filtered stratified data
#' @param omics_data Filtered omics data
#' @param drug_data Filtered drug data  
#' @param select_omics_type Type of omics data
#' @param merged_enabled Whether to create merged dataset
#' @param meta_enabled Whether to perform meta-analysis
#' @return Analysis results
#' @export
analyzeStratifiedData <- function(omics_data, drug_data, select_omics_type, 
                                 merged_enabled = TRUE, meta_enabled = TRUE) {
  
  # Check if we have data
  if (length(drug_data) == 0 || length(omics_data) == 0) {
    return(list())
  }
  
  # Perform analysis based on omics type
  if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
    # Continuous omics data
    pairs <- pairContinuousFeatures(omics_data, drug_data, merged = merged_enabled)
    
    if (!is.null(pairs)) {
      result <- list(
        plot = plotMultipleCorrelations(
          lapply(pairs, function(p) {
            if (!is.null(p$feature1) && !is.null(p$feature2)) {
              list(feature1 = p$feature1, feature2 = p$feature2)
            } else {
              NULL
            }
          })[!sapply(lapply(pairs, function(p) {
            if (!is.null(p$feature1) && !is.null(p$feature2)) {
              list(feature1 = p$feature1, feature2 = p$feature2)
            } else {
              NULL
            }
          }), is.null)],
          x_label = "Omic Feature",
          y_label = "Drug Response"
        ),
        data = pairs
      )
      
      # Perform meta-analysis if multiple studies
      if (length(pairs) > 1 && meta_enabled) {
        # Remove merged dataset for meta-analysis
        if ("merged_dataset" %in% names(pairs)) {
          meta_pairs <- pairs[names(pairs) != "merged_dataset"]
        } else {
          meta_pairs <- pairs
        }
        
        if (length(meta_pairs) > 0) {
          result$meta <- metaCalcConCon(meta_pairs)
        }
      }
    } else {
      result <- list()
    }
  } else {
    # Discrete omics data
    pairs <- pairDiscreteFeatures(omics_data, drug_data, merged = merged_enabled)
    
    if (!is.null(pairs)) {
      result <- list(
        plot = plotMultipleGroupComparisons(pairs, y_label = "Drug Response"),
        data = pairs
      )
      
      # Perform meta-analysis if multiple studies
      if (length(pairs) > 1 && meta_enabled) {
        # Remove merged dataset for meta-analysis
        if ("merged_dataset" %in% names(pairs)) {
          meta_pairs <- pairs[names(pairs) != "merged_dataset"]
        } else {
          meta_pairs <- pairs
        }
        
        if (length(meta_pairs) > 0) {
          result$meta <- metaCalcConDis(meta_pairs)
        }
      }
    } else {
      result <- list()
    }
  }
  
  return(result)
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x