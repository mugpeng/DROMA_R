# Meta-analysis Functions ----
# Note: Meta-analysis functions are now centralized in FuncMetaAnalysis module

# Utility Functions ----
# Note: Utility functions are now centralized in FuncDataPairing module

# Pairing Functions ----
# Note: Pairing functions are now centralized in FuncDataPairing module



# Plot Functions ----
# Note: Plot functions are now centralized in FuncMetaAnalysis module

# Helper Functions for Data Loading ----

#' Load feature data from DromaSet or MultiDromaSet
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature_type Type of feature to load
#' @param feature_name Name of feature to load
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples
#' @param is_continuous Whether the feature is continuous
#' @return List of feature data by dataset
#' @keywords internal
loadFeatureData <- function(dromaset_object, feature_type, feature_name,
                           data_type = "all", tumor_type = "all",
                           overlap_only = FALSE, is_continuous = TRUE) {

  # Validate inputs
  if (is.null(dromaset_object)) {
    stop("dromaset_object cannot be NULL")
  }
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("dromaset_object must be a DromaSet or MultiDromaSet object")
  }
  if (is.null(feature_type) || !is.character(feature_type) || length(feature_type) == 0) {
    stop("feature_type must be a non-empty character string")
  }
  if (is.null(feature_name) || !is.character(feature_name) || length(feature_name) == 0) {
    stop("feature_name must be a non-empty character string")
  }

  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    if (feature_type %in% c("drug")) {
      feature_data <- loadTreatmentResponseNormalized(dromaset_object,
                                           drugs = feature_name,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE)

      if (is.matrix(feature_data) && feature_name %in% rownames(feature_data)) {
        feature_vector <- as.numeric(feature_data[feature_name, ])
        names(feature_vector) <- colnames(feature_data)
        result <- list()
        result[[dromaset_object@name]] <- feature_vector[!is.na(feature_vector)]
        return(result)
      }
    } else {
      # Molecular profile data
      feature_data <- loadMolecularProfilesNormalized(dromaset_object,
                                           molecular_type = feature_type,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE)

      if (is_continuous) {
        # Continuous data
        if (is.matrix(feature_data) && feature_name %in% rownames(feature_data)) {
          feature_vector <- as.numeric(feature_data[feature_name, ])
          names(feature_vector) <- colnames(feature_data)
          result <- list()
          result[[dromaset_object@name]] <- feature_vector[!is.na(feature_vector)]
          return(result)
        }
      } else {
        # Discrete data
        if (is.data.frame(feature_data)) {
          if ("cells" %in% colnames(feature_data)) {
            sample_ids <- feature_data$cells[feature_data$genes == feature_name]
          } else if ("samples" %in% colnames(feature_data)) {
            sample_ids <- feature_data$samples[feature_data$genes == feature_name]
          } else {
            sample_ids <- character(0)
          }
          result <- list()
          result[[dromaset_object@name]] <- sample_ids
          return(result)
        }
      }
    }
  } else {
    # MultiDromaSet
    if (feature_type %in% c("drug")) {
      feature_data <- loadMultiProjectTreatmentResponseNormalized(dromaset_object,
                                                        drugs = feature_name,
                                                        overlap_only = overlap_only,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type)

      result <- lapply(feature_data, function(data_matrix) {
        if (is.matrix(data_matrix) && feature_name %in% rownames(data_matrix)) {
          data_vector <- as.numeric(data_matrix[feature_name, ])
          names(data_vector) <- colnames(data_matrix)
          return(data_vector[!is.na(data_vector)])
        }
        return(NULL)
      })
    } else {
      # Molecular profile data
      feature_data <- loadMultiProjectMolecularProfilesNormalized(dromaset_object,
                                                        molecular_type = feature_type,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type)

      if (is_continuous) {
        # Continuous data
        result <- lapply(feature_data, function(omics_matrix) {
          if (is.matrix(omics_matrix) && feature_name %in% rownames(omics_matrix)) {
            omics_vector <- as.numeric(omics_matrix[feature_name, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          }
          return(NULL)
        })
      } else {
        # Discrete data
        result <- lapply(feature_data, function(omics_df) {
          if (is.data.frame(omics_df)) {
            if ("cells" %in% colnames(omics_df)) {
              sample_ids <- omics_df$cells[omics_df$genes == feature_name]
            } else if ("samples" %in% colnames(omics_df)) {
              sample_ids <- omics_df$samples[omics_df$genes == feature_name]
            } else {
              sample_ids <- character(0)
            }
            return(sample_ids)
          }
          return(NULL)
        })
      }
    }

    # Remove NULL entries
    result <- result[!sapply(result, is.null)]
    return(result)
  }

  return(NULL)
}

#' Filter feature data to ensure minimum sample size
#' @param feature_data List of feature data by dataset
#' @param min_samples Minimum number of samples required per dataset
#' @return Filtered feature data
#' @keywords internal
filterFeatureData <- function(feature_data, min_samples = 3) {
  if (is.null(feature_data) || length(feature_data) == 0) {
    return(NULL)
  }

  filtered_data <- lapply(feature_data, function(dataset) {
    # Remove NA values
    dataset <- na.omit(dataset)
    # Check if the dataset has minimum samples after NA removal
    if (length(dataset) < min_samples) {
      return(NULL)
    } else {
      return(dataset)
    }
  })

  # Remove NULL entries
  filtered_data <- filtered_data[!sapply(filtered_data, is.null)]

  if (length(filtered_data) == 0) {
    return(NULL)
  }

  return(filtered_data)
}

#' Get sample metadata for discrete vs discrete analysis
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature1_type Type of first feature
#' @param feature2_type Type of second feature
#' @return Data frame with sample metadata
#' @keywords internal
getSampleMetadata <- function(dromaset_object, feature1_type, feature2_type) {
  # For DromaSet, create metadata from the object
  if (inherits(dromaset_object, "DromaSet")) {
    # Get all samples from the DromaSet
    all_samples <- c()

    # Try to get samples from molecular profiles
    available_profiles <- availableMolecularProfiles(dromaset_object)
    for (profile_type in available_profiles) {
      if (profile_type %in% c(feature1_type, feature2_type)) {
        profile_data <- loadMolecularProfilesNormalized(dromaset_object,
                                               molecular_type = profile_type,
                                               return_data = TRUE)
        if (is.matrix(profile_data)) {
          all_samples <- c(all_samples, colnames(profile_data))
        } else if (is.data.frame(profile_data)) {
          if ("cells" %in% colnames(profile_data)) {
            all_samples <- c(all_samples, unique(profile_data$cells))
          } else if ("samples" %in% colnames(profile_data)) {
            all_samples <- c(all_samples, unique(profile_data$samples))
          }
        }
      }
    }

    # Create metadata data frame
    if (length(all_samples) > 0) {
      all_samples <- unique(all_samples)
      return(data.frame(
        cells = all_samples,
        type = "sample",
        datasets = dromaset_object@name,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    # For MultiDromaSet, combine metadata from all DromaSets
    all_metadata <- data.frame()

    for (project_name in names(dromaset_object@DromaSets)) {
      dromaset <- dromaset_object@DromaSets[[project_name]]
      project_metadata <- getSampleMetadata(dromaset, feature1_type, feature2_type)
      if (!is.null(project_metadata) && nrow(project_metadata) > 0) {
        all_metadata <- rbind(all_metadata, project_metadata)
      }
    }

    return(all_metadata)
  }

  return(NULL)
}

# Main function----
#' Batch analysis to find significant features associated with a reference feature using DromaSet objects
#'
#' @description Performs batch analysis to identify features significantly associated with a reference feature using DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature1_type Type of the reference feature (e.g., "drug", "mRNA")
#' @param feature1_name Name of the reference feature
#' @param feature2_type Type of features to test against (e.g., "mRNA", "mutation_gene")
#' @param feature2_name Optional vector of specific feature names to test (default: NULL for all features)
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param cores Number of CPU cores to use for parallel processing
#' @param progress_callback Optional callback function for progress updates
#' @param test_top_n Integer, number of top features to test (default: NULL for all features)
#' @return Data frame with meta-analysis results (p_value, effect_size, name)
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' results <- batchFindSignificantFeatures(gCSI, "drug", "Paclitaxel", "mRNA")
#'
#' # Using MultiDromaSet
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA",
#'                                   feature2_name = c("TP53", "EGFR"))
#'
#' # Test only top 50 features for debugging
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA", test_top_n = 50)
#' }
batchFindSignificantFeatures <- function(dromaset_object,
                                      feature1_type,
                                      feature1_name,
                                      feature2_type,
                                      feature2_name = NULL,
                                      data_type = "all",
                                      tumor_type = "all",
                                      overlap_only = FALSE,
                                      cores = 1,
                                      progress_callback = NULL,
                                      test_top_n = NULL
) {
  # Validate input object
  if (is.null(dromaset_object)) {
    stop("dromaset_object cannot be NULL")
  }
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }

  # Validate inputs
  valid_feature_types <- c("mRNA", "cnv", "meth", "proteinrppa", "proteinms",
                           "drug", "mutation_gene", "mutation_site", "fusion")

  if(!all(c(feature1_type, feature2_type) %in% valid_feature_types)) {
    stop(paste0("The selected feature type doesn't exist. Please choose from: ",
                paste(valid_feature_types, collapse = ", ")))
  }

  valid_data_types <- c("all", "CellLine", "PDO", "PDC", "PDX")
  if(!data_type %in% valid_data_types) {
    stop(paste0("Invalid data_type. Please choose from: ",
                paste(valid_data_types, collapse = ", ")))
  }

  # Validate feature1_name parameter
  if (is.null(feature1_name) || !is.character(feature1_name) || length(feature1_name) == 0) {
    stop("feature1_name must be a non-empty character string")
  }

  # Validate feature2_name parameter if provided
  if (!is.null(feature2_name)) {
    if (!is.character(feature2_name) || length(feature2_name) == 0) {
      stop("feature2_name must be NULL or a non-empty character vector")
    }
    # Remove duplicates
    feature2_name <- unique(feature2_name)
  }

  # Validate cores parameter
  max_cores <- ifelse(requireNamespace("parallel", quietly = TRUE),
                    parallel::detectCores(),
                    1)
  if (!is.numeric(cores) || cores < 1 || cores > max_cores) {
    stop(paste0("cores must be a positive integer between 1 and ", max_cores))
  }

  # Validate progress callback
  if (!is.null(progress_callback) && !is.function(progress_callback)) {
    stop("progress_callback must be a function or NULL")
  }

  # Validate test_top_n parameter
  if (!is.null(test_top_n)) {
    if (!is.numeric(test_top_n) || length(test_top_n) != 1 || test_top_n < 1 || test_top_n != as.integer(test_top_n)) {
      stop("test_top_n must be a positive integer or NULL")
    }
  }

  # Track timing for progress updates
  start_time <- Sys.time()

  # Determine feature types
  continuous_types <- c("drug", "cnv", "proteinrppa",
                        "proteinms", "meth", "mRNA")
  is_continuous1 <- feature1_type %in% continuous_types
  is_continuous2 <- feature2_type %in% continuous_types

  # Get selected specific feature1 data using helper function
  selected_feas1 <- loadFeatureData(dromaset_object, feature1_type, feature1_name,
                                    data_type, tumor_type, overlap_only, is_continuous1)

  # Pre-load sample metadata for discrete vs discrete analysis
  samples_search <- NULL
  if (!is_continuous1 && !is_continuous2) {
    samples_search <- getSampleMetadata(dromaset_object, feature1_type, feature2_type)
  }

  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop("No data available for the selected feature 1.")
  }

  # Filter feature1 data to ensure minimum sample size
  selected_feas1 <- filterFeatureData(selected_feas1)

  # Check if any data remains after filtering
  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop(paste0("No sufficient data available for feature '", feature1_name,
                "' of type '", feature1_type, "'. Each dataset needs at least 3 samples. ",
                "Please try with a different feature."))
  }

  # Get list of features to test from the dromaset object
  feature2_list <- NULL

  # If feature2_name is provided, use those features directly
  if (!is.null(feature2_name)) {
    feature2_list <- feature2_name
  } else if (inherits(dromaset_object, "DromaSet")) {
    if (feature2_type %in% c("drug")) {
      # Get available drugs
      available_features <- availableTreatmentResponses(dromaset_object)
      if ("drug" %in% available_features) {
        # Load all drugs to get feature names
        all_drug_data <- loadTreatmentResponseNormalized(dromaset_object, return_data = TRUE)
        if (is.matrix(all_drug_data)) {
          feature2_list <- rownames(all_drug_data)
        }
      }
    } else {
      # Get available molecular profiles
      available_profiles <- availableMolecularProfiles(dromaset_object)
      if (feature2_type %in% available_profiles) {
        # Load all features to get feature names
        all_omics_data <- loadMolecularProfilesNormalized(dromaset_object,
                                               molecular_type = feature2_type,
                                               data_type = data_type,
                                               tumor_type = tumor_type,
                                               return_data = TRUE)
        if (is.matrix(all_omics_data)) {
          feature2_list <- rownames(all_omics_data)
        } else if (is.data.frame(all_omics_data) && "genes" %in% colnames(all_omics_data)) {
          feature2_list <- unique(all_omics_data$genes)
        }
      }
    }
  } else {
    # MultiDromaSet - get features from all projects that have the data (union)
    all_feature_lists <- list()
    
    for (project_name in names(dromaset_object@DromaSets)) {
      dromaset <- dromaset_object@DromaSets[[project_name]]

      if (feature2_type %in% c("drug")) {
        available_features <- availableTreatmentResponses(dromaset)
        if ("drug" %in% available_features) {
          all_drug_data <- loadTreatmentResponseNormalized(dromaset, return_data = TRUE)
          if (is.matrix(all_drug_data)) {
            all_feature_lists[[project_name]] <- rownames(all_drug_data)
          }
        }
      } else {
        available_profiles <- availableMolecularProfiles(dromaset)
        if (feature2_type %in% available_profiles) {
          all_omics_data <- loadMolecularProfilesNormalized(dromaset,
                                                 molecular_type = feature2_type,
                                                 data_type = data_type,
                                                 tumor_type = tumor_type,
                                                 return_data = TRUE)
          if (is.matrix(all_omics_data)) {
            all_feature_lists[[project_name]] <- rownames(all_omics_data)
          } else if (is.data.frame(all_omics_data) && "genes" %in% colnames(all_omics_data)) {
            all_feature_lists[[project_name]] <- unique(all_omics_data$genes)
          }
        }
      }
    }
    
    # Take union of all features across projects
    if (length(all_feature_lists) > 0) {
      feature2_list <- unique(unlist(all_feature_lists))
    }
  }

  if (is.null(feature2_list)) {
    stop("Could not find available features for feature2_type: ", feature2_type)
  }

  # If feature2_name was provided, validate that all features exist
  if (!is.null(feature2_name)) {
    # Get all available features to validate against
    all_available_features <- NULL
    if (inherits(dromaset_object, "DromaSet")) {
      if (feature2_type %in% c("drug")) {
        available_features <- availableTreatmentResponses(dromaset_object)
        if ("drug" %in% available_features) {
          all_drug_data <- loadTreatmentResponseNormalized(dromaset_object, return_data = TRUE)
          if (is.matrix(all_drug_data)) {
            all_available_features <- rownames(all_drug_data)
          }
        }
      } else {
        available_profiles <- availableMolecularProfiles(dromaset_object)
        if (feature2_type %in% available_profiles) {
          all_omics_data <- loadMolecularProfilesNormalized(dromaset_object,
                                                 molecular_type = feature2_type,
                                                 data_type = data_type,
                                                 tumor_type = tumor_type,
                                                 return_data = TRUE)
          if (is.matrix(all_omics_data)) {
            all_available_features <- rownames(all_omics_data)
          } else if (is.data.frame(all_omics_data) && "genes" %in% colnames(all_omics_data)) {
            all_available_features <- unique(all_omics_data$genes)
          }
        }
      }
    } else {
      # MultiDromaSet - get features from all projects that have the data (union)
      all_feature_lists <- list()

      for (project_name in names(dromaset_object@DromaSets)) {
        dromaset <- dromaset_object@DromaSets[[project_name]]

        if (feature2_type %in% c("drug")) {
          available_features <- availableTreatmentResponses(dromaset)
          if ("drug" %in% available_features) {
            all_drug_data <- loadTreatmentResponseNormalized(dromaset, return_data = TRUE)
            if (is.matrix(all_drug_data)) {
              all_feature_lists[[project_name]] <- rownames(all_drug_data)
            }
          }
        } else {
          available_profiles <- availableMolecularProfiles(dromaset)
          if (feature2_type %in% available_profiles) {
            all_omics_data <- loadMolecularProfilesNormalized(dromaset,
                                                   molecular_type = feature2_type,
                                                   data_type = data_type,
                                                   tumor_type = tumor_type,
                                                   return_data = TRUE)
            if (is.matrix(all_omics_data)) {
              all_feature_lists[[project_name]] <- rownames(all_omics_data)
            } else if (is.data.frame(all_omics_data) && "genes" %in% colnames(all_omics_data)) {
              all_feature_lists[[project_name]] <- unique(all_omics_data$genes)
            }
          }
        }
      }

      # Take union of all features across projects
      if (length(all_feature_lists) > 0) {
        all_available_features <- unique(unlist(all_feature_lists))
      }
    }

    # Check which requested features don't exist
    if (is.null(all_available_features)) {
      stop("No features found for feature2_type: ", feature2_type)
    }

    invalid_features <- feature2_name[!feature2_name %in% all_available_features]
    if (length(invalid_features) > 0) {
      warning("The following feature2_name features were not found and will be skipped: ",
              paste(invalid_features, collapse = ", "))
      # Keep only valid features
      feature2_list <- feature2_list[feature2_list %in% all_available_features]
    }
  }

  # Apply test_top_n filter
  if(!is.null(test_top_n) && is.numeric(test_top_n) && test_top_n > 0 && length(feature2_list) > test_top_n){
    feature2_list <- feature2_list[1:min(test_top_n, length(feature2_list))]
  }

  if(length(feature2_list) == 0) {
    stop("No features found for the selected feature 2 type or all specified feature2_name features are invalid.")
  }

  # Define the worker function
  worker_function <- function(x) {
    results <- tryCatch({
      # Validate input index
      if (!is.numeric(x) || x < 1 || x > length(feature2_list)) {
        stop(sprintf("Invalid feature index: %d", x))
      }

      feature2_name <- feature2_list[x]

      # Load feature2 data using helper function
      selected_feas2 <- loadFeatureData(dromaset_object, feature2_type, feature2_name,
                                       data_type, tumor_type, overlap_only, is_continuous2)

      if (is.null(selected_feas2) || length(selected_feas2) == 0) {
        return(NULL)
      }

      # Filter feature2 data to ensure minimum sample size
      selected_feas2 <- filterFeatureData(selected_feas2)

      # Skip if no sufficient data remains after filtering
      if(is.null(selected_feas2) || length(selected_feas2) == 0) {
        return(NULL)
      }

      # do statistics test based on four circumstances
      # con vs con ----
      if (is_continuous1 && is_continuous2) {
        selected_pair <- pairContinuousFeatures(selected_feas1, selected_feas2)
        cal_meta_re <- metaCalcConCon(selected_pair)
        # dis vs con ----
      } else if ((is_continuous1 && !is_continuous2) || (!is_continuous1 && is_continuous2)) {
        if (is_continuous1 && !is_continuous2){
          selected_pair <- pairDiscreteFeatures(selected_feas2, selected_feas1)
        } else {
          selected_pair <- pairDiscreteFeatures(selected_feas1, selected_feas2)
        }
        cal_meta_re <- metaCalcConDis(selected_pair)
        # dis vs dis ----
      } else {
        # Discrete vs discrete analysis
        # Get sample metadata needed for pairing
        if (is.null(samples_search)) {
          # Get sample metadata if not already available
          samples_search <- getSampleMetadata(dromaset_object, feature1_type, feature2_type)
          if (is.null(samples_search) || nrow(samples_search) == 0) {
            return(NULL)
          }
        }

        selected_pair <- pairDiscreteDiscrete(selected_feas1, selected_feas2,
                                             feature1_type, feature2_type,
                                             samples_search)
        cal_meta_re <- metaCalcDisDis(selected_pair)
      }

      if(is.null(cal_meta_re)) return(NULL)
      results <- data.frame(
        p_value = cal_meta_re[["pval.random"]],
        effect_size = cal_meta_re[["TE.random"]]
      #   N = length(cal_meta_re[["studlab"]])
      )
    }, error = function(e) {
      # Log detailed error but continue processing
      error_msg <- sprintf("Error processing feature %d (%s) of type %s: %s",
                          x, feature2_list[x], feature2_type, e$message)
      message(error_msg)
      NULL
    })

    # Update progress if callback provided
    if(!is.null(progress_callback)) {
      progress_callback(x, length(feature2_list),
                        difftime(Sys.time(), start_time, units = "secs"))
    }

    return(results)
  }
  message("Please be patient, it may take long time to run.")
  # Use parallel processing if cores > 1
  if (cores > 1) {
    # Check if snowfall is available
    if (!requireNamespace("snowfall", quietly = TRUE)) {
      warning("snowfall package not available. Running sequentially instead.")
      cores <- 1
    } else {
      # Initialize snowfall
      tryCatch({
        sfInit(parallel = TRUE, cpus = cores)

        # Load required packages on worker nodes
        sfLibrary(meta)
        sfLibrary(metafor)
        sfLibrary(effsize)
        sfLibrary(DROMA.Set)
        sfLibrary(DROMA.R)

        # Export required data and functions
        sfExport("dromaset_object", "selected_feas1", "feature2_list",
                 "is_continuous1", "is_continuous2",
                 "feature1_type", "feature2_type", "data_type", "tumor_type", "overlap_only",
                 "start_time", "samples_search")

        # Run parallel computation
        cal_re_list <- sfLapply(1:length(feature2_list), worker_function)

        # Clean up snowfall
        sfStop()
      }, error = function(e) {
        message("Parallel processing failed. Running sequentially instead.")
        message("Error details: ", e$message)
        sfStop()  # Clean up any existing parallel cluster
        # Run sequential computation
        cal_re_list <- lapply(1:length(feature2_list), worker_function)
      })
    }
  }

  if (cores == 1) {
    # Run sequential computation
    cal_re_list <- lapply(1:length(feature2_list), worker_function)
  }

  # Process results
  valid_results <- !sapply(cal_re_list, is.null)
  fea_names <- feature2_list[valid_results]
  cal_re_list <- cal_re_list[valid_results]

  if (length(cal_re_list) == 0) {
    stop("No valid results found for the selected features. Please try others.")
  }

  cal_re_df <- do.call(rbind, cal_re_list)
  cal_re_df$name <- fea_names

  # Log completion
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  message(sprintf("Analysis completed in %s", formatTime(as.numeric(total_time))))
  message(sprintf("Found %d significant associations out of %d features.",
                  sum(cal_re_df$p_value < 0.05), nrow(cal_re_df)))

  return(cal_re_df)
}

