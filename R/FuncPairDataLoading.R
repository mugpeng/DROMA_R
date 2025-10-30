# Data Loading Helper Functions ----

#' Load feature data from DromaSet or MultiDromaSet
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature_type Type of feature to load
#' @param select_features Name(s) of feature(s) to load. Can be a single feature name or a vector of feature names. If multiple features are provided for continuous data, returns a matrix with all features.
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples
#' @param is_continuous Whether the feature is continuous
#' @param return_all_samples For discrete features, whether to return all profiled samples in addition to feature-present samples (default: FALSE)
#' @param zscore Whether to apply z-score normalization (default: TRUE)
#' @return For continuous data: single feature returns list of named vectors by dataset; multiple features return list of matrices by dataset. For discrete data: single feature returns vector of sample IDs (or list with 'present' and 'all' if return_all_samples=TRUE); multiple features return named list of sample vectors per feature (or list with 'present' as named list and 'all' if return_all_samples=TRUE).
#' @keywords internal
loadFeatureData <- function(dromaset_object, feature_type, select_features,
                           data_type = "all", tumor_type = "all",
                           overlap_only = FALSE, is_continuous = TRUE,
                           return_all_samples = FALSE, zscore = TRUE) {

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
  if (is.null(select_features) || !is.character(select_features) || length(select_features) == 0) {
    stop("select_features must be a non-empty character string or character vector")
  }
  
  # Detect mode: single feature vs multiple features
  is_batch_mode <- length(select_features) > 1

  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    if (feature_type %in% c("drug")) {
      feature_data <- loadTreatmentResponse(dromaset_object,
                                           select_drugs = select_features,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE,
                                           zscore = zscore)

      if (is.matrix(feature_data) && nrow(feature_data) > 0) {
        result <- list()
        if (is_batch_mode) {
          # Return full matrix for batch mode
          result[[dromaset_object@name]] <- feature_data
        } else {
          # Return vector for single feature mode
          feature_vector <- as.numeric(feature_data[1, ])
          names(feature_vector) <- colnames(feature_data)
          result[[dromaset_object@name]] <- feature_vector[!is.na(feature_vector)]
        }
        return(result)
      }
    } else {
      # Molecular profile data
      feature_data <- loadMolecularProfiles(dromaset_object,
                                           feature_type = feature_type,
                                           select_features = select_features,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE,
                                           zscore = zscore)

      if (is_continuous) {
        # Continuous data
        if (is.matrix(feature_data) && nrow(feature_data) > 0) {
          result <- list()
          if (is_batch_mode) {
            # Return full matrix for batch mode (preload scenario)
            result[[dromaset_object@name]] <- feature_data
          } else {
            # Return vector for single feature mode (backward compatible)
            feature_vector <- as.numeric(feature_data[1, ])
            names(feature_vector) <- colnames(feature_data)
            result[[dromaset_object@name]] <- feature_vector[!is.na(feature_vector)]
          }
          return(result)
        }
      } else {
        # Discrete data - long dataframe format with samples and features columns
        if (is.data.frame(feature_data) && "samples" %in% colnames(feature_data)) {
          result <- list()
          
          if (is_batch_mode) {
            # Batch mode: organize by feature
            present_list <- lapply(select_features, function(feat) {
              feature_data$samples[feature_data$features == feat]
            })
            names(present_list) <- select_features
            
            if (return_all_samples) {
              all_profiled_samples <- listDROMASamples(dromaset_object@name,
                                                       feature_type = feature_type,
                                                       data_type = data_type,
                                                       tumor_type = tumor_type)
              result[[dromaset_object@name]] <- list(present = present_list, all = all_profiled_samples)
            } else {
              result[[dromaset_object@name]] <- present_list
            }
          } else {
            # Single mode: backward compatible
            present_samples <- feature_data$samples
            
            if (return_all_samples) {
              all_profiled_samples <- listDROMASamples(dromaset_object@name,
                                                       feature_type = feature_type,
                                                       data_type = data_type,
                                                       tumor_type = tumor_type)
              result[[dromaset_object@name]] <- list(present = present_samples, all = all_profiled_samples)
            } else {
              result[[dromaset_object@name]] <- present_samples
            }
          }
          return(result)
        }
      }
    }
  } else {
    # MultiDromaSet
    if (feature_type %in% c("drug")) {
      feature_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                        select_drugs = select_features,
                                                        overlap_only = overlap_only,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type,
                                                        zscore = zscore)

      result <- lapply(feature_data, function(data_matrix) {
        if (is.matrix(data_matrix) && nrow(data_matrix) > 0) {
          if (is_batch_mode) {
            # Return full matrix for batch mode
            return(data_matrix)
          } else {
            # Return vector for single feature mode
            data_vector <- as.numeric(data_matrix[1, ])
            names(data_vector) <- colnames(data_matrix)
            return(data_vector[!is.na(data_vector)])
          }
        }
        return(NULL)
      })
    } else {
      # Molecular profile data
      feature_data <- loadMultiProjectMolecularProfiles(dromaset_object,
                                                        feature_type = feature_type,
                                                        select_features = select_features,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type,
                                                        zscore = zscore)

      if (is_continuous) {
        # Continuous data
        result <- lapply(feature_data, function(omics_matrix) {
          if (is.matrix(omics_matrix) && nrow(omics_matrix) > 0) {
            if (is_batch_mode) {
              # Return full matrix for batch mode (preload scenario)
              return(omics_matrix)
            } else {
              # Return vector for single feature mode (backward compatible)
              omics_vector <- as.numeric(omics_matrix[1, ])
              names(omics_vector) <- colnames(omics_matrix)
              return(omics_vector[!is.na(omics_vector)])
            }
          }
          return(NULL)
        })
      } else {
        # Discrete data - long dataframe format with samples and features columns
        project_names <- names(feature_data)
        result <- lapply(seq_along(feature_data), function(i) {
          omics_df <- feature_data[[i]]
          projects <- project_names[i]
          if (is.data.frame(omics_df) && "samples" %in% colnames(omics_df)) {
            
            if (is_batch_mode) {
              # Batch mode: organize by feature
              present_list <- lapply(select_features, function(feat) {
                omics_df$samples[omics_df$features == feat]
              })
              names(present_list) <- select_features
              
              if (return_all_samples) {
                all_profiled_samples <- listDROMASamples(dromaset_object@DromaSets[[projects]]@name,
                                                         feature_type = feature_type,
                                                         data_type = data_type,
                                                         tumor_type = tumor_type)
                return(list(present = present_list, all = all_profiled_samples))
              } else {
                return(present_list)
              }
            } else {
              # Single mode: backward compatible
              present_samples <- omics_df$samples
              
              if (return_all_samples) {
                all_profiled_samples <- listDROMASamples(dromaset_object@DromaSets[[projects]]@name,
                                                         feature_type = feature_type,
                                                         data_type = data_type,
                                                         tumor_type = tumor_type)
                return(list(present = present_samples, all = all_profiled_samples))
              } else {
                return(present_samples)
              }
            }
          }
          return(NULL)
        })
        names(result) <- project_names
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
#' @param is_discrete_with_all For discrete data with 'present'/'all' format, whether to apply special handling (default: FALSE)
#' @return Filtered feature data
#' @keywords internal
filterFeatureData <- function(feature_data, min_samples = 3, is_discrete_with_all = FALSE) {
  if (is.null(feature_data) || length(feature_data) == 0) {
    return(NULL)
  }

  filtered_data <- lapply(feature_data, function(dataset) {
    # Handle discrete data with present/all format
    if (is_discrete_with_all && is.list(dataset) && "present" %in% names(dataset) && "all" %in% names(dataset)) {
      # For discrete data, check if we have minimum samples in the 'all' profiled samples
      if (length(dataset$all) < min_samples) {
        return(NULL)
      }
      return(dataset)
    } else {
      # Handle continuous data or simple discrete data
      # Remove NA values
      dataset <- na.omit(dataset)
      # Check if the dataset has minimum samples after NA removal
      if (length(dataset) < min_samples) {
        return(NULL)
      } else {
        return(dataset)
      }
    }
  })

  # Remove NULL entries
  filtered_data <- filtered_data[!sapply(filtered_data, is.null)]

  if (length(filtered_data) == 0) {
    return(NULL)
  }

  return(filtered_data)
}

