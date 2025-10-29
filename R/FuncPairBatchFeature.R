# Helper Functions for Data Loading ----

#' Process a single feature pair for batch analysis
#' @param feature_idx Index of feature2 in feature2_list
#' @param feature2_list Vector of feature2 names to test
#' @param dromaset_object DromaSet or MultiDromaSet object
#' @param selected_feas1 Feature1 data
#' @param feature1_type Type of feature1
#' @param feature2_type Type of feature2
#' @param is_continuous1 Whether feature1 is continuous
#' @param is_continuous2 Whether feature2 is continuous
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @param overlap_only Whether to use only overlapping samples
#' @param samples_search Sample metadata for discrete vs discrete analysis
#' @param sample_intersection_cache Cached sample intersections
#' @param progress_callback Optional progress callback function
#' @param start_time Start time for progress tracking
#' @return Data frame with p_value, effect_size, n_datasets or NULL
#' @keywords internal
processFeaturePair <- function(feature_idx, feature2_list, dromaset_object,
                               selected_feas1,
                               feature1_type, feature2_type,
                               is_continuous1, is_continuous2,
                               data_type, tumor_type, overlap_only,
                               samples_search, sample_intersection_cache,
                               progress_callback = NULL, start_time = NULL) {
  results <- tryCatch({
    # Validate input index
    if (!is.numeric(feature_idx) || feature_idx < 1 || feature_idx > length(feature2_list)) {
      stop(sprintf("Invalid feature index: %d", feature_idx))
    }

    feature2_name <- feature2_list[feature_idx]

    # Load feature2 data on-demand
    selected_feas2 <- loadFeatureData(dromaset_object, feature2_type, feature2_name,
                                     data_type, tumor_type, overlap_only, is_continuous2,
                                     return_all_samples = !is_continuous2)

    if (is.null(selected_feas2) || length(selected_feas2) == 0) {
      return(NULL)
    }

    # Filter feature2 data to ensure minimum sample size
    selected_feas2 <- filterFeatureData(selected_feas2, min_samples = 3, is_discrete_with_all = !is_continuous2)

    # Skip if no sufficient data remains after filtering
    if(is.null(selected_feas2) || length(selected_feas2) == 0) {
      return(NULL)
    }

    # do statistics test based on four circumstances
    # con vs con ----
    if (is_continuous1 && is_continuous2) {
      selected_pair <- pairContinuousFeatures(selected_feas1, selected_feas2,
                            merged = FALSE,
                            intersection_cache = sample_intersection_cache)
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
    
    # Extract p_value and effect_size
    p_val <- cal_meta_re[["pval.random"]]
    eff_size <- cal_meta_re[["TE.random"]]
    n_datasets <- length(cal_meta_re[["studlab"]])
    
    # Replace NA values: p_value -> 1, effect_size -> 0
    if(is.na(p_val)) p_val <- 1
    if(is.na(eff_size)) eff_size <- 0
    
    results <- data.frame(
      p_value = p_val,
      effect_size = eff_size,
      n_datasets = n_datasets
    )
  }, error = function(e) {
    # Silently continue processing on error
    NULL
  })

  # Update progress if callback provided
  if(!is.null(progress_callback) && !is.null(start_time)) {
    progress_callback(feature_idx, length(feature2_list),
                      difftime(Sys.time(), start_time, units = "secs"))
  }

  return(results)
}

#' Get and validate feature list from DromaSet object
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature_type Type of features to retrieve
#' @param feature_names Optional vector of specific feature names to validate (default: NULL for all features)
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @return Character vector of valid feature names
#' @keywords internal
getAndValidateFeatureList <- function(dromaset_object, feature_type, feature_names = NULL,
                                      data_type = "all", tumor_type = "all") {
  feature_list <- NULL
  
  # If specific feature names are provided, use them as the starting point
  if (!is.null(feature_names)) {
    feature_list <- unique(feature_names)
  } else {
    # Get feature list using listDROMAFeatures (optimized approach)
    if (inherits(dromaset_object, "DromaSet")) {
      # For DromaSet, use listDROMAFeatures
      tryCatch({
        feature_list <- listDROMAFeatures(dromaset_object@name, feature_type,
                        data_type = data_type, tumor_type = tumor_type)
      }, error = function(e) {
        feature_list <<- NULL
      })
    } else {
      # MultiDromaSet - get features from all projects
      all_feature_lists <- list()
      for (projects in names(dromaset_object@DromaSets)) {
        tryCatch({
          project_features <- listDROMAFeatures(projects, feature_type,
                          data_type = data_type, tumor_type = tumor_type)
          if (length(project_features) > 0) {
            all_feature_lists[[projects]] <- project_features
          }
        }, error = function(e) {
          # Silently continue on error
        })
      }
      # Take union of all features across projects
      if (length(all_feature_lists) > 0) {
        feature_list <- unique(unlist(all_feature_lists))
      }
    }
  }
  
  # Validate that we have a feature list
  if (is.null(feature_list) || length(feature_list) == 0) {
    stop("Could not retrieve feature list for feature_type: ", feature_type)
  }
  
  # If specific feature names were provided, validate they exist
  if (!is.null(feature_names)) {
    # Get all available features for validation
    all_available_features <- NULL
    
    if (inherits(dromaset_object, "DromaSet")) {
      tryCatch({
        all_available_features <- listDROMAFeatures(dromaset_object@name, feature_type,
                        data_type = data_type, tumor_type = tumor_type)
      }, error = function(e) {
        # Silently continue on error
      })
    } else {
      # MultiDromaSet
      all_feature_lists <- list()
      for (projects in names(dromaset_object@DromaSets)) {
        tryCatch({
          project_features <- listDROMAFeatures(projects, feature_type,
                          data_type = data_type, tumor_type = tumor_type)
          if (length(project_features) > 0) {
            all_feature_lists[[projects]] <- project_features
          }
        }, error = function(e) {
          # Silently continue on error
        })
      }
      if (length(all_feature_lists) > 0) {
        all_available_features <- unique(unlist(all_feature_lists))
      }
    }
    
    # Check for invalid features
    if (!is.null(all_available_features)) {
      invalid_features <- feature_names[!feature_names %in% all_available_features]
      if (length(invalid_features) > 0) {
        warning("The following features were not found and will be skipped: ",
                paste(invalid_features, collapse = ", "))
        # Keep only valid features
        feature_list <- feature_list[feature_list %in% all_available_features]
      }
      
      if (length(feature_list) == 0) {
        stop("None of the specified features are available")
      }
    }
  }
  
  return(feature_list)
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
        profile_data <- loadMolecularProfiles(dromaset_object,
                                               feature_type = profile_type,
                                               return_data = TRUE,
                                               zscore = TRUE)
        if (is.matrix(profile_data)) {
          all_samples <- c(all_samples, colnames(profile_data))
        } else if (is.data.frame(profile_data)) {
          if ("samples" %in% colnames(profile_data)) {
            all_samples <- c(all_samples, unique(profile_data$samples))
          }
        }
      }
    }

    # Create metadata data frame
    if (length(all_samples) > 0) {
      all_samples <- unique(all_samples)
      return(data.frame(
        samples = all_samples,
        type = "sample",
        datasets = dromaset_object@name,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    # For MultiDromaSet, combine metadata from all DromaSets
    all_metadata <- data.frame()

    for (projects in names(dromaset_object@DromaSets)) {
      dromaset <- dromaset_object@DromaSets[[projects]]
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
#'
#' @section Execution Flow:
#' The function follows this high-level workflow:
#' 1. Input validation and parameter setup
#' 2. Feature type determination (continuous vs discrete)
#' 3. Load and validate reference feature (feature1) data
#' 4. Get list of features to test (feature2) - uses listDROMAFeatures() for efficient database queries
#' 5. Parallel/sequential processing of feature pairs with on-demand data loading
#' 6. Statistical testing based on feature types
#' 7. Result aggregation and reporting
#'
#' @section Optimizations:
#' This function implements several optimizations for better performance and memory usage:
#' \itemize{
#'   \item{Feature List Retrieval: Uses listDROMAFeatures() from DROMA.Set package to efficiently 
#'         query available features directly from the database without loading full data matrices}
#'   \item{On-Demand Data Loading: Loads feature2 data only when needed for each feature pair,
#'         reducing memory footprint especially for large datasets with thousands of features}
#'   \item{Automatic Fallback: If database queries fail, automatically falls back to traditional
#'         data loading methods for compatibility}
#' }
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
#' @param verbose Logical, whether to display detailed messages from internal functions (default: FALSE)
#' @return Data frame with meta-analysis results (p_value, effect_size, n_datasets, name, q_value). NA values are replaced: p_value -> 1, effect_size -> 0. q_value is calculated using Benjamini-Hochberg method
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
                                      test_top_n = NULL,
                                      verbose = FALSE
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
  if (!verbose) {
    selected_feas1 <- suppressMessages(
      loadFeatureData(dromaset_object, feature1_type, feature1_name,
                      data_type, tumor_type, overlap_only, is_continuous1,
                      return_all_samples = !is_continuous1)
    )
  } else {
    selected_feas1 <- loadFeatureData(dromaset_object, feature1_type, feature1_name,
                                      data_type, tumor_type, overlap_only, is_continuous1,
                                      return_all_samples = !is_continuous1)
  }

  # Pre-load sample metadata for discrete vs discrete analysis
  samples_search <- NULL
  if (!is_continuous1 && !is_continuous2) {
    if (!verbose) {
      samples_search <- suppressMessages(
        getSampleMetadata(dromaset_object, feature1_type, feature2_type)
      )
    } else {
      samples_search <- getSampleMetadata(dromaset_object, feature1_type, feature2_type)
    }
  }

  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop("No data available for the selected feature 1.")
  }

  # Filter feature1 data to ensure minimum sample size
  if (!verbose) {
    selected_feas1 <- suppressMessages(
      filterFeatureData(selected_feas1, min_samples = 3, is_discrete_with_all = !is_continuous1)
    )
  } else {
    selected_feas1 <- filterFeatureData(selected_feas1, min_samples = 3, is_discrete_with_all = !is_continuous1)
  }

  # Check if any data remains after filtering
  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop(paste0("No sufficient data available for feature '", feature1_name,
                "' of type '", feature1_type, "'. Each dataset needs at least 3 samples. ",
                "Please try with a different feature."))
  }

  # Get and validate feature2 list using optimized helper function
  if (!verbose) {
    feature2_list <- suppressMessages(
      getAndValidateFeatureList(
        dromaset_object = dromaset_object,
        feature_type = feature2_type,
        feature_names = feature2_name,
        data_type = data_type,
        tumor_type = tumor_type
      )
    )
  } else {
    feature2_list <- getAndValidateFeatureList(
      dromaset_object = dromaset_object,
      feature_type = feature2_type,
      feature_names = feature2_name,
      data_type = data_type,
      tumor_type = tumor_type
    )
  }

  # Apply test_top_n filter
  if(!is.null(test_top_n) && is.numeric(test_top_n) && test_top_n > 0 && length(feature2_list) > test_top_n){
    feature2_list <- feature2_list[seq_len(min(test_top_n, length(feature2_list)))]
  }

  if(length(feature2_list) == 0) {
    stop("No features found for the selected feature 2 type or all specified feature2_name features are invalid.")
  }

  message("Starting batch analysis with on-demand data loading...")

  # Sample intersection cache - set to NULL for on-demand computation
  # Sample intersections will be computed as needed during pairing
  sample_intersection_cache <- NULL

  # Define the worker function as a wrapper for processFeaturePair
  worker_function <- function(x) {
    if (!verbose) {
      suppressMessages(
        processFeaturePair(
          feature_idx = x,
          feature2_list = feature2_list,
          dromaset_object = dromaset_object,
          selected_feas1 = selected_feas1,
          feature1_type = feature1_type,
          feature2_type = feature2_type,
          is_continuous1 = is_continuous1,
          is_continuous2 = is_continuous2,
          data_type = data_type,
          tumor_type = tumor_type,
          overlap_only = overlap_only,
          samples_search = samples_search,
          sample_intersection_cache = sample_intersection_cache,
          progress_callback = progress_callback,
          start_time = start_time
        )
      )
    } else {
      processFeaturePair(
        feature_idx = x,
        feature2_list = feature2_list,
        dromaset_object = dromaset_object,
        selected_feas1 = selected_feas1,
        feature1_type = feature1_type,
        feature2_type = feature2_type,
        is_continuous1 = is_continuous1,
        is_continuous2 = is_continuous2,
        data_type = data_type,
        tumor_type = tumor_type,
        overlap_only = overlap_only,
        samples_search = samples_search,
        sample_intersection_cache = sample_intersection_cache,
        progress_callback = progress_callback,
        start_time = start_time
      )
    }
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
        # Note: export dromaset_object for feature2 loading
        sfExport("selected_feas1", "feature2_list",
                 "is_continuous1", "is_continuous2",
                 "feature1_type", "feature2_type", "data_type", "tumor_type", "overlap_only",
                 "start_time", "samples_search", "sample_intersection_cache", "dromaset_object",
                 "verbose")

        # Run parallel computation
        cal_re_list <- sfLapply(seq_along(feature2_list), worker_function)

        # Clean up snowfall
        sfStop()
      }, error = function(e) {
        message("Parallel processing failed. Running sequentially instead.")
        message("Error details: ", e$message)
        sfStop()  # Clean up any existing parallel cluster
        # Run sequential computation
        cal_re_list <- lapply(seq_along(feature2_list), worker_function)
      })
    }
  }

  if (cores == 1) {
    # Run sequential computation
    cal_re_list <- lapply(seq_along(feature2_list), worker_function)
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
  cal_re_df$q_value <- p.adjust(cal_re_df$p_value, method = "BH")

  # Log completion
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  message(sprintf("Analysis completed in %s", formatTime(as.numeric(total_time))))
  message(sprintf("Found %d significant associations out of %d features.",
                  sum(cal_re_df$q_value < 0.01), nrow(cal_re_df)))

  return(cal_re_df)
}

