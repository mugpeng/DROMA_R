# Helper Functions for Data Loading ----

#' Create default progress callback function
#'
#' @description Creates a progress tracking function that displays processing status,
#'   elapsed time, and estimated time remaining
#' @param show_progress Logical, whether to show progress messages (default: TRUE)
#' @param update_interval Integer, minimum seconds between progress updates (default: 10)
#' @return A progress callback function
#' @keywords internal
createDefaultProgressCallback <- function(show_progress = TRUE, update_interval = 10) {
  if (!show_progress) {
    # Return a no-op function
    return(function(done, total, elapsed) {})
  }
  
  # Initialize tracking variables in closure
  last_update_time <- 0
  
  # Return the callback function
  function(done, total, elapsed) {
    # Only update if enough time has passed since last update
    if (elapsed - last_update_time < update_interval && done < total) {
      return(invisible(NULL))
    }
    
    # Update last update time
    last_update_time <<- elapsed
    
    # Calculate progress
    progress_pct <- (done / total) * 100
    
    # Estimate remaining time
    time_remaining <- estimateTimeRemaining(done, total, elapsed)
    
    # Format messages
    elapsed_str <- formatTime(elapsed)
    remaining_str <- if (is.finite(time_remaining)) {
      formatTime(time_remaining)
    } else {
      "calculating..."
    }
    
    # Display progress
    if (done < total) {
      message(sprintf(
        "Progress: %d/%d (%.1f%%) | Elapsed: %s | ETA: %s",
        done, total, progress_pct, elapsed_str, remaining_str
      ))
    } else {
      # Final message
      message(sprintf(
        "Progress: %d/%d (100%%) | Total time: %s",
        done, total, elapsed_str
      ))
    }
    
    invisible(NULL)
  }
}

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
#' @param preloaded_data Preloaded feature2 data (NULL for on-demand loading)
#' @return Data frame with p_value, effect_size, n_datasets or NULL
#' @keywords internal
processFeaturePair <- function(feature_idx, feature2_list, dromaset_object,
                               selected_feas1,
                               feature1_type, feature2_type,
                               is_continuous1, is_continuous2,
                               data_type, tumor_type, overlap_only,
                               samples_search, sample_intersection_cache,
                               preloaded_data = NULL) {
  results <- tryCatch({
    # Validate input index
    if (!is.numeric(feature_idx) || feature_idx < 1 || feature_idx > length(feature2_list)) {
      stop(sprintf("Invalid feature index: %d", feature_idx))
    }

    feature2_name <- feature2_list[feature_idx]

    # Load feature2 data: use preloaded data if available, otherwise load on-demand
    if (!is.null(preloaded_data)) {
      if (is_continuous2) {
        # Extract from preloaded continuous data
        # Preloaded_data is a list of matrices (features as rows, samples as columns)
        selected_feas2 <- lapply(preloaded_data, function(data_matrix) {
          if (is.matrix(data_matrix) && feature2_name %in% rownames(data_matrix)) {
            feature_vector <- as.numeric(data_matrix[feature2_name, ])
            names(feature_vector) <- colnames(data_matrix)
            return(feature_vector[!is.na(feature_vector)])
          }
          return(NULL)
        })
        selected_feas2 <- selected_feas2[!sapply(selected_feas2, is.null)]
      } else {
        # Extract from preloaded discrete data
        # Preloaded_data format: list(present = list(feat1 = vec1, feat2 = vec2, ...), all = vec) per dataset
        selected_feas2 <- lapply(preloaded_data, function(dataset) {
          if (is.list(dataset) && "present" %in% names(dataset)) {
            # Extract the specific feature from present list
            if (feature2_name %in% names(dataset$present)) {
              return(list(present = dataset$present[[feature2_name]], all = dataset$all))
            }
          }
          return(NULL)
        })
        selected_feas2 <- selected_feas2[!sapply(selected_feas2, is.null)]
      }
    } else {
      # Load feature2 data on-demand (continuous features only)
      if (!is_continuous2) {
        stop("Discrete features require preloaded data. This should not happen.")
      }
      selected_feas2 <- loadFeatureData(dromaset_object, feature2_type, feature2_name,
                                       data_type, tumor_type, overlap_only, is_continuous2,
                                       return_all_samples = !is_continuous2)
    }

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

#' Get significant features from meta-analysis results
#'
#' @description Filters meta-analysis results to identify significant features based on effect size, q-value or p-value, and optionally minimum number of datasets
#' @param meta_df Data frame containing meta-analysis results with columns: effect_size, q_value (or p_value), n_datasets, and name
#' @param es_t Effect size threshold (default: 0.4)
#' @param P_t Q-value or p-value threshold (default: 0.01)
#' @param n_datasets_t Minimum number of datasets threshold (default: NULL for no threshold)
#' @param use_p_value Logical, whether to use p_value column instead of q_value (default: FALSE)
#' @return Data frame with significant features only, including a "direction" column ("Up" or "Down")
#' @export
#' @examples
#' \dontrun{
#' # Get all significant features using q-values
#' sig_features <- getSignificantFeatures(meta_df)
#'
#' # Get all significant features using p-values
#' sig_features <- getSignificantFeatures(meta_df, use_p_value = TRUE)
#'
#' # With custom thresholds
#' sig_features <- getSignificantFeatures(meta_df, es_t = 0.5, P_t = 0.05)
#'
#' # Require features in at least 3 datasets
#' sig_features <- getSignificantFeatures(meta_df, n_datasets_t = 3)
#' }
getSignificantFeatures <- function(meta_df,
                                   es_t = 0.4,
                                   P_t = 0.01,
                                   n_datasets_t = NULL,
                                   use_p_value = FALSE) {
  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  required_cols <- c("effect_size", "name", "n_datasets")
  p_val_col <- if(use_p_value) "p_value" else "q_value"
  if(!all(c(required_cols, p_val_col) %in% colnames(meta_df))) {
    stop("meta_df must contain columns: effect_size, ", p_val_col, ", n_datasets, and name")
  }
  # Transfrom to data.frame
  meta_df <- as.data.frame(meta_df)
  # Filter based on thresholds
  if(is.null(n_datasets_t)) {
    sig_df <- meta_df[(abs(meta_df$effect_size) > es_t & meta_df[[p_val_col]] < P_t), ]
  } else {
    sig_df <- meta_df[(abs(meta_df$effect_size) > es_t &
                       meta_df[[p_val_col]] < P_t &
                       meta_df$n_datasets >= n_datasets_t), ]
  }
  
  # Add direction column
  if(nrow(sig_df) > 0) {
    sig_df$direction <- ifelse(sig_df$effect_size > 0, "Up", "Down")
    # Reorder columns to put direction after name
    col_order <- c("name", "direction", setdiff(names(sig_df), c("name", "direction")))
    sig_df <- sig_df[, col_order]
  }
  
  return(sig_df)
}

#' Get intersection of significant features across multiple analyses
#'
#' @description Identifies features that are significant across multiple meta-analysis results and merges their statistics
#' @param ... Named data frames from getSignificantFeatures(), or a named list of such data frames.
#'   Names will be used as suffixes for columns (e.g., effect_size_drug1, effect_size_drug2)
#' @return Data frame with intersecting features and their statistics from each analysis (columns suffixed by input names)
#' @export
#' @examples
#' \dontrun{
#' # Get significant features for two drugs
#' sig_drug1 <- getSignificantFeatures(meta_df1)
#' sig_drug2 <- getSignificantFeatures(meta_df2)
#'
#' # Find common significant features (using named arguments)
#' common <- getIntersectSignificantFeatures(drug1 = sig_drug1, drug2 = sig_drug2)
#'
#' # Or using a named list
#' sig_list <- list(Paclitaxel = sig_drug1, Gemcitabine = sig_drug2, Cisplatin = sig_drug3)
#' common <- getIntersectSignificantFeatures(sig_list)
#' }
getIntersectSignificantFeatures <- function(...) {
  # Get input data frames
  args <- list(...)
  
  # If first argument is a list, use it
  if (length(args) == 1 && is.list(args[[1]]) && !is.data.frame(args[[1]])) {
    df_list <- args[[1]]
  } else {
    df_list <- args
  }
  
  # Validate inputs
  if (length(df_list) < 2) {
    stop("At least 2 data frames are required for intersection")
  }
  
  if (is.null(names(df_list)) || any(names(df_list) == "")) {
    stop("All input data frames must be named. Use named arguments or a named list.")
  }
  
  # Validate each data frame has 'name' column
  for (i in seq_along(df_list)) {
    if (!is.data.frame(df_list[[i]])) {
      stop(paste0("Input '", names(df_list)[i], "' is not a data frame"))
    }
    if (!"name" %in% colnames(df_list[[i]])) {
      stop(paste0("Data frame '", names(df_list)[i], "' must contain a 'name' column"))
    }
  }
  
  # Find intersection of feature names
  common_names <- Reduce(intersect, lapply(df_list, function(df) df$name))
  
  if (length(common_names) == 0) {
    warning("No common features found across all inputs")
    return(data.frame(name = character(0)))
  }
  
  # Filter for consistent direction across all inputs
  if (all(sapply(df_list, function(df) "direction" %in% colnames(df)))) {
    consistent_names <- common_names[sapply(common_names, function(fname) {
      directions <- sapply(df_list, function(df) df$direction[df$name == fname])
      length(unique(directions)) == 1
    })]
    
    if (length(consistent_names) < length(common_names)) {
      n_filtered <- length(common_names) - length(consistent_names)
      message(sprintf("Filtered out %d feature(s) with inconsistent directions", n_filtered))
    }
    common_names <- consistent_names
    
    if (length(common_names) == 0) {
      warning("No features with consistent directions found across all inputs")
      return(data.frame(name = character(0)))
    }
  }
  
  # Build result data frame
  result <- data.frame(name = common_names, stringsAsFactors = FALSE)
  
  # Merge each data frame with suffixes
  for (i in seq_along(df_list)) {
    df_name <- names(df_list)[i]
    df <- df_list[[i]]
    
    # Filter to common names only
    df_filtered <- df[df$name %in% common_names, ]
    
    # Add suffix to all columns except 'name'
    colnames_to_rename <- setdiff(colnames(df_filtered), "name")
    for (col in colnames_to_rename) {
      colnames(df_filtered)[colnames(df_filtered) == col] <- paste0(col, "_", df_name)
    }
    
    # Merge with result
    result <- merge(result, df_filtered, by = "name", all.x = TRUE)
  }
  
  # Sort by name
  result <- result[order(result$name), ]
  rownames(result) <- NULL
  
  message(sprintf("Found %d common significant features across %d analyses", 
                  nrow(result), length(df_list)))
  
  return(result)
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
#' 5. Data loading strategy: discrete features use preload mode (required), continuous features can use preload or on-demand
#' 6. Parallel/sequential processing of feature pairs
#' 7. Statistical testing based on feature types
#' 8. Result aggregation and reporting
#'
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature1_type Type of the reference feature (e.g., "drug", "mRNA")
#' @param feature1_name Name of the reference feature
#' @param feature2_type Type of features to test against (e.g., "mRNA", "mutation_gene")
#' @param feature2_name Optional vector of specific feature names to test (default: NULL for all features)
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param cores Number of CPU cores to use for parallel processing (uses future + furrr with task chunking for efficiency)
#' @param show_progress Logical, whether to show progress updates (default: TRUE). 
#'   Progress updates include: current count, percentage, elapsed time, and estimated time remaining.
#'   Serial mode (cores=1): uses built-in progress tracker. Parallel mode (cores>1): uses progressr package if available.
#' @param progress_callback Optional custom callback function for progress updates (serial mode only). 
#'   If NULL (default), uses built-in progress tracking. Function signature: function(done, total, elapsed_seconds)
#' @param test_top_n Integer, number of top features to test (default: NULL for all features)
#' @param preloaded Logical or NULL, whether to preload all feature2 data at once (default: NULL for auto).
#'   Note: Discrete features (mutation_gene, mutation_site, fusion) ONLY support preload mode.
#'   For continuous features:
#'   NULL (auto): Preload when feature2 count > 1000, otherwise load on-demand
#'   TRUE: Always preload all feature2 data (faster for large datasets, uses more memory)
#'   FALSE: Load feature2 data on-demand (slower but memory efficient)
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
#'
#' # Force preload for faster processing (recommended for >1000 continuous features)
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA", 
#'                                         preloaded = TRUE, cores = 8)
#'
#' # Force on-demand loading for memory efficiency (continuous features only)
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA", 
#'                                         preloaded = FALSE)
#'
#' # Auto mode (default): preloads when >1000 continuous features
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA")
#'
#' # Discrete features always use preload mode (required)
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mutation_gene", 
#'                                         cores = 8)
#'
#' # With progress tracking (serial mode, cores = 1)
#' results <- batchFindSignificantFeatures(gCSI, "drug", "Paclitaxel", "mRNA",
#'                                         cores = 1, show_progress = TRUE)
#'
#' # Parallel with progressr (requires: install.packages("progressr"))
#' library(progressr)
#' handlers(global = TRUE)
#' handlers("progress")  # or "txtprogressbar"
#' results <- batchFindSignificantFeatures(gCSI, "drug", "Paclitaxel", "mRNA",
#'                                         cores = 8, show_progress = TRUE)
#'
#' # Custom progress callback (serial mode only)
#' my_progress <- function(done, total, elapsed) {
#'   message(sprintf("Processed %d/%d features in %.1f seconds", done, total, elapsed))
#' }
#' results <- batchFindSignificantFeatures(gCSI, "drug", "Paclitaxel", "mRNA",
#'                                         progress_callback = my_progress, cores = 1)
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
                                      show_progress = TRUE,
                                      progress_callback = NULL,
                                      test_top_n = NULL,
                                      preloaded = NULL,
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
                    parallel::detectCores() - 1,
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

  # Validate preloaded parameter
  if (!is.null(preloaded) && !is.logical(preloaded)) {
    stop("preloaded must be NULL, TRUE, or FALSE")
  }

  # Validate show_progress parameter
  if (!is.logical(show_progress) || length(show_progress) != 1) {
    stop("show_progress must be TRUE or FALSE")
  }

  # Setup progress tracking
  # Serial mode: use custom callback or default tracker
  if (cores == 1 && is.null(progress_callback)) {
    progress_callback <- createDefaultProgressCallback(
      show_progress = show_progress,
      update_interval = 10  # Update every 10 seconds
    )
  } else if (cores > 1 && !is.null(progress_callback)) {
    message("Note: Custom progress_callback only works in serial mode")
    progress_callback <- NULL
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
    selected_feas1 <- suppressWarnings(suppressMessages(
      loadFeatureData(dromaset_object, feature1_type, feature1_name,
                      data_type, tumor_type, overlap_only, is_continuous1,
                      return_all_samples = !is_continuous1)
    ))
  } else {
    selected_feas1 <- loadFeatureData(dromaset_object, feature1_type, feature1_name,
                                      data_type, tumor_type, overlap_only, is_continuous1,
                                      return_all_samples = !is_continuous1)
  }

  # Pre-load sample metadata for discrete vs discrete analysis
  samples_search <- NULL
  if (!is_continuous1 && !is_continuous2) {
    if (!verbose) {
      samples_search <- suppressWarnings(suppressMessages(
        getSampleMetadata(dromaset_object, feature1_type, feature2_type)
      ))
    } else {
      samples_search <- getSampleMetadata(dromaset_object, feature1_type, feature2_type)
    }
  }

  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop("No data available for the selected feature 1.")
  }

  # Filter feature1 data to ensure minimum sample size
  if (!verbose) {
    selected_feas1 <- suppressWarnings(suppressMessages(
      filterFeatureData(selected_feas1, min_samples = 3, is_discrete_with_all = !is_continuous1)
    ))
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
    feature2_list <- suppressWarnings(suppressMessages(
      getAndValidateFeatureList(
        dromaset_object = dromaset_object,
        feature_type = feature2_type,
        feature_names = feature2_name,
        data_type = data_type,
        tumor_type = tumor_type
      )
    ))
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

  # Determine whether to preload feature2 data based on preloaded parameter
  # Note: Discrete features ONLY support preload mode
  use_preload <- FALSE
  if (!is_continuous2) {
    # Discrete features always use preload mode
    use_preload <- TRUE
    if (!is.null(preloaded) && preloaded == FALSE) {
      message("Note: Discrete features only support preload mode. Switching to preload=TRUE.")
    }
  } else if (is.null(preloaded)) {
    # Auto mode for continuous: preload if more than 1000 features
    use_preload <- length(feature2_list) > 1000
  } else {
    # User-specified mode for continuous
    use_preload <- preloaded
  }

  # Preload all feature2 data if use_preload is TRUE
  preloaded_feas2 <- NULL
  if (use_preload) {
    message("Preloading all feature2 data (", length(feature2_list), " features)...")
    if (!verbose) {
      preloaded_feas2 <- suppressWarnings(suppressMessages(
        loadFeatureData(dromaset_object, feature2_type, feature2_list,
                       data_type, tumor_type, overlap_only, is_continuous2,
                       return_all_samples = !is_continuous2)
      ))
    } else {
      preloaded_feas2 <- loadFeatureData(dromaset_object, feature2_type, feature2_list,
                                        data_type, tumor_type, overlap_only, is_continuous2,
                                        return_all_samples = !is_continuous2)
    }
    
    # [OPTIMIZATION 3] Pre-filter feature2 data to avoid repeated filtering in worker
    if (!is.null(preloaded_feas2)) {
      preloaded_feas2 <- filterFeatureData(preloaded_feas2, min_samples = 3, 
                                          is_discrete_with_all = !is_continuous2)
    }
    message("Preloading completed. Starting batch analysis...")
  } else {
    message("Starting batch analysis with on-demand data loading...")
  }

  # [OPTIMIZATION 1] Pre-compute sample intersection cache for continuous vs continuous
  sample_intersection_cache <- NULL
  if (is_continuous1 && is_continuous2 && !is.null(preloaded_feas2)) {
    sample_intersection_cache <- list()
    for (i in seq_along(selected_feas1)) {
      for (j in seq_along(preloaded_feas2)) {
        cache_key <- paste(names(selected_feas1)[i], names(preloaded_feas2)[j], sep = "||")
        sample_intersection_cache[[cache_key]] <- intersect(
          names(selected_feas1[[i]]), 
          colnames(preloaded_feas2[[j]])
        )
      }
    }
  }

  # [OPTIMIZATION 2] Create worker function once, outside hot loop
  base_worker <- function(x) {
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
      preloaded_data = preloaded_feas2
    )
  }
  
  # Wrap with suppress only if needed
  worker_function <- if (!verbose) {
    function(x) {
      result <- suppressWarnings(suppressMessages(base_worker(x)))
      if (!is.null(progress_callback)) {
        progress_callback(x, length(feature2_list),
                         as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      }
      result
    }
  } else {
    function(x) {
      result <- base_worker(x)
      if (!is.null(progress_callback)) {
        progress_callback(x, length(feature2_list),
                         as.numeric(difftime(Sys.time(), start_time, units = "secs")))
      }
      result
    }
  }
  message("Please be patient, it may take long time to run.")
  
  # Use parallel processing if cores > 1
  if (cores > 1) {
    # Warn if using on-demand loading without preload
    if (is.null(preloaded_feas2)) {
      message("Note: Parallel processing without preloaded data may have database contention.")
      message("Consider using preloaded=TRUE for better performance with cores>1")
    }
    
    # Setup future plan: use multicore (memory-efficient) on Unix/Mac, multisession on Windows
    # Increase memory limit for both backends
    old_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)  # 2GB

    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = cores)
      message(sprintf("Using fork-based parallelization with %d cores (memory-efficient)", cores))
    } else {
      future::plan(future::multisession, workers = cores)
      message(sprintf("Using multisession parallelization with %d cores (Windows)", cores))
    }

    # Unified cleanup to ensure proper execution order
    on.exit({
      future::plan(future::sequential)
      options(future.globals.maxSize = old_size)
    }, add = TRUE)
    
    # [OPTIMIZATION 4] Adaptive task chunking based on dataset size
    n_features <- length(feature2_list)

    # Warn if over-parallelization (more cores than features/100)
    if (cores > n_features / 100) {
      recommended_cores <- max(1, floor(n_features / 1000))
      message(sprintf("Note: Using %d cores for %d features may cause overhead. Consider cores=%d for optimal performance.",
                     cores, n_features, recommended_cores))
    }

    chunk_size <- if (n_features < 100) {
      max(10, ceiling(n_features / cores))  # Small: minimal chunking
    } else if (n_features < 1000) {
      ceiling(n_features / (cores * 4))  # Medium: 4 tasks/core
    } else {
      # Large: 8 tasks/core with minimum chunk_size to avoid excessive scheduling overhead
      max(50, min(200, ceiling(n_features / (cores * 8))))
    }
    chunks <- split(seq_along(feature2_list), 
                   ceiling(seq_along(feature2_list) / chunk_size))
    
    # Process chunks in parallel with progressr support
    if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = n_features)
        chunk_results <- furrr::future_map(chunks, function(indices) {
          lapply(indices, function(idx) {
            result <- worker_function(idx)
            p()  # Update progress
            result
          })
        }, .options = furrr::furrr_options(seed = TRUE))
      })
    } else {
      if (show_progress) {
        message("Note: Install 'progressr' for progress tracking: install.packages('progressr')")
      }
      chunk_results <- furrr::future_map(chunks, function(indices) {
        lapply(indices, worker_function)
      }, .options = furrr::furrr_options(seed = TRUE))
    }
    
    # Flatten results
    cal_re_list <- unlist(chunk_results, recursive = FALSE)
  } else {
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
  
  # Report filtered features if any
  n_filtered <- length(feature2_list) - nrow(cal_re_df)
  if (n_filtered > 0) {
    message(sprintf("Note: %d feature(s) were filtered out due to insufficient data or failed processing.", n_filtered))
  }
  
  message(sprintf("Found %d significant associations out of %d valid features.",
                  sum(cal_re_df$q_value < 0.01), nrow(cal_re_df)))

  return(cal_re_df)
}

