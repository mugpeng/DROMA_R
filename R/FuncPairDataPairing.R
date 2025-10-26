# Data Pairing Functions Module ----

#' Pair continuous features from two datasets
#'
#' @description Creates paired datasets of two continuous features
#' @param dataset1 First dataset list with named vectors
#' @param dataset2 Second dataset list with named vectors
#' @param merged Logical, if TRUE, creates an additional merged dataset combining all pairs
#' @param intersection_cache Optional pre-computed sample intersection cache (NULL for auto-compute)
#' @return List of paired data with feature1 and feature2 values
#' @export
pairContinuousFeatures <- function(dataset1, dataset2, merged = FALSE, intersection_cache = NULL) {
  pair_list2 <- lapply(1:length(dataset1), function(x) {
    feat1_sel <- dataset1[[x]]
    pair_list <- lapply(1:length(dataset2), function(y) {
      feat2_sel <- dataset2[[y]]
      feat1_sel <- na.omit(feat1_sel)
      feat2_sel <- na.omit(feat2_sel)

      if(length(feat1_sel) == 0 | length(feat2_sel) == 0) {
        return(NULL)
      }

      # Use cached intersection if available, otherwise compute
      if (!is.null(intersection_cache)) {
        cache_key <- paste(names(dataset1)[x], names(dataset2)[y], sep = "||")
        intersected_cells <- intersection_cache[[cache_key]]
        if (is.null(intersected_cells)) {
          intersected_cells <- intersect(names(feat1_sel), names(feat2_sel))
        }
      } else {
        intersected_cells <- intersect(names(feat1_sel), names(feat2_sel))
      }
      
      if(length(intersected_cells) < 3) {
        return(NULL)
      }

      # Vectorized operation: direct indexing instead of match()
      feat1_sel <- feat1_sel[intersected_cells]
      feat2_sel <- feat2_sel[intersected_cells]

      list(feature1 = feat1_sel,
           feature2 = feat2_sel)
    })

    names(pair_list) <- paste0(names(dataset1)[x], "_",
                               names(dataset2))
    pair_list
  })

  pair_list2 <- unlist(pair_list2, recursive = FALSE)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]

  if(length(pair_list2) < 1) {
    stop("No valid pairs found. Please check your data.")
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

#' Pair discrete with continuous data
#'
#' @description Creates paired datasets of discrete and continuous features
#' @param discrete_dataset List with discrete feature data (samples with feature present)
#' @param continuous_dataset List with continuous feature data (numeric values)
#' @param merged Logical, if TRUE, creates an additional merged dataset
#' @return List of paired data with yes/no groups for continuous values
#' @export
pairDiscreteFeatures <- function(discrete_dataset, continuous_dataset, merged = FALSE) {
  pair_list2 <- lapply(1:length(discrete_dataset), function(x) {
    discrete_sel <- discrete_dataset[[x]]
    present_samples <- discrete_sel$present
    all_profiled_samples <- discrete_sel$all
    
    pair_list <- lapply(1:length(continuous_dataset), function(y) {
      continuous_sel <- continuous_dataset[[y]]
      
      # Get intersection of continuous data and all profiled samples
      valid_samples <- intersect(names(continuous_sel), all_profiled_samples)
      valid_samples <- valid_samples[!is.na(continuous_sel[valid_samples])]
      
      # yes_values: samples in present_samples
      yes_values <- continuous_sel[names(continuous_sel) %in% present_samples & names(continuous_sel) %in% valid_samples]
      yes_values <- na.omit(yes_values)
      
      # no_values: samples in valid_samples but NOT in present_samples
      no_values <- continuous_sel[names(continuous_sel) %in% valid_samples & !names(continuous_sel) %in% present_samples]
      no_values <- na.omit(no_values)

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
    stop("No valid pairs found. Please check your data.")
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

#' Pair discrete with discrete data
#'
#' @description Creates paired datasets of two discrete features with contingency tables
#' @param discrete_dataset1 First discrete feature data list (samples with feature present)
#' @param discrete_dataset2 Second discrete feature data list (samples with feature present)
#' @param feature1_type Type of first feature
#' @param feature2_type Type of second feature
#' @param samples_search Reference data for all samples
#' @return List of paired data with contingency tables
#' @export
pairDiscreteDiscrete <- function(discrete_dataset1, discrete_dataset2,
                                feature1_type, feature2_type,
                                samples_search) {
  # Create pairs list across all features in dataset 1
  pair_list3 <- lapply(1:length(discrete_dataset1), function(x) {
    feat1_sel <- discrete_dataset1[[x]]

    # For each feature in dataset 1, compare with all features in dataset 2
    pair_list <- lapply(1:length(discrete_dataset2), function(y) {
      feat2_sel <- discrete_dataset2[[y]]

      # Skip if either feature has no data
      if (length(feat1_sel) == 0 || length(feat2_sel) == 0) {
        return(NULL)
      }

      # Get relevant cell/sample universe
      cells_search_sel <- samples_search[samples_search$type %in% c(feature1_type, feature2_type) &
                                         samples_search$datasets %in% c(names(discrete_dataset1)[x], names(discrete_dataset2)[y]),]
      all_cells <- unique(cells_search_sel$cells)

      # Create 2x2 contingency table
      yes_yes <- length(intersect(feat1_sel, feat2_sel))
      yes_no <- length(feat1_sel) - yes_yes
      no_yes <- length(feat2_sel) - yes_yes
      no_no <- length(all_cells) - (yes_yes + yes_no + no_yes)

      # Skip if any cell count is too low or negative (invalid)
      if (any(c(yes_yes, yes_no, no_yes, no_no) < 3)) {
        return(NULL)
      }

      # Create contingency table
      cont_table <- matrix(
        c(yes_yes, yes_no,
          no_yes, no_no),
        nrow = 2,
        dimnames = list(
          Feature1 = c("Yes", "No"),
          Feature2 = c("Yes", "No")
        )
      )

      list(
        "cont_table" = cont_table
      )
    })

    names(pair_list) <- paste0(names(discrete_dataset1)[x], "_",
                               names(discrete_dataset2))
    pair_list
  })

  # Flatten list and remove NULL entries
  pair_list3 <- unlist(pair_list3, recursive = FALSE)
  pair_list3 <- pair_list3[!sapply(pair_list3, is.null)]

  return(pair_list3)
}

#' Format seconds into a human-readable time string
#'
#' @description Converts seconds into a more readable time format (hours, minutes, seconds)
#' @param seconds Number of seconds
#' @return Formatted time string
#' @export
formatTime <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%d seconds", round(seconds)))
  } else if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    remaining_seconds <- round(seconds %% 60)
    return(sprintf("%d minutes %d seconds", minutes, remaining_seconds))
  } else {
    hours <- floor(seconds / 3600)
    remaining_minutes <- floor((seconds %% 3600) / 60)
    return(sprintf("%d hours %d minutes", hours, remaining_minutes))
  }
}

#' Calculate estimated time remaining based on progress
#'
#' @description Estimates remaining time for a batch process based on current progress
#' @param done Number of items processed
#' @param total Total number of items
#' @param elapsed_time Time elapsed so far in seconds
#' @return Estimated time remaining in seconds
#' @export
estimateTimeRemaining <- function(done, total, elapsed_time) {
  if (done == 0) return(Inf)
  rate <- elapsed_time / done
  remaining <- total - done
  remaining * rate
}