# meta calculation ----
#' Calculate meta-analysis for continuous vs continuous features
#' @param selected_pair List of paired data
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcConCon <- function(selected_pair){
  if(length(selected_pair) < 1) return(NULL)
  # test pairs one by one
  cal_list <- lapply(1:length(selected_pair), function(y){
    fea1_sel <- selected_pair[[y]][[1]]
    fea2_sel <- selected_pair[[y]][[2]]
    # Check for minimum length
    if(length(fea1_sel) < 3 || length(fea2_sel) < 3) return(NULL)
    options(warn = -1)
    cor_re <- tryCatch(
      cor.test(fea1_sel, fea2_sel,
               method = "spearman"),
      error = function(x){return(NULL)}
    )
    if(is.null(cor_re)) return(NULL)
    data.frame(
      p = cor_re$p.value,
      effect = cor_re$estimate,
      N = length(fea2_sel)
    )
  })
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 1) return(NULL)
  cal_re <- do.call(rbind, cal_list)
  cal_re$se <- sqrt((1 - cal_re$effect^2) / (cal_re$N - 2))
  cal_re$z <- 0.5 * log((1 + cal_re$effect) / (1 - cal_re$effect))  # Fisher's z
  cal_re$se_z <- 1 / sqrt(cal_re$N - 3)

  cal_meta_re <- tryCatch(
    suppressWarnings({metagen(TE = z, seTE = se_z, data = cal_re, sm = "Z",
                              control = list(maxiter = 2000,
                                             stepadj = 0.1,
                                             threshold = 0.000001)
    )}),
    error = function(x){return(NULL)}
  )
  cal_meta_re
}

#' Calculate meta-analysis for continuous vs discrete features
#' @param selected_pair List of paired data
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcConDis <- function(selected_pair){
  options(warn = -1)
  if(length(selected_pair) < 1) return(NULL)
  cal_list <- lapply(1:length(selected_pair), function(y){
    yes_drugs <- selected_pair[[y]][[1]]
    no_drugs <- selected_pair[[y]][[2]]

    # Check for minimum length
    if(length(yes_drugs) < 3 || length(no_drugs) < 3) return(NULL)
    wilcox_re <- tryCatch(
      wilcox.test(no_drugs, yes_drugs),
      error = function(x){return(NULL)}
    )
    if(is.null(wilcox_re)) return(NULL)

    cliff_delta <- tryCatch(
      cliff.delta(no_drugs, yes_drugs),
      error = function(x){return(NULL)}
    )
    if(is.null(cliff_delta)) return(NULL)

    data.frame(
      p = wilcox_re$p.value,
      effect = cliff_delta$estimate,
      N = length(yes_drugs) + length(no_drugs),
      n1 = length(yes_drugs),
      n2 = length(no_drugs)
    )
  })
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 1) return(NULL)
  cal_re <- do.call(rbind, cal_list)
  # Calculate standard error for Cliff's Delta
  cal_re$se <- sqrt((1 - cal_re$effect^2) * (cal_re$n1 + cal_re$n2 + 1) /
                      (12 * cal_re$n1 * cal_re$n2))
  cal_meta_re <- tryCatch(
    suppressWarnings({meta_result <- metagen(TE = effect,
                                             seTE = se,
                                             data = cal_re,
                                             control = list(maxiter = 2000,
                                                            stepadj = 0.1,
                                                            threshold = 0.000001),
                                             sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
    )
    }),
    error = function(x){return(NULL)}
  )
  cal_meta_re
}

#' Calculate meta-analysis for discrete vs discrete features
#' @param selected_pair List of paired data
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcDisDis <- function(selected_pair) {
  # Check if we have enough pairs for meta-analysis
  if(length(selected_pair) < 1) return(NULL)

  # Calculate statistics for each pair
  cal_list <- lapply(1:length(selected_pair), function(y) {
    cont_table <- selected_pair[[y]]$cont_table

    # Skip if any cell has too few observations (e.g., < 3)
    if(any(cont_table < 3)) return(NULL)

    # Calculate odds ratio and its standard error
    tryCatch({
      # Extract values from contingency table
      a <- cont_table[1,1] # yes-yes
      b <- cont_table[1,2] # yes-no
      c <- cont_table[2,1] # no-yes
      d <- cont_table[2,2] # no-no

      # Calculate log odds ratio and its standard error
      log_or <- log((a * d)/(b * c))
      se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)

      # Calculate Fisher's exact test p-value
      fisher_test <- fisher.test(cont_table)

      data.frame(
        log_or = log_or,
        se = se_log_or,
        p = fisher_test$p.value,
        N = sum(cont_table)
      )
    }, error = function(x) NULL)
  })

  # Remove NULL results and combine
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if(length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)

  # Perform meta-analysis using random effects model
  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = log_or,
              seTE = se,
              data = cal_re,
              sm = "OR", # Specify odds ratio as summary measure
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )
    }),
    error = function(x) NULL
  )

  cal_meta_re
}

# Utility Functions ----

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

# Pairing Functions ----

#' Pair continuous with continuous data
#'
#' @description Creates paired datasets of two continuous features
#' @param myOmics First feature data list
#' @param myDrugs Second feature data list
#' @return List of paired data
#' @export
pairContinuousFeatures <- function(myOmics, myDrugs){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      omic_sel <- na.omit(omic_sel); drug_sel <- na.omit(drug_sel)
      if(length(na.omit(omic_sel)) == 0 | length(na.omit(drug_sel)) == 0){ return(NULL) }
      intersected_cells <- intersect(names(omic_sel), names(drug_sel))
      if(length(intersected_cells) < 3){ return(NULL) }
      omic_sel <- omic_sel[match(intersected_cells, names(omic_sel))]
      drug_sel <- drug_sel[match(intersected_cells, names(drug_sel))]
      list(omic_sel, drug_sel)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  pair_list2
}

#' Pair discrete with continuous data
#'
#' @description Creates paired datasets of discrete and continuous features
#' @param myOmics Discrete feature data list (samples with feature present)
#' @param myDrugs Continuous feature data list (numeric values)
#' @return List of paired data with yes/no groups
#' @export
pairDiscreteFeatures <- function(myOmics, myDrugs){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      yes_drugs <- na.omit(drug_sel[names(drug_sel) %in% omic_sel])
      no_drugs <- na.omit(drug_sel[!names(drug_sel) %in% omic_sel])
      if(length(yes_drugs) < 3 | length(no_drugs) < 3){
        return(NULL)
      }
      list(yes = yes_drugs,
           no = no_drugs)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  return(pair_list2)
}

#' Pair discrete with discrete data
#'
#' @description Creates paired datasets of two discrete features with contingency tables
#' @param my_feas1 First discrete feature data list (samples with feature present)
#' @param my_feas2 Second discrete feature data list (samples with feature present)
#' @param feature1_type Type of first feature
#' @param feature2_type Type of second feature
#' @param samples_search Reference data for all cells
#' @return List of paired data with contingency tables
#' @export
pairDiscreteDiscrete <- function(my_feas1, my_feas2,
                                feature1_type, feature2_type,
                                samples_search) {
  # Create pairs list across all features in dataset 1
  pair_list3 <- lapply(1:length(my_feas1), function(x) {
    fea1_sel <- my_feas1[[x]]
    # For each feature in dataset 1, compare with all features in dataset 2
    pair_list <- lapply(1:length(my_feas2), function(y) {
      fea2_sel <- my_feas2[[y]]
      # Skip if either feature has no data
      if (length(fea1_sel) == 0 || length(fea2_sel) == 0) {
        return(NULL)
      }

      # Get relevant cell/sample universe
      cells_search_sel <- samples_search[samples_search$type %in% c(feature1_type, feature2_type) &
                                         samples_search$datasets %in% c(names(my_feas1)[x], names(my_feas2)[y]),]
      all_cells <- unique(cells_search_sel$cells)

      # Create 2x2 contingency table
      yes_yes <- length(intersect(fea1_sel, fea2_sel))
      yes_no <- length(fea1_sel) - yes_yes
      no_yes <- length(fea2_sel) - yes_yes
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

    names(pair_list) <- paste0(names(my_feas1)[x], "_",
                               names(my_feas2))
    pair_list
  })

  # Flatten list and remove NULL entries
  pair_list3 <- unlist(pair_list3, recursive = FALSE)
  pair_list3 <- pair_list3[!sapply(pair_list3, is.null)]

  return(pair_list3)
}

# Plot ----
#' Create a volcano plot from meta-analysis results
#'
#' @param meta_df Data frame containing meta-analysis results with columns:
#'                effect_size, p_value, and name
#' @param es_t Effect size threshold to consider significant
#' @param P_t P-value threshold to consider significant
#' @param label Whether to add labels to top points (TRUE/FALSE)
#' @param top_label_each Number of top points in each direction to label
#' @param label_size Size of text labels
#' @param point_size Size of points
#' @param point_alpha Alpha transparency of points
#' @param title Plot title (NULL for no title)
#' @param p_adj_method Method for p-value adjustment ("none", "BH", "bonferroni")
#' @param custom_colors Custom color vector for Up, NS, Down (NULL for defaults)
#' @return ggplot object with volcano plot
#' @export
plotMetaVolcano <- function(meta_df,
                            es_t = .4,
                            P_t = .001,
                            label = TRUE,
                            top_label_each = 5,
                            label_size = 5,
                            point_size = 2.5,
                            point_alpha = 0.6,
                            title = NULL,
                            p_adj_method = "none",
                            custom_colors = NULL) {

  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  if(!all(c("effect_size", "p_value", "name") %in% colnames(meta_df))) {
    stop("meta_df must contain columns: effect_size, p_value, and name")
  }

  # Handle p-value adjustment if requested
  if(p_adj_method != "none") {
    meta_df$p_value <- p.adjust(meta_df$p_value, method = p_adj_method)
  }

  # Default colors
  if(is.null(custom_colors)) {
    custom_colors <- c("Down" = "#44bce4", "NS" = "grey", "Up" = "#fc7474")
  }

  # Group the points based on thresholds
  meta_df$group <- dplyr::case_when(
    meta_df$effect_size > es_t & meta_df$p_value < P_t ~ "Up",
    meta_df$effect_size < -es_t & meta_df$p_value < P_t ~ "Down",
    TRUE ~ "NS"
  )

  # Count significant findings
  sig_counts <- table(meta_df$group)
  sig_text <- paste0(
    "Up: ", sum(meta_df$group == "Up"), ", ",
    "Down: ", sum(meta_df$group == "Down"), ", ",
    "Total: ", nrow(meta_df)
  )

  # Basic volcano plot
  p <- ggplot(data = meta_df,
              aes(x = effect_size,
                  y = -log10(p_value))) +
    geom_point(size = point_size, alpha = point_alpha,
               aes(color = group)) +
    theme_bw() +
    theme(
      legend.position = "none",
      title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, colour = "black"),
      axis.text = element_text(size = 15, color = "black"),
      legend.title = element_text(size = 15, colour = "black"),
      legend.text = element_text(size = 15),
      text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black")
    ) +
    ylab("-log10(Pvalue)") +
    xlab("Effect Size") +
    scale_color_manual(values = custom_colors) +
    geom_vline(xintercept = c(-es_t, es_t), lty = 4, col = "black", lwd = 0.5) +
    geom_hline(yintercept = -log10(P_t), lty = 4, col = "black", lwd = 0.5) +
    annotate("text", x = min(meta_df$effect_size, na.rm = TRUE) * 0.8,
             y = max(-log10(meta_df$p_value), na.rm = TRUE) * 0.9,
             label = sig_text, hjust = 0, size = 5)

  # Add title if provided
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }

  # Add labels if requested
  if(label) {
    meta_df2 <- meta_df[meta_df$group != "NS",]

    # Skip labeling if there are no significant points
    if(nrow(meta_df2) > 0) {
      # Get top points to label
      low_indices <- head(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      high_indices <- tail(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      forlabel_names <- c(meta_df2$name[low_indices], meta_df2$name[high_indices])
      forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]

      p <- p +
        geom_point(size = point_size + 0.5, shape = 1, data = forlabel_df) +
        ggrepel::geom_text_repel(
          data = forlabel_df,
          aes(label = name),
          size = label_size,
          color = "black",
          box.padding = 0.5,
          point.padding = 0.3,
          force = 5,
          max.overlaps = 20
        )
    }
  }

  p
}

# Main function----
#' Batch analysis to find significant features associated with a reference feature using DromaSet objects
#'
#' @description Performs batch analysis to identify features significantly associated with a reference feature using DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature1_type Type of the reference feature (e.g., "drug", "mRNA")
#' @param feature1_name Name of the reference feature
#' @param feature2_type Type of features to test against (e.g., "mRNA", "mutation_gene")
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param cores Number of CPU cores to use for parallel processing
#' @param progress_callback Optional callback function for progress updates
#' @param test_top_100 Logical, whether to test only top 100 features (for debugging)
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
#' results <- batchFindSignificantFeatures(multi_set, "drug", "Paclitaxel", "mRNA")
#' }
batchFindSignificantFeatures <- function(dromaset_object,
                                      feature1_type,
                                      feature1_name,
                                      feature2_type,
                                      data_type = "all",
                                      tumor_type = "all",
                                      overlap_only = FALSE,
                                      cores = 1,
                                      progress_callback = NULL,
                                      test_top_100 = FALSE
) {
  # Validate input object
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

  # Track timing for progress updates
  start_time <- Sys.time()

  # Determine feature types
  continuous_types <- c("drug", "cnv", "proteinrppa",
                        "proteinms", "meth", "mRNA")
  is_continuous1 <- feature1_type %in% continuous_types
  is_continuous2 <- feature2_type %in% continuous_types

  # Get selected specific feature1 data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    if (feature1_type %in% c("drug")) {
      feature1_data <- loadTreatmentResponseNormalized(dromaset_object,
                                           drugs = feature1_name,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE)

      if (is.matrix(feature1_data) && feature1_name %in% rownames(feature1_data)) {
        feature1_vector <- as.numeric(feature1_data[feature1_name, ])
        names(feature1_vector) <- colnames(feature1_data)
        selected_feas1 <- list()
        selected_feas1[[dromaset_object@name]] <- feature1_vector[!is.na(feature1_vector)]
      } else {
        stop("Feature1 not found in treatment response data")
      }
    } else {
      # Molecular profile data
      feature1_data <- loadMolecularProfilesNormalized(dromaset_object,
                                           molecular_type = feature1_type,
                                           data_type = data_type,
                                           tumor_type = tumor_type,
                                           return_data = TRUE)

      if (is_continuous1) {
        # Continuous data
        if (is.matrix(feature1_data) && feature1_name %in% rownames(feature1_data)) {
          feature1_vector <- as.numeric(feature1_data[feature1_name, ])
          names(feature1_vector) <- colnames(feature1_data)
          selected_feas1 <- list()
          selected_feas1[[dromaset_object@name]] <- feature1_vector[!is.na(feature1_vector)]
        } else {
          stop("Feature1 not found in molecular profile data")
        }
      } else {
        # Discrete data
        if (is.data.frame(feature1_data)) {
          if ("cells" %in% colnames(feature1_data)) {
            sample_ids <- feature1_data$cells[feature1_data$genes == feature1_name]
          } else if ("samples" %in% colnames(feature1_data)) {
            sample_ids <- feature1_data$samples[feature1_data$genes == feature1_name]
          } else {
            sample_ids <- character(0)
          }
          selected_feas1 <- list()
          selected_feas1[[dromaset_object@name]] <- sample_ids
        } else {
          stop("Feature1 not found in molecular profile data")
        }
      }
    }
  } else {
    # MultiDromaSet
    if (feature1_type %in% c("drug")) {
      feature1_data <- loadMultiProjectTreatmentResponseNormalized(dromaset_object,
                                                        drugs = feature1_name,
                                                        overlap_only = overlap_only,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type)

      selected_feas1 <- lapply(feature1_data, function(drug_matrix) {
        if (is.matrix(drug_matrix) && feature1_name %in% rownames(drug_matrix)) {
          drug_vector <- as.numeric(drug_matrix[feature1_name, ])
          names(drug_vector) <- colnames(drug_matrix)
          return(drug_vector[!is.na(drug_vector)])
        }
        return(NULL)
      })
    } else {
      # Molecular profile data
      feature1_data <- loadMultiProjectMolecularProfilesNormalized(dromaset_object,
                                                        molecular_type = feature1_type,
                                                        data_type = data_type,
                                                        tumor_type = tumor_type)

      if (is_continuous1) {
        # Continuous data
        selected_feas1 <- lapply(feature1_data, function(omics_matrix) {
          if (is.matrix(omics_matrix) && feature1_name %in% rownames(omics_matrix)) {
            omics_vector <- as.numeric(omics_matrix[feature1_name, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          }
          return(NULL)
        })
      } else {
        # Discrete data
        selected_feas1 <- lapply(feature1_data, function(omics_df) {
          if (is.data.frame(omics_df)) {
            if ("cells" %in% colnames(omics_df)) {
              sample_ids <- omics_df$cells[omics_df$genes == feature1_name]
            } else if ("samples" %in% colnames(omics_df)) {
              sample_ids <- omics_df$samples[omics_df$genes == feature1_name]
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
    selected_feas1 <- selected_feas1[!sapply(selected_feas1, is.null)]
  }

  if(is.null(selected_feas1) || length(selected_feas1) == 0) {
    stop("No data available for the selected feature 1.")
  }

  # Filter selected_feas1 to ensure each dataset has at least 3 samples
  selected_feas1 <- lapply(selected_feas1, function(dataset) {
    # Remove NA values
    dataset <- na.omit(dataset)
    # Check if the dataset has at least 3 samples after NA removal
    if (length(dataset) < 3) {
      return(NULL)
    } else {
      return(dataset)
    }
  })

  # Remove NULL entries
  selected_feas1 <- selected_feas1[!sapply(selected_feas1, is.null)]

  # Check if any data remains after filtering
  if(length(selected_feas1) == 0) {
    stop(paste0("No sufficient data available for feature '", feature1_name,
                "' of type '", feature1_type, "'. Each dataset needs at least 3 samples. ",
                "Please try with a different feature."))
  }

  # Get list of features to test from the dromaset object
  feature2_list <- NULL
  if (inherits(dromaset_object, "DromaSet")) {
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
    # MultiDromaSet - get features from first project that has the data
    for (project_name in names(dromaset_object@DromaSets)) {
      dromaset <- dromaset_object@DromaSets[[project_name]]

      if (feature2_type %in% c("drug")) {
        available_features <- availableTreatmentResponses(dromaset)
        if ("drug" %in% available_features) {
          all_drug_data <- loadTreatmentResponseNormalized(dromaset, return_data = TRUE)
          if (is.matrix(all_drug_data)) {
            feature2_list <- rownames(all_drug_data)
            break
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
            feature2_list <- rownames(all_omics_data)
            break
          } else if (is.data.frame(all_omics_data) && "genes" %in% colnames(all_omics_data)) {
            feature2_list <- unique(all_omics_data$genes)
            break
          }
        }
      }
    }
  }

  if (is.null(feature2_list)) {
    stop("Could not find available features for feature2_type: ", feature2_type)
  }

  # Apply test_top_100 filter
  if(test_top_100 && length(feature2_list) > 100){
    feature2_list <- feature2_list[1:100]
  }

  if(length(feature2_list) == 0) {
    stop("No features found for the selected feature 2 type.")
  }

  # Define the worker function
  worker_function <- function(x) {
    results <- tryCatch({
      feature2_name <- feature2_list[x]

      # Get feature2 data using the same logic as feature1
      if (inherits(dromaset_object, "DromaSet")) {
        # Single DromaSet
        if (feature2_type %in% c("drug")) {
          feature2_data <- loadTreatmentResponseNormalized(dromaset_object,
                                               drugs = feature2_name,
                                               data_type = data_type,
                                               tumor_type = tumor_type,
                                               return_data = TRUE)

          if (is.matrix(feature2_data) && feature2_name %in% rownames(feature2_data)) {
            feature2_vector <- as.numeric(feature2_data[feature2_name, ])
            names(feature2_vector) <- colnames(feature2_data)
            selected_feas2 <- list()
            selected_feas2[[dromaset_object@name]] <- feature2_vector[!is.na(feature2_vector)]
          } else {
            return(NULL)
          }
        } else {
          # Molecular profile data
          feature2_data <- loadMolecularProfilesNormalized(dromaset_object,
                                               molecular_type = feature2_type,
                                               data_type = data_type,
                                               tumor_type = tumor_type,
                                               return_data = TRUE)

          if (is_continuous2) {
            # Continuous data
            if (is.matrix(feature2_data) && feature2_name %in% rownames(feature2_data)) {
              feature2_vector <- as.numeric(feature2_data[feature2_name, ])
              names(feature2_vector) <- colnames(feature2_data)
              selected_feas2 <- list()
              selected_feas2[[dromaset_object@name]] <- feature2_vector[!is.na(feature2_vector)]
            } else {
              return(NULL)
            }
          } else {
            # Discrete data
            if (is.data.frame(feature2_data)) {
              if ("cells" %in% colnames(feature2_data)) {
                sample_ids <- feature2_data$cells[feature2_data$genes == feature2_name]
              } else if ("samples" %in% colnames(feature2_data)) {
                sample_ids <- feature2_data$samples[feature2_data$genes == feature2_name]
              } else {
                sample_ids <- character(0)
              }
              selected_feas2 <- list()
              selected_feas2[[dromaset_object@name]] <- sample_ids
            } else {
              return(NULL)
            }
          }
        }
      } else {
        # MultiDromaSet
        if (feature2_type %in% c("drug")) {
          feature2_data <- loadMultiProjectTreatmentResponseNormalized(dromaset_object,
                                                            drugs = feature2_name,
                                                            overlap_only = overlap_only,
                                                            data_type = data_type,
                                                            tumor_type = tumor_type)

          selected_feas2 <- lapply(feature2_data, function(drug_matrix) {
            if (is.matrix(drug_matrix) && feature2_name %in% rownames(drug_matrix)) {
              drug_vector <- as.numeric(drug_matrix[feature2_name, ])
              names(drug_vector) <- colnames(drug_matrix)
              return(drug_vector[!is.na(drug_vector)])
            }
            return(NULL)
          })
        } else {
          # Molecular profile data
          feature2_data <- loadMultiProjectMolecularProfilesNormalized(dromaset_object,
                                                            molecular_type = feature2_type,
                                                            data_type = data_type,
                                                            tumor_type = tumor_type)

          if (is_continuous2) {
            # Continuous data
            selected_feas2 <- lapply(feature2_data, function(omics_matrix) {
              if (is.matrix(omics_matrix) && feature2_name %in% rownames(omics_matrix)) {
                omics_vector <- as.numeric(omics_matrix[feature2_name, ])
                names(omics_vector) <- colnames(omics_matrix)
                return(omics_vector[!is.na(omics_vector)])
              }
              return(NULL)
            })
          } else {
            # Discrete data
            selected_feas2 <- lapply(feature2_data, function(omics_df) {
              if (is.data.frame(omics_df)) {
                if ("cells" %in% colnames(omics_df)) {
                  sample_ids <- omics_df$cells[omics_df$genes == feature2_name]
                } else if ("samples" %in% colnames(omics_df)) {
                  sample_ids <- omics_df$samples[omics_df$genes == feature2_name]
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
        selected_feas2 <- selected_feas2[!sapply(selected_feas2, is.null)]
      }

      if (is.null(selected_feas2) || length(selected_feas2) == 0) return(NULL)

      # Filter selected_feas2 to ensure each dataset has at least 3 samples
      selected_feas2 <- lapply(selected_feas2, function(dataset) {
        # Remove NA values
        dataset <- na.omit(dataset)
        # Check if the dataset has at least 3 samples after NA removal
        if (length(dataset) < 3) {
          return(NULL)
        } else {
          return(dataset)
        }
      })

      # Remove NULL entries
      selected_feas2 <- selected_feas2[!sapply(selected_feas2, is.null)]

      # Skip if no sufficient data remains after filtering
      if(length(selected_feas2) == 0) return(NULL)

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
        # For discrete vs discrete, we need sample metadata
        # This is a limitation - we'll skip this for now or implement a workaround
        return(NULL)
      }

      if(is.null(cal_meta_re)) return(NULL)
      results <- data.frame(
        p_value = cal_meta_re[["pval.random"]],
        effect_size = cal_meta_re[["TE.random"]]
      #   N = length(cal_meta_re[["studlab"]])
      )
    }, error = function(e) {
      # Log error but continue processing
      message(sprintf("Error processing feature %d (%s): %s", x, feature2_list[x], e$message))
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
    # Initialize snowfall
    sfInit(parallel = TRUE, cpus = cores)

    # Export required data and functions
    sfExport("dromaset_object", "selected_feas1", "feature2_list",
             "is_continuous1", "is_continuous2",
             "feature1_type", "feature2_type", "data_type", "tumor_type", "overlap_only",
             "start_time")

    # Export functions
    # sfExport("loadTreatmentResponseNormalized", "loadMolecularProfilesNormalized",
    #          "loadMultiProjectTreatmentResponseNormalized", "loadMultiProjectMolecularProfilesNormalized",
    #          "pairContinuousFeatures", "pairDiscreteFeatures",
    #          "metaCalcConCon", "metaCalcConDis")

    # Load required packages on worker nodes
    sfLibrary(meta)
    sfLibrary(metafor)
    sfLibrary(effsize)
    sfLibrary(DROMA.Set)
    sfLibrary(Droma.R)
    # Run parallel computation
    cal_re_list <- sfLapply(1:length(feature2_list), worker_function)

    # Clean up snowfall
    sfStop()
  } else {
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

