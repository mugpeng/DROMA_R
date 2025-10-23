# Drug-Omic Pair Analysis Functions ----
# These functions leverage the modular FuncMetaAnalysis, FuncDataPairing, and FuncVisualization modules

#' Pair continuous omics and drug data
#'
#' @description Creates paired datasets of omics and drug response data from multiple sources
#' @param myOmics List of omics data vectors
#' @param myDrugs List of drug response data vectors
#' @param merged Logical, if TRUE, creates an additional merged dataset combining all pairs
#' @return A list of paired omics-drug datasets
#' @export
pairDrugOmic <- function(myOmics, myDrugs, merged = FALSE) {
  # Convert to standard format for pairContinuousFeatures
  pairs <- pairContinuousFeatures(myOmics, myDrugs, merged = merged)

  # Convert feature1/feature2 to omic/drug format for backward compatibility
  if (!is.null(pairs)) {
    pairs <- lapply(pairs, function(pair) {
      if (!is.null(pair$feature1) && !is.null(pair$feature2)) {
        list(omic = pair$feature1, drug = pair$feature2)
      } else {
        pair
      }
    })
  }

  return(pairs)
}

#' Analyze continuous drug-omics pairs
#'
#' @description Performs meta-analysis on continuous drug-omic pairs using Spearman correlation
#' @param myPairs List of paired omics-drug datasets from pairDrugOmic
#' @return Meta-analysis results object or NULL if analysis couldn't be performed
#' @export
analyzeContinuousDrugOmic <- function(myPairs) {
  # Convert omic/drug format to feature1/feature2 format for metaCalcConCon
  converted_pairs <- lapply(myPairs, function(pair) {
    if (!is.null(pair$omic) && !is.null(pair$drug)) {
      list(feature1 = pair$omic, feature2 = pair$drug)
    } else {
      NULL
    }
  })

  # Remove NULL entries
  converted_pairs <- converted_pairs[!sapply(converted_pairs, is.null)]

  if (length(converted_pairs) == 0) return(NULL)

  # Use the centralized meta-analysis function
  metaCalcConCon(converted_pairs)
}

#' Pair discrete omics and drug data
#'
#' @description Creates paired datasets of discrete omics and drug response data
#' @param myOmics List of discrete omics data (samples with feature present)
#' @param myDrugs List of drug response data vectors
#' @param merged Logical, if TRUE, creates an additional merged dataset combining all pairs
#' @return A list of paired drug-omics datasets with yes/no groups
#' @export
pairDiscreteDrugOmic <- function(myOmics, myDrugs, merged = FALSE) {
  # Use the centralized pairing function
  pairDiscreteFeatures(myOmics, myDrugs, merged = merged)
}

#' Analyze discrete drug-omics pairs
#'
#' @description Performs meta-analysis on discrete drug-omic pairs using Wilcoxon test and Cliff's Delta
#' @param myPairs List of paired drug-omics datasets from pairDiscreteDrugOmic
#' @return Meta-analysis results object or NULL if analysis couldn't be performed
#' @export
analyzeDiscreteDrugOmic <- function(myPairs) {
  # Use the centralized meta-analysis function
  metaCalcConDis(myPairs)
}

#' Plot correlation between continuous drug and omic data
#'
#' @description Creates a scatter plot with correlation statistics for a single drug-omic pair
#' @param omic_values Vector of omic feature values
#' @param drug_values Vector of drug response values
#' @param study_name Name of the study for plot title
#' @return A ggplot2 object with scatter plot and correlation statistics
#' @export
plotContinuousDrugOmic <- function(omic_values, drug_values, study_name) {
  # Use the centralized plotting function
  plotCorrelation(omic_values, drug_values,
                  x_label = "Omic Feature",
                  y_label = "Drug Response",
                  title = study_name,
                  method = "spearman")
}

#' Plot all continuous drug-omic pairs
#'
#' @description Creates and combines plots for all continuous drug-omic pairs
#' @param pairs_list List of paired drug-omic datasets
#' @return A combined plot with all drug-omic correlations or NULL if no valid pairs
#' @export
plotAllContinuousDrugOmic <- function(pairs_list) {
  # Convert omic/drug format to feature1/feature2 format for plotting
  converted_pairs <- lapply(pairs_list, function(pair) {
    if (!is.null(pair$omic) && !is.null(pair$drug)) {
      list(feature1 = pair$omic, feature2 = pair$drug)
    } else {
      NULL
    }
  })

  # Remove NULL entries
  converted_pairs <- converted_pairs[!sapply(converted_pairs, is.null)]

  if (length(converted_pairs) == 0) return(NULL)

  # Skip merged dataset for individual plots
  if ("merged_dataset" %in% names(converted_pairs)) {
    converted_pairs <- converted_pairs[names(converted_pairs) != "merged_dataset"]
  }

  # Use the centralized plotting function
  plotMultipleCorrelations(converted_pairs,
                            x_label = "Omic Feature",
                            y_label = "Drug Response")
}

#' Plot drug response by discrete omic feature
#'
#' @description Creates a boxplot comparing drug response between samples with and without a discrete omic feature
#' @param yes_values Drug response values for samples with the feature
#' @param no_values Drug response values for samples without the feature
#' @param study_name Name of the study for plot title
#' @return A ggplot2 object with boxplot and statistical test
#' @export
plotDiscreteDrugOmic <- function(yes_values, no_values, study_name) {
  # Use the centralized plotting function
  plotGroupComparison(yes_values, no_values,
                       group_labels = c("Without Feature", "With Feature"),
                       title = study_name,
                       y_label = "Drug Response")
}

#' Plot all discrete drug-omic pairs
#'
#' @description Creates and combines plots for all discrete drug-omic pairs
#' @param pairs_list List of paired drug-omic datasets with yes/no groups
#' @return A combined plot with all discrete drug-omic comparisons or NULL if no valid pairs
#' @export
plotAllDiscreteDrugOmic <- function(pairs_list) {
  # Skip merged dataset for individual plots
  if ("merged_dataset" %in% names(pairs_list)) {
    plot_list <- pairs_list[names(pairs_list) != "merged_dataset"]
  } else {
    plot_list <- pairs_list
  }

  # Use the centralized plotting function
  plotMultipleGroupComparisons(plot_list,
                               y_label = "Drug Response")
}

#' Create a forest plot for meta-analysis results
#'
#' @description Creates a standardized forest plot for visualizing meta-analysis results
#' @param meta_obj Meta-analysis object from metagen() function
#' @param xlab Label for x-axis
#' @param show_common Logical, whether to show common effect model
#' @return A forest plot visualization
#' @export
createForestPlot <- function(meta_obj,
                             xlab = "Effect Size (95% CI)",
                             show_common = FALSE) {
  # Use the centralized function from FuncMetaAnalysis
  createForestPlot(meta_obj, xlab, show_common)
}

#' Create plot with common axis labels
#'
#' @description Creates a plot with common axis labels for multiple subplots
#' @param p A patchwork object containing multiple plots
#' @param x_title Common x-axis title
#' @param y_title Common y-axis title
#' @return A function that generates the plot when called
#' @export
createPlotWithCommonAxes <- function(p, x_title = "Common X-Axis Title",
                                     y_title = "Common Y-Axis Title") {
  # Use the centralized function from FuncVisualization
  FuncVisualization::createPlotWithCommonAxes(p, x_title, y_title)
}

# Main Analysis Function ----

#' Analyze drug-omic pair associations using DromaSet objects
#'
#' @description Main function for analyzing associations between a drug and an omic feature using DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_omics_type Type of omics data to analyze (e.g., "mRNA", "mutation_gene")
#' @param select_omics Name of the specific omics feature
#' @param select_drugs Name of the drug to analyze
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param merged_enabled Logical, whether to create a merged dataset from all studies
#' @param meta_enabled Logical, whether to perform meta-analysis
#' @return A list containing plot (individual study plots), merged_plot (merged dataset plot if merged_enabled=TRUE), meta-analysis results, and data
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' result <- analyzeDrugOmicPair(gCSI, "mRNA", "ABCB1", "Paclitaxel")
#'
#' # Using MultiDromaSet
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")
#' result <- analyzeDrugOmicPair(multi_set, "mRNA", "ABCB1", "Paclitaxel")
#' }
analyzeDrugOmicPair <- function(dromaset_object, select_omics_type, select_omics,
                                select_drugs,
                                data_type = "all", tumor_type = "all",
                                overlap_only = FALSE,
                                merged_enabled = TRUE,
                                meta_enabled = TRUE){

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }

  # Load and extract drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    myDrugs <- loadTreatmentResponse(dromaset_object,
                                    drugs = select_drugs,
                                    data_type = data_type,
                                    tumor_type = tumor_type,
                                    return_data = TRUE,
                                    zscore = TRUE)

    # Convert matrix to list format for compatibility
    if (is.matrix(myDrugs) && select_drugs %in% rownames(myDrugs)) {
      drug_vector <- as.numeric(myDrugs[select_drugs, ])
      names(drug_vector) <- colnames(myDrugs)
      myDrugs <- list()
      myDrugs[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    }

  } else {
    # MultiDromaSet
    myDrugs <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                drugs = select_drugs,
                                                overlap_only = overlap_only,
                                                data_type = data_type,
                                                tumor_type = tumor_type,
                                                zscore = TRUE)

    # Extract specific drug from each project
    myDrugs <- lapply(myDrugs, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    myDrugs <- myDrugs[!sapply(myDrugs, is.null)]
  }

  # Load and extract omics data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    if (select_omics_type %in% c("drug", "drug_raw")) {
      # This is actually another drug, use treatment response
      myOmics <- loadTreatmentResponse(dromaset_object,
                                      drugs = select_omics,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = TRUE)

      if (is.matrix(myOmics) && select_omics %in% rownames(myOmics)) {
        omics_vector <- as.numeric(myOmics[select_omics, ])
        names(omics_vector) <- colnames(myOmics)
        myOmics <- list()
        myOmics[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
      }
    } else {
      # Molecular profile data
      myOmics <- loadMolecularProfiles(dromaset_object,
                                      molecular_type = select_omics_type,
                                      features = select_omics,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = TRUE)

      # Handle different data types
      if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
        # Continuous data
        if (is.matrix(myOmics) && select_omics %in% rownames(myOmics)) {
          omics_vector <- as.numeric(myOmics[select_omics, ])
          names(omics_vector) <- colnames(myOmics)
          myOmics <- list()
          myOmics[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
        }
      } else {
        # Discrete data (mutations, fusions) - now in matrix format like continuous data
        if (is.matrix(myOmics) && select_omics %in% rownames(myOmics)) {
          # Get the row for the selected gene
          gene_row <- myOmics[select_omics, ]
          # Get sample IDs where the feature is present (value != 0)
          # When extracting a row from matrix, it becomes a named vector
          present_samples <- names(gene_row)[gene_row != 0]
          myOmics <- list()
          myOmics[[dromaset_object@name]] <- present_samples
        } else {
          myOmics <- list()
          myOmics[[dromaset_object@name]] <- character(0)
        }
      }
    }

  } else {
    # MultiDromaSet
    if (select_omics_type %in% c("drug", "drug_raw")) {
      # This is actually another drug, use treatment response
      myOmics <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                  drugs = select_omics,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  zscore = TRUE)

      # Extract specific drug from each project
      myOmics <- lapply(myOmics, function(drug_matrix) {
        if (is.matrix(drug_matrix) && select_omics %in% rownames(drug_matrix)) {
          drug_vector <- as.numeric(drug_matrix[select_omics, ])
          names(drug_vector) <- colnames(drug_matrix)
          return(drug_vector[!is.na(drug_vector)])
        }
        return(NULL)
      })
    } else {
      # Molecular profile data
      myOmics <- loadMultiProjectMolecularProfiles(dromaset_object,
                                                  molecular_type = select_omics_type,
                                                  features = select_omics,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  zscore = TRUE)

      # Handle different data types
      if (select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
        # Continuous data
        myOmics <- lapply(myOmics, function(omics_matrix) {
          if (is.matrix(omics_matrix) && select_omics %in% rownames(omics_matrix)) {
            omics_vector <- as.numeric(omics_matrix[select_omics, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          }
          return(NULL)
        })
      } else {
        # Discrete data (mutations, fusions) - now in matrix format like continuous data
        myOmics <- lapply(myOmics, function(omics_matrix) {
          if (is.matrix(omics_matrix) && select_omics %in% rownames(omics_matrix)) {
            # Get the row for the selected gene
            gene_row <- omics_matrix[select_omics, ]
            # Get sample IDs where the feature is present (value != 0)
            # When extracting a row from matrix, it becomes a named vector
            present_samples <- names(gene_row)[gene_row != 0]
            return(present_samples)
          }
          return(NULL)
        })
      }
    }

    # Remove NULL entries
    myOmics <- myOmics[!sapply(myOmics, is.null)]
  }

  # Check if we have data
  if (length(myDrugs) == 0 || length(myOmics) == 0) {
    stop("No data found for the specified drug-omics pair")
  }

  # Initialize result list
  result <- list()

  # Handle continuous omics data
  if(select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")){
    # Pair data
    myPairs <- pairDrugOmic(myOmics, myDrugs, merged = merged_enabled)

    # Separate individual pairs from merged pairs for analysis
    individual_pairs <- myPairs
    merged_pair <- NULL

    # Extract merged dataset if it exists
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      # Remove merged dataset from individual pairs for meta-analysis
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }

    # Create plots for individual studies
    result$plot <- plotAllContinuousDrugOmic(individual_pairs)

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotContinuousDrugOmic(merged_pair$omic,
                                                   merged_pair$drug,
                                                   "Merged Dataset")
    }

    # Perform meta-analysis on individual pairs only (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 0){
      meta_result <- analyzeContinuousDrugOmic(individual_pairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }

    # Store data
    result$data <- myPairs
  } else {
    # Pair data
    myPairs <- pairDiscreteDrugOmic(myOmics, myDrugs, merged = merged_enabled)

    # Separate individual pairs from merged pairs for analysis
    individual_pairs <- myPairs
    merged_pair <- NULL

    # Extract merged dataset if it exists
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      # Remove merged dataset from individual pairs for meta-analysis
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }

    # Create plots for individual studies
    result$plot <- plotAllDiscreteDrugOmic(individual_pairs)

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotDiscreteDrugOmic(merged_pair$yes, merged_pair$no, "Merged Dataset")
    }

    # Perform meta-analysis on individual pairs only (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 0){
      meta_result <- analyzeDiscreteDrugOmic(individual_pairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }

    # Store data
    result$data <- myPairs
  }

  # Return results if there's a plot or merged plot
  if (is.null(result$plot) && is.null(result$merged_plot)) {
    return(list())
  } else {
    return(result)
  }
}

