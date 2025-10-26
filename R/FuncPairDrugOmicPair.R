# Drug-Omic Pair Analysis Functions ----

#' Analyze drug-omic pair associations using DromaSet objects
#'
#' @description Main function for analyzing associations between a drug and an omic feature using DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param feature_type Type of omics data to analyze (e.g., "mRNA", "mutation_gene")
#' @param select_features Name of the specific omics feature
#' @param select_drugs Name of the drug to analyze
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param merged_enabled Logical, whether to create a merged dataset from all studies
#' @param meta_enabled Logical, whether to perform meta-analysis
#' @param zscore Logical, whether to apply z-score normalization to treatment response and molecular profiles (default: TRUE). If FALSE, merged_enabled should be set to FALSE to avoid combining non-normalized data from different studies.
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
analyzeDrugOmicPair <- function(dromaset_object, feature_type, select_features,
                                select_drugs,
                                data_type = "all", tumor_type = "all",
                                overlap_only = FALSE,
                                merged_enabled = TRUE,
                                meta_enabled = TRUE,
                                zscore = TRUE){

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  # Validate omics type
  supported_omics_types <- c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", 
                             "mutation_gene", "mutation_site", "fusion")
  if (!feature_type %in% supported_omics_types) {
    warning(sprintf("'%s' is not a supported omics type. Supported types are: %s", 
                    feature_type, paste(supported_omics_types, collapse = ", ")))
  }
  
  # Warning if zscore is FALSE but merged_enabled is TRUE
  if (!zscore && merged_enabled) {
    warning("Without z-score normalization (zscore=FALSE), merging data from different studies may not be appropriate. Consider setting merged_enabled=FALSE.")
  }

  # Load and extract drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    myDrugs <- loadTreatmentResponse(dromaset_object,
                                    select_drugs = select_drugs,
                                    data_type = data_type,
                                    tumor_type = tumor_type,
                                    return_data = TRUE,
                                    zscore = zscore)

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
                                                select_drugs = select_drugs,
                                                overlap_only = overlap_only,
                                                data_type = data_type,
                                                tumor_type = tumor_type,
                                                zscore = zscore)

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
    # Single DromaSet - Molecular profile data
    myOmics <- loadMolecularProfiles(dromaset_object,
                                    feature_type = feature_type,
                                    select_features = select_features,
                                    data_type = data_type,
                                    tumor_type = tumor_type,
                                    return_data = TRUE,
                                    zscore = zscore)

    # Handle different data types
    if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
      # Continuous data
      if (is.matrix(myOmics) && select_features %in% rownames(myOmics)) {
        omics_vector <- as.numeric(myOmics[select_features, ])
        names(omics_vector) <- colnames(myOmics)
        myOmics <- list()
        myOmics[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
      }
    } else {
      # Discrete data (mutations, fusions) - long dataframe format with samples and features columns
      if (is.data.frame(myOmics) && "samples" %in% colnames(myOmics) && "features" %in% colnames(myOmics)) {
        # Extract samples where the feature matches select_features
        present_samples <- myOmics$samples[myOmics$features == select_features]
        # Get all profiled samples for this omics type
        all_profiled_samples <- listDROMASamples(dromaset_object@name,
                                                 feature_type = feature_type,
                                                 data_type = data_type,
                                                 tumor_type = tumor_type)
        myOmics <- list()
        myOmics[[dromaset_object@name]] <- list(present = present_samples, all = all_profiled_samples)
      } else {
        myOmics <- list()
        myOmics[[dromaset_object@name]] <- list(present = character(0), all = character(0))
      }
    }

  } else {
    # MultiDromaSet - Molecular profile data
    myOmics <- loadMultiProjectMolecularProfiles(dromaset_object,
                                                feature_type = feature_type,
                                                select_features = select_features,
                                                overlap_only = overlap_only,
                                                data_type = data_type,
                                                tumor_type = tumor_type,
                                                zscore = zscore)

    # Handle different data types
    if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
      # Continuous data
      myOmics <- lapply(myOmics, function(omics_matrix) {
        if (is.matrix(omics_matrix) && select_features %in% rownames(omics_matrix)) {
          omics_vector <- as.numeric(omics_matrix[select_features, ])
          names(omics_vector) <- colnames(omics_matrix)
          return(omics_vector[!is.na(omics_vector)])
        }
        return(NULL)
      })
    } else {
      # Discrete data (mutations, fusions) - long dataframe format with samples and features columns
      project_names <- names(myOmics)
      myOmics <- lapply(seq_along(myOmics), function(i) {
        omics_df <- myOmics[[i]]
        projects <- project_names[i]
        if (is.data.frame(omics_df) && "samples" %in% colnames(omics_df) && "features" %in% colnames(omics_df)) {
          # Extract samples where the feature matches select_features
          present_samples <- omics_df$samples[omics_df$features == select_features]
          # Get all profiled samples for this omics type from the specific project
          all_profiled_samples <- listDROMASamples(dromaset_object@DromaSets[[projects]]@name,
                                                   feature_type = feature_type,
                                                   data_type = data_type,
                                                   tumor_type = tumor_type)
          return(list(present = present_samples, all = all_profiled_samples))
        }
        return(NULL)
      })
      names(myOmics) <- project_names
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
  if(feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")){
    # Pair data using pairContinuousFeatures
    myPairs <- pairContinuousFeatures(myOmics, myDrugs, merged = merged_enabled)

    # Separate individual pairs from merged pairs for analysis
    individual_pairs <- myPairs
    merged_pair <- NULL

    # Extract merged dataset if it exists
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      # Remove merged dataset from individual pairs for meta-analysis
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }

    # Create plots for individual studies using plotMultipleCorrelations
    if (length(individual_pairs) > 0) {
      multi_plot <- plotMultipleCorrelations(individual_pairs,
                                              x_label = paste0(select_features, " (", feature_type, ")"),
                                              y_label = "Drug Response")
      # Add common axis labels
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                              x_title = paste(feature_type, "expression"),
                                              y_title = "drug sensitivity(Area Above Curve)")
    }

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotCorrelation(merged_pair$feature1,
                                           merged_pair$feature2,
                                           x_label = paste(feature_type, "expression"),
                                           y_label = "drug sensitivity(Area Above Curve)",
                                           title = paste(feature_type, ":", select_features, "vs", select_drugs),
                                           method = "spearman")
    }

    # Perform meta-analysis on individual pairs only (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 0){
      meta_result <- metaCalcConCon(individual_pairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }

    # Store data
    result$data <- myPairs
  } else {
    # Pair data using pairDiscreteFeatures
    myPairs <- pairDiscreteFeatures(myOmics, myDrugs, merged = merged_enabled)

    # Separate individual pairs from merged pairs for analysis
    individual_pairs <- myPairs
    merged_pair <- NULL

    # Extract merged dataset if it exists
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      # Remove merged dataset from individual pairs for meta-analysis
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }

    # Create plots for individual studies using plotMultipleGroupComparisons
    if (length(individual_pairs) > 0) {
      multi_plot <- plotMultipleGroupComparisons(individual_pairs,
                                                x_label = paste0(select_features, " (", feature_type, ")"),
                                                y_label = "Drug Response")
      # Add common axis label
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                              y_title = "drug sensitivity(Area Above Curve)")
    }

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotGroupComparison(merged_pair$yes, merged_pair$no,
                                                group_labels = c(paste("Without", select_features), paste("With", select_features)),
                                                title = paste(feature_type, ":", select_features, "vs", select_drugs),
                                                y_label = "drug sensitivity(Area Above Curve)")
    }

    # Perform meta-analysis on individual pairs only (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 0){
      meta_result <- metaCalcConDis(individual_pairs)
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

