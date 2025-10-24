# Drug-Omic Pair Analysis Functions ----
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
analyzeDrugOmicPair <- function(dromaset_object, select_omics_type, select_omics,
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
  
  # Warning if zscore is FALSE but merged_enabled is TRUE
  if (!zscore && merged_enabled) {
    warning("Without z-score normalization (zscore=FALSE), merging data from different studies may not be appropriate. Consider setting merged_enabled=FALSE.")
  }

  # Load and extract drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    myDrugs <- loadTreatmentResponse(dromaset_object,
                                    drugs = select_drugs,
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
                                                drugs = select_drugs,
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
    # Single DromaSet
    if (select_omics_type %in% c("drug", "drug_raw")) {
      # This is actually another drug, use treatment response
      myOmics <- loadTreatmentResponse(dromaset_object,
                                      drugs = select_omics,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = zscore)

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
                                      zscore = zscore)

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
                                                  zscore = zscore)

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
                                                  zscore = zscore)

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
      result$plot <- plotMultipleCorrelations(individual_pairs,
                                              x_label = "Omic Feature",
                                              y_label = "Drug Response")
    }

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotCorrelation(merged_pair$feature1,
                                           merged_pair$feature2,
                                           x_label = "Omic Feature",
                                           y_label = "Drug Response",
                                           title = "Merged Dataset",
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
      result$plot <- plotMultipleGroupComparisons(individual_pairs,
                                                  y_label = "Drug Response")
    }

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotGroupComparison(merged_pair$yes, merged_pair$no,
                                                group_labels = c("Without Feature", "With Feature"),
                                                title = "Merged Dataset",
                                                y_label = "Drug Response")
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

