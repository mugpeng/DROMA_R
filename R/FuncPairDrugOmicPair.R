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
#' @param data_type_anno Optional character string to add as annotation in plot titles (e.g., "Cell lines"). If provided, it will be appended to titles as "(annotation)"
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
                                zscore = TRUE,
                                data_type_anno = NULL){

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
  
  # Create title suffix for data type annotation
  title_suffix <- if (!is.null(data_type_anno) && nchar(data_type_anno) > 0) {
    paste0(" (", data_type_anno, ")")
  } else {
    ""
  }
  
  # Determine if omics feature is continuous
  is_continuous_omics <- feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")

  # Load drug data using helper function
  myDrugs <- loadFeatureData(dromaset_object, "drug", select_drugs,
                            data_type = data_type, tumor_type = tumor_type,
                            overlap_only = overlap_only, is_continuous = TRUE,
                            zscore = zscore)

  # Load omics data using helper function
  # For discrete features (mutations, fusions), we need the full format with present/all samples
  myOmics <- loadFeatureData(dromaset_object, feature_type, select_features,
                            data_type = data_type, tumor_type = tumor_type,
                            overlap_only = overlap_only, is_continuous = is_continuous_omics,
                            return_all_samples = !is_continuous_omics, zscore = zscore)

  # Filter data to ensure minimum sample size
  myDrugs <- filterFeatureData(myDrugs, min_samples = 3, is_discrete_with_all = FALSE)
  myOmics <- filterFeatureData(myOmics, min_samples = 3, is_discrete_with_all = !is_continuous_omics)

  # Check if we have sufficient data after filtering
  if (is.null(myDrugs) || length(myDrugs) == 0 || is.null(myOmics) || length(myOmics) == 0) {
    stop("No sufficient data found for the specified drug-omics pair. Each dataset needs at least 3 samples.")
  }

  # Initialize result list
  result <- list()

  # Handle continuous omics data
  if(is_continuous_omics){
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

    # Create plots for individual studies
    if (length(individual_pairs) == 1) {
      # When only one pair, use single plot method (like merged_plot)
      single_pair <- individual_pairs[[1]]
      result$plot <- plotCorrelation(single_pair$feature1,
                                     single_pair$feature2,
                                     x_label = paste(feature_type, "expression"),
                                     y_label = "drug sensitivity(Area Above Curve)",
                                     title = paste0(paste(feature_type, ":", select_features, "vs", select_drugs), title_suffix),
                                     method = "spearman")
    } else if (length(individual_pairs) > 1) {
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
                                           title = paste0(paste(feature_type, ":", select_features, "vs", select_drugs), title_suffix),
                                           method = "spearman")
    }

    # Perform meta-analysis only when there are multiple pairs (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 1){
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

    # Create plots for individual studies
    if (length(individual_pairs) == 1) {
      # When only one pair, use single plot method (like merged_plot)
      single_pair <- individual_pairs[[1]]
      result$plot <- plotGroupComparison(single_pair$no, single_pair$yes,
                                        group_labels = c(paste("Without", select_features), paste("With", select_features)),
                                        title = paste0(paste(feature_type, ":", select_features, "vs", select_drugs), title_suffix),
                                        y_label = "drug sensitivity(Area Above Curve)")
    } else if (length(individual_pairs) > 1) {
      multi_plot <- plotMultipleGroupComparisons(individual_pairs,
                                                group_labels = c(paste("Without", select_features), paste("With", select_features)),
                                                x_label = paste0(select_features, " (", feature_type, ")"),
                                                y_label = "Drug Response")
      # Add common axis label
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                              y_title = "drug sensitivity(Area Above Curve)")
    }

    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotGroupComparison(merged_pair$no, merged_pair$yes,
                                                group_labels = c(paste("Without", select_features), paste("With", select_features)),
                                                title = paste0(paste(feature_type, ":", select_features, "vs", select_drugs), title_suffix),
                                                y_label = "drug sensitivity(Area Above Curve)")
    }

    # Perform meta-analysis only when there are multiple pairs (exclude merged data)
    if(meta_enabled && length(individual_pairs) > 1){
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

