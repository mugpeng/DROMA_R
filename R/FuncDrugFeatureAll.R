# All Drugs Analysis Functions ----

#' Get All Drugs Sensitivity Data (Internal)
#'
#' @description Internal function to get drug sensitivity data for all drugs from DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param include_annotations Logical indicating whether to include sample annotations
#' @param sample_annotations Optional dataframe containing sample annotations
#' @param db_path Optional path to SQLite database for loading sample annotations
#' @param summarize Method to summarize drug sensitivity per project ("median", "mean", "none"). Default: "median"
#' @param summarize_by What to summarize by: "drug" (per drug per project) or "sample" (per sample per project). Default: "drug"
#' @return A dataframe with drug sensitivity data for all drugs
#' @keywords internal
getAllDrugSensitivityData <- function(dromaset_object,
                                      data_type = "all",
                                      tumor_type = "all",
                                      overlap_only = FALSE,
                                      include_annotations = TRUE,
                                      sample_annotations = NULL,
                                      db_path = NULL,
                                      summarize = "median",
                                      summarize_by = "drug") {
  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  # Validate summarize parameter
  summarize <- match.arg(summarize, c("median", "mean", "none"))
  summarize_by <- match.arg(summarize_by, c("drug", "sample"))
  
  # Load drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet - Load normalized data
    drug_data <- loadTreatmentResponse(dromaset_object,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = TRUE)
    
    # Load raw data
    raw_drug_data <- loadTreatmentResponse(dromaset_object,
                                          data_type = data_type,
                                          tumor_type = tumor_type,
                                          zscore = FALSE,
                                          return_data = TRUE)
    
    # Get drug names
    drug_names <- rownames(drug_data)
    
    # Process each drug - optimized version
    all_data_list <- lapply(drug_names, function(drug) {
      # Extract drug vectors
      drug_vector <- drug_data[drug, ]
      raw_drug_vector <- raw_drug_data[drug, ]
      
      # Find valid samples (non-NA in both)
      valid_idx <- !is.na(drug_vector) & !is.na(raw_drug_vector)
      
      if (!any(valid_idx)) return(NULL)
      
      data.frame(
        SampleID = colnames(drug_data)[valid_idx],
        zscore_value = as.numeric(drug_vector[valid_idx]),
        raw_value = as.numeric(raw_drug_vector[valid_idx]),
        ProjectID = dromaset_object@name,
        DrugName = drug,
        stringsAsFactors = FALSE
      )
    })
    
    # Remove NULL entries
    all_data_list <- all_data_list[!sapply(all_data_list, is.null)]
    
  } else {
    # MultiDromaSet - Load normalized data
    drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  zscore = TRUE)
    
    # Load raw data
    raw_drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                      overlap_only = overlap_only,
                                                      data_type = data_type,
                                                      tumor_type = tumor_type,
                                                      zscore = FALSE)
    
    # Get all drug names across projects
    drug_names <- unique(unlist(lapply(drug_data, rownames)))
    
    # Process each drug across projects - optimized version
    all_data_list <- lapply(drug_names, function(drug) {
      project_data_list <- lapply(names(drug_data), function(project) {
        drug_matrix <- drug_data[[project]]
        raw_matrix <- raw_drug_data[[project]]
        
        if (!(drug %in% rownames(drug_matrix) && drug %in% rownames(raw_matrix))) {
          return(NULL)
        }
        
        # Extract drug vectors
        drug_vector <- drug_matrix[drug, ]
        raw_drug_vector <- raw_matrix[drug, ]
        
        # Find valid samples (non-NA in both)
        valid_idx <- !is.na(drug_vector) & !is.na(raw_drug_vector)
        
        if (!any(valid_idx)) return(NULL)
        
        data.frame(
          SampleID = colnames(drug_matrix)[valid_idx],
          zscore_value = as.numeric(drug_vector[valid_idx]),
          raw_value = as.numeric(raw_drug_vector[valid_idx]),
          ProjectID = project,
          DrugName = drug,
          stringsAsFactors = FALSE
        )
      })
      
      # Remove NULL entries and combine
      project_data_list <- project_data_list[!sapply(project_data_list, is.null)]
      if (length(project_data_list) > 0) {
        do.call(rbind, project_data_list)
      } else {
        NULL
      }
    })
    
    # Remove NULL entries
    all_data_list <- all_data_list[!sapply(all_data_list, is.null)]
  }
  
  # Combine all drugs into a single dataframe
  if (length(all_data_list) == 0) {
    return(NULL)
  }
  
  combined_data <- do.call(rbind, all_data_list)
  combined_data <- na.omit(combined_data)
  
  # Apply summarization if requested
  if (summarize != "none") {
    summary_func <- switch(summarize,
                          "median" = median,
                          "mean" = mean)
    
    if (summarize_by == "drug") {
      # Summarize per drug per project
      combined_data <- aggregate(
        cbind(zscore_value, raw_value) ~ ProjectID + DrugName,
        data = combined_data,
        FUN = summary_func
      )
      # Create pseudo SampleID for summarized data
      combined_data$SampleID <- paste0(combined_data$ProjectID, "_", combined_data$DrugName, "_", summarize)
    } else {
      # Summarize per sample per project
      combined_data <- aggregate(
        cbind(zscore_value, raw_value) ~ ProjectID + SampleID,
        data = combined_data,
        FUN = summary_func
      )
      # Create pseudo DrugName for summarized data
      combined_data$DrugName <- paste0("AllDrugs_", summarize)
    }
  }
  
  # Add annotations if requested
  if (include_annotations && summarize == "none") {
    combined_data <- annotateDrugData(combined_data, sample_annotations, db_path)
  }

  return(combined_data)
}

#' Plot Continuous Variable Comparison for All Drugs
#'
#' @description Creates a scatter plot with correlation statistics combining all drugs against a continuous variable
#' @param all_drug_data Return value from getAllDrugSensitivityData function
#' @param cont_column Name of the continuous column to use for x-axis
#' @param value_column Name of the value column to use for y-axis (default: "zscore_value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with all drugs combined
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet with median summarization by drug (default)
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' all_data <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "drug")
#' plot <- plotContinuousComparisonAll(all_data, "Age")
#' 
#' # Summarize by sample instead of by drug
#' all_data_sample <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "sample")
#' 
#' # Without summarization (all samples, all drugs)
#' all_data_full <- getAllDrugSensitivityData(gCSI, summarize = "none")
#' }
plotContinuousComparisonAll <- function(all_drug_data,
                                       cont_column,
                                       value_column = "zscore_value",
                                       value_label = "Drug Sensitivity",
                                       title = NULL) {
  if (is.null(all_drug_data) || !is.data.frame(all_drug_data)) {
    stop("Invalid all_drug_data parameter. Should be return value from getAllDrugSensitivityData()")
  }

  # Check if continuous column exists
  if (!cont_column %in% colnames(all_drug_data)) {
    stop("Column '", cont_column, "' not found in data")
  }

  # Use combined data for all drugs
  combined_data <- all_drug_data

  # Remove rows with NA in the continuous column
  combined_data <- combined_data[!is.na(combined_data[[cont_column]]), ]

  if (nrow(combined_data) == 0) {
    stop("No valid data after removing missing values")
  }

  # Create plot with all drugs combined
  p <- plotContinuousComparison(
    data = combined_data,
    cont_column = cont_column,
    value_column = value_column,
    value_label = value_label,
    title = title
  )

  return(p)
}

#' Plot Continuous Variable as Groups for All Drugs
#'
#' @description Creates boxplots for binned groups of a continuous variable combining all drugs
#' @param all_drug_data Return value from getAllDrugSensitivityData function
#' @param cont_column Name of the continuous column to bin into groups
#' @param value_column Name of the value column to use for y-axis (default: "zscore_value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param num_bins Number of bins to create from the continuous variable (default: 4)
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with all drugs combined
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet with median summarization by drug (default)
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' all_data <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "drug")
#' plot <- plotContinuousGroupsAll(all_data, "Age", num_bins = 4)
#' 
#' # Summarize by sample to see average drug response per sample
#' all_data_sample <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "sample")
#' plot <- plotContinuousGroupsAll(all_data_sample, "Age", num_bins = 4)
#' }
plotContinuousGroupsAll <- function(all_drug_data,
                                   cont_column,
                                   value_column = "zscore_value",
                                   value_label = "Drug Sensitivity",
                                   num_bins = 4,
                                   title = NULL) {
  if (is.null(all_drug_data) || !is.data.frame(all_drug_data)) {
    stop("Invalid all_drug_data parameter. Should be return value from getAllDrugSensitivityData()")
  }

  # Check if continuous column exists
  if (!cont_column %in% colnames(all_drug_data)) {
    stop("Column '", cont_column, "' not found in data")
  }

  # Use combined data for all drugs
  combined_data <- all_drug_data

  # Remove rows with NA in the continuous column
  combined_data <- combined_data[!is.na(combined_data[[cont_column]]), ]

  if (nrow(combined_data) == 0) {
    stop("No valid data after removing missing values")
  }

  # Create plot with all drugs combined
  p <- plotContinuousGroups(
    data = combined_data,
    cont_column = cont_column,
    value_column = value_column,
    value_label = value_label,
    num_bins = num_bins,
    title = title
  )

  return(p)
}

#' Plot Category Comparison for All Drugs
#'
#' @description Creates boxplots comparing drug sensitivity values across categories combining all drugs
#' @param all_drug_data Return value from getAllDrugSensitivityData function
#' @param category_column Name of the categorical column to use for grouping
#' @param value_column Name of the value column to use for y-axis (default: "zscore_value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with all drugs combined
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet with median summarization by drug (default)
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' all_data <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "drug")
#' plot <- plotCategoryComparisonAll(all_data, "TumorType")
#' 
#' # Summarize by sample to compare average drug response per sample across tumor types
#' all_data_sample <- getAllDrugSensitivityData(gCSI, summarize = "median", summarize_by = "sample")
#' plot <- plotCategoryComparisonAll(all_data_sample, "TumorType")
#' }
plotCategoryComparisonAll <- function(all_drug_data,
                                     category_column,
                                     value_column = "zscore_value",
                                     value_label = "Drug Sensitivity",
                                     title = NULL) {
  if (is.null(all_drug_data) || !is.data.frame(all_drug_data)) {
    stop("Invalid all_drug_data parameter. Should be return value from getAllDrugSensitivityData()")
  }

  # Check if category column exists
  if (!category_column %in% colnames(all_drug_data)) {
    stop("Column '", category_column, "' not found in data")
  }

  # Use combined data for all drugs
  combined_data <- all_drug_data

  # Remove rows with NA in the category column
  combined_data <- combined_data[!is.na(combined_data[[category_column]]), ]

  if (nrow(combined_data) == 0) {
    stop("No valid data after removing missing values")
  }

  # Create plot with all drugs combined
  p <- plotCategoryComparison(
    data = combined_data,
    category_column = category_column,
    value_column = value_column,
    value_label = value_label,
    title = title
  )

  return(p)
}

