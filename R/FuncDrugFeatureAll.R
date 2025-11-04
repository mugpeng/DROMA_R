# All Drugs Analysis Functions ----

#' Get All Drugs Sensitivity Data
#'
#' @description function to get drug sensitivity data for all drugs from DromaSet or MultiDromaSet objects
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
#' @export
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