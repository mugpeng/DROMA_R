# Data Processing Functions ----

#' Process Drug Sensitivity Data using DromaSet objects
#'
#' @description Creates a combined dataframe with raw and normalized drug sensitivity values from DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE).
#'   TRUE: Use only sample types present in all projects (recommended for meta-analysis)
#'   FALSE: Use all available samples from each project (may increase power but introduce bias)
#' @return A dataframe with combined raw and normalized drug sensitivity values
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' drug_data <- processDrugData(gCSI, "Paclitaxel")
#'
#' # Using MultiDromaSet with overlapping samples (recommended for meta-analysis)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))
#' drug_data <- processDrugData(multi_set, "Paclitaxel", overlap_only = FALSE)
#'
#' # Using MultiDromaSet with all samples (maximum sample size)
#' drug_data_all <- processDrugData(multi_set, "Paclitaxel", overlap_only = FALSE)
#' }
processDrugData <- function(dromaset_object, select_drugs, data_type = "all", tumor_type = "all", overlap_only = FALSE) {
  if (select_drugs == "") {
    stop("Please select a drug.")
  }

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }

  # Load drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet - Load normalized data
    drug_data <- loadTreatmentResponse(dromaset_object,
                                      select_drugs = select_drugs,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = TRUE)

    # Convert matrix to list format for compatibility
    if (is.matrix(drug_data) && select_drugs %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[select_drugs, ])
      names(drug_vector) <- colnames(drug_data)
      drug_data_list <- list()
      drug_data_list[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    } else {
      stop("Drug not found in treatment response data")
    }

    # Load raw data (non-normalized)
    raw_drug_data <- loadTreatmentResponse(dromaset_object,
                                          select_drugs = select_drugs,
                                          data_type = data_type,
                                          tumor_type = tumor_type,
                                          zscore = FALSE,
                                          return_data = TRUE)

    # Convert raw matrix to list format for compatibility
    if (is.matrix(raw_drug_data) && select_drugs %in% rownames(raw_drug_data)) {
      raw_drug_vector <- as.numeric(raw_drug_data[select_drugs, ])
      names(raw_drug_vector) <- colnames(raw_drug_data)
      raw_data_list <- list()
      raw_data_list[[dromaset_object@name]] <- raw_drug_vector[!is.na(raw_drug_vector)]
    } else {
      stop("Drug not found in raw treatment response data")
    }

  } else {
    # MultiDromaSet - Load normalized data
    drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                  select_drugs = select_drugs,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  zscore = TRUE)

    # Extract specific drug from each project
    drug_data_list <- lapply(drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    drug_data_list <- drug_data_list[!sapply(drug_data_list, is.null)]

    # Load raw data (non-normalized)
    raw_drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                      select_drugs = select_drugs,
                                                      overlap_only = overlap_only,
                                                      data_type = data_type,
                                                      tumor_type = tumor_type,
                                                      zscore = FALSE)

    # Extract specific drug from each project for raw data
    raw_data_list <- lapply(raw_drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    raw_data_list <- raw_data_list[!sapply(raw_data_list, is.null)]
  }

  # Get all study names (using union of keys from both lists)
  all_studies <- union(names(drug_data_list), names(raw_data_list))

  # Create an empty list to hold combined data for each study
  combined_data <- list()

  # For each study, process and combine the data
  for (study in all_studies) {
    if (study %in% names(drug_data_list) && study %in% names(raw_data_list)) {
      # Get values for this study
      drug_values <- drug_data_list[[study]]
      raw_values <- raw_data_list[[study]]

      # Get common sample IDs
      common_samples <- intersect(names(drug_values), names(raw_values))

      if (length(common_samples) > 0) {
        # Create dataframe with both values
        study_df <- data.frame(
          SampleID = common_samples,
          zscore_value = as.numeric(drug_values[common_samples]),
          raw_value = as.numeric(raw_values[common_samples]),
          ProjectID = study,
          stringsAsFactors = FALSE
        )

        combined_data[[study]] <- study_df
      }
    }
  }

  # Combine all studies into a single dataframe
  if (length(combined_data) > 0) {
    result_df <- bind_rows(combined_data)
    result_df <- na.omit(result_df)

    return(result_df)
  } else {
    return(NULL)
  }
}

#' Annotate Drug Sensitivity Data
#'
#' @description Adds sample annotations to drug sensitivity data
#' @details This function merges drug sensitivity data with sample annotations.
#'   It can obtain sample annotations from three sources, in order of preference:
#'   1. From the provided `sample_annotations` parameter
#'   2. From the global environment variable `sample_anno`
#'   3. From the SQLite database specified by `db_path`
#'
#' @param drug_data Dataframe containing drug sensitivity data from processDrugData
#' @param sample_annotations Dataframe containing sample annotations (default: uses global sample_anno)
#' @param db_path Optional path to SQLite database for loading sample annotations if not provided and not in global environment
#' @return A dataframe with drug sensitivity data and sample annotations
#' @export
#' @examples
#' \dontrun{
#' # Using existing sample_anno in global environment
#' drug_data <- processDrugData(gCSI, "Paclitaxel")
#' annotated_data <- annotateDrugData(drug_data)
#'
#' # Using provided sample annotations
#' my_annotations <- data.frame(SampleID = c("sample1", "sample2"),
#'                             TumorType = c("BRCA", "LUAD"))
#' annotated_data <- annotateDrugData(drug_data, sample_annotations = my_annotations)
#'
#' # Loading directly from database
#' annotated_data <- annotateDrugData(drug_data, db_path = "path/to/droma.sqlite")
#' }
annotateDrugData <- function(drug_data, sample_annotations = NULL, db_path = NULL) {
  if (is.null(drug_data)) {
    return(NULL)
  }

  # Use provided sample annotations or global sample_anno
  annotations <- sample_annotations
  if (is.null(annotations)) {
    if (!is.null(db_path) && file.exists(db_path)) {
      # Try to load sample annotations from the database
      tryCatch({
        # Load required package
        if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) {
          warning("DBI and/or RSQLite packages not available. Cannot load sample annotations from database.")
          return(drug_data)
        }

        # Connect to database
        con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
        on.exit(DBI::dbDisconnect(con), add = TRUE)

        # Query sample annotations
        annotations <- DBI::dbGetQuery(con, "SELECT * FROM sample_anno")
      }, error = function(e) {
        warning(paste("Failed to load sample annotations from database:", e$message,
                      "\nReturning non-annotated data."))
        return(drug_data)
      })
    } else {
      warning("sample_anno object not found and no valid db_path provided. Returning non-annotated data.")
      return(drug_data)
    }
  }

  # Merge drug data with annotations
  merged_df <- left_join(drug_data,
                         annotations,
                         by = c("SampleID", "ProjectID"))

  # Clean up
  if ("ProjectRawName" %in% colnames(merged_df)) {
    merged_df$ProjectRawName <- NULL
  }

  merged_df <- unique(merged_df)
  return(merged_df)
}

#' Format Drug Data Table
#'
#' @description Creates a formatted datatable for drug sensitivity data
#' @param drug_data Dataframe containing drug sensitivity data
#' @param caption Caption text for the table
#' @return A formatted DT::datatable object
#' @export
formatDrugTable <- function(drug_data, caption = "Drug sensitivity data - showing both raw and Z-score normalized values") {
  if (is.null(drug_data) || nrow(drug_data) == 0) {
    return(NULL)
  }

  DT::datatable(
    drug_data,
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
      htmltools::strong(caption)
    ),
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons',
    rownames = FALSE,
    filter = 'top'
  ) %>%
  # Format numeric columns to 3 decimal places
  DT::formatRound(columns = c('zscore_value', 'raw_value'), digits = 3)
}

#' Get Drug Sensitivity Data using DromaSet objects
#'
#' @description Wrapper function that processes and returns drug sensitivity data from DromaSet or MultiDromaSet objects
#' @details This function combines the functionality of `processDrugData()` and `annotateDrugData()`.
#'   When `include_annotations = TRUE`, it will add sample annotations to the drug data.
#'   Sample annotations can be provided directly via the `sample_annotations` parameter,
#'   or loaded from the global `sample_anno` variable, or retrieved from the SQLite
#'   database specified by `db_path`.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE).
#'   TRUE: Use only sample types present in all projects (recommended for meta-analysis)
#'   FALSE: Use all available samples from each project (may increase power but introduce bias)
#' @param include_annotations Logical indicating whether to include sample annotations
#' @param sample_annotations Optional dataframe containing sample annotations
#' @param db_path Optional path to SQLite database for loading sample annotations if not provided and not in global environment
#' @return A dataframe with drug sensitivity data
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel")
#'
#' # Using MultiDromaSet with overlapping samples (recommended)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))
#' drug_data <- getDrugSensitivityData(multi_set, "Paclitaxel",
#'                                    include_annotations = TRUE)
#'
#' # Using MultiDromaSet with all available samples
#' drug_data_all <- getDrugSensitivityData(multi_set, "Paclitaxel",
#'                                        overlap_only = FALSE,
#'                                        include_annotations = TRUE)
#'
#' # Using database path to load sample annotations
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel",
#'                                    include_annotations = TRUE,
#'                                    db_path = "path/to/droma.sqlite")
#'
#' # Using custom sample annotations
#' my_annotations <- data.frame(SampleID = c("sample1", "sample2"),
#'                             TumorType = c("BRCA", "LUAD"))
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel",
#'                                    include_annotations = TRUE,
#'                                    sample_annotations = my_annotations)
#' }
getDrugSensitivityData <- function(dromaset_object,
                                   select_drugs,
                                   data_type = "all",
                                   tumor_type = "all",
                                   overlap_only = FALSE,
                                   include_annotations = TRUE,
                                   sample_annotations = NULL,
                                   db_path = NULL) {
  # Process drug data
  drug_data <- processDrugData(dromaset_object, select_drugs, data_type, tumor_type, overlap_only)

  # Add annotations if requested
  if (include_annotations) {
    drug_data <- annotateDrugData(drug_data, sample_annotations, db_path)
  }

  return(drug_data)
}

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