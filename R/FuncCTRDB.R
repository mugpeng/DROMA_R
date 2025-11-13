# Clinical Trial Database (CTRDB) Analysis Functions ----

#' Get patient expression data from CTRDB database
#'
#' @description Retrieves expression data for a specific patient from CTRDB database.
#' Automatically detects and applies log transformation if needed based on data range.
#' @param patient_id Patient identifier (will be sanitized for database table name)
#' @param connection Database connection object
#' @param auto_log Logical, whether to automatically detect and apply log transformation (default: TRUE)
#' @param log_threshold Maximum value threshold to determine if data is already log-transformed (default: 50)
#' @return Expression matrix with genes as rows and samples as columns, or NULL if retrieval fails
#' @export
#' @examples
#' \dontrun{
#' # Connect to CTRDB database
#' con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#'
#' # Get expression data for a patient
#' expr_data <- getPatientExpressionData("Patient-001", con)
#' }
getPatientExpressionData <- function(patient_id,
                                    connection,
                                    auto_log = TRUE,
                                    log_threshold = 50) {

  # Validate inputs
  if (missing(patient_id) || is.null(patient_id) || patient_id == "") {
    stop("patient_id must be specified")
  }
  if (is.null(connection)) {
    stop("Database connection must be provided")
  }

  # Replace "-" with "_" for database table name
  table_name_db <- gsub("-", "_", patient_id)

  # Query expression matrix
  expr_data <- tryCatch({
    expr_query <- paste0("SELECT * FROM ", table_name_db)
    full_expr <- DBI::dbGetQuery(connection, expr_query)

    # Check if data was retrieved
    if (is.null(full_expr) || nrow(full_expr) == 0) {
      warning("No data found in table ", table_name_db, " for patient ", patient_id)
      return(NULL)
    }

    # Convert to matrix with feature_id as rownames
    if (!"feature_id" %in% colnames(full_expr)) {
      stop("feature_id column not found in expression data for patient ", patient_id)
    }

    rownames(full_expr) <- full_expr$feature_id
    full_expr$feature_id <- NULL

    expr_matrix <- as.matrix(full_expr)

    # Check if log transformation is needed
    if (auto_log && nrow(expr_matrix) > 0 && ncol(expr_matrix) > 0) {
      # Calculate statistics to determine if data is log-transformed
      max_value <- max(expr_matrix, na.rm = TRUE)
      min_value <- min(expr_matrix, na.rm = TRUE)
      median_value <- median(expr_matrix, na.rm = TRUE)

      # If min value is negative, skip log transformation
      if (min_value < 0) {
        message("Expression data for patient ", patient_id,
                " contains negative values (min = ", round(min_value, 2),
                "). Skipping log transformation.")
      } else if (max_value > log_threshold || median_value > 20) {
        # If max value exceeds threshold, likely not log-transformed
        # Also check if median is suspiciously high (> 20 typically indicates non-log data)
        message("Expression data for patient ", patient_id,
                " appears to be non-log-transformed (max = ", round(max_value, 2),
                ", median = ", round(median_value, 2),
                "). Applying log2(x + 1) transformation...")
        expr_matrix <- log2(expr_matrix + 1)
      } else {
        message("Expression data for patient ", patient_id,
                " appears to be log-transformed (max = ", round(max_value, 2),
                ", median = ", round(median_value, 2), ")")
      }
    }

    return(expr_matrix)

  }, error = function(e) {
    warning("Failed to retrieve expression data for patient ", patient_id, ": ", e$message)
    return(NULL)
  })

  return(expr_data)
}

#' @title Analyze Clinical Drug Response with Omics Data from CTRDB
#' @description Performs analysis of clinical drug response associations with omics features using CTRDB data
#' @param select_omics Character string specifying the omics feature name
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param meta_enabled Logical, whether to perform meta-analysis across datasets
#' @param zscore Logical, whether to apply z-score normalization to omics expression data (default: TRUE). If FALSE, merged_enabled should be set to FALSE to avoid combining non-normalized data from different patients.
#' @param merged_enabled Logical, whether to create a merged dataset from all patients
#' @return A list containing plot (individual study plots), merged_plot (merged dataset plot if merged_enabled=TRUE), meta-analysis results, and data
#' @export
#' @examples
#' \dontrun{
#' # Connect to CTRDB database first
#' con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#'
#' # Analyze clinical drug response
#' result <- analyzeClinicalDrugResponse("EGFR", "Erlotinib")
#'
#' # View individual patient plots
#' result$plot
#'
#' # View merged dataset plot
#' result$merged_plot
#'
#' # View meta-analysis results
#' result$meta
#' }
analyzeClinicalDrugResponse <- function(select_omics,
                                      select_drugs,
                                      data_type = "all",
                                      tumor_type = "all",
                                      connection = NULL,
                                      meta_enabled = TRUE,
                                      zscore = TRUE,
                                      merged_enabled = TRUE) {

  # Validate inputs
  if (missing(select_omics) || is.null(select_omics) || select_omics == "") {
    stop("select_omics must be specified")
  }

  if (missing(select_drugs) || is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
  }

  # Warning if zscore is FALSE but merged_enabled is TRUE
  if (!zscore && merged_enabled) {
    warning("Without z-score normalization (zscore=FALSE), merging data from different patients may not be appropriate. Consider setting merged_enabled=FALSE.")
  }

  # Get CTRDB connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("ctrdb_connection", envir = .GlobalEnv)) {
      stop("No CTRDB database connection found. Connect first with connectCTRDatabase()")
    }
    connection <- get("ctrdb_connection", envir = .GlobalEnv)
  }

  # Get sample annotations from DROMA
  sample_anno <- tryCatch({
    # Use DROMA.Set function to get sample annotations
    getDROMAAnnotation("sample")
  }, error = function(e) {
    stop("Failed to load sample annotations from DROMA: ", e$message)
  })

  if (is.null(sample_anno) || nrow(sample_anno) == 0) {
    stop("No sample annotations found in database")
  }

  # Check if CliUsedDrug column exists
  if (!"CliUsedDrug" %in% colnames(sample_anno)) {
    stop("CliUsedDrug column not found in sample annotations")
  }

  # Find samples with the selected drug
  drug_samples <- sample_anno[!is.na(sample_anno$CliUsedDrug) &
                             grepl(select_drugs, sample_anno$CliUsedDrug, ignore.case = TRUE), ]

  if (nrow(drug_samples) == 0) {
    stop("No samples found with drug: ", select_drugs)
  }

  # Filter for CTRDB project
  if ("ProjectID" %in% colnames(drug_samples)) {
    ctrdb_samples <- drug_samples[drug_samples$ProjectID == "CTRDB", ]
  } else {
    stop("ProjectID column not found in sample annotations")
  }

  if (nrow(ctrdb_samples) == 0) {
    stop("No CTRDB samples found with drug: ", select_drugs)
  }

  # Apply data type filter if specified
  if (!is.null(data_type) && data_type != "all") {
    if ("DataType" %in% colnames(ctrdb_samples)) {
      ctrdb_samples <- ctrdb_samples[ctrdb_samples$DataType == data_type, ]
      if (nrow(ctrdb_samples) == 0) {
        stop("No samples found with data type: ", data_type)
      }
    } else {
      warning("DataType column not found. Skipping data type filter.")
    }
  }

  # Apply tumor type filter if specified
  if (!is.null(tumor_type) && tumor_type != "all") {
    if ("TumorType" %in% colnames(ctrdb_samples)) {
      # Handle case-insensitive matching
      tumor_pattern <- paste0("^", tumor_type, "$", collapse = "|")
      ctrdb_samples <- ctrdb_samples[grepl(tumor_pattern, ctrdb_samples$TumorType, ignore.case = TRUE), ]
      if (nrow(ctrdb_samples) == 0) {
        stop("No samples found with tumor type: ", tumor_type)
      }
    } else {
      warning("TumorType column not found. Skipping tumor type filter.")
    }
  }

  # Get unique PatientIDs
  if (!"PatientID" %in% colnames(ctrdb_samples)) {
    stop("PatientID column not found in sample annotations")
  }

  unique_datasets <- unique(ctrdb_samples$PatientID)
  unique_datasets <- unique_datasets[!is.na(unique_datasets)]

  if (length(unique_datasets) == 0) {
    stop("No valid PatientIDs found")
  }

  message("Found ", length(unique_datasets), " datasets with drug ", select_drugs)

  # Process each patient
  patient_data_list <- list()
  valid_datasets <- character()

  for (patient_id in unique_datasets) {
    tryCatch({
      # Get metadata for this patient
      patient_samples <- ctrdb_samples[ctrdb_samples$PatientID == patient_id, ]

      # Check required columns
      if (!"SampleID" %in% colnames(patient_samples) ||
          !"Response" %in% colnames(patient_samples)) {
        warning("Missing SampleID or Response columns for patient: ", patient_id)
        next
      }

      # Get expression data using the reusable function
      expr_data <- getPatientExpressionData(patient_id, connection, auto_log = TRUE)

      if (is.null(expr_data) || nrow(expr_data) == 0) {
        warning("No expression data found for patient: ", patient_id)
        next
      }

      # Match samples between metadata and expression data
      common_samples <- intersect(patient_samples$SampleID, colnames(expr_data))

      if (length(common_samples) < 3) {
        warning("Insufficient samples for patient: ", patient_id, " (", length(common_samples), " samples)")
        next
      }

      # Filter data to common samples
      patient_meta <- patient_samples[patient_samples$SampleID %in% common_samples, ]
      expr_values <- as.numeric(expr_data[select_omics, common_samples])
      names(expr_values) <- common_samples

      # Apply z-score normalization if enabled
      if (zscore) {
        mean_expr <- mean(expr_values, na.rm = TRUE)
        sd_expr <- sd(expr_values, na.rm = TRUE)
        if (sd_expr > 0) {
          expr_values <- (expr_values - mean_expr) / sd_expr
        }
      }

      # Group by response
      response_groups <- split(expr_values, patient_meta$Response[match(names(expr_values), patient_meta$SampleID)])

      # Check if we have both response groups
      if (!"Response" %in% names(response_groups) || !"Non_response" %in% names(response_groups)) {
        warning("Missing response groups for patient: ", patient_id)
        next
      }

      response_values <- response_groups$Response
      non_response_values <- response_groups$Non_response

      # Check minimum group sizes
      if (length(response_values) < 2 || length(non_response_values) < 2) {
        warning("Insufficient samples in response groups for patient: ", patient_id)
        next
      }

      # Store data for this patient
      patient_data_list[[patient_id]] <- list(
        response = response_values,
        non_response = non_response_values,
        metadata = patient_meta
      )

      valid_datasets <- c(valid_datasets, patient_id)

    }, error = function(e) {
      warning("Error processing patient ", patient_id, ": ", e$message)
    })
  }

  if (length(patient_data_list) == 0) {
    stop("No valid patient data found for analysis")
  }

  message("Successfully processed ", length(patient_data_list), " datasets")

  # Initialize result list
  result <- list()

  # Transform patient data to match expected format for plotting functions
  transformed_pairs <- lapply(patient_data_list, function(patient_data) {
    list(
      yes = patient_data$response,      # Response samples
      no = patient_data$non_response    # Non-response samples
    )
  })
  names(transformed_pairs) <- names(patient_data_list)

  # Create plots for individual patients/studies
  if (length(transformed_pairs) > 0) {
    if (length(transformed_pairs) == 1) {
      # Single plot case
      result$plot <- plotGroupComparison(
        transformed_pairs[[1]]$no,
        transformed_pairs[[1]]$yes,
        group_labels = c("Non_response", "Response"),
        title = names(transformed_pairs)[1],
        y_label = if(zscore) paste0(select_omics, " Expression (z-score)") else paste(select_omics, "Expression")
      )
    } else if (length(transformed_pairs) > 1) {
      # Multiple plots case
      multi_plot <- plotMultipleGroupComparisons(
        transformed_pairs,
        group_labels = c("Non_response", "Response"),
        x_label = select_omics,
        y_label = if(zscore) "Expression (z-score)" else "Expression"
      )
      # Add common axis label
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                              y_title = if(zscore) paste0(select_omics, " Expression (z-score)") else paste(select_omics, "Expression"))
    }
  }

  # Create merged plot if we have multiple datasets and merged_enabled is TRUE
  if (length(transformed_pairs) > 1 && merged_enabled) {
    # Combine all data for merged analysis
    all_response <- unlist(lapply(transformed_pairs, function(x) x$yes))
    all_non_response <- unlist(lapply(transformed_pairs, function(x) x$no))

    # Apply overall z-score to merged data if zscore is enabled
    if (zscore) {
      all_values <- c(all_response, all_non_response)
      mean_all <- mean(all_values, na.rm = TRUE)
      sd_all <- sd(all_values, na.rm = TRUE)
      if (sd_all > 0) {
        all_response <- (all_response - mean_all) / sd_all
        all_non_response <- (all_non_response - mean_all) / sd_all
      }
    }

    result$merged_plot <- plotGroupComparison(
      all_non_response,
      all_response,
      group_labels = c("Non_response", "Response"),
      title = paste("Merged:", select_omics, "vs", select_drugs),
      y_label = if(zscore) paste0(select_omics, " Expression (z-score)") else paste(select_omics, "Expression")
    )
  }

  # Perform meta-analysis if enabled and we have multiple datasets
  if (meta_enabled && length(transformed_pairs) > 1) {
    meta_result <- analyzeClinicalMeta(patient_data_list)
    if (!is.null(meta_result)) {
      result$meta <- meta_result
    }
  }

  # Store data
  result$data <- patient_data_list

  # Return results if there's a plot or merged plot
  if (is.null(result$plot) && is.null(result$merged_plot)) {
    return(list())
  } else {
    return(result)
  }
}

#' Analyze clinical meta-analysis across datasets
#'
#' @description Performs meta-analysis on clinical drug response data across multiple datasets
#' @param patient_data_list List of patient data from analyzeClinicalDrugResponse
#' @return Meta-analysis results object or NULL if analysis couldn't be performed
#' @export
analyzeClinicalMeta <- function(patient_data_list) {
  # Transform patient data to match expected format for metaCalcConDis
  # Response = yes, Non_response = no
  transformed_pairs <- lapply(patient_data_list, function(patient_data) {
    list(
      yes = patient_data$response,      # Response samples
      no = patient_data$non_response    # Non-response samples
    )
  })
  names(transformed_pairs) <- names(patient_data_list)

  # Use existing meta-analysis function from FuncPairMetaAnalysis
  meta_result <- tryCatch({
    metaCalcConDis(transformed_pairs)
  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(NULL)
  })

  return(meta_result)
}
