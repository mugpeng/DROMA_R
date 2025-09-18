# Clinical Trial Database (CTRDB) Analysis Functions ----

#' Analyze clinical drug response with omics data from CTRDB
#'
#' @description Performs analysis of clinical drug response associations with omics features using CTRDB data
#' @param select_omics Character string specifying the omics feature name
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param meta_enabled Logical, whether to perform meta-analysis across datasets
#' @return A list containing plot, meta-analysis results, and data
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
#' # View meta-analysis forest plot
#' result$forest_plot
#' }
analyzeClinicalDrugResponse <- function(select_omics,
                                      select_drugs,
                                      data_type = "all",
                                      tumor_type = "all",
                                      connection = NULL,
                                      meta_enabled = TRUE) {

  # Validate inputs
  if (missing(select_omics) || is.null(select_omics) || select_omics == "") {
    stop("select_omics must be specified")
  }

  if (missing(select_drugs) || is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
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
    DROMA.Set::getDROMAAnnotation("sample")
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

      # Get expression data directly from CTRDB
      table_name_db <- gsub("-", "_", patient_id)  # Replace "-" with "_" for database table name

      expr_data <- tryCatch({
        # Query expression matrix directly
        expr_query <- paste0("SELECT * FROM ", table_name_db)
        full_expr <- DBI::dbGetQuery(connection, expr_query)

        # Convert to matrix with feature_id as rownames
        rownames(full_expr) <- full_expr$feature_id
        full_expr$feature_id <- NULL

        # Return matrix
        as.matrix(full_expr)
      }, error = function(e) {
        warning("Failed to retrieve expression data for patient ", patient_id, ": ", e$message)
        return(NULL)
      })

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

  # Create plots for each patient
  result$plot <- plotAllClinicaldatasets(patient_data_list, select_omics, select_drugs)

  # Perform meta-analysis if enabled and we have multiple datasets
  if (meta_enabled && length(patient_data_list) > 1) {
    meta_result <- analyzeClinicalMeta(patient_data_list)
    if (!is.null(meta_result)) {
      result$meta <- meta_result
      result$forest_plot <- createClinicalForestPlot(meta_result,
                                                   xlab = paste("Effect Size (", select_omics, " vs ", select_drugs, " Response)"))
    }
  }

  # Store data
  result$data <- patient_data_list
  result$summary <- data.frame(
    PatientID = names(patient_data_list),
    N_Response = sapply(patient_data_list, function(x) length(x$response)),
    N_Non_Response = sapply(patient_data_list, function(x) length(x$non_response)),
    stringsAsFactors = FALSE
  )

  return(result)
}

#' Analyze clinical meta-analysis across datasets
#'
#' @description Performs meta-analysis on clinical drug response data across multiple datasets
#' @param patient_data_list List of patient data from analyzeClinicalDrugResponse
#' @return Meta-analysis results object or NULL if analysis couldn't be performed
#' @export
analyzeClinicalMeta <- function(patient_data_list) {
  options(warn = -1)

  # Initialize list to store test results
  test_list <- list()
  valid_datasets <- character()

  # Analyze each patient
  for (patient_id in names(patient_data_list)) {
    tryCatch({
      patient_data <- patient_data_list[[patient_id]]
      response_values <- patient_data$response
      non_response_values <- patient_data$non_response

      # Check for valid data
      if (length(response_values) < 2 || length(non_response_values) < 2) next

      # Perform statistical test (Wilcoxon test)
      wilcox_test <- wilcox.test(non_response_values, response_values)

      # Calculate effect size (Cliff's Delta)
      cliff_delta <- tryCatch({
        effsize::cliff.delta(non_response_values, response_values)
      }, error = function(e) {
        # Fallback effect size calculation
        list(estimate = (median(response_values, na.rm = TRUE) - median(non_response_values, na.rm = TRUE)) /
                       (mad(c(response_values, non_response_values), na.rm = TRUE) + 1e-10))
      })

      # Store results
      test_list[[patient_id]] <- data.frame(
        p = wilcox_test$p.value,
        effect = cliff_delta$estimate,
        N = length(response_values) + length(non_response_values),
        n1 = length(response_values),
        n2 = length(non_response_values)
      )
      valid_datasets <- c(valid_datasets, patient_id)

    }, error = function(e) {
      warning("Error analyzing patient ", patient_id, ": ", e$message)
    })
  }

  # Return NULL if no tests could be performed
  if (length(test_list) < 1) return(NULL)

  # Prepare data for meta-analysis
  meta_df <- do.call(rbind, test_list)
  meta_df$study <- valid_datasets

  # Calculate standard error for Cliff's Delta
  meta_df$se <- sqrt((1 - meta_df$effect^2) * (meta_df$n1 + meta_df$n2 + 1) /
                     (12 * meta_df$n1 * meta_df$n2))

  # Perform meta-analysis
  tryCatch({
    # Only perform meta-analysis if we have at least 2 studies
    if (nrow(meta_df) >= 2) {
      meta_result <- meta::metagen(TE = effect,
                                   seTE = se,
                                   data = meta_df,
                                   sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
                                   studlab = study)
      return(meta_result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(NULL)
  })
}

#' Plot clinical drug response for all datasets
#'
#' @description Creates and combines plots for all clinical datasets showing drug response differences
#' @param patient_data_list List of patient data
#' @param select_omics Name of the omics feature
#' @param select_drugs Name of the drug
#' @return A combined plot with all patient comparisons or NULL if no valid data
#' @export
plotAllClinicaldatasets <- function(patient_data_list, select_omics, select_drugs) {

  # Initialize list to store plots
  p_list <- list()

  # Create plot for each patient
  for (patient_id in names(patient_data_list)) {
    tryCatch({
      patient_data <- patient_data_list[[patient_id]]
      response_values <- patient_data$response
      non_response_values <- patient_data$non_response

      # Ensure adequate data for plotting
      if (length(response_values) < 2 || length(non_response_values) < 2) next

      # Create plot and add to list
      p_list[[patient_id]] <- plotClinicalPatient(response_values, non_response_values, patient_id)

    }, error = function(e) {
      warning("Error creating plot for patient ", patient_id, ": ", e$message)
    })
  }

  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]

  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (length(p_list) <= 3) {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = length(p_list))
      } else {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = 3)
      }

      # Add overall title
      combined_plot <- combined_plot +
        patchwork::plot_annotation(
          title = paste(select_omics, "Expression vs", select_drugs, "Response (Clinical Data)"),
          theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
        )

      return(combined_plot)
    } else {
      warning("patchwork package not available. Returning list of individual plots.")
      return(p_list)
    }
  } else {
    return(NULL)
  }
}

#' Plot clinical drug response for a single patient
#'
#' @description Creates a boxplot comparing omics expression between response and non-response samples
#' @param response_values Expression values for response samples
#' @param non_response_values Expression values for non-response samples
#' @param patient_id Patient identifier for plot title
#' @return A ggplot2 object with boxplot and statistical test
#' @export
plotClinicalPatient <- function(response_values, non_response_values, patient_id) {

  # Combine data into dataframe
  box_df <- data.frame(
    expression = c(non_response_values, response_values),
    response = rep(c("Non-Response", "Response"),
                   times = c(length(non_response_values), length(response_values)))
  )

  # Create boxplot with statistical test
  p <- ggplot2::ggplot(box_df, ggplot2::aes(x = response, y = expression, fill = response)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    ggplot2::scale_fill_manual(values = c("Non-Response" = "#BEBADAFF", "Response" = "#FB8072FF")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      title = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::ylab("Expression Level") +
    ggplot2::ggtitle(paste("Patient:", patient_id))

  # Add statistical test if ggpubr is available
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    p <- p + ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      size = 3.5,
      label.x = 1.5,
      label.y = max(box_df$expression) * 1.1
    )
  }

  return(p)
}

#' Create a forest plot for clinical meta-analysis results
#'
#' @description Creates a standardized forest plot for visualizing clinical meta-analysis results
#' @param meta_obj Meta-analysis object from metagen() function
#' @param xlab Label for x-axis
#' @param show_common Logical, whether to show common effect model
#' @return A forest plot visualization
#' @export
createClinicalForestPlot <- function(meta_obj,
                                   xlab = "Effect Size (95% CI)",
                                   show_common = FALSE) {

  # Validate input
  if (!inherits(meta_obj, "meta") && !inherits(meta_obj, "metagen")) {
    stop("Input must be a meta-analysis object from the 'meta' package")
  }

  # Format p-value text for random effects model
  p_val <- meta_obj$pval.random
  p_text <- if(p_val < 0.001) {
    paste("Random-Effects Model (p =", format(p_val, scientific = TRUE, digits = 3), ")")
  } else {
    paste("Random-Effects Model (p =", round(p_val, 3), ")")
  }

  # Create forest plot with clinical-specific styling
  tryCatch({
    meta::forest(meta_obj,
                 xlab = xlab,
                 slab = paste("Patient", meta_obj$studlab),
                 print.pval.common = show_common,
                 boxsize = 0.3,
                 lineheight = "auto",
                 print.pval.Q = TRUE,
                 print.I2 = TRUE,
                 print.tau2 = FALSE,
                 common = show_common,
                 text.random = p_text,
                 col.diamond = "red",
                 col.diamond.lines = "red",
                 colgap.forest.left = "2cm",
                 colgap.forest.right = "2cm")
  }, error = function(e) {
    warning("Could not create forest plot: ", e$message)
    return(NULL)
  })
}

#' Get clinical drug response summary statistics
#'
#' @description Provides summary statistics for clinical drug response analysis
#' @param result_obj Result object from analyzeClinicalDrugResponse
#' @return Data frame with summary statistics
#' @export
getClinicalSummary <- function(result_obj) {

  if (is.null(result_obj$data) || length(result_obj$data) == 0) {
    return(NULL)
  }

  # Calculate summary statistics for each patient
  summary_stats <- data.frame(
    PatientID = character(),
    N_Response = integer(),
    N_Non_Response = integer(),
    Mean_Response = numeric(),
    Mean_Non_Response = numeric(),
    P_Value = numeric(),
    Effect_Size = numeric(),
    stringsAsFactors = FALSE
  )

  for (patient_id in names(result_obj$data)) {
    patient_data <- result_obj$data[[patient_id]]

    # Calculate basic statistics
    n_resp <- length(patient_data$response)
    n_non_resp <- length(patient_data$non_response)
    mean_resp <- mean(patient_data$response, na.rm = TRUE)
    mean_non_resp <- mean(patient_data$non_response, na.rm = TRUE)

    # Perform statistical test
    test_result <- tryCatch({
      wilcox.test(patient_data$non_response, patient_data$response)
    }, error = function(e) list(p.value = NA))

    # Calculate effect size
    effect_size <- tryCatch({
      effsize::cliff.delta(patient_data$non_response, patient_data$response)$estimate
    }, error = function(e) NA)

    # Add to summary
    summary_stats <- rbind(summary_stats, data.frame(
      PatientID = patient_id,
      N_Response = n_resp,
      N_Non_Response = n_non_resp,
      Mean_Response = round(mean_resp, 3),
      Mean_Non_Response = round(mean_non_resp, 3),
      P_Value = round(test_result$p.value, 4),
      Effect_Size = round(effect_size, 3),
      stringsAsFactors = FALSE
    ))
  }

  # Add meta-analysis results if available
  if (!is.null(result_obj$meta)) {
    meta_row <- data.frame(
      PatientID = "Meta-Analysis",
      N_Response = sum(summary_stats$N_Response),
      N_Non_Response = sum(summary_stats$N_Non_Response),
      Mean_Response = NA,
      Mean_Non_Response = NA,
      P_Value = round(result_obj$meta$pval.random, 4),
      Effect_Size = round(result_obj$meta$TE.random, 3),
      stringsAsFactors = FALSE
    )

    summary_stats <- rbind(summary_stats, meta_row)
  }

  return(summary_stats)
}
