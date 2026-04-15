# Clinical Analysis and Prioritization Functions ----

#' Compute ROC coordinates and AUC for response prediction from expression
#'
#' @description
#' Builds an empirical ROC curve using expression as a continuous predictor
#' for binary response (response vs non-response). With \code{positive_direction = "auto"},
#' the score sign is chosen so that higher values imply response when the response
#' group has higher median expression, and conversely.
#'
#' @param no_values Numeric vector, non-response group (negative class, label 0).
#' @param yes_values Numeric vector, response group (positive class, label 1).
#' @param positive_direction One of \code{"auto"}, \code{"high"}, \code{"low"}.
#'   \code{"high"} means higher expression associates with response; \code{"low"} negates scores.
#' @return List with \code{df} (\code{data.frame} with \code{FPR}, \code{TPR}) and
#'   numeric \code{auc}, or \code{list(df = NULL, auc = NA_real_)} if insufficient data.
#' @keywords internal
.computeClinicalResponseRoc <- function(no_values,
                                        yes_values,
                                        positive_direction = c("auto", "high", "low")) {
  positive_direction <- match.arg(positive_direction)
  no_values <- as.numeric(no_values)
  yes_values <- as.numeric(yes_values)
  no_values <- no_values[!is.na(no_values)]
  yes_values <- yes_values[!is.na(yes_values)]
  if (length(no_values) < 1L || length(yes_values) < 1L) {
    return(list(df = NULL, auc = NA_real_))
  }
  score <- c(no_values, yes_values)
  label <- c(rep(0L, length(no_values)), rep(1L, length(yes_values)))
  if (positive_direction == "auto") {
    if (stats::median(yes_values) < stats::median(no_values)) {
      score <- -score
    }
  } else if (positive_direction == "low") {
    score <- -score
  }
  ord <- order(score, decreasing = TRUE)
  y <- label[ord]
  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  if (n_pos < 1L || n_neg < 1L) {
    return(list(df = NULL, auc = NA_real_))
  }
  cum_tp <- cumsum(y == 1L)
  cum_fp <- cumsum(y == 0L)
  tpr <- c(0, cum_tp / n_pos)
  fpr <- c(0, cum_fp / n_neg)
  roc_df <- data.frame(FPR = fpr, TPR = tpr)
  ord_f <- order(roc_df$FPR, roc_df$TPR)
  rf <- roc_df$FPR[ord_f]
  rt <- roc_df$TPR[ord_f]
  auc <- sum(diff(rf) * (head(rt, -1L) + tail(rt, -1L)) / 2)
  list(df = roc_df, auc = auc)
}

#' Analyze clinical meta-analysis across datasets
#'
#' @description Performs meta-analysis on clinical drug response data across multiple datasets
#' @param patient_data_list List of patient data from analyzeClinicalDrugResponse
#' @return Meta-analysis results object or NULL if analysis couldn't be performed
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
#' @param roc_plot Logical, whether to add an ROC / AUC plot (\code{roc_plot}) using the same
#'   pooled expression scale as the merged comparison when multiple cohorts exist and
#'   \code{merged_enabled} is TRUE; otherwise pooled per-cohort values without a second global z-score.
#' @return A list containing plot (individual study plots), merged_plot (merged dataset plot if merged_enabled=TRUE), optional \code{roc_plot} and \code{roc_auc}, meta-analysis results, and data
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
#' # ROC / AUC (expression classifying responders)
#' result$roc_plot
#' result$roc_auc
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
                                      merged_enabled = TRUE,
                                      roc_plot = TRUE) {

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

      # Get expression data using the function from DROMA.Set package (CTRDB_SQLManager.R)
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
        group_labels = c("Non response", "Response"),
        title = names(transformed_pairs)[1],
        y_label = if(zscore) paste0(select_omics, " Expression (z-score)") else paste(select_omics, "Expression")
      )
    } else if (length(transformed_pairs) > 1) {
      # Multiple plots case
      multi_plot <- plotMultipleGroupComparisons(
        transformed_pairs,
        group_labels = c("Non response", "Response"),
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
      group_labels = c("Non response", "Response"),
      title = paste("Merged:", select_omics, "vs", select_drugs),
      y_label = if(zscore) paste0(select_omics, " Expression (z-score)") else paste(select_omics, "Expression")
    )
  }

  # ROC / AUC: same pooling / scaling as merged plot when multiple cohorts and merge enabled
  if (isTRUE(roc_plot) && length(transformed_pairs) > 0L) {
    if (length(transformed_pairs) == 1L) {
      roc_no <- transformed_pairs[[1]]$no
      roc_yes <- transformed_pairs[[1]]$yes
    } else {
      roc_yes <- unlist(lapply(transformed_pairs, function(x) x$yes), use.names = FALSE)
      roc_no <- unlist(lapply(transformed_pairs, function(x) x$no), use.names = FALSE)
      if (isTRUE(zscore) && isTRUE(merged_enabled)) {
        all_roc <- c(roc_yes, roc_no)
        mean_roc <- mean(all_roc, na.rm = TRUE)
        sd_roc <- stats::sd(all_roc, na.rm = TRUE)
        if (!is.na(sd_roc) && sd_roc > 0) {
          roc_yes <- (roc_yes - mean_roc) / sd_roc
          roc_no <- (roc_no - mean_roc) / sd_roc
        }
      }
    }
    roc_gg <- plotClinicalDrugResponseRoc(
      roc_no,
      roc_yes,
      title = "Clinical readout",
      subtitle = sprintf("%s vs %s — response classification", select_omics, select_drugs)
    )
    if (!is.null(roc_gg)) {
      result$roc_plot <- roc_gg
      result$roc_auc <- attr(roc_gg, "roc_auc")
    } else {
      result$roc_auc <- NA_real_
    }
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

#' Batch analysis of clinical drug response for multiple omics features
#'
#' @description For each omics feature, runs analyzeClinicalDrugResponse against a drug,
#'   extracts meta-analysis results, and returns a data frame in the same format as batchFindSignificantFeatures.
#' @param select_omics Character vector of omics feature names (e.g., gene symbols)
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param cores Number of CPU cores (default: 1). Parallelization not implemented; kept for API consistency.
#' @return Data frame with p_value, effect_size, n_datasets, name, q_value. For multi-dataset cases, meta-analysis is used; for single-dataset cases, Wilcoxon test and Cliff's Delta are used (consistent with metaCalcConDis).
#' @export
#' @examples
#' \dontrun{
#' connectCTRDatabase("path/to/ctrdb.sqlite")
#' results <- batchFindClinicalSigResponse(
#'   select_omics = c("EPB41L4B", "FAM83H", "PLEKHA1"),
#'   select_drugs = "Bortezomib"
#' )
#' }
batchFindClinicalSigResponse <- function(select_omics,
                                        select_drugs,
                                        data_type = "all",
                                        tumor_type = "all",
                                        cores = 1) {
  if (is.null(select_omics) || length(select_omics) == 0) {
    stop("select_omics must be a non-empty character vector")
  }
  if (is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
  }

  # Match analyzeClinicalDrugResponse() logic, but preload patient-level
  # metadata and expression matrices once for all genes.
  if (!exists("ctrdb_connection", envir = .GlobalEnv)) {
    stop("No CTRDB database connection found. Connect first with connectCTRDatabase()")
  }
  connection <- get("ctrdb_connection", envir = .GlobalEnv)

  sample_anno <- tryCatch({
    getDROMAAnnotation("sample")
  }, error = function(e) {
    stop("Failed to load sample annotations from DROMA: ", e$message)
  })

  if (is.null(sample_anno) || nrow(sample_anno) == 0) {
    stop("No sample annotations found in database")
  }
  if (!"CliUsedDrug" %in% colnames(sample_anno)) {
    stop("CliUsedDrug column not found in sample annotations")
  }

  drug_samples <- sample_anno[!is.na(sample_anno$CliUsedDrug) &
                               grepl(select_drugs, sample_anno$CliUsedDrug, ignore.case = TRUE), ]
  if (nrow(drug_samples) == 0) {
    stop("No samples found with drug: ", select_drugs)
  }

  if ("ProjectID" %in% colnames(drug_samples)) {
    ctrdb_samples <- drug_samples[drug_samples$ProjectID == "CTRDB", ]
  } else {
    stop("ProjectID column not found in sample annotations")
  }
  if (nrow(ctrdb_samples) == 0) {
    stop("No CTRDB samples found with drug: ", select_drugs)
  }

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

  if (!is.null(tumor_type) && tumor_type != "all") {
    if ("TumorType" %in% colnames(ctrdb_samples)) {
      tumor_pattern <- paste0("^", tumor_type, "$", collapse = "|")
      ctrdb_samples <- ctrdb_samples[grepl(tumor_pattern, ctrdb_samples$TumorType, ignore.case = TRUE), ]
      if (nrow(ctrdb_samples) == 0) {
        stop("No samples found with tumor type: ", tumor_type)
      }
    } else {
      warning("TumorType column not found. Skipping tumor type filter.")
    }
  }

  if (!"PatientID" %in% colnames(ctrdb_samples)) {
    stop("PatientID column not found in sample annotations")
  }

  unique_datasets <- unique(ctrdb_samples$PatientID)
  unique_datasets <- unique_datasets[!is.na(unique_datasets)]
  if (length(unique_datasets) == 0) {
    stop("No valid PatientIDs found")
  }

  message("Found ", length(unique_datasets), " datasets with drug ", select_drugs)

  patient_context <- list()
  for (patient_id in unique_datasets) {
    tryCatch({
      patient_samples <- ctrdb_samples[ctrdb_samples$PatientID == patient_id, ]

      if (!"SampleID" %in% colnames(patient_samples) ||
          !"Response" %in% colnames(patient_samples)) {
        warning("Missing SampleID or Response columns for patient: ", patient_id)
        next
      }

      expr_data <- getPatientExpressionData(patient_id, connection, auto_log = TRUE)
      if (is.null(expr_data) || nrow(expr_data) == 0) {
        warning("No expression data found for patient: ", patient_id)
        next
      }

      common_samples <- intersect(patient_samples$SampleID, colnames(expr_data))
      if (length(common_samples) < 3) {
        warning("Insufficient samples for patient: ", patient_id, " (", length(common_samples), " samples)")
        next
      }

      patient_meta <- patient_samples[patient_samples$SampleID %in% common_samples, ]

      patient_context[[patient_id]] <- list(
        expr_data = expr_data,
        patient_meta = patient_meta,
        common_samples = common_samples
      )
    }, error = function(e) {
      warning("Error processing patient ", patient_id, ": ", e$message)
    })
  }

  if (length(patient_context) == 0) {
    return(data.frame(p_value = numeric(0), effect_size = numeric(0),
                      n_datasets = integer(0), name = character(0), q_value = numeric(0)))
  }

  message("Successfully processed ", length(patient_context), " datasets")

  cal_re_list <- lapply(select_omics, function(gene) {
    patient_data_list <- list()

    for (patient_id in names(patient_context)) {
      ctx <- patient_context[[patient_id]]

      gene_result <- tryCatch({
        expr_values <- as.numeric(ctx$expr_data[gene, ctx$common_samples])
        names(expr_values) <- ctx$common_samples

        mean_expr <- mean(expr_values, na.rm = TRUE)
        sd_expr <- sd(expr_values, na.rm = TRUE)
        if (sd_expr > 0) {
          expr_values <- (expr_values - mean_expr) / sd_expr
        }

        response_groups <- split(
          expr_values,
          ctx$patient_meta$Response[match(names(expr_values), ctx$patient_meta$SampleID)]
        )

        if (!"Response" %in% names(response_groups) || !"Non_response" %in% names(response_groups)) {
          return(NULL)
        }

        response_values <- response_groups$Response
        non_response_values <- response_groups$Non_response

        if (length(response_values) < 2 || length(non_response_values) < 2) {
          return(NULL)
        }

        list(
          response = response_values,
          non_response = non_response_values,
          metadata = ctx$patient_meta
        )
      }, error = function(e) NULL)

      if (!is.null(gene_result)) {
        patient_data_list[[patient_id]] <- gene_result
      }
    }

    if (length(patient_data_list) == 0) return(NULL)

    if (length(patient_data_list) > 1) {
      meta <- analyzeClinicalMeta(patient_data_list)
      if (is.null(meta)) return(NULL)
      p_val <- meta[["pval.random"]]
      eff_size <- meta[["TE.random"]]
      n_datasets <- length(meta[["studlab"]])
    } else {
      dat <- patient_data_list[[1]]
      resp <- dat$response
      non_resp <- dat$non_response
      if (length(resp) < 2 || length(non_resp) < 2) return(NULL)
      wilcox_re <- tryCatch(
        suppressWarnings(stats::wilcox.test(resp, non_resp)),
        error = function(e) NULL
      )
      cliff_re <- tryCatch(
        suppressWarnings(effsize::cliff.delta(resp, non_resp)),
        error = function(e) NULL
      )
      if (is.null(wilcox_re) || is.null(cliff_re)) return(NULL)
      p_val <- wilcox_re$p.value
      eff_size <- cliff_re$estimate
      n_datasets <- 1L
    }

    if (is.na(p_val)) p_val <- 1
    if (is.na(eff_size)) eff_size <- 0
    data.frame(p_value = p_val, effect_size = eff_size, n_datasets = n_datasets)
  })

  valid <- !sapply(cal_re_list, is.null)
  fea_names <- select_omics[valid]
  cal_re_list <- cal_re_list[valid]

  if (length(cal_re_list) == 0) {
    return(data.frame(p_value = numeric(0), effect_size = numeric(0),
                      n_datasets = integer(0), name = character(0), q_value = numeric(0)))
  }

  cal_re_df <- do.call(rbind, cal_re_list)
  cal_re_df$name <- fea_names
  cal_re_df$q_value <- p.adjust(cal_re_df$p_value, method = "BH")
  cal_re_df
}

#' Select significant in vitro candidates from one or more meta-analysis tables
#'
#' @description Filters one or more batch/meta-analysis result tables and returns
#'   significant candidate features. Multiple tables can be combined by
#'   intersection or union.
#' @param ... Named meta-analysis result tables.
#' @param result_list Optional named list of meta-analysis result tables.
#' @param combine One of "single", "intersect", or "union".
#' @param es_t Effect size threshold passed to getSignificantFeatures().
#' @param P_t P-value or q-value threshold passed to getSignificantFeatures().
#' @param use_p_value Whether to use p-values instead of q-values.
#' @param n_datasets_t Optional minimum number of datasets.
#' @param direction_cols Optional direction columns passed to getIntersectSignificantFeatures().
#' @return A data frame of significant candidate features with aggregated summary columns.
#' @export
getInVitroCandidateFeatures <- function(...,
                                        result_list = NULL,
                                        combine = c("single", "intersect", "union"),
                                        es_t = 0.1,
                                        P_t = 0.05,
                                        use_p_value = FALSE,
                                        n_datasets_t = NULL,
                                        direction_cols = NULL) {
  combine <- match.arg(combine)

  df_list <- if (!is.null(result_list)) {
    result_list
  } else {
    list(...)
  }

  if (length(df_list) == 1 && is.list(df_list[[1]]) && !is.data.frame(df_list[[1]])) {
    df_list <- df_list[[1]]
  }

  if (length(df_list) == 0) {
    stop("At least one meta-analysis result table must be provided")
  }

  if (is.null(names(df_list)) || any(names(df_list) == "")) {
    names(df_list) <- paste0("source", seq_along(df_list))
  }

  sig_list <- lapply(df_list, function(df) {
    getSignificantFeatures(
      meta_df = df,
      es_t = es_t,
      P_t = P_t,
      use_p_value = use_p_value,
      n_datasets_t = n_datasets_t
    )
  })
  sig_list <- sig_list[vapply(sig_list, nrow, integer(1)) > 0]

  if (length(sig_list) == 0) {
    return(data.frame(name = character(0)))
  }

  if (length(sig_list) == 1 || combine == "single") {
    out <- sig_list[[1]]
    if (!"direction" %in% colnames(out) && "effect_size" %in% colnames(out)) {
      out$direction <- ifelse(out$effect_size >= 0, "Up", "Down")
    }
    out$sources <- names(sig_list)[1]
    out$n_sources <- 1L
    return(out)
  }

  if (combine == "intersect") {
    out <- do.call(
      getIntersectSignificantFeatures,
      c(sig_list, list(direction_cols = direction_cols))
    )
    if (nrow(out) == 0) {
      return(out)
    }

    effect_cols <- grep("^effect_size_", colnames(out), value = TRUE)
    p_cols <- grep("^p_value_", colnames(out), value = TRUE)
    q_cols <- grep("^q_value_", colnames(out), value = TRUE)
    n_cols <- grep("^n_datasets_", colnames(out), value = TRUE)
    dir_cols <- grep("^direction_", colnames(out), value = TRUE)

    if (length(effect_cols) > 0) out$effect_size <- rowMeans(out[, effect_cols, drop = FALSE], na.rm = TRUE)
    if (length(p_cols) > 0) out$p_value <- apply(out[, p_cols, drop = FALSE], 1, max, na.rm = TRUE)
    if (length(q_cols) > 0) out$q_value <- apply(out[, q_cols, drop = FALSE], 1, max, na.rm = TRUE)
    if (length(n_cols) > 0) out$n_datasets <- rowSums(out[, n_cols, drop = FALSE], na.rm = TRUE)
    if (length(dir_cols) > 0) out$direction <- out[[dir_cols[1]]]
    out$sources <- paste(names(sig_list), collapse = ";")
    out$n_sources <- length(sig_list)
    return(out)
  }

  union_df <- do.call(rbind, lapply(names(sig_list), function(nm) {
    df <- sig_list[[nm]]
    df$source <- nm
    if (!"direction" %in% colnames(df) && "effect_size" %in% colnames(df)) {
      df$direction <- ifelse(df$effect_size >= 0, "Up", "Down")
    }
    df
  }))

  split_list <- split(union_df, union_df$name)
  out <- do.call(rbind, lapply(split_list, function(df) {
    eff <- if ("effect_size" %in% names(df)) mean(df$effect_size, na.rm = TRUE) else NA_real_
    p_val <- if ("p_value" %in% names(df)) min(df$p_value, na.rm = TRUE) else NA_real_
    q_val <- if ("q_value" %in% names(df)) min(df$q_value, na.rm = TRUE) else NA_real_
    n_ds <- if ("n_datasets" %in% names(df)) max(df$n_datasets, na.rm = TRUE) else NA_real_
    direction <- if ("direction" %in% names(df)) {
      names(sort(table(df$direction), decreasing = TRUE))[1]
    } else if (!is.na(eff)) {
      ifelse(eff >= 0, "Up", "Down")
    } else {
      NA_character_
    }

    data.frame(
      name = df$name[1],
      effect_size = eff,
      p_value = p_val,
      q_value = q_val,
      n_datasets = n_ds,
      direction = direction,
      sources = paste(unique(df$source), collapse = ";"),
      n_sources = length(unique(df$source)),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out[order(out$p_value, decreasing = FALSE), , drop = FALSE]
}

#' Extract pathway-level support for candidate genes
#'
#' @description Summarizes pathway support for candidate genes using either
#'   enrichment results (with a geneID column) or significant pathway-level
#'   meta-analysis results plus a gene set definition list.
#' @param candidates Optional candidate names or a data frame containing a name column.
#' @param enrich_result Optional enrichResult object or data frame from runGO()/runKEGG().
#' @param pathway_meta Optional pathway-level result table (for example GSVA meta results).
#' @param gene_sets Optional named list mapping pathway names to genes.
#' @param P_t Significance threshold.
#' @param use_p_value Whether to use p-values instead of q-values/p.adjust.
#' @param es_t Optional effect size threshold for pathway_meta.
#' @param top_n Optional number of top pathways to retain.
#' @return A list with detail, summary, and pathway_results data frames.
#' @export
extractPathwayGeneSupport <- function(candidates = NULL,
                                      enrich_result = NULL,
                                      pathway_meta = NULL,
                                      gene_sets = NULL,
                                      P_t = 0.05,
                                      use_p_value = FALSE,
                                      es_t = NULL,
                                      top_n = NULL) {
  if (is.null(enrich_result) && (is.null(pathway_meta) || is.null(gene_sets))) {
    stop("Provide either enrich_result, or pathway_meta together with gene_sets")
  }

  candidate_names <- if (is.null(candidates)) {
    NULL
  } else if (is.data.frame(candidates)) {
    candidates$name
  } else {
    as.character(candidates)
  }

  if (!is.null(enrich_result)) {
    pathway_df <- as.data.frame(enrich_result)
    if (nrow(pathway_df) == 0) {
      detail_df <- data.frame(name = character(0), pathway = character(0), pathway_p = numeric(0))
      summary_df <- data.frame(name = character(0))
      return(list(detail = detail_df, summary = summary_df, pathway_results = pathway_df))
    }

    p_col <- if (use_p_value && "pvalue" %in% colnames(pathway_df)) {
      "pvalue"
    } else if ("p.adjust" %in% colnames(pathway_df)) {
      "p.adjust"
    } else if ("qvalue" %in% colnames(pathway_df)) {
      "qvalue"
    } else if ("q_value" %in% colnames(pathway_df)) {
      "q_value"
    } else if ("p_value" %in% colnames(pathway_df)) {
      "p_value"
    } else {
      stop("Could not identify p-value column in enrich_result")
    }

    term_col <- if ("Description" %in% colnames(pathway_df)) {
      "Description"
    } else if ("ID" %in% colnames(pathway_df)) {
      "ID"
    } else if ("name" %in% colnames(pathway_df)) {
      "name"
    } else {
      stop("Could not identify pathway term column in enrich_result")
    }

    gene_col <- if ("geneID" %in% colnames(pathway_df)) {
      "geneID"
    } else {
      stop("enrich_result must contain a geneID column")
    }

    pathway_df <- pathway_df[pathway_df[[p_col]] <= P_t, , drop = FALSE]
    pathway_df <- pathway_df[order(pathway_df[[p_col]], decreasing = FALSE), , drop = FALSE]
    if (!is.null(top_n) && nrow(pathway_df) > top_n) {
      pathway_df <- pathway_df[seq_len(top_n), , drop = FALSE]
    }

    detail_list <- lapply(seq_len(nrow(pathway_df)), function(i) {
      genes <- unique(trimws(unlist(strsplit(pathway_df[[gene_col]][i], "[/;]"))))
      data.frame(
        name = genes,
        pathway = pathway_df[[term_col]][i],
        pathway_p = pathway_df[[p_col]][i],
        stringsAsFactors = FALSE
      )
    })
    detail_df <- do.call(rbind, detail_list)
  } else {
    pathway_df <- as.data.frame(pathway_meta)
    p_col <- if (use_p_value && "p_value" %in% colnames(pathway_df)) {
      "p_value"
    } else if ("q_value" %in% colnames(pathway_df)) {
      "q_value"
    } else if ("p_value" %in% colnames(pathway_df)) {
      "p_value"
    } else {
      stop("Could not identify p-value column in pathway_meta")
    }

    term_col <- if ("name" %in% colnames(pathway_df)) {
      "name"
    } else if ("Description" %in% colnames(pathway_df)) {
      "Description"
    } else {
      stop("Could not identify pathway term column in pathway_meta")
    }

    pathway_df <- pathway_df[pathway_df[[p_col]] <= P_t, , drop = FALSE]
    if (!is.null(es_t) && "effect_size" %in% colnames(pathway_df)) {
      pathway_df <- pathway_df[abs(pathway_df$effect_size) >= es_t, , drop = FALSE]
    }
    pathway_df <- pathway_df[order(pathway_df[[p_col]], decreasing = FALSE), , drop = FALSE]
    if (!is.null(top_n) && nrow(pathway_df) > top_n) {
      pathway_df <- pathway_df[seq_len(top_n), , drop = FALSE]
    }

    valid_terms <- intersect(pathway_df[[term_col]], names(gene_sets))
    pathway_df <- pathway_df[pathway_df[[term_col]] %in% valid_terms, , drop = FALSE]
    detail_list <- lapply(seq_len(nrow(pathway_df)), function(i) {
      genes <- unique(as.character(gene_sets[[pathway_df[[term_col]][i]]]))
      data.frame(
        name = genes,
        pathway = pathway_df[[term_col]][i],
        pathway_p = pathway_df[[p_col]][i],
        stringsAsFactors = FALSE
      )
    })
    detail_df <- if (length(detail_list) > 0) do.call(rbind, detail_list) else data.frame(name = character(0), pathway = character(0), pathway_p = numeric(0))
  }

  if (nrow(detail_df) == 0) {
    summary_df <- if (is.null(candidate_names)) data.frame(name = character(0)) else data.frame(name = candidate_names, mechanism_supported = FALSE, n_supported_pathways = 0L, top_pathway = NA_character_, supporting_pathways = NA_character_)
    return(list(detail = detail_df, summary = summary_df, pathway_results = pathway_df))
  }

  if (!is.null(candidate_names)) {
    detail_df <- detail_df[detail_df$name %in% candidate_names, , drop = FALSE]
  }

  split_detail <- split(detail_df, detail_df$name)
  summary_df <- do.call(rbind, lapply(split_detail, function(df) {
    df <- df[order(df$pathway_p, decreasing = FALSE), , drop = FALSE]
    data.frame(
      name = df$name[1],
      mechanism_supported = TRUE,
      n_supported_pathways = length(unique(df$pathway)),
      top_pathway = df$pathway[1],
      supporting_pathways = paste(unique(df$pathway)[seq_len(min(5, length(unique(df$pathway))))], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  rownames(summary_df) <- NULL

  if (!is.null(candidate_names)) {
    summary_df <- merge(
      data.frame(name = candidate_names, stringsAsFactors = FALSE),
      summary_df,
      by = "name",
      all.x = TRUE
    )
    summary_df$mechanism_supported[is.na(summary_df$mechanism_supported)] <- FALSE
    summary_df$n_supported_pathways[is.na(summary_df$n_supported_pathways)] <- 0L
  }

  list(detail = detail_df, summary = summary_df, pathway_results = pathway_df)
}

#' Prioritize candidates by integrating in vitro, clinical, and pathway evidence
#'
#' @description Merges feature-level, clinical, and mechanism-level evidence into
#'   a single ranking table with an interpretable priority score.
#' @param invitro_candidates In vitro candidate table.
#' @param clinical_candidates Optional clinical candidate table.
#' @param pathway_support Optional pathway support summary table, or the summary
#'   element returned by extractPathwayGeneSupport().
#' @param direction_cols Named vector specifying the direction columns for
#'   in vitro and clinical tables.
#' @param multi_dataset_threshold Minimum in vitro dataset count for multi-dataset support.
#' @param weights Named numeric vector for scoring evidence types.
#' @return A ranked data frame of integrated candidates.
#' @export
prioritizeIntegratedCandidates <- function(invitro_candidates,
                                           clinical_candidates = NULL,
                                           pathway_support = NULL,
                                           direction_cols = c(invitro = "direction", clinical = "direction"),
                                           multi_dataset_threshold = 2L,
                                           weights = c(invitro = 2, clinical = 2, pathway = 1, multi = 1, direction = 1)) {
  if (is.null(invitro_candidates) || !"name" %in% colnames(invitro_candidates)) {
    stop("invitro_candidates must be a data frame containing a 'name' column")
  }

  defaults <- c(invitro = 2, clinical = 2, pathway = 1, multi = 1, direction = 1)
  defaults[names(weights)] <- weights
  weights <- defaults

  clinical_df <- if (is.null(clinical_candidates)) {
    data.frame(name = character(0))
  } else {
    clinical_candidates
  }

  pathway_df <- if (is.null(pathway_support)) {
    data.frame(name = character(0))
  } else if (is.list(pathway_support) && "summary" %in% names(pathway_support)) {
    pathway_support$summary
  } else {
    pathway_support
  }

  base_names <- unique(c(
    invitro_candidates$name,
    if ("name" %in% colnames(clinical_df)) clinical_df$name else character(0),
    if ("name" %in% colnames(pathway_df)) pathway_df$name else character(0)
  ))
  out <- data.frame(name = base_names, stringsAsFactors = FALSE)

  iv_dir_col <- if ("invitro" %in% names(direction_cols)) direction_cols[["invitro"]] else "direction"
  cl_dir_col <- if ("clinical" %in% names(direction_cols)) direction_cols[["clinical"]] else "direction"

  iv <- data.frame(
    name = invitro_candidates$name,
    effect_size_invitro = if ("effect_size" %in% colnames(invitro_candidates)) invitro_candidates$effect_size else NA_real_,
    p_value_invitro = if ("p_value" %in% colnames(invitro_candidates)) invitro_candidates$p_value else NA_real_,
    q_value_invitro = if ("q_value" %in% colnames(invitro_candidates)) invitro_candidates$q_value else NA_real_,
    n_datasets_invitro = if ("n_datasets" %in% colnames(invitro_candidates)) invitro_candidates$n_datasets else NA_real_,
    direction_invitro = if (iv_dir_col %in% colnames(invitro_candidates)) invitro_candidates[[iv_dir_col]] else if ("effect_size" %in% colnames(invitro_candidates)) ifelse(invitro_candidates$effect_size >= 0, "Up", "Down") else NA_character_,
    sources_invitro = if ("sources" %in% colnames(invitro_candidates)) invitro_candidates$sources else NA_character_,
    stringsAsFactors = FALSE
  )

  out <- merge(out, iv, by = "name", all.x = TRUE)

  if ("name" %in% colnames(clinical_df)) {
    cl <- data.frame(
      name = clinical_df$name,
      effect_size_clinical = if ("effect_size" %in% colnames(clinical_df)) clinical_df$effect_size else NA_real_,
      p_value_clinical = if ("p_value" %in% colnames(clinical_df)) clinical_df$p_value else NA_real_,
      q_value_clinical = if ("q_value" %in% colnames(clinical_df)) clinical_df$q_value else NA_real_,
      n_datasets_clinical = if ("n_datasets" %in% colnames(clinical_df)) clinical_df$n_datasets else NA_real_,
      direction_clinical = if (cl_dir_col %in% colnames(clinical_df)) clinical_df[[cl_dir_col]] else if ("effect_size" %in% colnames(clinical_df)) ifelse(clinical_df$effect_size >= 0, "Up", "Down") else NA_character_,
      stringsAsFactors = FALSE
    )
    out <- merge(out, cl, by = "name", all.x = TRUE)
  }

  if ("name" %in% colnames(pathway_df)) {
    keep_cols <- intersect(c("name", "mechanism_supported", "n_supported_pathways", "top_pathway", "supporting_pathways"), colnames(pathway_df))
    out <- merge(out, pathway_df[, keep_cols, drop = FALSE], by = "name", all.x = TRUE)
  }

  if (!"mechanism_supported" %in% colnames(out)) out$mechanism_supported <- FALSE
  if (!"n_supported_pathways" %in% colnames(out)) out$n_supported_pathways <- 0L
  out$mechanism_supported[is.na(out$mechanism_supported)] <- FALSE
  out$n_supported_pathways[is.na(out$n_supported_pathways)] <- 0L

  out$invitro_supported <- !is.na(out$p_value_invitro)
  out$clinical_supported <- !is.na(out$p_value_clinical)
  out$multi_dataset_supported <- !is.na(out$n_datasets_invitro) & out$n_datasets_invitro >= multi_dataset_threshold
  out$direction_consistent <- !is.na(out$direction_invitro) & !is.na(out$direction_clinical) &
    out$direction_invitro == out$direction_clinical

  out$priority_score <- weights["invitro"] * as.numeric(out$invitro_supported) +
    weights["clinical"] * as.numeric(out$clinical_supported) +
    weights["pathway"] * as.numeric(out$mechanism_supported) +
    weights["multi"] * as.numeric(out$multi_dataset_supported) +
    weights["direction"] * as.numeric(out$direction_consistent)

  out$evidence_class <- ifelse(
    out$invitro_supported & out$clinical_supported & out$mechanism_supported,
    "Clinical + mechanism + in vitro",
    ifelse(
      out$invitro_supported & out$clinical_supported,
      "Clinical + in vitro",
      ifelse(
        out$invitro_supported & out$mechanism_supported,
        "Mechanism + in vitro",
        ifelse(out$invitro_supported, "In vitro only", "Other")
      )
    )
  )

  out$priority_tier <- ifelse(
    out$priority_score >= 6, "High",
    ifelse(out$priority_score >= 4, "Medium", "Low")
  )

  out <- out[order(
    -out$priority_score,
    -out$n_supported_pathways,
    -abs(ifelse(is.na(out$effect_size_clinical), 0, out$effect_size_clinical)),
    -abs(ifelse(is.na(out$effect_size_invitro), 0, out$effect_size_invitro)),
    out$name
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}
