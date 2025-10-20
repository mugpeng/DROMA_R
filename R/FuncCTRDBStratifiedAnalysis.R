# Stratified Analysis Functions for CTRDB Data ----

# Required packages: GSVA, ggplot2, ggpubr, patchwork, meta, dplyr

# Core principle: Perform stratified analysis on CTRDB clinical data
# to identify drug response signatures and apply them across different drugs

#' Stratified analysis for CTRDB clinical data
#'
#' @description Performs stratified analysis on CTRDB data to identify drug response signatures
#' and analyze the correlation between specific omics features and the signature scores.
#' Each patient's B drug signature is calculated using their own top upregulated genes.
#' @param drug_b_name Name of drug B used for signature generation
#' @param drug_a_name Name of drug A for signature application
#' @param select_omics Character string specifying the omics feature name to analyze in drug A
#' @param connection CTRDB database connection object
#' @param top_n_genes Number of top genes to select from each dataset (default: 100)
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param min_response_size Minimum samples required in response group (default: 3)
#' @param min_non_response_size Minimum samples required in non-response group (default: 3)
#' @param meta_enabled Logical, whether to perform meta-analysis (default: TRUE)
#' @return A list containing analysis results, plots, and signature genes
#' @export
#' @examples
#' \dontrun{
#' # Connect to CTRDB database first
#' con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#'
#' # Perform stratified analysis with specific omics feature
#' result <- analyzeStratifiedCTRDB(
#'   drug_b_name = "Cisplatin",
#'   drug_a_name = "Paclitaxel",
#'   select_omics = "EGFR",
#'   connection = con
#' )
#'
#' # View signature genes
#' result$signature_genes
#'
#' # View comprehensive correlation comparison plot
#' result$correlation_results$comprehensive_plot
#'
#' # View correlation data
#' result$correlation_results$meta_data
#'
#' # View meta-analysis forest plot
#' result$meta_analysis$forest_plot
#' }
analyzeStratifiedCTRDB <- function(drug_b_name,
                                  drug_a_name,
                                  select_omics,
                                  connection,
                                  top_n_genes = 100,
                                  data_type = "all",
                                  tumor_type = "all",
                                  min_response_size = 3,
                                  min_non_response_size = 3,
                                  meta_enabled = TRUE) {

  # Validate inputs
  if (missing(drug_b_name) || is.null(drug_b_name) || drug_b_name == "") {
    stop("drug_b_name must be specified")
  }
  if (missing(drug_a_name) || is.null(drug_a_name) || drug_a_name == "") {
    stop("drug_a_name must be specified")
  }
  if (missing(select_omics) || is.null(select_omics) || select_omics == "") {
    stop("select_omics must be specified")
  }
  if (is.null(connection)) {
    stop("CTRDB database connection must be provided")
  }

  # Initialize result list
  result <- list()

  # Step 1: Get patient data for Drug B and identify upregulated genes
  cat("Step 1: Analyzing Drug B (", drug_b_name, ") for signature generation...\n", sep = "")
  drug_b_result <- analyzeDrugBForSignature(
    drug_name = drug_b_name,
    connection = connection,
    top_n_genes = top_n_genes,
    data_type = data_type,
    tumor_type = tumor_type,
    min_response_size = min_response_size,
    min_non_response_size = min_non_response_size
  )

  if (is.null(drug_b_result) || length(drug_b_result$patient_data) == 0) {
    stop("No valid data found for Drug B: ", drug_b_name)
  }

  result$drug_b_analysis <- drug_b_result

  # Step 2: Select final signature genes based on intersection and logFC ranking
  cat("Step 2: Selecting final signature genes...\n")
  final_signature_genes <- selectFinalSignatureGenes(drug_b_result$top_genes_per_patient)
  result$signature_genes <- final_signature_genes

  # Step 3: Calculate B drug response scores for each patient in Drug B
  cat("Step 3: Calculating B drug response scores...\n")
  drug_b_scores <- calculateDrugBScores(
    patient_data_list = drug_b_result$patient_data
  )
  result$drug_b_scores <- drug_b_scores

  # Step 4: Analyze scores in Drug B datasets (NR vs Response)
  cat("Step 4: Analyzing B drug scores in NR vs Response groups...\n")
  drug_b_score_analysis <- analyzeDrugBScores(drug_b_scores, drug_b_name)
  result$drug_b_score_analysis <- drug_b_score_analysis

  # Step 5: Apply B drug signature to Drug A datasets
  cat("Step 5: Applying B drug signature to Drug A (", drug_a_name, ") datasets...\n", sep = "")
  drug_a_result <- analyzeDrugAWithSignature(
    drug_name = drug_a_name,
    connection = connection,
    signature_genes = final_signature_genes,
    data_type = data_type,
    tumor_type = tumor_type
  )

  if (is.null(drug_a_result) || length(drug_a_result$patient_data) == 0) {
    warning("No valid data found for Drug A: ", drug_a_name)
    result$drug_a_analysis <- NULL
    result$drug_a_scores <- NULL
    result$correlation_results <- NULL
    result$meta_analysis <- NULL
  } else {
    result$drug_a_analysis <- drug_a_result

    # Step 6: Calculate A drug response scores using B drug signature
    cat("Step 6: Calculating A drug response scores using B drug signature...\n")
    drug_a_scores <- calculateDrugAScores(
      patient_data_list = drug_a_result$patient_data,
      signature_genes = final_signature_genes
    )
    result$drug_a_scores <- drug_a_scores

    # Step 7: Analyze correlation between specific omics feature and A drug scores in Drug A
    cat("Step 7: Analyzing correlations for", select_omics, "in Drug A datasets...\n")
    correlation_results <- analyzeCorrelationsInDrugA(
      drug_a_data = drug_a_result$patient_data,
      drug_a_scores = drug_a_scores,
      select_omics = select_omics,
      drug_a_name = drug_a_name
    )
    result$correlation_results <- correlation_results

    # Step 8: Perform meta-analysis if enabled
    if (meta_enabled && !is.null(correlation_results$meta_data) &&
        nrow(correlation_results$meta_data) > 1) {
      cat("Step 8: Performing meta-analysis...\n")
      meta_result <- performMetaAnalysisForCorrelations(correlation_results$meta_data)
      result$meta_analysis <- meta_result
    }
  }

  return(result)
}

#' Analyze Drug B to identify response signature
#'
#' @description Analyzes Drug B response data to identify upregulated genes in responders
#' @param drug_name Name of the drug
#' @param connection Database connection
#' @param top_n_genes Number of top genes to select
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @param min_response_size Minimum response samples
#' @param min_non_response_size Minimum non-response samples
#' @return List containing patient data and top genes per patient
#' @keywords internal
analyzeDrugBForSignature <- function(drug_name, connection, top_n_genes,
                                   data_type, tumor_type, min_response_size,
                                   min_non_response_size) {

  # Get sample annotations from DROMA
  sample_anno <- tryCatch({
    DROMA.Set::getDROMAAnnotation("sample")
  }, error = function(e) {
    stop("Failed to load sample annotations from DROMA: ", e$message)
  })

  # Find samples with the selected drug
  drug_samples <- sample_anno[!is.na(sample_anno$CliUsedDrug) &
                             grepl(drug_name, sample_anno$CliUsedDrug, ignore.case = TRUE), ]

  # Filter for CTRDB project
  if ("ProjectID" %in% colnames(drug_samples)) {
    ctrdb_samples <- drug_samples[drug_samples$ProjectID == "CTRDB", ]
  } else {
    stop("ProjectID column not found in sample annotations")
  }

  # Apply filters
  if (!is.null(data_type) && data_type != "all" && "DataType" %in% colnames(ctrdb_samples)) {
    ctrdb_samples <- ctrdb_samples[ctrdb_samples$DataType == data_type, ]
  }
  if (!is.null(tumor_type) && tumor_type != "all" && "TumorType" %in% colnames(ctrdb_samples)) {
    tumor_pattern <- paste0("^", tumor_type, "$", collapse = "|")
    ctrdb_samples <- ctrdb_samples[grepl(tumor_pattern, ctrdb_samples$TumorType, ignore.case = TRUE), ]
  }

  # Get unique PatientIDs
  unique_datasets <- unique(ctrdb_samples$PatientID)
  unique_datasets <- unique_datasets[!is.na(unique_datasets)]

  # Process each patient
  patient_data_list <- list()
  top_genes_per_patient <- list()

  for (patient_id in unique_datasets) {
    tryCatch({
      # Get metadata for this patient
      patient_samples <- ctrdb_samples[ctrdb_samples$PatientID == patient_id, ]

      # Get expression data
      table_name_db <- gsub("-", "_", patient_id)
      expr_data <- tryCatch({
        expr_query <- paste0("SELECT * FROM ", table_name_db)
        full_expr <- DBI::dbGetQuery(connection, expr_query)
        rownames(full_expr) <- full_expr$feature_id
        full_expr$feature_id <- NULL
        as.matrix(full_expr)
      }, error = function(e) {
        warning("Failed to retrieve expression data for patient ", patient_id)
        return(NULL)
      })

      if (is.null(expr_data) || nrow(expr_data) == 0) next

      # Match samples
      common_samples <- intersect(patient_samples$SampleID, colnames(expr_data))
      if (length(common_samples) < (min_response_size + min_non_response_size)) next

      # Filter data
      patient_meta <- patient_samples[patient_samples$SampleID %in% common_samples, ]
      expr_data_subset <- expr_data[, common_samples]

      # Group by response
      response_groups <- split(colnames(expr_data_subset),
                              patient_meta$Response[match(colnames(expr_data_subset),
                                                         patient_meta$SampleID)])

      # Check groups
      if (!"Response" %in% names(response_groups) || !"Non_response" %in% names(response_groups)) next

      response_samples <- response_groups$Response
      non_response_samples <- response_groups$Non_response

      if (length(response_samples) < min_response_size ||
          length(non_response_samples) < min_non_response_size) next

      # Perform differential expression analysis
      response_expr <- expr_data_subset[, response_samples, drop = FALSE]
      non_response_expr <- expr_data_subset[, non_response_samples, drop = FALSE]

      # Calculate log2 fold change (mean difference)
      mean_response <- rowMeans(response_expr, na.rm = TRUE)
      mean_non_response <- rowMeans(non_response_expr, na.rm = TRUE)
      logFC <- mean_response - mean_non_response

      # Perform Wilcoxon test for each gene
      p_values <- apply(expr_data_subset, 1, function(gene_expr) {
        wilcox.test(gene_expr[response_samples],
                   gene_expr[non_response_samples])$p.value
      })

      # Adjust p-values
      p_adj <- p.adjust(p_values, method = "BH")

      # Create results data frame
      de_results <- data.frame(
        gene = rownames(expr_data_subset),
        logFC = logFC,
        p_value = p_values,
        p_adj = p_adj,
        stringsAsFactors = FALSE
      )

      # Sort by logFC (descending) and get top upregulated genes
      de_results <- de_results[order(de_results$logFC, decreasing = TRUE), ]
      n_top <- min(top_n_genes, nrow(de_results))
      top_genes <- if (n_top > 0) de_results$gene[seq_len(n_top)] else character(0)

      # Store data
      patient_data_list[[patient_id]] <- list(
        response_expr = response_expr,
        non_response_expr = non_response_expr,
        response_samples = response_samples,
        non_response_samples = non_response_samples,
        metadata = patient_meta,
        top_genes = top_genes
      )

      top_genes_per_patient[[patient_id]] <- top_genes

    }, error = function(e) {
      warning("Error processing patient ", patient_id, ": ", e$message)
    })
  }

  return(list(
    patient_data = patient_data_list,
    top_genes_per_patient = top_genes_per_patient
  ))
}

#' Select final signature genes based on intersection and logFC ranking
#'
#' @description Selects final signature genes by prioritizing genes that appear
#' in multiple datasets and have high logFC values
#' @param top_genes_per_patient List of top genes for each patient
#' @return Vector of final signature genes
#' @keywords internal
selectFinalSignatureGenes <- function(top_genes_per_patient) {

  if (length(top_genes_per_patient) == 0) {
    stop("No top genes data available")
  }

  # Count gene occurrences across all datasets
  all_genes <- unlist(top_genes_per_patient)
  gene_counts <- table(all_genes)

  # Create ranking data frame
  gene_ranking <- data.frame(
    gene = names(gene_counts),
    count = as.integer(gene_counts),
    stringsAsFactors = FALSE
  )

  # Calculate average rank by logFC for each gene
  avg_logfc_rank <- numeric(nrow(gene_ranking))
  names(avg_logfc_rank) <- gene_ranking$gene

  for (patient_id in names(top_genes_per_patient)) {
    genes <- top_genes_per_patient[[patient_id]]
    for (i in seq_along(genes)) {
      gene <- genes[i]
      if (gene %in% names(avg_logfc_rank)) {
        avg_logfc_rank[gene] <- avg_logfc_rank[gene] + i
      }
    }
  }

  # Average the ranks
  avg_logfc_rank <- avg_logfc_rank / length(top_genes_per_patient)

  gene_ranking$avg_logfc_rank <- avg_logfc_rank[gene_ranking$gene]

  # Sort by count (descending) then by avg_logfc_rank (descending)
  gene_ranking <- gene_ranking[order(-gene_ranking$count, -gene_ranking$avg_logfc_rank), ]

  # Select top 100 genes
  n_final <- min(100, nrow(gene_ranking))
  final_genes <- if (n_final > 0) gene_ranking$gene[seq_len(n_final)] else character(0)

  message("Selected ", length(final_genes), " signature genes")
  if (length(final_genes) > 0) {
    message("Genes appear in ", round(mean(gene_ranking$count[seq_along(final_genes)]), 1),
            " datasets on average")
  }

  return(final_genes)
}

#' Calculate B drug response scores using ssGSEA
#'
#' @description Calculates drug response scores using ssGSEA on each patient's top genes
#' @param patient_data_list List of patient expression data with top_genes included
#' @return List of scores for each patient
#' @keywords internal
calculateDrugBScores <- function(patient_data_list) {

  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for ssGSEA calculation")
  }

  scores_list <- list()

  for (patient_id in names(patient_data_list)) {
    patient_data <- patient_data_list[[patient_id]]

    # Check if top_genes are available
    if (!"top_genes" %in% names(patient_data)) {
      warning("No top_genes found for patient ", patient_id)
      next
    }

    # Get patient's top genes
    patient_top_genes <- patient_data$top_genes

    # Combine response and non-response expression data
    all_expr <- cbind(patient_data$response_expr, patient_data$non_response_expr)

    if (nrow(all_expr) == 0) {
      warning("No top genes found in expression data for patient ", patient_id)
      next
    }

    # Create gene set list
    gene_sets <- list(signature = patient_top_genes)

    # Calculate ssGSEA scores
    tryCatch({
      # Create ssGSEA parameter object for new GSVA API
      ssgsea_param <- GSVA::ssgseaParam(
        exprData = as.matrix(all_expr),
        geneSets = gene_sets,
        alpha = 0.25,
        normalize = TRUE
      )

      # Run ssGSEA with parameter object
      gsva_result <- GSVA::gsva(ssgsea_param)

      # Extract scores
      scores <- as.numeric(gsva_result["signature", ])
      names(scores) <- colnames(gsva_result)

      scores_list[[patient_id]] <- list(
        scores = scores,
        response_samples = patient_data$response_samples,
        non_response_samples = patient_data$non_response_samples
      )

    }, error = function(e) {
      warning("Failed to calculate ssGSEA scores for patient ", patient_id, ": ", e$message)
    })
  }

  return(scores_list)
}

#' Calculate A drug response scores using B drug signature
#'
#' @description Calculates drug A response scores using the signature genes from drug B
#' @param patient_data_list List of patient expression data from Drug A
#' @param signature_genes Final signature genes from Drug B analysis
#' @return List of scores for each patient in Drug A
#' @keywords internal
calculateDrugAScores <- function(patient_data_list, signature_genes) {

  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for ssGSEA calculation")
  }

  scores_list <- list()

  for (patient_id in names(patient_data_list)) {
    patient_data <- patient_data_list[[patient_id]]

    # Get all expression data for this patient
    all_expr <- patient_data$expr_data

    if (nrow(all_expr) == 0) {
      warning("No expression data found for patient ", patient_id)
      next
    }

    # Filter signature genes present in expression data
    available_genes <- intersect(signature_genes, rownames(all_expr))

    if (length(available_genes) < 10) {
      warning("Too few signature genes (", length(available_genes),
              ") found in expression data for patient ", patient_id)
      next
    }

    # Create gene set list
    gene_sets <- list(signature = available_genes)

    # Calculate ssGSEA scores
    tryCatch({
      # Create ssGSEA parameter object for new GSVA API
      ssgsea_param <- GSVA::ssgseaParam(
        exprData = as.matrix(all_expr),
        geneSets = gene_sets,
        alpha = 0.25,
        normalize = TRUE
      )

      # Run ssGSEA with parameter object
      gsva_result <- GSVA::gsva(ssgsea_param)

      # Extract scores
      scores <- as.numeric(gsva_result["signature", ])
      names(scores) <- colnames(gsva_result)

      scores_list[[patient_id]] <- list(
        scores = scores,
        response_samples = patient_data$response_samples,
        non_response_samples = patient_data$non_response_samples
      )

    }, error = function(e) {
      warning("Failed to calculate ssGSEA scores for patient ", patient_id, ": ", e$message)
    })
  }

  return(scores_list)
}

#' Analyze B drug scores in NR vs Response groups
#'
#' @description Compares B drug response scores between non-responders and responders
#' @param drug_b_scores List of B drug scores
#' @param drug_name Name of drug B
#' @return List containing analysis results and plots
#' @keywords internal
analyzeDrugBScores <- function(drug_b_scores, drug_name) {

  results <- list()
  plots <- list()
  summary_data <- data.frame()

  for (patient_id in names(drug_b_scores)) {
    score_data <- drug_b_scores[[patient_id]]

    # Extract scores for each group
    response_scores <- score_data$scores[score_data$response_samples]
    non_response_scores <- score_data$scores[score_data$non_response_samples]

    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(non_response_scores, response_scores)

    # Store results
    patient_result <- data.frame(
      PatientID = patient_id,
      N_Response = length(response_scores),
      N_Non_Response = length(non_response_scores),
      Mean_Response = mean(response_scores, na.rm = TRUE),
      Mean_Non_Response = mean(non_response_scores, na.rm = TRUE),
      P_Value = wilcox_result$p.value,
      stringsAsFactors = FALSE
    )

    summary_data <- rbind(summary_data, patient_result)

    # Create plot
    plot_data <- data.frame(
      score = c(non_response_scores, response_scores),
      group = rep(c("Non-response", "Response"),
                  times = c(length(non_response_scores), length(response_scores)))
    )

    p <- ggpubr::ggboxplot(plot_data, x = "group", y = "score",
                           fill = "group", palette = c("#BEBADAFF", "#FB8072FF"),
                           add = "jitter", add.params = list(alpha = 0.5)) +
      ggpubr::stat_compare_means(method = "wilcox.test",
                                 label.x = 0.8,
                                 label.y = max(plot_data$score, na.rm = TRUE) * 0.95) +
      ggplot2::labs(title = paste("B Drug Response Scores -", patient_id),
           x = "", y = "ssGSEA Score") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    plots[[patient_id]] <- p
  }

  # Combine plots
  if (length(plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = paste("B Drug (", drug_name, ") Response Scores Analysis", sep = ""),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
      )
    results$combined_plot <- combined_plot
  }

  results$plots <- plots
  results$summary <- summary_data

  return(results)
}

#' Analyze Drug A with B drug signature
#'
#' @description Loads and processes Drug A data for correlation analysis
#' @param drug_name Name of drug A
#' @param connection Database connection
#' @param signature_genes B drug signature genes
#' @param data_type Filter by data type
#' @param tumor_type Filter by tumor type
#' @return List containing patient data for Drug A
#' @keywords internal
analyzeDrugAWithSignature <- function(drug_name, connection, signature_genes,
                                   data_type, tumor_type) {

  # Get sample annotations
  sample_anno <- tryCatch({
    DROMA.Set::getDROMAAnnotation("sample")
  }, error = function(e) {
    stop("Failed to load sample annotations: ", e$message)
  })

  # Find samples with drug A
  drug_samples <- sample_anno[!is.na(sample_anno$CliUsedDrug) &
                             grepl(drug_name, sample_anno$CliUsedDrug, ignore.case = TRUE), ]

  # Filter for CTRDB project
  if ("ProjectID" %in% colnames(drug_samples)) {
    ctrdb_samples <- drug_samples[drug_samples$ProjectID == "CTRDB", ]
  } else {
    stop("ProjectID column not found")
  }

  # Apply filters
  if (!is.null(data_type) && data_type != "all" && "DataType" %in% colnames(ctrdb_samples)) {
    ctrdb_samples <- ctrdb_samples[ctrdb_samples$DataType == data_type, ]
  }
  if (!is.null(tumor_type) && tumor_type != "all" && "TumorType" %in% colnames(ctrdb_samples)) {
    tumor_pattern <- paste0("^", tumor_type, "$", collapse = "|")
    ctrdb_samples <- ctrdb_samples[grepl(tumor_pattern, ctrdb_samples$TumorType, ignore.case = TRUE), ]
  }

  # Get unique PatientIDs
  unique_datasets <- unique(ctrdb_samples$PatientID)
  unique_datasets <- unique_datasets[!is.na(unique_datasets)]

  # Process each patient
  patient_data_list <- list()

  for (patient_id in unique_datasets) {
    tryCatch({
      # Get metadata
      patient_samples <- ctrdb_samples[ctrdb_samples$PatientID == patient_id, ]

      # Get expression data
      table_name_db <- gsub("-", "_", patient_id)
      expr_data <- tryCatch({
        expr_query <- paste0("SELECT * FROM ", table_name_db)
        full_expr <- DBI::dbGetQuery(connection, expr_query)
        rownames(full_expr) <- full_expr$feature_id
        full_expr$feature_id <- NULL
        as.matrix(full_expr)
      }, error = function(e) {
        warning("Failed to retrieve expression data for patient ", patient_id)
        return(NULL)
      })

      if (is.null(expr_data) || nrow(expr_data) == 0) next

      # Match samples
      common_samples <- intersect(patient_samples$SampleID, colnames(expr_data))
      if (length(common_samples) < 3) next

      # Filter data
      patient_meta <- patient_samples[patient_samples$SampleID %in% common_samples, ]
      expr_data_subset <- expr_data[, common_samples]

      # Group by response
      response_groups <- split(colnames(expr_data_subset),
                              patient_meta$Response[match(colnames(expr_data_subset),
                                                         patient_meta$SampleID)])

      if (!"Response" %in% names(response_groups) || !"Non_response" %in% names(response_groups)) next

      # Store data
      patient_data_list[[patient_id]] <- list(
        expr_data = expr_data_subset,
        response_samples = response_groups$Response,
        non_response_samples = response_groups$Non_response,
        metadata = patient_meta
      )

    }, error = function(e) {
      warning("Error processing patient ", patient_id, ": ", e$message)
    })
  }

  return(list(patient_data = patient_data_list))
}

#' Analyze correlations in Drug A datasets
#'
#' @description Analyzes correlation between specific omics feature expression and A drug scores
#' @param drug_a_data Drug A patient data
#' @param drug_a_scores Drug A scores calculated using B drug signature
#' @param select_omics Name of the omics feature to analyze
#' @param drug_a_name Name of drug A
#' @return List containing correlation results and plots
#' @keywords internal
analyzeCorrelationsInDrugA <- function(drug_a_data, drug_a_scores, select_omics, drug_a_name) {

  # Data frame to store all correlation results for meta-analysis
  meta_data <- data.frame()

  # List to store individual scatter plots
  scatter_plots <- list()

  for (patient_id in names(drug_a_data)) {
    patient_data <- drug_a_data[[patient_id]]

    # Check if we have scores for this patient
    if (!patient_id %in% names(drug_a_scores)) {
      warning("No scores found for patient ", patient_id)
      next
    }

    score_data <- drug_a_scores[[patient_id]]

    # Get expression data for the selected omics feature
    if (!select_omics %in% rownames(patient_data$expr_data)) {
      warning("Feature ", select_omics, " not found in patient ", patient_id)
      next
    }

    # Calculate correlations for each response group
    for (group in c("Response", "Non_response")) {
      group_samples <- if (group == "Response") {
        score_data$response_samples
      } else {
        score_data$non_response_samples
      }

      if (length(group_samples) < 3) next

      # Get expression and score values for this group
      group_expr <- as.numeric(patient_data$expr_data[select_omics, group_samples])
      group_scores <- score_data$scores[group_samples]

      # Remove NA values
      valid_idx <- !is.na(group_expr) & !is.na(group_scores)
      group_expr <- group_expr[valid_idx]
      group_scores <- group_scores[valid_idx]

      if (length(group_expr) < 3) {
        warning("Insufficient observations for correlation in ", patient_id, " ", group)
        next
      }

      # Calculate Spearman correlation
      cor_test <- tryCatch({
        cor.test(group_expr, group_scores, method = "spearman")
      }, error = function(e) {
        warning("Correlation test failed for ", patient_id, " ", group, ": ", e$message)
        return(NULL)
      })

      if (is.null(cor_test)) next

      # Store results
      result_row <- data.frame(
        patient_id = patient_id,
        group = group,
        cor = as.numeric(cor_test$estimate),
        p_value = cor_test$p.value,
        n_samples = length(group_expr),
        stringsAsFactors = FALSE
      )

      meta_data <- rbind(meta_data, result_row)

      # Create scatter plot
      plot_data <- data.frame(
        expression = group_expr,
        score = group_scores
      )

      cor_text <- paste0("rho = ", round(cor_test$estimate, 3),
                        "\np = ", format.pval(cor_test$p.value, digits = 2))

      p_scatter <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!rlang::sym("expression"), y = !!rlang::sym("score"))) +
        ggplot2::geom_point(alpha = 0.6, size = 2.5, color = ifelse(group == "Response", "#FB8072", "#BEBADA")) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 0.8) +
        ggplot2::labs(
          title = paste(patient_id, "-", group),
          x = paste(select_omics, "Expression"),
          y = "Drug Response Score"
        ) +
        ggplot2::annotate("text", x = Inf, y = Inf,
                 label = cor_text,
                 hjust = 1.1, vjust = 1.1,
                 size = 3.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.title = ggplot2::element_text(size = 9)
        )

      scatter_plots[[paste0(patient_id, "_", group)]] <- p_scatter
    }
  }

  # Create comprehensive comparison plot across all datasets
  comprehensive_plot <- NULL
  if (nrow(meta_data) > 0) {
    # Prepare data for comprehensive comparison
    plot_data <- meta_data
    plot_data$group <- factor(plot_data$group, levels = c("Response", "Non_response"))

    # Create comprehensive boxplot/violin plot
    comprehensive_plot <- ggpubr::ggviolin(
      plot_data,
      x = "group",
      y = "cor",
      fill = "group",
      palette = c("Response" = "#FB8072FF", "Non_response" = "#BEBADAFF"),
      add = "boxplot",
      add.params = list(fill = "white", width = 0.1)
    ) +
      ggplot2::geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
      ggpubr::stat_compare_means(
        method = "wilcox.test",
        label.x = 1.3,
        label.y = max(plot_data$cor, na.rm = TRUE) * 1.05,
        size = 4
      ) +
      ggplot2::labs(
        title = paste("Correlation of", select_omics, "with", drug_a_name, "Response Signature"),
        subtitle = paste("Across", length(unique(plot_data$patient_id)), "datasets"),
        x = "",
        y = "Spearman Correlation Coefficient"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 11),
        legend.position = "none"
      )
  }

  # Combine scatter plots (optional visualization)
  combined_scatter_plot <- NULL
  if (length(scatter_plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined_scatter_plot <- patchwork::wrap_plots(
      scatter_plots,
      ncol = min(4, ceiling(sqrt(length(scatter_plots))))
    ) +
      patchwork::plot_annotation(
        title = paste("Individual Dataset Correlations:", select_omics, "vs Drug Response Score"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
      )
  }

  return(list(
    meta_data = meta_data,
    scatter_plots = scatter_plots,
    comprehensive_plot = comprehensive_plot,
    combined_scatter_plot = combined_scatter_plot
  ))
}

#' Perform meta-analysis for correlation differences
#'
#' @description Performs meta-analysis on correlation differences between groups
#' @param meta_data Data frame with correlation results from all studies and groups
#' @return Meta-analysis result
#' @keywords internal
performMetaAnalysisForCorrelations <- function(meta_data) {

  if (nrow(meta_data) == 0) {
    warning("No data available for meta-analysis")
    return(NULL)
  }

  # Reshape data to get Response vs Non_response for each patient
  response_data <- meta_data[meta_data$group == "Response", ]
  non_response_data <- meta_data[meta_data$group == "Non_response", ]

  # Match by patient_id
  common_patients <- intersect(response_data$patient_id, non_response_data$patient_id)

  if (length(common_patients) < 2) {
    warning("Insufficient studies with both Response and Non_response data for meta-analysis")
    return(NULL)
  }

  # Create paired data frame
  meta_df <- data.frame()
  for (patient_id in common_patients) {
    resp_row <- response_data[response_data$patient_id == patient_id, ]
    non_resp_row <- non_response_data[non_response_data$patient_id == patient_id, ]

    meta_df <- rbind(meta_df, data.frame(
      study = patient_id,
      cor_resp = resp_row$cor,
      cor_non_resp = non_resp_row$cor,
      diff_cor = resp_row$cor - non_resp_row$cor,
      n_resp = resp_row$n_samples,
      n_non_resp = non_resp_row$n_samples,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate standard error for correlation difference
  # Using Fisher's z-transformation approximation
  meta_df$se_diff <- sqrt(
    (1 - meta_df$cor_resp^2)^2 / pmax(meta_df$n_resp - 3, 1) +
    (1 - meta_df$cor_non_resp^2)^2 / pmax(meta_df$n_non_resp - 3, 1)
  )

  # Perform meta-analysis
  tryCatch({
    meta_result <- meta::metagen(
      TE = meta_df$diff_cor,
      seTE = meta_df$se_diff,
      data = meta_df,
      sm = "SMD",  # Standardized Mean Difference
      studlab = meta_df$study,
      comb.fixed = FALSE,
      comb.random = TRUE
    )

    # Create forest plot
    forest_plot <- meta::forest(
      meta_result,
      xlab = "Difference in Correlation (Response - Non-response)",
      slab = "study",
      boxsize = 0.2,
      lineheight = "auto",
      print.pval.Q = TRUE,
      print.I2 = TRUE,
      print.tau2 = TRUE
    )

    return(list(
      meta_result = meta_result,
      forest_plot = forest_plot,
      meta_df = meta_df
    ))

  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(NULL)
  })
}

#' Get summary of stratified CTRDB analysis
#'
#' @description Provides a comprehensive summary of the stratified analysis results
#' @param result_obj Result object from analyzeStratifiedCTRDB
#' @return A summary list with key statistics and findings
#' @export
getStratifiedCTRDBSummary <- function(result_obj) {

  if (is.null(result_obj)) {
    stop("Result object is NULL")
  }

  summary <- list()

  # Drug B analysis summary
  if (!is.null(result_obj$drug_b_analysis)) {
    summary$drug_b_summary <- list(
      n_patients = length(result_obj$drug_b_analysis$patient_data),
      signature_genes = length(result_obj$signature_genes),
      top_genes_per_patient = sapply(result_obj$drug_b_analysis$top_genes_per_patient, length)
    )
  }

  # Score analysis summary
  if (!is.null(result_obj$drug_b_score_analysis)) {
    summary$score_summary <- result_obj$drug_b_score_analysis$summary
  }

  # Drug A analysis summary
  if (!is.null(result_obj$drug_a_analysis)) {
    summary$drug_a_summary <- list(
      n_patients = length(result_obj$drug_a_analysis$patient_data)
    )
  }

  # Drug A scores summary
  if (!is.null(result_obj$drug_a_scores)) {
    summary$drug_a_scores_summary <- list(
      n_scored_patients = length(result_obj$drug_a_scores)
    )
  }

  # Correlation results summary
  if (!is.null(result_obj$correlation_results$meta_data)) {
    meta_data <- result_obj$correlation_results$meta_data

    response_data <- meta_data[meta_data$group == "Response", ]
    non_response_data <- meta_data[meta_data$group == "Non_response", ]

    summary$correlation_summary <- list(
      n_datasets = length(unique(meta_data$patient_id)),
      n_response_correlations = nrow(response_data),
      n_non_response_correlations = nrow(non_response_data),
      mean_response_cor = mean(response_data$cor, na.rm = TRUE),
      mean_non_response_cor = mean(non_response_data$cor, na.rm = TRUE),
      median_response_cor = median(response_data$cor, na.rm = TRUE),
      median_non_response_cor = median(non_response_data$cor, na.rm = TRUE)
    )
  }

  # Meta-analysis summary
  if (!is.null(result_obj$meta_analysis)) {
    meta_result <- result_obj$meta_analysis$meta_result
    summary$meta_summary <- list(
      n_studies = meta_result$k,
      overall_effect = meta_result$TE.random,
      ci_lower = meta_result$lower.random,
      ci_upper = meta_result$upper.random,
      p_value = meta_result$pval.random,
      i_squared = meta_result$I2,
      heterogeneity_p = meta_result$pval.Q
    )
  }

  return(summary)
}

#' @examples
#' # Example usage of FuncCTRDBStratifiedAnalysis
#' \dontrun{
#' # Load required libraries
#' library(GSVA)
#' library(ggplot2)
#' library(ggpubr)
#' library(patchwork)
#' library(meta)
#' library(DROMA.Set)
#'
#' # Connect to CTRDB database
#' # ctrdb_con <- connectCTRDatabase("path/to/ctrdb.sqlite")
#' # assign("ctrdb_connection", ctrdb_con, envir = .GlobalEnv)
#'
#' # Example 1: Basic stratified analysis with specific omics feature
#' result <- analyzeStratifiedCTRDB(
#'   drug_b_name = "Cisplatin",      # Drug for signature generation
#'   drug_a_name = "Paclitaxel",     # Drug for signature application
#'   select_omics = "EGFR",          # Omics feature to analyze
#'   connection = ctrdb_con,         # Database connection
#'   top_n_genes = 100,              # Number of top genes from each dataset
#'   data_type = "all",              # Data type filter
#'   tumor_type = "all",             # Tumor type filter
#'   meta_enabled = TRUE             # Enable meta-analysis
#' )
#'
#' # View results
#' # 1. Signature genes (from Drug B analysis)
#' print(result$signature_genes)
#'
#' # 2. B drug score analysis plots
#' if (!is.null(result$drug_b_score_analysis$combined_plot)) {
#'   print(result$drug_b_score_analysis$combined_plot)
#' }
#'
#' # 3. Comprehensive correlation comparison plot across all datasets
#' if (!is.null(result$correlation_results$comprehensive_plot)) {
#'   print(result$correlation_results$comprehensive_plot)
#' }
#'
#' # 4. Individual scatter plots for each patient/group
#' if (!is.null(result$correlation_results$combined_scatter_plot)) {
#'   print(result$correlation_results$combined_scatter_plot)
#' }
#'
#' # 5. Meta-analysis forest plot
#' if (!is.null(result$meta_analysis$forest_plot)) {
#'   print(result$meta_analysis$forest_plot)
#' }
#'
#' # 6. View correlation meta data
#' print(result$correlation_results$meta_data)
#'
#' # Get comprehensive summary
#' summary <- getStratifiedCTRDBSummary(result)
#' print(summary)
#'
#' # Example 2: Analyzing different omics features
#' # Compare correlation patterns for different genes
#' genes_to_test <- c("EGFR", "ERBB2", "MET", "VEGFA")
#' results_list <- list()
#'
#' for (gene in genes_to_test) {
#'   cat("\nAnalyzing", gene, "...\n")
#'   result_gene <- analyzeStratifiedCTRDB(
#'     drug_b_name = "Cisplatin",
#'     drug_a_name = "Carboplatin",
#'     select_omics = gene,
#'     connection = ctrdb_con,
#'     data_type = "all",
#'     tumor_type = "Lung",
#'     top_n_genes = 50
#'   )
#'   results_list[[gene]] <- result_gene
#' }
#'
#' # Compare correlation strengths across genes
#' cor_comparison <- data.frame(
#'   Gene = character(),
#'   Dataset = character(),
#'   Group = character(),
#'   Correlation = numeric(),
#'   P_Value = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (gene in names(results_list)) {
#'   if (!is.null(results_list[[gene]]$correlation_results$meta_data)) {
#'     meta_data <- results_list[[gene]]$correlation_results$meta_data
#'     meta_data$Gene <- gene
#'     cor_comparison <- rbind(cor_comparison, data.frame(
#'       Gene = meta_data$Gene,
#'       Dataset = meta_data$patient_id,
#'       Group = meta_data$group,
#'       Correlation = meta_data$cor,
#'       P_Value = meta_data$p_value,
#'       stringsAsFactors = FALSE
#'     ))
#'   }
#' }
#'
#' print(cor_comparison)
#'
#' # Example 3: Using specific data type and tumor type
#' result_filtered <- analyzeStratifiedCTRDB(
#'   drug_b_name = "5-FU",
#'   drug_a_name = "Oxaliplatin",
#'   select_omics = "TYMS",          # Thymidylate synthetase
#'   connection = ctrdb_con,
#'   data_type = "PDX",              # Only PDX models
#'   tumor_type = "Colorectal",      # Only colorectal cancer
#'   min_response_size = 5,          # Require at least 5 responders
#'   min_non_response_size = 5       # Require at least 5 non-responders
#' )
#'
#' # Example 4: Extract and analyze correlation results
#' if (!is.null(result$correlation_results$meta_data)) {
#'   # Use the meta_data directly for summary
#'   cor_summary <- result$correlation_results$meta_data
#'
#'   print("Correlation Summary:")
#'   print(cor_summary)
#'
#'   # Calculate average correlations by group
#'   response_cor <- cor_summary[cor_summary$group == "Response", "cor"]
#'   non_response_cor <- cor_summary[cor_summary$group == "Non_response", "cor"]
#'
#'   cat("\nAverage correlation in Response group:",
#'       round(mean(response_cor, na.rm = TRUE), 3), "\n")
#'   cat("Average correlation in Non-response group:",
#'       round(mean(non_response_cor, na.rm = TRUE), 3), "\n")
#'
#'   # Statistical test between groups
#'   if (length(response_cor) > 0 && length(non_response_cor) > 0) {
#'     test_result <- wilcox.test(response_cor, non_response_cor)
#'     cat("Wilcoxon test p-value:", format.pval(test_result$p.value), "\n")
#'   }
#' }
#'
#' # Example 5: Analyze signature genes across datasets
#' if (!is.null(result$drug_b_analysis$top_genes_per_patient)) {
#'   # Count how many times each gene appears in top 100
#'   all_top_genes <- unlist(result$drug_b_analysis$top_genes_per_patient)
#'   gene_counts <- sort(table(all_top_genes), decreasing = TRUE)
#'
#'   # Print genes that appear in multiple datasets
#'   print("Genes appearing in multiple datasets:")
#'   print(head(gene_counts[gene_counts > 1]))
#'
#'   # Create heatmap of gene presence across datasets
#'   gene_matrix <- table(
#'     Gene = all_top_genes,
#'     Dataset = rep(names(result$drug_b_analysis$top_genes_per_patient),
#'                  times = sapply(result$drug_b_analysis$top_genes_per_patient, length))
#'   )
#'
#'   # Plot heatmap (requires pheatmap or similar)
#'   # pheatmap::pheatmap(gene_matrix[1:20, ],
#'   #                    cluster_rows = TRUE,
#'   #                    cluster_cols = TRUE,
#'   #                    main = "Top 20 Signature Genes Across Datasets")
#' }
#' }
