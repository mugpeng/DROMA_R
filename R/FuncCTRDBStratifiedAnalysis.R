# Stratified Analysis Functions for CTRDB Data ----

# Required packages: GSVA, ggplot2, ggpubr, patchwork, meta, dplyr

# Core principle: Perform stratified analysis on CTRDB clinical data
# to identify drug response signatures and apply them across different drugs

#' Stratified analysis for CTRDB clinical data
#'
#' @description Performs stratified analysis on CTRDB data to identify drug response signatures
#' and analyze the correlation between specific omics features and the signature scores
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
#' # View correlation plots for EGFR
#' result$correlation_plots
#'
#' # View meta-analysis forest plot
#' result$meta_forest_plot
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
    patient_data_list = drug_b_result$patient_data,
    signature_genes = final_signature_genes
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
    result$correlation_results <- NULL
    result$meta_analysis <- NULL
  } else {
    result$drug_a_analysis <- drug_a_result

    # Step 6: Analyze correlation between specific omics feature and B drug scores in Drug A
    cat("Step 6: Analyzing correlations for", select_omics, "in Drug A datasets...\n")
    correlation_results <- analyzeCorrelationsInDrugA(
      drug_a_data = drug_a_result$patient_data,
      select_omics = select_omics,
      drug_b_scores = drug_b_scores
    )
    result$correlation_results <- correlation_results

    # Step 7: Perform meta-analysis if enabled
    if (meta_enabled && length(correlation_results$correlation_data) > 1) {
      cat("Step 7: Performing meta-analysis...\n")
      meta_result <- performMetaAnalysisForCorrelations(correlation_results$correlation_data)
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
      top_genes <- de_results$gene[1:min(top_n_genes, nrow(de_results))]

      # Store data
      patient_data_list[[patient_id]] <- list(
        response_expr = response_expr,
        non_response_expr = non_response_expr,
        response_samples = response_samples,
        non_response_samples = non_response_samples,
        metadata = patient_meta
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
  final_genes <- gene_ranking$gene[1:min(100, nrow(gene_ranking))]

  message("Selected ", length(final_genes), " signature genes")
  message("Genes appear in ", round(mean(gene_ranking$count[1:length(final_genes)]), 1),
          " datasets on average")

  return(final_genes)
}

#' Calculate B drug response scores using ssGSEA
#'
#' @description Calculates drug response scores using ssGSEA on signature genes
#' @param patient_data_list List of patient expression data
#' @param signature_genes Vector of signature genes
#' @return List of scores for each patient
#' @keywords internal
calculateDrugBScores <- function(patient_data_list, signature_genes) {

  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required for ssGSEA calculation")
  }

  scores_list <- list()

  for (patient_id in names(patient_data_list)) {
    patient_data <- patient_data_list[[patient_id]]

    # Combine response and non-response expression data
    all_expr <- cbind(patient_data$response_expr, patient_data$non_response_expr)

    # Filter to signature genes
    signature_expr <- all_expr[rownames(all_expr) %in% signature_genes, ]

    if (nrow(signature_expr) == 0) {
      warning("No signature genes found in patient ", patient_id)
      next
    }

    # Create gene set list
    gene_sets <- list(signature = signature_genes)

    # Calculate ssGSEA scores
    tryCatch({
      # Create ssGSEA parameter object for new GSVA API
      ssgsea_param <- GSVA::ssgseaParam(
        exprData = as.matrix(signature_expr),
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
      labs(title = paste("B Drug Response Scores -", patient_id),
           x = "", y = "ssGSEA Score") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    plots[[patient_id]] <- p
  }

  # Combine plots
  if (length(plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = paste("B Drug (", drug_name, ") Response Scores Analysis", sep = ""),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
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
#' @description Analyzes correlation between specific omics feature expression and B drug scores
#' @param drug_a_data Drug A patient data
#' @param select_omics Name of the omics feature to analyze
#' @param drug_b_scores B drug scores for reference
#' @return List containing correlation results and plots
#' @keywords internal
analyzeCorrelationsInDrugA <- function(drug_a_data, select_omics, drug_b_scores) {

  correlation_results <- list()
  correlation_data <- list()
  plots <- list()

  # Calculate reference B drug score once (average across all datasets)
  all_b_scores <- unlist(lapply(drug_b_scores, function(x) x$scores))
  all_b_scores <- all_b_scores[!is.na(all_b_scores)]

  if (length(all_b_scores) == 0) {
    ref_b_score <- 0
    warning("No valid B drug scores found across all datasets")
  } else {
    ref_b_score <- mean(all_b_scores)
  }

  for (patient_id in names(drug_a_data)) {
    patient_data <- drug_a_data[[patient_id]]

    # Get expression data for the selected omics feature
    if (!select_omics %in% rownames(patient_data$expr_data)) {
      warning("Feature ", select_omics, " not found in patient ", patient_id)
      next
    }

    omics_expr <- patient_data$expr_data[select_omics, , drop = FALSE]

    # Calculate correlations for each response group
    results_list <- list()

    for (group in c("Response", "Non_response")) {
      group_samples <- if (group == "Response") {
        patient_data$response_samples
      } else {
        patient_data$non_response_samples
      }

      if (length(group_samples) < 3) next

      # Get expression values for this group
      group_expr <- as.numeric(omics_expr[, group_samples])
      names(group_expr) <- group_samples

      # Remove NA values
      group_expr <- group_expr[!is.na(group_expr)]

      # Check if we have enough finite observations
      if (length(group_expr) < 3) {
        warning("Insufficient finite observations for correlation in ", patient_id, " ", group)
        next
      }

      # Since we're correlating with a constant, correlation is always 0
      cor_estimate <- 0
      cor_pvalue <- 1

      # Additional check: if variance of group_expr is 0, correlation is undefined
      if (var(group_expr) == 0) {
        warning("Zero variance in expression values for ", select_omics, " in ", patient_id, " ", group)
      }

      results_list[[group]] <- data.frame(
        feature = select_omics,
        cor = cor_estimate,
        p_value = cor_pvalue,
        n_samples = length(group_expr),
        stringsAsFactors = FALSE
      )
    }

    correlation_results[[patient_id]] <- results_list

    # Create plots
    if ("Response" %in% names(results_list) && "Non_response" %in% names(results_list)) {
      resp_data <- results_list$Response
      non_resp_data <- results_list$Non_response

      # Create scatter plots for each group
      # Get expression values without NA
      resp_expr <- as.numeric(patient_data$expr_data[select_omics, patient_data$response_samples])
      resp_expr <- resp_expr[!is.na(resp_expr)]
      names(resp_expr) <- patient_data$response_samples[!is.na(as.numeric(patient_data$expr_data[select_omics, patient_data$response_samples]))]

      non_resp_expr <- as.numeric(patient_data$expr_data[select_omics, patient_data$non_response_samples])
      non_resp_expr <- non_resp_expr[!is.na(non_resp_expr)]
      names(non_resp_expr) <- patient_data$non_response_samples[!is.na(as.numeric(patient_data$expr_data[select_omics, patient_data$non_response_samples]))]

      # Create stratified score vectors with proper names
      if (length(resp_expr) > 0) {
        resp_scores <- rep(ref_b_score, length(resp_expr))
        names(resp_scores) <- names(resp_expr)

        plots[[paste0(patient_id, "_Response")]] <- plotCorrelationWithStratifiedScore(
          gene_expr = resp_expr,
          stratified_score = resp_scores,
          patient_id = paste(patient_id, "Response"),
          group = "Response"
        )
      }

      if (length(non_resp_expr) > 0) {
        non_resp_scores <- rep(ref_b_score, length(non_resp_expr))
        names(non_resp_scores) <- names(non_resp_expr)

        plots[[paste0(patient_id, "_Non_response")]] <- plotCorrelationWithStratifiedScore(
          gene_expr = non_resp_expr,
          stratified_score = non_resp_scores,
          patient_id = paste(patient_id, "Non_response"),
          group = "Non_response"
        )
      }

      # Create comparison boxplot
      plot_df <- data.frame(
        cor = c(resp_data$cor, non_resp_data$cor),
        group = c("Response", "Non_response")
      )

      p_compare <- ggpubr::ggboxplot(plot_df, x = "group", y = "cor",
                                    fill = "group", palette = c("#FB8072FF", "#BEBADAFF"),
                                    add = "jitter", add.params = list(width = 0.1, alpha = 0.5)) +
        ggpubr::stat_compare_means(method = "wilcox.test",
                                   label.x = 0.8,
                                   label.y = ifelse(is.na(max(plot_df$cor)), 0, max(plot_df$cor) * 0.9)) +
        labs(title = paste(select_omics, "Correlation with B Drug Score -", patient_id),
             x = "", y = "Spearman Correlation") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))

      plots[[paste0(patient_id, "_comparison")]] <- p_compare

      # Prepare data for meta-analysis
      correlation_data[[patient_id]] <- data.frame(
        study = patient_id,
        cor_resp = resp_data$cor,
        cor_non_resp = non_resp_data$cor,
        diff_cor = resp_data$cor - non_resp_data$cor,
        p_value = wilcox.test(resp_data$cor, non_resp_data$cor)$p.value,
        n_resp = resp_data$n_samples,
        n_non_resp = non_resp_data$n_samples
      )
    }
  }

  # Combine comparison plots
  comparison_plots <- plots[grepl("_comparison$", names(plots))]
  if (length(comparison_plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(comparison_plots, ncol = min(3, length(comparison_plots)))
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = paste("Correlation of", select_omics, "with B Drug Response Signature"),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
      )
  } else {
    combined_plot <- NULL
  }

  return(list(
    correlation_results = correlation_results,
    correlation_data = correlation_data,
    plots = plots,
    combined_plot = combined_plot,
    comparison_plots = comparison_plots
  ))
}

#' Perform meta-analysis for correlation differences
#'
#' @description Performs meta-analysis on correlation differences between groups
#' @param correlation_data List of correlation data from each study
#' @return Meta-analysis result
#' @keywords internal
performMetaAnalysisForCorrelations <- function(correlation_data) {

  if (length(correlation_data) < 2) {
    warning("Insufficient studies for meta-analysis")
    return(NULL)
  }

  # Combine data
  meta_df <- do.call(rbind, correlation_data)

  # Calculate standard error for correlation difference
  # Using Fisher's z-transformation approximation
  meta_df$se_diff <- sqrt(
    (1 - meta_df$cor_resp^2)^2 / (meta_df$n_resp - 3) +
    (1 - meta_df$cor_non_resp^2)^2 / (meta_df$n_non_resp - 3)
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
      forest_plot = forest_plot
    ))

  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(NULL)
  })
}

#' Plot correlation with stratified score for single patient
#'
#' @description Creates a scatter plot showing correlation between gene expression
#' and stratified score for a single patient
#' @param gene_expr Gene expression values
#' @param stratified_score Stratified score values
#' @param patient_id Patient identifier
#' @param group Response group ("Response" or "Non_response")
#' @return A ggplot2 scatter plot
#' @export
plotCorrelationWithStratifiedScore <- function(gene_expr, stratified_score,
                                              patient_id, group) {

  # Ensure data alignment
  common_samples <- intersect(names(gene_expr), names(stratified_score))
  if (length(common_samples) == 0) {
    stop("No common samples between gene expression and stratified score")
  }

  plot_data <- data.frame(
    expression = gene_expr[common_samples],
    score = stratified_score[common_samples]
  )

  # Calculate correlation
  cor_test <- cor.test(plot_data$expression, plot_data$score, method = "spearman")
  cor_text <- paste0("rho = ", round(cor_test$estimate, 3),
                     "\np = ", format.pval(cor_test$p.value))

  # Create plot
  ggplot2::ggplot(plot_data, aes(x = expression, y = score)) +
    ggplot2::geom_point(alpha = 0.6, size = 3) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue") +
    ggplot2::labs(title = paste(patient_id, "-", group),
         x = "Gene Expression",
         y = "Stratified Score") +
    ggplot2::annotate("text", x = Inf, y = Inf,
             label = cor_text,
             hjust = 1.1, vjust = 1.1,
             size = 4, family = "mono") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
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
#' # 3. Correlation analysis plots for EGFR
#' if (!is.null(result$correlation_results$combined_plot)) {
#'   print(result$correlation_results$combined_plot)
#' }
#'
#' # 4. Individual scatter plots for each patient/group
#' # View scatter plot for first patient's response group
#' patient_names <- names(result$drug_a_analysis$patient_data)
#' if (length(patient_names) > 0) {
#'   resp_plot_name <- paste0(patient_names[1], "_Response")
#'   if (resp_plot_name %in% names(result$correlation_results$plots)) {
#'     print(result$correlation_results$plots[[resp_plot_name]])
#'   }
#' }
#'
#' # 5. Meta-analysis forest plot
#' if (!is.null(result$meta_analysis$forest_plot)) {
#'   print(result$meta_analysis$forest_plot)
#' }
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
#'   Response_Cor = numeric(),
#'   NonResponse_Cor = numeric(),
#'   stringsAsFactors = FALSE
#' )
#'
#' for (gene in names(results_list)) {
#'   if (!is.null(results_list[[gene]]$correlation_results)) {
#'     for (dataset in names(results_list[[gene]]$correlation_results$correlation_results)) {
#'       corr_data <- results_list[[gene]]$correlation_results$correlation_results[[dataset]]
#'       if ("Response" %in% names(corr_data) && "Non_response" %in% names(corr_data)) {
#'         cor_comparison <- rbind(cor_comparison, data.frame(
#'           Gene = gene,
#'           Dataset = dataset,
#'           Response_Cor = corr_data$Response$cor,
#'           NonResponse_Cor = corr_data$Non_response$cor
#'         ))
#'       }
#'     }
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
#' if (!is.null(result$correlation_results$correlation_results)) {
#'   # Create summary table of correlations
#'   cor_summary <- data.frame()
#'
#'   for (patient_id in names(result$correlation_results$correlation_results)) {
#'     patient_corr <- result$correlation_results$correlation_results[[patient_id]]
#'
#'     if ("Response" %in% names(patient_corr) && "Non_response" %in% names(patient_corr)) {
#'       cor_summary <- rbind(cor_summary, data.frame(
#'         Patient = patient_id,
#'         Response_Cor = patient_corr$Response$cor,
#'         Response_P = patient_corr$Response$p_value,
#'         NonResponse_Cor = patient_corr$Non_response$cor,
#'         NonResponse_P = patient_corr$Non_response$p_value,
#'         stringsAsFactors = FALSE
#'       ))
#'     }
#'   }
#'
#'   print("Correlation Summary:")
#'   print(cor_summary)
#'
#'   # Calculate average correlations
#'   cat("\nAverage correlation in Response group:",
#'       round(mean(cor_summary$Response_Cor, na.rm = TRUE), 3), "\n")
#'   cat("Average correlation in Non-response group:",
#'       round(mean(cor_summary$NonResponse_Cor, na.rm = TRUE), 3), "\n")
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