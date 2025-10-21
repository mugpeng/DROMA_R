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
#' # View forest plot with meta-analyzed correlations
#' result$correlation_results$forest_plot
#'
#' # View combined scatter plots
#' result$correlation_results$combined_scatter_plot
#'
#' # View correlation data
#' result$correlation_results$meta_data
#' }
analyzeStratifiedCTRDB <- function(drug_b_name,
                                  drug_a_name,
                                  select_omics,
                                  connection,
                                  top_n_genes = 100,
                                  data_type = "all",
                                  tumor_type = "all",
                                  min_response_size = 3,
                                  min_non_response_size = 3) {

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
      drug_a_name = drug_a_name,
      drug_b_name = drug_b_name
    )
    result$correlation_results <- correlation_results
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

      # Get expression data using the reusable function from DROMA package
      expr_data <- getPatientExpressionData(patient_id, connection, auto_log = FALSE)

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

      # Get expression data using the reusable function from DROMA package
      expr_data <- getPatientExpressionData(patient_id, connection, auto_log = FALSE)

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

#' Create forest plot for correlation meta-analysis
#'
#' @description Creates a forest plot comparing meta-analyzed correlations between Response and Non_response groups
#' @param meta_data Data frame with correlation results (patient_id, group, cor, p_value, n_samples)
#' @param select_omics Name of the omics feature
#' @param drug_a_name Name of drug A
#' @param drug_b_name Name of drug B
#' @param p_value_digits Number of decimal places for p-values (default: 3)
#' @param add_difference_row Whether to add difference row (default: TRUE)
#' @return ggplot object or NULL
#' @keywords internal
createCorrelationForestPlot <- function(meta_data, select_omics, drug_a_name,
                                       drug_b_name = NULL, p_value_digits = 3,
                                       add_difference_row = TRUE) {

  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for plotting")
    return(NULL)
  }

  if (!requireNamespace("meta", quietly = TRUE)) {
    warning("meta package is required for meta-analysis")
    return(NULL)
  }

  # Separate data by group
  response_data <- meta_data[meta_data$group == "Response", ]
  non_response_data <- meta_data[meta_data$group == "Non_response", ]

  # Check if we have data for both groups
  if (nrow(response_data) < 2 || nrow(non_response_data) < 2) {
    warning("Insufficient data for meta-analysis in one or both groups")
    return(NULL)
  }

  # Perform meta-analysis for Response group
  response_meta <- performCorrelationMetaAnalysis(
    correlations = response_data$cor,
    sample_sizes = response_data$n_samples,
    study_labels = response_data$patient_id
  )

  # Perform meta-analysis for Non_response group
  non_response_meta <- performCorrelationMetaAnalysis(
    correlations = non_response_data$cor,
    sample_sizes = non_response_data$n_samples,
    study_labels = non_response_data$patient_id
  )

  # Check if meta-analyses were successful
  if (is.null(response_meta) || is.null(non_response_meta)) {
    warning("Meta-analysis failed for one or both groups")
    return(NULL)
  }

  # Extract meta-analysis results
  response_cor <- response_meta$TE.random
  response_lower <- response_meta$lower.random
  response_upper <- response_meta$upper.random
  response_p <- response_meta$pval.random
  response_se <- response_meta$seTE.random

  non_response_cor <- non_response_meta$TE.random
  non_response_lower <- non_response_meta$lower.random
  non_response_upper <- non_response_meta$upper.random
  non_response_p <- non_response_meta$pval.random
  non_response_se <- non_response_meta$seTE.random

  # Calculate difference and p-value
  if (add_difference_row && !is.na(response_cor) && !is.na(non_response_cor)) {
    z_diff <- (response_cor - non_response_cor) / sqrt(response_se^2 + non_response_se^2)
    p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)
  } else {
    p_diff <- NA
  }

  # Format p-value for display
  format_pvalue <- function(p, digits = 3) {
    if (is.na(p) || is.null(p)) return("NA")
    if (p < 0.001) return("<0.001")
    if (p < 0.01) return(sprintf("<%.2f", p))
    if (p < 0.05) return(sprintf("<%.3f", p))
    return(sprintf("%.3f", round(p, digits)))
  }

  # Create plot data
  groups <- c("Response", "Non_response")
  correlations <- c(response_cor, non_response_cor)
  ci_lowers <- c(response_lower, non_response_lower)
  ci_uppers <- c(response_upper, non_response_upper)
  p_values <- c(response_p, non_response_p)

  # Add difference row if requested
  if (add_difference_row && !is.na(p_diff)) {
    groups <- c(groups, "Difference")
    correlations <- c(correlations, response_cor - non_response_cor)
    ci_lowers <- c(ci_lowers, NA)
    ci_uppers <- c(ci_uppers, NA)
    p_values <- c(p_values, p_diff)
  }

  # Create group labels with p-values using newline
  p_labels <- sapply(p_values, format_pvalue, digits = p_value_digits)
  group_labels <- paste0(groups, "\np=", p_labels)

  plot_data <- data.frame(
    Group = group_labels,
    Group_Name = groups,
    Correlation = correlations,
    CI_lower = ci_lowers,
    CI_upper = ci_uppers,
    P_value = p_values,
    stringsAsFactors = FALSE
  )

  # Create title
  if (!is.null(drug_b_name)) {
    plot_title <- paste0(select_omics, " vs ", drug_a_name, " (using ", drug_b_name, " signature)")
  } else {
    plot_title <- paste0(select_omics, " vs ", drug_a_name, " Response Signature")
  }

  # Create forest plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "Group", y = "Correlation")) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      data = plot_data[!is.na(plot_data$CI_lower), ],
      ggplot2::aes_string(ymin = "CI_lower", ymax = "CI_upper"),
      width = 0.2
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::coord_flip(ylim = range(c(plot_data$Correlation, plot_data$CI_lower, plot_data$CI_upper), na.rm = TRUE) * c(1.1, 1.1)) +
    ggplot2::labs(
      title = plot_title,
      subtitle = "Meta-analyzed Correlation Coefficients with 95% CI",
      y = "Correlation Coefficient (Spearman's rho)",
      x = ""
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
      axis.text.y = ggplot2::element_text(hjust = 1, face = "bold", lineheight = 0.8, size = 11),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      panel.grid.major.y = ggplot2::element_blank()
    )

  return(p)
}

#' Perform meta-analysis for correlations using Fisher's z transformation
#'
#' @description Performs meta-analysis on correlation coefficients using Fisher's z transformation
#' @param correlations Vector of correlation coefficients
#' @param sample_sizes Vector of sample sizes
#' @param study_labels Vector of study labels
#' @return Meta-analysis result from meta::metagen or NULL
#' @keywords internal
performCorrelationMetaAnalysis <- function(correlations, sample_sizes, study_labels) {

  if (!requireNamespace("meta", quietly = TRUE)) {
    warning("meta package is required for meta-analysis")
    return(NULL)
  }

  # Fisher's z transformation
  z_values <- 0.5 * log((1 + correlations) / (1 - correlations))

  # Standard error of z (using n-3 for correlation)
  se_z <- 1 / sqrt(pmax(sample_sizes - 3, 1))

  # Perform meta-analysis
  tryCatch({
    meta_result <- meta::metagen(
      TE = z_values,
      seTE = se_z,
      data = data.frame(z = z_values, se = se_z),
      sm = "COR",
      studlab = study_labels,
      comb.fixed = FALSE,
      comb.random = TRUE,
      method.tau = "REML",
      hakn = TRUE
    )

    return(meta_result)

  }, error = function(e) {
    warning("Meta-analysis failed: ", e$message)
    return(NULL)
  })
}

#' Analyze correlations in Drug A datasets
#'
#' @description Analyzes correlation between specific omics feature expression and A drug scores
#' @param drug_a_data Drug A patient data
#' @param drug_a_scores Drug A scores calculated using B drug signature
#' @param select_omics Name of the omics feature to analyze
#' @param drug_a_name Name of drug A
#' @param drug_b_name Name of drug B (for plot titles)
#' @return List containing correlation results, forest plot, and combined scatter plot
#' @keywords internal
analyzeCorrelationsInDrugA <- function(drug_a_data, drug_a_scores, select_omics, drug_a_name, drug_b_name = NULL) {

  # Data frame to store all correlation results for meta-analysis
  meta_data <- data.frame()

  # List to temporarily store scatter plots (not returned)
  temp_scatter_plots <- list()

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

      # Create scatter plot (store temporarily, not in final output)
      plot_data <- data.frame(
        expression = group_expr,
        score = group_scores
      )

      cor_text <- paste0("rho = ", round(cor_test$estimate, 3),
                        "\np = ", format.pval(cor_test$p.value, digits = 2))

      p_scatter <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "expression", y = "score")) +
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

      temp_scatter_plots[[paste0(patient_id, "_", group)]] <- p_scatter
    }
  }

  # Create forest plot with meta-analysis for each group
  forest_plot <- NULL
  if (nrow(meta_data) > 0) {
    forest_plot <- createCorrelationForestPlot(
      meta_data = meta_data,
      select_omics = select_omics,
      drug_a_name = drug_a_name,
      drug_b_name = drug_b_name
    )
  }

  # Combine scatter plots (keep combined version only)
  combined_scatter_plot <- NULL
  if (length(temp_scatter_plots) > 0 && requireNamespace("patchwork", quietly = TRUE)) {
    combined_scatter_plot <- patchwork::wrap_plots(
      temp_scatter_plots,
      ncol = min(4, ceiling(sqrt(length(temp_scatter_plots))))
    ) +
      patchwork::plot_annotation(
        title = paste("Individual Dataset Correlations:", select_omics, "vs Drug Response Score"),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
      )
  }

  return(list(
    meta_data = meta_data,
    forest_plot = forest_plot,
    combined_scatter_plot = combined_scatter_plot
  ))
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
#'   tumor_type = "all"              # Tumor type filter
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
#' # 3. Forest plot showing meta-analyzed correlations for Response and Non_response groups
#' if (!is.null(result$correlation_results$forest_plot)) {
#'   print(result$correlation_results$forest_plot)
#' }
#'
#' # 4. Combined scatter plots for individual datasets
#' if (!is.null(result$correlation_results$combined_scatter_plot)) {
#'   print(result$correlation_results$combined_scatter_plot)
#' }
#'
#' # 5. View correlation meta data
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
