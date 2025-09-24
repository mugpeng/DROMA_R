# CTRDB Stratified Analysis Example
# Analysis of clinical response signatures across different drugs

# This example demonstrates how to perform stratified analysis on CTRDB clinical data
# to identify drug response signatures and apply them across different drugs.
# The core principle is to:
# 1. Generate response signatures from Drug B (e.g., Cisplatin)
# 2. Apply these signatures to Drug A (e.g., Paclitaxel)
# 3. Analyze correlation between specific omics features and signature scores
#
# This workflow is particularly useful for:
# - Identifying biomarkers that predict cross-drug response
# - Understanding molecular mechanisms underlying drug response
# - Discovering context-dependent associations between genes and drug efficacy

# Load required libraries
library(DROMA.Set)
library(DROMA.R)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(meta)
library(dplyr)

# Note: This example assumes a CTRDB database connection is available
# Uncomment and modify the following lines to connect to your CTRDB database:
#
# ctrdb_path <- "Data/ctrdb.sqlite"
# ctrdb_con <- connectCTRDatabase(ctrdb_path)
# assign("ctrdb_connection", ctrdb_con, envir = .GlobalEnv)

# Check if CTRDB connection exists
if (!exists("ctrdb_connection")) {
  warning("CTRDB database connection not found. Please connect to CTRDB database first.")
  warning("Example will use simulated data for demonstration purposes.")
  use_simulated_data <- TRUE
} else {
  use_simulated_data <- FALSE
}

# Example 1: Basic stratified analysis with specific omics feature
# ----------------------------------------------------------------
# Generate response signature from Cisplatin and apply to Paclitaxel
# Analyze correlation between EGFR expression and signature scores

cat("Example 1: Basic stratified analysis with EGFR\n")
cat("Signature from Cisplatin, Application to Paclitaxel\n\n")

if (!use_simulated_data) {
  # Perform stratified analysis with real data
  result_egfr <- tryCatch({
    analyzeStratifiedCTRDB(
      drug_b_name = "Cisplatin",      # Drug for signature generation
      drug_a_name = "Paclitaxel",     # Drug for signature application
      select_omics = "EGFR",          # Omics feature to analyze
      connection = ctrdb_connection,  # Database connection
      top_n_genes = 100,              # Number of top genes from each dataset
      data_type = "all",              # Data type filter
      tumor_type = "all",             # Tumor type filter
      min_response_size = 3,         # Minimum response samples
      min_non_response_size = 3,      # Minimum non-response samples
      meta_enabled = TRUE             # Enable meta-analysis
    )
  }, error = function(e) {
    cat("Error in real data analysis:", e$message, "\n")
    return(NULL)
  })
} else {
  # Create simulated result for demonstration
  result_egfr <- createSimulatedCTRDBResult()
  cat("Using simulated data for demonstration\n")
}

# Understanding the result structure:
# ----------------------------------
# result$drug_b_analysis     - Analysis results for Drug B (signature generation)
# result$signature_genes     - Final signature genes selected
# result$drug_b_scores       - B drug response scores for each patient
# result$drug_b_score_analysis - Score analysis results and plots
# result$drug_a_analysis     - Drug A patient data
# result$correlation_results - Correlation analysis results
# result$meta_analysis       - Meta-analysis results

# Display key results
if (!is.null(result_egfr)) {
  # Signature genes summary
  cat("Signature Generation Summary:\n")
  cat("- Drug B (signature source):", result_egfr$drug_b_analysis$drug_name %||% "Cisplatin", "\n")
  cat("- Signature genes selected:", length(result_egfr$signature_genes), "\n")

  if (!is.null(result_egfr$drug_b_analysis$top_genes_per_patient)) {
    n_patients <- length(result_egfr$drug_b_analysis$top_genes_per_patient)
    cat("- Patients analyzed for signature:", n_patients, "\n")

    # Show genes appearing in multiple datasets
    all_genes <- unlist(result_egfr$drug_b_analysis$top_genes_per_patient)
    gene_counts <- table(all_genes)
    multi_dataset_genes <- gene_counts[gene_counts > 1]
    if (length(multi_dataset_genes) > 0) {
      cat("- Genes in multiple datasets:", length(multi_dataset_genes), "\n")
      cat("  Top multi-dataset genes:", names(head(sort(multi_dataset_genes, decreasing = TRUE), 5)), "\n")
    }
  }

  # B drug score analysis
  if (!is.null(result_egfr$drug_b_score_analysis$summary)) {
    cat("\nB Drug Score Analysis:\n")
    score_summary <- result_egfr$drug_b_score_analysis$summary
    cat("- Patients with score analysis:", nrow(score_summary), "\n")

    # Show significant score differences
    sig_patients <- score_summary[score_summary$P_Value < 0.05, ]
    if (nrow(sig_patients) > 0) {
      cat("- Patients with significant score differences:", nrow(sig_patients), "\n")
      for (i in seq_len(min(3, nrow(sig_patients)))) {
        patient <- sig_patients[i, ]
        cat("  *", patient$PatientID, ": p =", round(patient$P_Value, 4), "\n")
      }
    }
  }

  # Correlation analysis results
  if (!is.null(result_egfr$correlation_results$correlation_results)) {
    cat("\nCorrelation Analysis for EGFR:\n")
    corr_results <- result_egfr$correlation_results$correlation_results

    # Create summary table
    cor_summary <- data.frame()
    for (patient_id in names(corr_results)) {
      patient_corr <- corr_results[[patient_id]]
      if ("Response" %in% names(patient_corr) && "Non_response" %in% names(patient_corr)) {
        cor_summary <- rbind(cor_summary, data.frame(
          Patient = patient_id,
          Response_Cor = patient_corr$Response$cor,
          Response_P = patient_corr$Response$p_value,
          NonResponse_Cor = patient_corr$Non_response$cor,
          NonResponse_P = patient_corr$Non_response$p_value,
          stringsAsFactors = FALSE
        ))
      }
    }

    if (nrow(cor_summary) > 0) {
      cat("- Patients with correlation data:", nrow(cor_summary), "\n")
      cat("- Average correlation (Response):", round(mean(cor_summary$Response_Cor, na.rm = TRUE), 3), "\n")
      cat("- Average correlation (Non-response):", round(mean(cor_summary$NonResponse_Cor, na.rm = TRUE), 3), "\n")

      # Show patients with strongest correlation differences
      cor_summary$cor_diff <- abs(cor_summary$Response_Cor - cor_summary$NonResponse_Cor)
      cor_summary <- cor_summary[order(-cor_summary$cor_diff), ]
      cat("- Largest correlation differences:\n")
      for (i in seq_len(min(3, nrow(cor_summary)))) {
        patient <- cor_summary[i, ]
        cat("  *", patient$PatientID, ": diff =", round(patient$cor_diff, 3), "\n")
      }
    }
  }

  # Meta-analysis results
  if (!is.null(result_egfr$meta_analysis)) {
    cat("\nMeta-Analysis Results:\n")
    meta_result <- result_egfr$meta_analysis$meta_result
    cat("- Studies included:", meta_result$k, "\n")
    cat("- Overall effect:", round(meta_result$TE.random, 3), "\n")
    cat("- 95% CI: [", round(meta_result$lower.random, 3), ",", round(meta_result$upper.random, 3), "]\n")
    cat("- P-value:", format.pval(meta_result$pval.random), "\n")
    cat("- Heterogeneity I²:", round(meta_result$I2, 1), "%\n")
  }
}

# Example 2: Multiple omics features comparison
# ---------------------------------------------
# Compare correlation patterns for different genes with the same drug pair

cat("\nExample 2: Multiple omics features comparison\n\n")

genes_to_test <- c("EGFR", "ERBB2", "MET", "VEGFA", "KRAS")
multi_gene_results <- list()

for (gene in genes_to_test) {
  cat("Analyzing", gene, "...\n")

  if (!use_simulated_data) {
    result_gene <- tryCatch({
      analyzeStratifiedCTRDB(
        drug_b_name = "Cisplatin",
        drug_a_name = "Paclitaxel",
        select_omics = gene,
        connection = ctrdb_connection,
        top_n_genes = 100,
        data_type = "all",
        tumor_type = "all",
        meta_enabled = TRUE
      )
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
      return(NULL)
    })
  } else {
    # Create varied simulated results
    result_gene <- createSimulatedCTRDBResult(gene_name = gene)
  }

  if (!is.null(result_gene)) {
    multi_gene_results[[gene]] <- result_gene

    # Extract key metrics
    if (!is.null(result_gene$correlation_results$correlation_results)) {
      cor_results <- result_gene$correlation_results$correlation_results

      # Calculate average correlations
      resp_cor <- numeric()
      non_resp_cor <- numeric()

      for (patient_id in names(cor_results)) {
        patient_corr <- cor_results[[patient_id]]
        if ("Response" %in% names(patient_corr) && "Non_response" %in% names(patient_corr)) {
          resp_cor <- c(resp_cor, patient_corr$Response$cor)
          non_resp_cor <- c(non_resp_cor, patient_corr$Non_response$cor)
        }
      }

      if (length(resp_cor) > 0 && length(non_resp_cor) > 0) {
        avg_resp <- mean(resp_cor, na.rm = TRUE)
        avg_non_resp <- mean(non_resp_cor, na.rm = TRUE)
        cat("  Response: r =", round(avg_resp, 3), "\n")
        cat("  Non-response: r =", round(avg_non_resp, 3), "\n")
        cat("  Difference:", round(abs(avg_resp - avg_non_resp), 3), "\n\n")
      }
    }
  }
}

# Create comparison table
if (length(multi_gene_results) > 0) {
  comparison_table <- data.frame(
    Gene = character(),
    N_Patients = numeric(),
    Avg_Response_Cor = numeric(),
    Avg_NonResponse_Cor = numeric(),
    Abs_Difference = numeric(),
    stringsAsFactors = FALSE
  )

  for (gene in names(multi_gene_results)) {
    result <- multi_gene_results[[gene]]

    if (!is.null(result$correlation_results$correlation_results)) {
      cor_results <- result$correlation_results$correlation_results

      resp_cor <- numeric()
      non_resp_cor <- numeric()
      n_patients <- 0

      for (patient_id in names(cor_results)) {
        patient_corr <- cor_results[[patient_id]]
        if ("Response" %in% names(patient_corr) && "Non_response" %in% names(patient_corr)) {
          resp_cor <- c(resp_cor, patient_corr$Response$cor)
          non_resp_cor <- c(non_resp_cor, patient_corr$Non_response$cor)
          n_patients <- n_patients + 1
        }
      }

      if (length(resp_cor) > 0 && length(non_resp_cor) > 0) {
        comparison_table <- rbind(comparison_table, data.frame(
          Gene = gene,
          N_Patients = n_patients,
          Avg_Response_Cor = round(mean(resp_cor, na.rm = TRUE), 3),
          Avg_NonResponse_Cor = round(mean(non_resp_cor, na.rm = TRUE), 3),
          Abs_Difference = round(abs(mean(resp_cor, na.rm = TRUE) - mean(non_resp_cor, na.rm = TRUE)), 3)
        ))
      }
    }
  }

  if (nrow(comparison_table) > 0) {
    # Sort by absolute difference
    comparison_table <- comparison_table[order(-comparison_table$Abs_Difference), ]

    cat("Multi-Gene Analysis Summary:\n")
    print(comparison_table, row.names = FALSE)

    # Save comparison table
    write.csv(comparison_table, "ctrdb_multi_gene_comparison.csv", row.names = FALSE)
    cat("\nSaved comparison table to ctrdb_multi_gene_comparison.csv\n")
  }
}

# Example 3: Different drug pairs analysis
# ---------------------------------------
# Compare signature generation and application across different drug combinations

cat("\nExample 3: Different drug pairs analysis\n\n")

drug_pairs <- list(
  list(drug_b = "Cisplatin", drug_a = "Paclitaxel"),
  list(drug_b = "Carboplatin", drug_a = "Gemcitabine"),
  list(drug_b = "5-FU", drug_a = "Oxaliplatin")
)

drug_pair_results <- list()

for (pair in drug_pairs) {
  cat(sprintf("Analyzing %s -> %s...\n", pair$drug_b, pair$drug_a))

  if (!use_simulated_data) {
    result_pair <- tryCatch({
      analyzeStratifiedCTRDB(
        drug_b_name = pair$drug_b,
        drug_a_name = pair$drug_a,
        select_omics = "EGFR",
        connection = ctrdb_connection,
        top_n_genes = 100,
        data_type = "all",
        tumor_type = "all",
        meta_enabled = TRUE
      )
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
      return(NULL)
    })
  } else {
    # Create varied simulated results for different drug pairs
    result_pair <- createSimulatedCTRDBResult(
      drug_b_name = pair$drug_b,
      drug_a_name = pair$drug_a
    )
  }

  if (!is.null(result_pair)) {
    drug_pair_results[[paste(pair$drug_b, pair$drug_a, sep = "_")]] <- result_pair

    # Report signature gene overlap
    if (!is.null(result_pair$signature_genes)) {
      cat("  Signature genes:", length(result_pair$signature_genes), "\n")

      # Check overlap with known drug response genes
      known_genes <- c("ERCC1", "ERCC2", "XPA", "XPC", "BRCA1", "BRCA2", "TP53")
      overlap <- intersect(result_pair$signature_genes, known_genes)
      if (length(overlap) > 0) {
        cat("  Known drug response genes in signature:", paste(overlap, collapse = ", "), "\n")
      }
    }

    # Report correlation strength
    if (!is.null(result_pair$correlation_results$correlation_results)) {
      cor_results <- result_pair$correlation_results$correlation_results

      resp_cor <- numeric()
      non_resp_cor <- numeric()

      for (patient_id in names(cor_results)) {
        patient_corr <- cor_results[[patient_id]]
        if ("Response" %in% names(patient_corr) && "Non_response" %in% names(patient_corr)) {
          resp_cor <- c(resp_cor, patient_corr$Response$cor)
          non_resp_cor <- c(non_resp_cor, patient_corr$Non_response$cor)
        }
      }

      if (length(resp_cor) > 0 && length(non_resp_cor) > 0) {
        avg_resp <- mean(resp_cor, na.rm = TRUE)
        avg_non_resp <- mean(non_resp_cor, na.rm = TRUE)
        cat("  Avg correlation difference:", round(abs(avg_resp - avg_non_resp), 3), "\n\n")
      }
    }
  }
}

# Example 4: Visualization of results
# ------------------------------------
# Create comprehensive visualizations for stratified analysis

cat("Example 4: Creating visualizations\n\n")

if (!is.null(result_egfr)) {
  # 1. B drug score analysis plots
  if (!is.null(result_egfr$drug_b_score_analysis$combined_plot)) {
    ggsave("ctrdb_b_drug_scores.png", result_egfr$drug_b_score_analysis$combined_plot,
           width = 12, height = 8, dpi = 300)
    cat("Saved B drug score analysis plot to ctrdb_b_drug_scores.png\n")
  }

  # 2. Correlation analysis plots
  if (!is.null(result_egfr$correlation_results$combined_plot)) {
    ggsave("ctrdb_correlation_analysis.png", result_egfr$correlation_results$combined_plot,
           width = 12, height = 8, dpi = 300)
    cat("Saved correlation analysis plot to ctrdb_correlation_analysis.png\n")
  }

  # 3. Meta-analysis forest plot
  if (!is.null(result_egfr$meta_analysis$forest_plot)) {
    ggsave("ctrdb_meta_forest.png", result_egfr$meta_analysis$forest_plot,
           width = 10, height = 8, dpi = 300)
    cat("Saved meta-analysis forest plot to ctrdb_meta_forest.png\n")
  }

  # 4. Individual scatter plots
  if (!is.null(result_egfr$correlation_results$plots)) {
    # Save individual scatter plots
    scatter_plots <- result_egfr$correlation_results$plots
    scatter_plots <- scatter_plots[!grepl("_comparison$", names(scatter_plots))]

    if (length(scatter_plots) > 0) {
      # Combine scatter plots
      if (requireNamespace("patchwork", quietly = TRUE)) {
        combined_scatter <- patchwork::wrap_plots(scatter_plots, ncol = min(3, length(scatter_plots)))
        combined_scatter <- combined_scatter +
          patchwork::plot_annotation(
            title = "EGFR Expression vs B Drug Signature Scores",
            theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
          )

        ggsave("ctrdb_scatter_plots.png", combined_scatter,
               width = 15, height = 10, dpi = 300)
        cat("Saved combined scatter plots to ctrdb_scatter_plots.png\n")
      }
    }
  }
}

# Example 5: Signature gene analysis
# ---------------------------------
# Analyze signature genes across different datasets

cat("\nExample 5: Signature gene analysis\n")

if (!is.null(result_egfr$drug_b_analysis$top_genes_per_patient)) {
  # Analyze gene frequency across datasets
  all_top_genes <- unlist(result_egfr$drug_b_analysis$top_genes_per_patient)
  gene_counts <- sort(table(all_top_genes), decreasing = TRUE)

  cat("Top 20 most frequent signature genes:\n")
  print(head(gene_counts, 20))

  # Create gene presence matrix
  gene_matrix <- table(
    Gene = all_top_genes,
    Dataset = rep(names(result_egfr$drug_b_analysis$top_genes_per_patient),
                 times = sapply(result_egfr$drug_b_analysis$top_genes_per_patient, length))
  )

  # Function significance analysis
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cat("\nPerforming gene ontology enrichment analysis...\n")

    # Convert gene symbols to Entrez IDs
    signature_genes <- result_egfr$signature_genes
    entrez_ids <- clusterProfiler::bitr(signature_genes,
                                       fromType = "SYMBOL",
                                       toType = "ENTREZID",
                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db)

    if (nrow(entrez_ids) > 0) {
      # Perform GO enrichment
      go_result <- clusterProfiler::enrichGO(
        gene = entrez_ids$ENTREZID,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = "BP",  # Biological Process
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )

      if (!is.null(go_result) && nrow(go_result) > 0) {
        cat("Top enriched GO terms:\n")
        print(head(go_result[, c("Description", "GeneRatio", "p.adjust", "qvalue")], 10))

        # Save GO results
        write.csv(as.data.frame(go_result), "ctrdb_signature_go_enrichment.csv", row.names = FALSE)
        cat("Saved GO enrichment results to ctrdb_signature_go_enrichment.csv\n")
      }
    }
  }
}

# Example 6: Summary and statistics
# ---------------------------------
# Generate comprehensive summary statistics

cat("\nExample 6: Comprehensive analysis summary\n\n")

if (!is.null(result_egfr)) {
  # Get structured summary
  summary <- getStratifiedCTRDBSummary(result_egfr)

  cat("Analysis Summary:\n")
  cat("================\n")

  # Drug B analysis
  if (!is.null(summary$drug_b_summary)) {
    cat("\nDrug B (Signature Generation):\n")
    cat("- Patients analyzed:", summary$drug_b_summary$n_patients, "\n")
    cat("- Signature genes:", summary$drug_b_summary$signature_genes, "\n")
    cat("- Genes per patient (range):",
        min(summary$drug_b_summary$top_genes_per_patient), "-",
        max(summary$drug_b_summary$top_genes_per_patient), "\n")
  }

  # Score analysis
  if (!is.null(summary$score_summary)) {
    cat("\nB Drug Score Analysis:\n")
    score_summary <- summary$score_summary
    cat("- Total analyses:", nrow(score_summary), "\n")

    sig_analyses <- score_summary[score_summary$P_Value < 0.05, ]
    cat("- Significant score differences:", nrow(sig_analyses), "\n")

    if (nrow(sig_analyses) > 0) {
      cat("- Average score difference (Non-response - Response):",
          round(mean(sig_analyses$Mean_Non_Response - sig_analyses$Mean_Response, na.rm = TRUE), 3), "\n")
    }
  }

  # Drug A analysis
  if (!is.null(summary$drug_a_summary)) {
    cat("\nDrug A (Signature Application):\n")
    cat("- Patients analyzed:", summary$drug_a_summary$n_patients, "\n")
  }

  # Meta-analysis
  if (!is.null(summary$meta_summary)) {
    cat("\nMeta-Analysis:\n")
    meta_summary <- summary$meta_summary
    cat("- Studies combined:", meta_summary$n_studies, "\n")
    cat("- Overall effect size:", round(meta_summary$overall_effect, 3), "\n")
    cat("- 95% Confidence interval: [", round(meta_summary$ci_lower, 3),
        ", ", round(meta_summary$ci_upper, 3), "]\n")
    cat("- P-value:", format.pval(meta_summary$p_value), "\n")
    cat("- Heterogeneity (I²):", round(meta_summary$i_squared, 1), "%\n")
    cat("- Heterogeneity p-value:", format.pval(meta_summary$heterogeneity_p), "\n")
  }

  # Save summary
  if (!is.null(summary)) {
    # Convert to data frame for saving
    summary_df <- data.frame(
      Metric = character(),
      Value = character(),
      stringsAsFactors = FALSE
    )

    # Add all summary metrics
    if (!is.null(summary$drug_b_summary)) {
      summary_df <- rbind(summary_df, data.frame(
        Metric = c("Drug B Patients", "Signature Genes"),
        Value = c(as.character(summary$drug_b_summary$n_patients),
                  as.character(summary$drug_b_summary$signature_genes))
      ))
    }

    if (!is.null(summary$meta_summary)) {
      summary_df <- rbind(summary_df, data.frame(
        Metric = c("Meta Studies", "Overall Effect", "P-value", "I²"),
        Value = c(as.character(summary$meta_summary$n_studies),
                  as.character(round(summary$meta_summary$overall_effect, 3)),
                  format.pval(summary$meta_summary$p_value),
                  paste0(round(summary$meta_summary$i_squared, 1), "%"))
      ))
    }

    write.csv(summary_df, "ctrdb_analysis_summary.csv", row.names = FALSE)
    cat("\nSaved analysis summary to ctrdb_analysis_summary.csv\n")
  }
}

# Summary
# -------
# This example demonstrates the power of stratified analysis for CTRDB clinical data:
#
# 1. Cross-drug signature identification: Generate response signatures from one drug
#    and apply them to understand mechanisms in another drug
# 2. Context-dependent biomarker discovery: Identify biomarkers that show different
#    correlation patterns in responders vs non-responders
# 3. Meta-analysis across patients: Combine evidence from multiple patients to
#    increase statistical power and robustness
# 4. Comprehensive visualization: Create plots for score analysis, correlations,
#    and meta-analysis results
# 5. Biological interpretation: Perform pathway enrichment analysis on signature genes
#
# Key methodological features:
# - Signature generation using differential expression analysis
# - Score calculation using ssGSEA (single-sample Gene Set Enrichment Analysis)
# - Correlation analysis between omics features and signature scores
# - Meta-analysis using random-effects models
# - Comprehensive result structure with organized outputs
#
# Potential applications:
# - Biomarker discovery for combination therapies
# - Understanding resistance mechanisms
# - Patient stratification for clinical trials
# - Drug repurposing studies

cat("\nCTRDB stratified analysis examples completed!\n")
cat("For more information on the functions, see ?analyzeStratifiedCTRDB\n")

# Helper function for creating simulated results (for demonstration only)
createSimulatedCTRDBResult <- function(gene_name = "EGFR",
                                       drug_b_name = "Cisplatin",
                                       drug_a_name = "Paclitaxel") {

  # Simulate signature genes
  set.seed(Sys.getpid())
  signature_genes <- sample(c("ERCC1", "ERCC2", "XPA", "XPC", "BRCA1", "BRCA2",
                            "TP53", "EGFR", "ERBB2", "MET", "VEGFA", "KRAS",
                            "TYMS", "MTHFR", "DPYD", "UGT1A1", "GSTP1"), 100, replace = TRUE)
  signature_genes <- unique(signature_genes)[1:min(100, length(unique(signature_genes)))]

  # Simulate top genes per patient
  n_patients <- sample(3:8, 1)
  top_genes_per_patient <- list()
  for (i in 1:n_patients) {
    patient_genes <- sample(signature_genes, sample(50:100, 1), replace = FALSE)
    top_genes_per_patient[[paste0("Patient_", i)]] <- patient_genes
  }

  # Simulate drug B scores
  drug_b_scores <- list()
  for (patient_id in names(top_genes_per_patient)) {
    n_samples <- sample(10:30, 1)
    scores <- rnorm(n_samples, mean = 0, sd = 1)
    response_samples <- sample(1:n_samples, floor(n_samples/2))
    non_response_samples <- setdiff(1:n_samples, response_samples)

    drug_b_scores[[patient_id]] <- list(
      scores = scores,
      response_samples = response_samples,
      non_response_samples = non_response_samples
    )
  }

  # Simulate drug B score analysis
  score_summary <- data.frame()
  for (patient_id in names(drug_b_scores)) {
    score_data <- drug_b_scores[[patient_id]]
    resp_scores <- score_data$scores[score_data$response_samples]
    non_resp_scores <- score_data$scores[score_data$non_response_samples]

    score_summary <- rbind(score_summary, data.frame(
      PatientID = patient_id,
      N_Response = length(resp_scores),
      N_Non_Response = length(non_resp_scores),
      Mean_Response = mean(resp_scores),
      Mean_Non_Response = mean(non_resp_scores),
      P_Value = wilcox.test(non_resp_scores, resp_scores)$p.value
    ))
  }

  # Simulate correlation results
  correlation_results <- list()
  correlation_data <- list()

  for (patient_id in names(top_genes_per_patient)) {
    n_samples <- sample(15:25, 1)

    # Simulate gene expression
    gene_expr <- rnorm(n_samples, mean = 10, sd = 2)

    # Simulate correlations
    resp_cor <- rnorm(1, mean = 0.3, sd = 0.2)
    non_resp_cor <- rnorm(1, mean = -0.1, sd = 0.2)

    correlation_results[[patient_id]] <- list(
      Response = data.frame(
        feature = gene_name,
        cor = resp_cor,
        p_value = runif(1, 0.001, 0.1),
        n_samples = floor(n_samples/2)
      ),
      Non_response = data.frame(
        feature = gene_name,
        cor = non_resp_cor,
        p_value = runif(1, 0.001, 0.1),
        n_samples = ceil(n_samples/2)
      )
    )

    correlation_data[[patient_id]] <- data.frame(
      study = patient_id,
      cor_resp = resp_cor,
      cor_non_resp = non_resp_cor,
      diff_cor = resp_cor - non_resp_cor,
      p_value = runif(1, 0.01, 0.5),
      n_resp = floor(n_samples/2),
      n_non_resp = ceil(n_samples/2)
    )
  }

  # Simulate meta-analysis
  if (length(correlation_data) >= 2) {
    meta_df <- do.call(rbind, correlation_data)
    meta_df$se_diff <- sqrt(
      (1 - meta_df$cor_resp^2)^2 / (meta_df$n_resp - 3) +
      (1 - meta_df$cor_non_resp^2)^2 / (meta_df$n_non_resp - 3)
    )

    meta_result <- list(
      k = nrow(meta_df),
      TE.random = mean(meta_df$diff_cor),
      lower.random = mean(meta_df$diff_cor) - 1.96 * sd(meta_df$diff_cor),
      upper.random = mean(meta_df$diff_cor) + 1.96 * sd(meta_df$diff_cor),
      pval.random = 0.05,
      I2 = 50,
      pval.Q = 0.1
    )
  } else {
    meta_result <- NULL
  }

  # Return comprehensive result structure
  list(
    drug_b_analysis = list(
      drug_name = drug_b_name,
      patient_data = list(),  # Simplified for simulation
      top_genes_per_patient = top_genes_per_patient
    ),
    signature_genes = signature_genes,
    drug_b_scores = drug_b_scores,
    drug_b_score_analysis = list(
      summary = score_summary,
      plots = list(),
      combined_plot = NULL
    ),
    drug_a_analysis = list(
      drug_name = drug_a_name,
      patient_data = list()  # Simplified for simulation
    ),
    correlation_results = list(
      correlation_results = correlation_results,
      correlation_data = correlation_data,
      plots = list(),
      combined_plot = NULL
    ),
    meta_analysis = if (!is.null(meta_result)) {
      list(
        meta_result = meta_result,
        forest_plot = NULL
      )
    } else {
      NULL
    }
  )
}

# Helper function for safe null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x