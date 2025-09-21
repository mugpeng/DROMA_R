# Stratified Drug-Omics Analysis Example
# Analysis of bortezomib response stratified by cisplatin sensitivity

# This example demonstrates how to perform stratified analysis to uncover
# context-dependent associations between molecular features and drug response.
# The core principle is that biological associations may differ between
# resistant and sensitive subgroups.

# Load required libraries
library(DROMA.Set)
library(DROMA.R)
library(ggplot2)
library(patchwork)

# Set database path
db_path <- "path/to/droma.sqlite"

# Create MultiDromaSet from multiple projects
# -------------------------------------------------
# This combines data from multiple sources to increase statistical power
# while accounting for study-specific effects through meta-analysis

cat("Creating MultiDromaSet from multiple projects...\n")
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE", "GDSC"),
  db_path = db_path
)

# Example 1: Stratified analysis of ERCC1 and bortezomib
# ------------------------------------------------------
# ERCC1 is involved in DNA repair and may affect response to
# DNA-damaging agents. Stratifying by cisplatin response may
# reveal context-dependent associations with bortezomib.

cat("\nExample 1: ERCC1 expression vs bortezomib response\n")
cat("Stratified by cisplatin sensitivity\n\n")

# Perform stratified analysis
result_ercc1 <- analyzeStratifiedDrugOmic(
  dromaset_object = multi_set,
  stratification_drug = "cisplatin",      # Drug for stratification
  strata_quantile = 0.33,                 # Use tertiles
  select_omics_type = "mRNA",             # Gene expression
  select_omics = "ERCC1",                 # Target gene
  select_drugs = "bortezomib",            # Target drug
  data_type = "CellLine",                 # Filter by data type
  tumor_type = "all",                     # All tumor types
  overlap_only = TRUE,                    # Use overlapping samples
  merged_enabled = TRUE,                  # Create merged dataset
  meta_enabled = TRUE                     # Perform meta-analysis
)

# Display stratification information
cat("Stratification Summary:\n")
cat("- Stratification drug:", result_ercc1$stratification_info$drug, "\n")
cat("- Quantile threshold:", result_ercc1$stratification_info$quantile, "\n")
cat("- Lower threshold:", round(result_ercc1$stratification_info$lower_threshold, 2), "\n")
cat("- Upper threshold:", round(result_ercc1$stratification_info$upper_threshold, 2), "\n")

# View sensitive group results
if (!is.null(result_ercc1$sensitive$meta)) {
  cat("\nSensitive Group (n =", result_ercc1$stratification_info$n_sensitive, "):\n")
  cat("Correlation with bortezomib:", round(result_ercc1$sensitive$meta$TE.random, 3),
      "(p =", round(result_ercc1$sensitive$meta$pval.random, 4), ")\n")
}

# View resistant group results
if (!is.null(result_ercc1$resistant$meta)) {
  cat("\nResistant Group (n =", result_ercc1$stratification_info$n_resistant, "):\n")
  cat("Correlation with bortezomib:", round(result_ercc1$resistant$meta$TE.random, 3),
      "(p =", round(result_ercc1$resistant$meta$pval.random, 4), ")\n")
}

# Display comparison
if (!is.null(result_ercc1$comparison$correlation_comparison)) {
  cat("\nBetween-Group Comparison:\n")
  print(result_ercc1$comparison$correlation_comparison)
}

# Example 2: Multiple genes analysis
# ----------------------------------
# Analyze multiple DNA repair genes in the context of cisplatin stratification

dna_repair_genes <- c("ERCC1", "ERCC2", "XPA", "XPC", "BRCA1")
stratified_results <- list()

cat("\nExample 2: Multiple DNA repair genes analysis\n\n")

for (gene in dna_repair_genes) {
  cat("Analyzing", gene, "...\n")

  result <- tryCatch({
    analyzeStratifiedDrugOmic(
      dromaset_object = multi_set,
      stratification_drug = "cisplatin",
      strata_quantile = 0.33,
      select_omics_type = "mRNA",
      select_omics = gene,
      select_drugs = "bortezomib",
      data_type = "CellLine",
      tumor_type = "all",
      overlap_only = TRUE,
      merged_enabled = FALSE,  # Skip merged for speed
      meta_enabled = TRUE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(result)) {
    stratified_results[[gene]] <- result

    # Extract key results
    if (!is.null(result$sensitive$meta) && !is.null(result$resistant$meta)) {
      sens_cor <- result$sensitive$meta$TE.random
      sens_p <- result$sensitive$meta$pval.random
      res_cor <- result$resistant$meta$TE.random
      res_p <- result$resistant$meta$pval.random

      cat("  Sensitive: r =", round(sens_cor, 3), "(p =", round(sens_p, 4), ")\n")
      cat("  Resistant: r =", round(res_cor, 3), "(p =", round(res_p, 4), ")\n")
      cat("  Difference:", round(abs(sens_cor - res_cor), 3), "\n\n")
    }
  }
}

# Example 3: Different stratification quantiles
# ----------------------------------------------
# Explore how different stratification thresholds affect results

cat("Example 3: Effect of stratification quantile\n")
cat("Analyzing ERCC1 with different quantile thresholds\n\n")

quantiles_to_test <- c(0.25, 0.33, 0.40)  # Quartiles, tertiles, and custom
quantile_results <- list()

for (q in quantiles_to_test) {
  cat("Testing quantile:", q, "\n")

  result <- tryCatch({
    analyzeStratifiedDrugOmic(
      dromaset_object = multi_set,
      stratification_drug = "cisplatin",
      strata_quantile = q,
      select_omics_type = "mRNA",
      select_omics = "ERCC1",
      select_drugs = "bortezomib",
      data_type = "CellLine",
      tumor_type = "all",
      overlap_only = TRUE,
      merged_enabled = FALSE,
      meta_enabled = TRUE
    )
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(result)) {
    quantile_results[[paste0("q", q)]] <- result

    # Report sample sizes
    info <- result$stratification_info
    cat("  Samples per stratum:", info$n_sensitive, "/", info$n_resistant, "\n")

    # Report correlations if available
    if (!is.null(result$sensitive$meta) && !is.null(result$resistant$meta)) {
      cat("  Sensitive r:", round(result$sensitive$meta$TE.random, 3), "\n")
      cat("  Resistant r:", round(result$resistant$meta$TE.random, 3), "\n\n")
    }
  }
}

# Example 4: Visualization of stratified results
# ----------------------------------------------
# Create comprehensive visualizations for stratified analysis

cat("Example 4: Creating visualizations\n\n")

if (!is.null(result_ercc1$sensitive$plot) && !is.null(result_ercc1$resistant$plot)) {
  # Combine plots from both strata
  combined_plot <- (result_ercc1$sensitive$plot / result_ercc1$resistant$plot) +
    plot_annotation(
      title = "ERCC1 Expression vs Bortezomib Response\nStratified by Cisplatin Sensitivity",
      subtitle = c("Sensitive Group", "Resistant Group")
    )

  # Save the plot
  ggsave("stratified_ercc1_bortezomib.png", combined_plot,
         width = 10, height = 12, dpi = 300)

  cat("Saved combined plot to stratified_ercc1_bortezomib.png\n")
}

# Create forest plots for meta-analysis results
if (!is.null(result_ercc1$sensitive$meta) && !is.null(result_ercc1$resistant$meta)) {
  # Forest plot for sensitive group
  forest_sensitive <- createForestPlot(
    result_ercc1$sensitive$meta,
    xlab = "Correlation with Bortezomib IC50 (Sensitive Group)"
  )

  # Forest plot for resistant group
  forest_resistant <- createForestPlot(
    result_ercc1$resistant$meta,
    xlab = "Correlation with Bortezomib IC50 (Resistant Group)"
  )

  # Combine forest plots
  forest_combined <- forest_sensitive + forest_resistant +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Meta-analysis Results by Cisplatin Sensitivity Stratum"
    )

  ggsave("stratified_forest_plots.png", forest_combined,
         width = 12, height = 8, dpi = 300)

  cat("Saved forest plots to stratified_forest_plots.png\n")
}

# Example 5: Batch stratified analysis
# ------------------------------------
# Perform stratified analysis for multiple drug-omics pairs

cat("\nExample 5: Batch stratified analysis\n")

# Define analysis parameters
analysis_pairs <- expand.grid(
  omics_gene = c("ERCC1", "BRCA1", "TP53"),
  target_drug = c("bortezomib", "carfilzomib"),
  stringsAsFactors = FALSE
)

batch_results <- list()

for (i in 1:nrow(analysis_pairs)) {
  gene <- analysis_pairs$omics_gene[i]
  drug <- analysis_pairs$target_drug[i]

  cat(sprintf("Analyzing %s vs %s... ", gene, drug))

  result <- tryCatch({
    analyzeStratifiedDrugOmic(
      dromaset_object = multi_set,
      stratification_drug = "cisplatin",
      strata_quantile = 0.33,
      select_omics_type = "mRNA",
      select_omics = gene,
      select_drugs = drug,
      data_type = "CellLine",
      tumor_type = "all",
      overlap_only = TRUE,
      merged_enabled = FALSE,
      meta_enabled = TRUE
    )
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(result)) {
    batch_results[[paste(gene, drug, sep = "_")]] <- result
    cat("Done\n")
  }
}

# Create summary table of batch results
if (length(batch_results) > 0) {
  summary_table <- data.frame(
    Gene = character(),
    Drug = character(),
    Sens_Cor = numeric(),
    Sens_P = numeric(),
    Res_Cor = numeric(),
    Res_P = numeric(),
    Abs_Diff = numeric(),
    stringsAsFactors = FALSE
  )

  for (result_name in names(batch_results)) {
    result <- batch_results[[result_name]]

    if (!is.null(result$sensitive$meta) && !is.null(result$resistant$meta)) {
      parts <- strsplit(result_name, "_")[[1]]
      gene <- parts[1]
      drug <- paste(parts[-1], collapse = "_")

      sens_cor <- result$sensitive$meta$TE.random
      sens_p <- result$sensitive$meta$pval.random
      res_cor <- result$resistant$meta$TE.random
      res_p <- result$resistant$meta$pval.random

      summary_table <- rbind(summary_table, data.frame(
        Gene = gene,
        Drug = drug,
        Sens_Cor = round(sens_cor, 3),
        Sens_P = format.pval(sens_p, digits = 2),
        Res_Cor = round(res_cor, 3),
        Res_P = format.pval(res_p, digits = 2),
        Abs_Diff = round(abs(sens_cor - res_cor), 3),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Sort by absolute difference
  summary_table <- summary_table[order(-summary_table$Abs_Diff), ]

  cat("\nBatch Analysis Summary (sorted by stratum difference):\n")
  print(summary_table, row.names = FALSE)

  # Save summary table
  write.csv(summary_table, "stratified_analysis_summary.csv", row.names = FALSE)
  cat("\nSaved summary table to stratified_analysis_summary.csv\n")
}

# Summary
# -------
# This example demonstrates the power of stratified analysis for uncovering
# context-dependent drug-omics associations. Key insights:
#
# 1. Stratification can reveal associations masked in pooled analysis
# 2. The choice of stratification threshold affects results
# 3. Different molecular features may show varying degrees of context-dependence
# 4. Meta-analysis across studies increases robustness of findings

cat("\nStratified analysis examples completed!\n")