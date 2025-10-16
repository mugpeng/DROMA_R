# Stratified Drug-Omics Analysis Example
# Analysis of Bortezomib response stratified by Cisplatin sensitivity

# This example demonstrates how to perform stratified analysis to uncover
# context-dependent associations between molecular features and drug response.
# The core principle is that biological associations may differ between
# resistant and sensitive subgroups.
#
# Updated for FuncStratifiedAnalysis.R v3.0
# Key changes:
# - Removed deprecated functions (loadStratifiedDataByProject)
# - Removed unused stratified pairing functions
# - Removed backward compatibility wrapper (loadAndFilterStratifiedData)
# - Implemented createStratifiedComparisonPlot for better visualization
# - Cleaner, more focused API with essential functions only

# Load required libraries
library(DROMA.Set)
library(DROMA.R)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(meta)


# Set database path
db_path <- "../Data/droma.sqlite"
connectDROMADatabase(db_path)

# Create MultiDromaSet from multiple projects
# -------------------------------------------------
# This combines data from multiple sources to increase statistical power
# while accounting for study-specific effects through meta-analysis

cat("Creating MultiDromaSet from multiple projects...\n")
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE", "GDSC"),
  db_path = db_path
)

project_anno <- DROMA.Set::listDROMAProjects()
cellline_names <- project_anno[project_anno$dataset_type %in% "CellLine",]$project_name

cellline_sets <- createMultiDromaSetFromDatabase(db_path = db_path,
                                                 project_names = cellline_names)
cellline_sets2 <- removeDromaSetFromMulti(cellline_sets, "NCI60")
cellline_sets <- cellline_sets2

# Example 1: Stratified analysis of ERCC1 and Bortezomib
# ------------------------------------------------------
# ERCC1 is involved in DNA repair and may affect response to
# DNA-damaging agents. Stratifying by Cisplatin response may
# reveal context-dependent associations with Bortezomib.

cat("\nExample 1: ERCC1 expression vs Bortezomib response\n")
cat("Stratified by Cisplatin sensitivity\n\n")

# Perform stratified analysis
result_ercc1 <- analyzeStratifiedDrugOmic(
  dromaset_object = cellline_sets,
  stratification_drug = "Cisplatin",      # Drug for stratification
  strata_quantile = 0.33,                 # Use tertiles
  select_omics_type = "mRNA",             # Gene expression
  select_omics = "ERCC1",                 # Target gene
  select_drugs = "Bortezomib",            # Target drug
  data_type = "all",                 # Filter by data type
  tumor_type = "all",                     # All tumor types
  overlap_only = TRUE,                    # Use overlapping samples
  merged_enabled = TRUE,                  # Create merged dataset
  meta_enabled = TRUE                     # Perform meta-analysis
)

# Understanding the result structure:
# ----------------------------------
# result$sensitive        - Analysis results for sensitive samples
# result$resistant        - Analysis results for resistant samples
# result$comparison       - Between-group statistical comparisons
# result$stratification_info - Project-specific sample assignments:
#   - Each project contains $sensitive and $resistant sample names
#   - Each project contains $thresholds with lower/upper cutoffs
# result$analysis_parameters - All analysis settings and parameters

# Display stratification information
cat("Stratification Summary:\n")
cat("- Stratification drug:", result_ercc1$analysis_parameters$stratification_drug, "\n")
cat("- Quantile threshold:", result_ercc1$analysis_parameters$strata_quantile, "\n")

# Display project-specific stratification details
for (project_name in names(result_ercc1$stratification_info)) {
  project_info <- result_ercc1$stratification_info[[project_name]]
  cat("- Project", project_name, ":\n")
  cat("  * Sensitive samples:", length(project_info$sensitive), "\n")
  cat("  * Resistant samples:", length(project_info$resistant), "\n")
  cat("  * Lower threshold:", round(project_info$thresholds["lower"], 2), "\n")
  cat("  * Upper threshold:", round(project_info$thresholds["upper"], 2), "\n")
}

# Calculate total sample sizes across projects
total_sensitive <- sum(sapply(result_ercc1$stratification_info, function(x) length(x$sensitive)))
total_resistant <- sum(sapply(result_ercc1$stratification_info, function(x) length(x$resistant)))

# View sensitive group results
if (!is.null(result_ercc1$sensitive$meta)) {
  cat("\nSensitive Group (n =", total_sensitive, "):\n")
  cat("Correlation with Bortezomib:", round(result_ercc1$sensitive$meta$TE.random, 3),
      "(p =", round(result_ercc1$sensitive$meta$pval.random, 4), ")\n")
}

# View resistant group results
if (!is.null(result_ercc1$resistant$meta)) {
  cat("\nResistant Group (n =", total_resistant, "):\n")
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
# Analyze multiple DNA repair genes in the context of Cisplatin stratification

dna_repair_genes <- c("ERCC1", "ERCC2", "XPA", "XPC", "BRCA1")
stratified_results <- list()

cat("\nExample 2: Multiple DNA repair genes analysis\n\n")

for (gene in dna_repair_genes) {
  cat("Analyzing", gene, "...\n")

  result <- tryCatch({
    analyzeStratifiedDrugOmic(
      dromaset_object = multi_set,
      stratification_drug = "Cisplatin",
      strata_quantile = 0.33,
      select_omics_type = "mRNA",
      select_omics = gene,
      select_drugs = "Bortezomib",
      data_type = "all",
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
      stratification_drug = "Cisplatin",
      strata_quantile = q,
      select_omics_type = "mRNA",
      select_omics = "ERCC1",
      select_drugs = "Bortezomib",
      data_type = "all",
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

    # Report sample sizes across projects
    total_sens <- sum(sapply(result$stratification_info, function(x) length(x$sensitive)))
    total_res <- sum(sapply(result$stratification_info, function(x) length(x$resistant)))
    cat("  Samples per stratum:", total_sens, "/", total_res, "\n")

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
# The new createStratifiedComparisonPlot function provides integrated visualization

cat("Example 4: Creating visualizations with createStratifiedComparisonPlot\n\n")

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

# Create comparison plots using the new createStratifiedComparisonPlot function
if (!is.null(result_ercc1$comparison$comparison_plot)) {
  # The comparison plot can be:
  # 1. A forest plot showing correlations/effect sizes for both strata
  # 2. Side-by-side individual plots from each stratum

  # Save the comparison plot
  ggsave("stratified_comparison_plot.png", result_ercc1$comparison$comparison_plot,
         width = 10, height = 8, dpi = 300)

  cat("Saved stratified comparison plot to stratified_comparison_plot.png\n")

  # If the result is a list of individual plots (when patchwork is not available)
  if (is.list(result_ercc1$comparison$comparison_plot)) {
    cat("Individual plots saved as list elements\n")
  }
}

# Example 5: Direct use of createStratifiedComparisonPlot
# ------------------------------------------------------
# This example shows how to directly call createStratifiedComparisonPlot
# with results from stratified analysis

cat("\nExample 5: Direct use of createStratifiedComparisonPlot function\n")

# First, let's create a simple stratified analysis result
cat("Performing stratified analysis for BRCA1 vs Cisplatin...\n")

result_brca1 <- tryCatch({
  analyzeStratifiedDrugOmic(
    dromaset_object = multi_set,
    stratification_drug = "Cisplatin",
    strata_quantile = 0.33,
    select_omics_type = "mRNA",
    select_omics = "BRCA1",
    select_drugs = "Cisplatin",
    data_type = "all",
    tumor_type = "all",
    overlap_only = TRUE,
    meta_enabled = TRUE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result_brca1)) {
  # Now directly use createStratifiedComparisonPlot
  cat("Creating comparison plot directly using createStratifiedComparisonPlot...\n")

  comparison_plot <- createStratifiedComparisonPlot(
    sensitive_result = result_brca1$sensitive,
    resistant_result = result_brca1$resistant,
    select_omics = "BRCA1",
    select_drugs = "Cisplatin",
    select_omics_type = "mRNA"
  )

  if (!is.null(comparison_plot)) {
    # Save the direct comparison plot
    if (is.list(comparison_plot)) {
      # If it returns a list (e.g., when patchwork is not available)
      cat("Comparison plot returned as a list with elements:\n")
      print(names(comparison_plot))

      # Save individual plots if available
      if ("forest_plot" %in% names(comparison_plot)) {
        ggsave("brca1_cisplatin_forest_plot.png",
               comparison_plot$forest_plot,
               width = 8, height = 6, dpi = 300)
        cat("Saved forest plot to brca1_cisplatin_forest_plot.png\n")
      }
    } else {
      # Single plot object
      ggsave("brca1_cisplatin_direct_comparison.png",
             comparison_plot,
             width = 10, height = 8, dpi = 300)
      cat("Saved direct comparison plot to brca1_cisplatin_direct_comparison.png\n")
    }
  }
}

# Example 5b: Using createStratifiedComparisonPlot with discrete data
# ----------------------------------------------------------------
cat("\nExample 5b: createStratifiedComparisonPlot with discrete data\n")

# Analyze mutation data (discrete)
result_tp53_mut <- tryCatch({
  analyzeStratifiedDrugOmic(
    dromaset_object = multi_set,
    stratification_drug = "Cisplatin",
    strata_quantile = 0.33,
    select_omics_type = "mutation",
    select_omics = "TP53",
    select_drugs = "Bortezomib",
    data_type = "CellLine",
    tumor_type = "all",
    overlap_only = TRUE,
    meta_enabled = TRUE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result_tp53_mut)) {
  # Create comparison plot for discrete data
  discrete_comparison <- createStratifiedComparisonPlot(
    sensitive_result = result_tp53_mut$sensitive,
    resistant_result = result_tp53_mut$resistant,
    select_omics = "TP53",
    select_drugs = "Bortezomib",
    select_omics_type = "mutation"
  )

  if (!is.null(discrete_comparison)) {
    ggsave("tp53_mutation_discrete_comparison.png",
           discrete_comparison,
           width = 10, height = 8, dpi = 300)
    cat("Saved discrete data comparison plot to tp53_mutation_discrete_comparison.png\n")
  }
}

# Example 6: Batch stratified analysis
# ------------------------------------
# Perform stratified analysis for multiple drug-omics pairs

cat("\nExample 6: Batch stratified analysis\n")

# Define analysis parameters
analysis_pairs <- expand.grid(
  omics_gene = c("ERCC1", "BRCA1", "TP53"),
  target_drug = c("Bortezomib", "Carfilzomib"),
  stringsAsFactors = FALSE
)

batch_results <- list()

for (i in seq_len(nrow(analysis_pairs))) {
  gene <- analysis_pairs$omics_gene[i]
  drug <- analysis_pairs$target_drug[i]

  cat(sprintf("Analyzing %s vs %s... ", gene, drug))

  result <- tryCatch({
    analyzeStratifiedDrugOmic(
      dromaset_object = multi_set,
      stratification_drug = "Cisplatin",
      strata_quantile = 0.33,
      select_omics_type = "mRNA",
      select_omics = gene,
      select_drugs = drug,
      data_type = "all",
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
# 5. The updated API provides project-specific stratification information
# 6. Result structure is more organized with separate analysis parameters
# 7. The createStratifiedComparisonPlot function provides integrated visualization

# Key methodological principles:
# - Follow pairing logic: filter samples within each analysis operation
# - Use quantile-based stratification for robust sample assignment
# - Leverage meta-analysis across projects for increased statistical power
# - Always validate sufficient sample sizes in each stratum before analysis
# - Use createStratifiedComparisonPlot for effective visualization of strata differences
# - The cleaned API removes deprecated functions for better maintainability

# Functions removed in v3.0:
# - loadStratifiedDataByProject (deprecated)
# - pairContinuousFeaturesStratified (unused)
# - pairDiscreteFeaturesStratified (unused)
# - loadAndFilterStratifiedData (backward compatibility wrapper)

# Functions enhanced in v3.0:
# - createStratifiedComparisonPlot: Now fully implemented with support for:
#   * Forest plots for correlation/effect size comparisons
#   * Side-by-side individual stratum plots
#   * Both continuous and discrete data types
#   * Automatic plot combination using patchwork (when available)

cat("\nStratified analysis examples completed!\n")
cat("For more information on the updated API, see ?analyzeStratifiedDrugOmic\n")
