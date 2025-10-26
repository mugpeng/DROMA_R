#!/usr/bin/env Rscript

# Example script for DROMA.R package: Batch Feature Analysis with DromaSet Objects
# This example demonstrates how to perform batch analysis to find features associated with drug response using DromaSet objects

# Load required libraries
library(meta)
library(metafor)
library(effsize)
library(ggpubr)
library(dplyr)
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions

# Setup: Create DromaSet Objects----

# Note: Replace with your actual database path
db_path <- "../Data/droma.sqlite"

# Connect to DROMA database
connectDROMADatabase(db_path)

# Create DromaSet objects
cat("Creating DromaSet objects...\n")
gCSI <- createDromaSetFromDatabase("gCSI", db_path)

# Create a MultiDromaSet for cross-project analysis
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE"),
  db_path = db_path
)

cat("DromaSet objects created successfully!\n\n")

# Example 1: Find genes associated with Paclitaxel response using single project ----

cat("Example 1: Finding genes associated with Paclitaxel response using single DromaSet...\n")
# Note: Setting test_top_n=100 to limit analysis time for the example
# In real usage, you might want to use test_top_n=NULL to test all features
batch_results_single <- batchFindSignificantFeatures(
  gCSI,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "cnv",
  data_type = "all",
  tumor_type = "all",
  test_top_n = 900  # Use only top 100 features for this example
)

# Display top features
if (!is.null(batch_results_single) && nrow(batch_results_single) > 0) {
  cat("Top genes associated with Paclitaxel response (single project):\n")
  # Sort by p-value
  sorted_results <- batch_results_single[order(batch_results_single$p_value),]
  print(head(sorted_results, 10))

  # Create a volcano plot of results
  cat("\nCreating volcano plot for single project results...\n")
  volcano_plot_single <- plotMetaVolcano(batch_results_single, es_t = 0.3, P_t = 0.05)
  print(volcano_plot_single)
} else {
  cat("No significant results found for single project analysis\n")
}

# Example 2: Find genes associated with Paclitaxel response using multiple projects ----

cat("\nExample 2: Finding genes associated with Paclitaxel response using MultiDromaSet...\n")
batch_results_multi <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "cnv",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  test_top_n = 100  # Use only top 100 features for this example
)

# Display top features
if (!is.null(batch_results_multi) && nrow(batch_results_multi) > 0) {
  cat("Top genes associated with Paclitaxel response (multi-project meta-analysis):\n")
  # Sort by p-value
  sorted_results_multi <- batch_results_multi[order(batch_results_multi$p_value),]
  print(head(sorted_results_multi, 10))

  # Create a volcano plot of results
  cat("\nCreating volcano plot for multi-project results...\n")
  volcano_plot_multi <- plotMetaVolcano(batch_results_multi, es_t = 0.3, P_t = 0.05,
                                       title = "Multi-Project Meta-Analysis")
  print(volcano_plot_multi)
} else {
  cat("No significant results found for multi-project analysis\n")
}

######################################
# Example 3: Find mutations associated with Paclitaxel response
######################################

cat("\nExample 3: Finding mutations associated with Paclitaxel response...\n")
mutation_results <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mutation_gene",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  test_top_n = 100  # Use only top 100 features for this example
)

# Display top mutations
if (!is.null(mutation_results) && nrow(mutation_results) > 0) {
  cat("Top mutations associated with Paclitaxel response:\n")
  # Sort by p-value
  sorted_mutations <- mutation_results[order(mutation_results$p_value),]
  print(head(sorted_mutations, 10))

  # Create a volcano plot for mutations
  cat("\nCreating volcano plot for mutation results...\n")
  volcano_plot_mut <- plotMetaVolcano(mutation_results, es_t = 0.3, P_t = 0.05,
                                     title = "Mutations Associated with Paclitaxel Response")
  print(volcano_plot_mut)
} else {
  cat("No significant mutation associations found\n")
}

######################################
# Example 4: Compare results between single and multi-project analyses
######################################

if (!is.null(batch_results_single) && !is.null(batch_results_multi) &&
    nrow(batch_results_single) > 0 && nrow(batch_results_multi) > 0) {

  cat("\nExample 4: Comparing single vs multi-project results...\n")

  # Find common genes
  common_genes <- intersect(batch_results_single$name, batch_results_multi$name)

  if (length(common_genes) > 0) {
    cat("Found", length(common_genes), "genes analyzed in both approaches\n")

    # Compare p-values for common genes
    single_pvals <- batch_results_single$p_value[match(common_genes, batch_results_single$name)]
    multi_pvals <- batch_results_multi$p_value[match(common_genes, batch_results_multi$name)]

    # Create comparison data frame
    comparison_df <- data.frame(
      gene = common_genes,
      single_project_pval = single_pvals,
      multi_project_pval = multi_pvals,
      stringsAsFactors = FALSE
    )

    # Show genes with improved significance in multi-project analysis
    improved_genes <- comparison_df[comparison_df$multi_project_pval < comparison_df$single_project_pval, ]
    if (nrow(improved_genes) > 0) {
      cat("Genes with improved significance in multi-project analysis:\n")
      print(head(improved_genes[order(improved_genes$multi_project_pval), ], 5))
    }
  } else {
    cat("No common genes found between single and multi-project analyses\n")
  }
}

######################################
# Example 5: Test specific features using feature2_name parameter
######################################

cat("\nExample 5: Testing specific genes instead of all features...\n")
cat("Using feature2_name parameter to test only specific genes of interest\n")

# Define specific genes to test
genes_of_interest <- c("TP53", "EGFR", "KRAS", "BRCA1", "BRCA2")

specific_gene_results <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mRNA",
  # feature2_name = genes_of_interest,  # Test only these specific genes
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

# Display results for specific genes
if (!is.null(specific_gene_results) && nrow(specific_gene_results) > 0) {
  cat("Results for specific genes of interest:\n")
  sorted_specific <- specific_gene_results[order(specific_gene_results$p_value),]
  print(sorted_specific)
} else {
  cat("No results found for the specified genes\n")
}

######################################
# Example 6: Parallel processing example
######################################

cat("\nExample 6: Using parallel processing for faster analysis...\n")
cat("Note: This example shows how to use multiple cores for batch analysis\n")

# Example with parallel processing (uncomment to use)
# parallel_results <- batchFindSignificantFeatures(
#   multi_set,
#   feature1_type = "drug",
#   feature1_name = "Paclitaxel",
#   feature2_type = "mRNA",
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   cores = 4,  # Use 4 CPU cores
#   test_top_n = 100
# )

cat("To use parallel processing, set cores > 1 in batchFindSignificantFeatures()\n")

cat("\nBatch feature analysis examples completed!\n")
cat("Key takeaways:\n")
cat("1. Use batchFindSignificantFeatures() with DromaSet objects for batch analysis\n")
cat("2. MultiDromaSet enables meta-analysis across multiple projects\n")
cat("3. Use plotMetaVolcano() to visualize batch analysis results\n")
cat("4. Set test_top_n to limit features for faster testing, or NULL for comprehensive analysis\n")
cat("5. Use parallel processing (cores > 1) for large-scale analyses\n")
cat("6. Compare single vs multi-project results to assess consistency\n")
cat("7. Use feature2_name parameter to test specific features instead of all features\n")
