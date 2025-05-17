#!/usr/bin/env Rscript

# Example script for DROMA package: Batch Feature Analysis
# This example demonstrates how to perform batch analysis to find features associated with drug response

# Load required libraries
library(meta)
library(metafor)
library(effsize)
library(ggplot2)
library(dplyr)
library(DROMA)

# Setup DROMA environment - load search vectors and sample annotations
setupDROMA()

# Load required data types
loadDROMA("drug")   # Load drug data and apply normalization
loadDROMA("mRNA")   # Load mRNA data and apply normalization
loadDROMA("mut")    # Load mutation data

# Example 1: Find genes associated with Paclitaxel response
cat("Finding genes associated with Paclitaxel response...\n")
# Note: Setting test_top_100=TRUE to limit analysis time for the example
# In real usage, you might want to use test_top_100=FALSE to test all features
batch_results <- batchFindSignificantFeatures(
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mRNA",
  data_type = "all",
  tumor_type = "all",
  test_top_100 = TRUE  # Use only top 100 features for this example
)

# Display top features
cat("Top genes associated with Paclitaxel response:\n")
# Sort by p-value
sorted_results <- batch_results[order(batch_results$p_value),]
print(head(sorted_results, 10))

# Example 2: Create a volcano plot of results
cat("\nCreating volcano plot...\n")
volcano_plot <- plotMetaVolcano(batch_results, es_t = 0.3, P_t = 0.05)

# Example 3: Find mutations associated with Paclitaxel response
cat("\nFinding mutations associated with Paclitaxel response...\n")
mutation_results <- batchFindSignificantFeatures(
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mutation_gene",
  data_type = "all",
  tumor_type = "all",
  test_top_100 = TRUE  # Use only top 100 features for this example
)

# Display top mutations
cat("Top mutations associated with Paclitaxel response:\n")
# Sort by p-value
sorted_mutations <- mutation_results[order(mutation_results$p_value),]
print(head(sorted_mutations, 10))
