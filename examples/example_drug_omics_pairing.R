#!/usr/bin/env Rscript

# Example script for DROMA package: Drug-Omics Pairing Analysis
# This example demonstrates how to analyze associations between drug sensitivity and omics features

# Load required libraries
library(DROMA)
library(ggplot2)
library(metafor)
library(meta)
library(effsize)
library(ggpubr)
library(patchwork)

# Setup DROMA environment - load search vectors and sample annotations
setupDROMA()

# Load required data types
loadDROMA("drug")   # Load drug data and apply normalization
loadDROMA("mRNA")   # Load mRNA data and apply normalization
loadDROMA("mut")    # Load mutation data

# Now we can use the search vectors to find available drugs
cat("Available drugs (first 5):\n")
print(head(drugs_search$name, 5))

# And available omics features
cat("\nAvailable omics types:\n")
print(unique(omics_search$type))

# Example 1: Analyze association between Paclitaxel and ABCB1 gene expression
cat("\nAnalyzing association between Paclitaxel and ABCB1 gene expression...\n")
result <- analyzeDrugOmicPair(
  select_omics_type = "mRNA",
  select_omics = "A1CF",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display results
cat("Meta-analysis results:\n")
createForestPlot(result$meta)
result$plot

# Example 2: Analyze association between Paclitaxel and TP53 mutation status
cat("\nAnalyzing association between Paclitaxel and TP53 mutation status...\n")
result_mut <- analyzeDrugOmicPair(
  select_omics_type = "mutation_gene",
  select_omics = "TP53",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display results
cat("Meta-analysis results for mutation association:\n")
print(names(result_mut$data))
createForestPlot(result_mut$meta)
result_mut$plot

# Example 3: Get drug sensitivity data for visualization
cat("\nRetrieving drug sensitivity data for Paclitaxel...\n")
drug_data <- getDrugSensitivityData(
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display summary of retrieved data
cat("Number of samples with Paclitaxel sensitivity data:", nrow(drug_data), "\n")
cat("Studies with Paclitaxel data:\n")
print(table(drug_data$study))
