#!/usr/bin/env Rscript

# Example script for DROMA package: Drug Feature Analysis
# This example demonstrates drug data processing and visualization functions

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DT)
library(htmltools)
library(gridExtra)
library(DROMA)

# Setup DROMA environment - load search vectors and sample annotations
setupDROMA()

# Load drug data
loadDROMA("drug")
loadDROMA("mRNA") # For visualization examples with mRNA data

######################################
# Example 1: Basic Drug Data Processing
######################################

# Get processed drug data for Paclitaxel
cat("Example 1: Processing data for Paclitaxel\n")
paclitaxel_data <- processDrugData(
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display summary information
cat("Retrieved Paclitaxel data for", nrow(paclitaxel_data), "samples\n")
cat("Data columns:", paste(colnames(paclitaxel_data), collapse = ", "), "\n")
cat("Number of studies:", length(unique(paclitaxel_data$study)), "\n")
cat("Studies:", paste(unique(paclitaxel_data$study), collapse = ", "), "\n\n")

######################################
# Example 2: Annotating Drug Data
######################################

# Annotate the drug data with sample information
cat("Example 2: Annotating Paclitaxel data with sample information\n")
annotated_data <- annotateDrugData(paclitaxel_data)

# Display summaries with annotations
if ("TumorType" %in% colnames(annotated_data)) {
  tumor_counts <- table(annotated_data$TumorType)
  cat("Tumor Types in the dataset:\n")
  print(tumor_counts)
}

if ("ModelType" %in% colnames(annotated_data)) {
  model_counts <- table(annotated_data$ModelType)
  cat("\nModel Types in the dataset:\n")
  print(model_counts)
}

######################################
# Example 3: Formatters and Wrappers
######################################

# Format drug data as a datatable (for interactive use)
cat("\nExample 3: Using the drug data table formatter\n")
cat("In an interactive R session, this would display a formatted datatable\n")

# Using the all-in-one wrapper to get drug sensitivity data
cat("Using the getDrugSensitivityData wrapper function:\n")
drug_data <- getDrugSensitivityData(
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)
cat("Retrieved", nrow(drug_data), "samples with annotations\n\n")

######################################
# Example 4: Visualization Functions
######################################

# Load gene expression data for ABCB1 (drug resistance gene)
cat("Example 4: Visualizing relationships between drug sensitivity and other features\n")
abcb1_data <- selectFeatures(
  select_feas_type = "mRNA",
  select_feas = "ABCB1",
  data_type = "all",
  tumor_type = "all"
)

# Find a dataset where we have both ABCB1 expression and Paclitaxel sensitivity
common_studies <- intersect(names(abcb1_data), unique(drug_data$study))

if (length(common_studies) > 0) {
  # Use the first common study
  study_name <- common_studies[1]
  cat("Creating visualizations for study:", study_name, "\n")

  # Get samples from this study
  study_samples <- drug_data$sampleid[drug_data$study == study_name]
  gene_expr <- abcb1_data[[study_name]]

  # Find common samples
  common_samples <- intersect(names(gene_expr), study_samples)

  if (length(common_samples) >= 5) { # Need enough samples for meaningful visualization
    # Create data for plotting
    plot_data <- data.frame(
      ABCB1_expression = gene_expr[common_samples],
      drug_response = drug_data$zscore_value[match(common_samples, drug_data$sampleid)],
      stringsAsFactors = FALSE
    )

    cat("Created plot data with", nrow(plot_data), "samples having both measurements\n")

    # Demonstrate continuous comparison plot
    cat("In an interactive R session, this would display:\n")
    cat("1. A scatter plot with correlation statistics\n")
    cat("2. A boxplot with binned expression groups\n")
    cat("3. A comparison plot combining both visualizations\n\n")

    # Example code for interactive sessions:
    # p1 <- plotContinuousComparison(plot_data, "ABCB1_expression", "drug_response", "Paclitaxel IC50")
    # p2 <- plotContinuousGroups(plot_data, "ABCB1_expression", "drug_response", "Paclitaxel IC50", num_bins = 3)
    # combined_plot <- createDrugComparisonPlot(plot_data, "ABCB1_expression", "drug_response", "Paclitaxel IC50")
  } else {
    cat("Not enough common samples for visualization\n")
  }
} else {
  cat("No common studies with both ABCB1 expression and Paclitaxel sensitivity\n")
}

######################################
# Example 5: Categorical Comparisons
######################################

# Compare drug sensitivity across tumor types
cat("Example 5: Comparing drug sensitivity across categories\n")
# Filter to ensure we have enough samples per category
if ("TumorType" %in% colnames(drug_data)) {
  # Count samples per tumor type
  tumor_counts <- table(drug_data$TumorType)
  valid_tumors <- names(tumor_counts)[tumor_counts >= 3]

  if (length(valid_tumors) >= 2) {
    # Filter data to include only valid tumor types
    filtered_data <- drug_data[drug_data$TumorType %in% valid_tumors,]
    cat("Comparing Paclitaxel sensitivity across", length(valid_tumors), "tumor types\n")
    cat("In an interactive R session, this would display a categorical comparison boxplot\n")

    # Example code for interactive sessions:
    # tumor_plot <- plotCategoryComparison(filtered_data, "TumorType", "zscore_value", "Paclitaxel IC50")
  } else {
    cat("Not enough tumor types with sufficient samples for comparison\n")
  }
}
