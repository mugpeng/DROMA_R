#!/usr/bin/env Rscript

# Example script for DROMA package: Feature Selection with selectFeatures
# This example demonstrates different ways to retrieve omics and drug data

# Load required libraries
library(DROMA)

# Setup DROMA environment - load search vectors and sample annotations
setupDROMA()

######################################
# Basic Feature Selection Examples
######################################

# Load necessary data
loadDROMA("drug")
loadDROMA("mRNA")
loadDROMA("mut")

# Example 1: Select a drug feature (Paclitaxel) across all data types and tumor types
cat("Example 1: Selecting drug data for Paclitaxel across all datasets\n")
paclitaxel_data <- selectFeatures(
  select_feas_type = "drug",
  select_feas = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Print summary of retrieved data
cat("Retrieved Paclitaxel data from", length(paclitaxel_data), "datasets:\n")
for (dataset_name in names(paclitaxel_data)) {
  cat(sprintf("- %s: %d samples\n", dataset_name, length(paclitaxel_data[[dataset_name]])))
}

# Example 2: Select an mRNA feature (ABCB1) across all data types and tumor types
cat("\nExample 2: Selecting mRNA data for ABCB1 across all datasets\n")
abcb1_data <- selectFeatures(
  select_feas_type = "mRNA",
  select_feas = "ABCB1",
  data_type = "all",
  tumor_type = "all"
)

# Print summary of retrieved data
cat("Retrieved ABCB1 data from", length(abcb1_data), "datasets:\n")
for (dataset_name in names(abcb1_data)) {
  cat(sprintf("- %s: %d samples\n", dataset_name, length(abcb1_data[[dataset_name]])))
}

# Example 3: Select a mutation feature (TP53) across all data types and tumor types
cat("\nExample 3: Selecting mutation data for TP53 across all datasets\n")
tp53_data <- selectFeatures(
  select_feas_type = "mutation_gene",
  select_feas = "TP53",
  data_type = "all",
  tumor_type = "all"
)

# Print summary of retrieved data
cat("Retrieved TP53 mutation data from", length(tp53_data), "datasets:\n")
for (dataset_name in names(tp53_data)) {
  cat(sprintf("- %s: %d samples\n", dataset_name, length(tp53_data[[dataset_name]])))
}

######################################
# Filtering by Data Type and Tumor Type
######################################

# Example 4: Get Paclitaxel data only for cell lines
cat("\nExample 4: Selecting Paclitaxel data only for cell lines\n")
paclitaxel_cellline_data <- selectFeatures(
  select_feas_type = "drug",
  select_feas = "Paclitaxel",
  data_type = "CellLine",
  tumor_type = "all"
)

# Print summary
cat("Retrieved Paclitaxel cell line data from", length(paclitaxel_cellline_data), "datasets:\n")
for (dataset_name in names(paclitaxel_cellline_data)) {
  cat(sprintf("- %s: %d samples\n", dataset_name, length(paclitaxel_cellline_data[[dataset_name]])))
}

# Example 5: Get ABCB1 data only for a specific tumor type (e.g., breast cancer)
# Note: This example assumes breast cancer data is available in the loaded datasets
cat("\nExample 5: Selecting ABCB1 data for breast cancer\n")
abcb1_breast_data <- selectFeatures(
  select_feas_type = "mRNA",
  select_feas = "ABCB1",
  data_type = "all",
  tumor_type = "breast cancer"
)

# Print summary (note this will show an empty list if breast cancer data is not available)
cat("Retrieved ABCB1 breast cancer data from", length(abcb1_breast_data), "datasets\n")
if (length(abcb1_breast_data) > 0) {
  for (dataset_name in names(abcb1_breast_data)) {
    cat(sprintf("- %s: %d samples\n", dataset_name, length(abcb1_breast_data[[dataset_name]])))
  }
} else {
  cat("No breast cancer data available in the loaded datasets\n")
}

######################################
# Example of Error Handling
######################################

# Example 6: Try to select a non-existent feature (should generate error)
cat("\nExample 6: Attempting to select a non-existent feature\n")
tryCatch({
  nonexistent_data <- selectFeatures(
    select_feas_type = "drug",
    select_feas = "NonExistentDrug",
    data_type = "all",
    tumor_type = "all"
  )
}, error = function(e) {
  cat("Error caught (as expected):", conditionMessage(e), "\n")
})

# Example 7: Try to use an invalid feature type (should generate error)
cat("\nExample 7: Attempting to use an invalid feature type\n")
tryCatch({
  invalid_type_data <- selectFeatures(
    select_feas_type = "invalid_type",
    select_feas = "Paclitaxel",
    data_type = "all",
    tumor_type = "all"
  )
}, error = function(e) {
  cat("Error caught (as expected):", conditionMessage(e), "\n")
})
