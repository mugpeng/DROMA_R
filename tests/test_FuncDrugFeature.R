#!/usr/bin/env Rscript

# Test script for FuncDrugFeature.R functions
# Tests processDrugData, annotateDrugData, and visualization functions

# Load required packages
library(testthat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DT)
library(htmltools)
library(gridExtra)
library(DROMA)
# Load necessary data files
setupDROMA()
loadDROMA("drug")
loadDROMA("mRNA")

# When running tests in the package environment, we should use the package functions directly
# rather than sourcing the files

context("Drug Feature Functions")

# Setup real data for tests
# Define example features to use
select_drugs <- "Paclitaxel"  # Example drug

# Test for processDrugData function
test_that("processDrugData correctly processes drug data", {
  skip_if_not_installed("DROMA")

  # Process drug data
  drug_data <- processDrugData(select_drugs, data_type = "all", tumor_type = "all")

  # Verify results
  expect_true(is.data.frame(drug_data))
  expect_true(all(c("sampleid", "zscore_value", "raw_value", "study") %in% colnames(drug_data)))
  expect_true(nrow(drug_data) > 0)

  # Test error when no drug name provided
  expect_error(processDrugData(""), "Please select a drug.")
})

# Test for annotateDrugData function
test_that("annotateDrugData correctly annotates drug data", {
  skip_if_not_installed("DROMA")

  # Process drug data
  drug_data <- processDrugData(select_drugs, data_type = "all", tumor_type = "all")

  # Annotate drug data
  annotated_data <- annotateDrugData(drug_data)

  # Verify results
  expect_true(is.data.frame(annotated_data))
  expect_true(all(c("TumorType", "Gender", "SimpleEthnicity",
                    "Age") %in% colnames(annotated_data)))
  expect_true(nrow(annotated_data) > 0)

  # Test handling of NULL input
  expect_null(annotateDrugData(NULL))
})

# Test for formatDrugTable function
test_that("formatDrugTable creates a formatted datatable", {
  skip_if_not_installed("DROMA")

  # Process and annotate drug data
  drug_data <- processDrugData(select_drugs, data_type = "all", tumor_type = "all")

  # Create formatted table
  table <- formatDrugTable(drug_data)

  # Verify table was created
  expect_true(inherits(table, "datatables"))

  # Test with NULL input
  expect_null(formatDrugTable(NULL))
  expect_null(formatDrugTable(data.frame()))
})

# Test for getDrugSensitivityData function
test_that("getDrugSensitivityData correctly retrieves drug data", {
  skip_if_not_installed("DROMA")

  # Get drug data
  drug_data <- getDrugSensitivityData(
    drug_name = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Verify results
  expect_true(is.data.frame(drug_data))
  expect_true(all(c("sampleid", "zscore_value", "raw_value", "study") %in% colnames(drug_data)))
  expect_true(nrow(drug_data) > 0)

  # Test without annotations
  drug_data_no_anno <- getDrugSensitivityData(
    drug_name = select_drugs,
    include_annotations = FALSE
  )

  expect_true(is.data.frame(drug_data_no_anno))
  expect_false("TumorType" %in% colnames(drug_data_no_anno))
})

# Test visualization functions with real data
test_that("plotContinuousComparison creates a scatter plot", {
  skip_if_not_installed("DROMA")
  skip_if_not_installed("ggpubr")

  # Get drug data
  drug_data <- getDrugSensitivityData(
    drug_name = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Get a subset of the data for one study
  study_data <- drug_data[drug_data$study %in% "gdsc1", ]

  # Get expression data for ABCB1 gene
  myOmics <- selectFeatures("mRNA", "ABCB1", data_type = "all", tumor_type = "all")

  # Get expression values for the first dataset
  first_dataset <- names(myOmics)[1]
  expr_values <- myOmics[[first_dataset]]

  # Find common samples between drug and expression data
  common_samples <- intersect(study_data$sampleid, names(expr_values))

  # Skip if not enough common samples
  skip_if(length(common_samples) < 5, "Not enough common samples between drug and expression data")

  # Create data frame for plotting
  plot_data <- data.frame(
    value = study_data$zscore_value[study_data$sampleid %in% common_samples],
    expression = expr_values[common_samples]
  )

  # Create plot
  plot <- plotContinuousComparison(plot_data, "expression", "value", "Drug IC50")

  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotContinuousGroups creates a boxplot with binned groups", {
  skip_if_not_installed("DROMA")
  skip_if_not_installed("ggpubr")

  # Get drug data
  drug_data <- getDrugSensitivityData(
    drug_name = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Get a subset of the data for one study
  study_data <- drug_data[drug_data$study %in% "gdsc1", ]

  # Get expression data for ABCB1 gene
  myOmics <- selectFeatures("mRNA", "ABCB1", data_type = "all", tumor_type = "all")

  # Get expression values for the first dataset
  first_dataset <- names(myOmics)[1]
  expr_values <- myOmics[[first_dataset]]

  # Find common samples between drug and expression data
  common_samples <- intersect(study_data$sampleid, names(expr_values))

  # Skip if not enough common samples
  skip_if(length(common_samples) < 5, "Not enough common samples between drug and expression data")

  # Create data frame for plotting
  plot_data <- data.frame(
    value = study_data$zscore_value[study_data$sampleid %in% common_samples],
    expression = expr_values[common_samples]
  )

  # Create plot
  plot <- plotContinuousGroups(plot_data, "expression", "value", "Drug IC50", num_bins = 2)

  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotCategoryComparison creates a categorical boxplot", {
  skip_if_not_installed("DROMA")
  skip_if_not_installed("ggpubr")

  # Get drug data
  drug_data <- getDrugSensitivityData(
    drug_name = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Ensure we have categorical data (TumorType)
  skip_if(!"TumorType" %in% colnames(drug_data), "Drug data doesn't have TumorType column")

  # Count samples per tumor type
  tumor_counts <- table(drug_data$TumorType)
  valid_tumors <- names(tumor_counts)[tumor_counts >= 3]

  # Skip if not enough tumor types with sufficient samples
  skip_if(length(valid_tumors) < 2, "Not enough tumor types with sufficient samples")

  # Filter data to include only valid tumor types
  plot_data <- drug_data[drug_data$TumorType %in% valid_tumors, ]

  # Create plot
  plot <- plotCategoryComparison(plot_data, "TumorType", "zscore_value", "Drug IC50")

  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("createDrugComparisonPlot creates appropriate plot by variable type", {
  skip_if_not_installed("DROMA")
  skip_if_not_installed("ggpubr")

  # Get drug data
  drug_data <- getDrugSensitivityData(
    drug_name = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Test with categorical variable (TumorType)
  if ("TumorType" %in% colnames(drug_data)) {
    # Count samples per tumor type
    tumor_counts <- table(drug_data$TumorType)
    valid_tumors <- names(tumor_counts)[tumor_counts >= 3]

    # Skip if not enough tumor types with sufficient samples
    if (length(valid_tumors) >= 2) {
      # Filter data to include only valid tumor types
      plot_data_cat <- drug_data[drug_data$TumorType %in% valid_tumors, ]

      # Create plot
      plot_cat <- createDrugComparisonPlot(plot_data_cat, "TumorType", "zscore_value", "Drug IC50")

      # Verify plot is created
      expect_true(inherits(plot_cat, "ggplot"))
    }
  }

  # Test with continuous variable (expression data)
  # Get expression data for ABCB1 gene
  myOmics <- selectFeatures("mRNA", "ABCB1", data_type = "all", tumor_type = "all")

  # Get expression values for the first dataset
  first_dataset <- names(myOmics)[1]
  expr_values <- myOmics[[first_dataset]]

  # Get a subset of the drug data for one study
  study_data <- drug_data[drug_data$study %in% "gdsc1", ]

  # Find common samples between drug and expression data
  common_samples <- intersect(study_data$sampleid, names(expr_values))

  # Skip if not enough common samples
  skip_if(length(common_samples) < 5, "Not enough common samples between drug and expression data")

  # Create data frame for plotting
  plot_data_cont <- data.frame(
    value = study_data$zscore_value[study_data$sampleid %in% common_samples],
    expression = expr_values[common_samples]
  )

  # Create plot without grouped boxplot
  plot_cont <- createDrugComparisonPlot(plot_data_cont, "expression", "value", "Drug IC50", show_groups_boxplot = FALSE)

  # Verify plot is created
  expect_true(inherits(plot_cont, "ggplot"))

  # Create plot with grouped boxplot
  plot_cont_group <- createDrugComparisonPlot(plot_data_cont, "expression", "value", "Drug IC50", show_groups_boxplot = TRUE)

  # Verify plot is created
  expect_true(inherits(plot_cont_group, "grob") || inherits(plot_cont_group, "gtable"))
})
