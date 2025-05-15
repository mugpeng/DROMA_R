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

# Source the function files
source("R/FuncDrugFeature.R")
source("R/FuncGetData.R") # For selectFeatures function

context("Drug Feature Functions")

# Mock the selectFeatures function for testing
mock_selectFeatures <- function(select_feas_type, select_feas, data_type, tumor_type) {
  if (select_feas_type == "drug") {
    return(list(
      dataset1 = c("sample1" = 0.1, "sample2" = 0.7, "sample3" = -0.3, "sample4" = 0.2),
      dataset2 = c("sample2" = 0.4, "sample3" = -0.2, "sample5" = 0.6, "sample6" = -0.4)
    ))
  } else if (select_feas_type == "drug_raw") {
    return(list(
      dataset1 = c("sample1" = 10, "sample2" = 70, "sample3" = 30, "sample4" = 20),
      dataset2 = c("sample2" = 40, "sample3" = 20, "sample5" = 60, "sample6" = 40)
    ))
  }
}

# Test for processDrugData function
test_that("processDrugData correctly processes drug data", {
  # This test requires mocking selectFeatures
  # We'll temporarily replace the function for this test
  original_selectFeatures <- selectFeatures
  
  # Setup the test environment
  temp_env <- environment()
  assign("selectFeatures", mock_selectFeatures, envir = temp_env)
  
  # Call processDrugData within this environment
  local({
    # Process drug data
    drug_data <- processDrugData("Paclitaxel", data_type = "all", tumor_type = "all")
    
    # Verify results
    expect_true(is.data.frame(drug_data))
    expect_true(all(c("sampleid", "zscore_value", "raw_value", "study") %in% colnames(drug_data)))
    
    # Check specific values
    dataset1_rows <- drug_data[drug_data$study == "dataset1", ]
    expect_equal(nrow(dataset1_rows), 4)
    
    # Check for proper merging of normalized and raw values
    sample2_row <- drug_data[drug_data$sampleid == "sample2" & drug_data$study == "dataset1", ]
    expect_equal(sample2_row$zscore_value, 0.7)
    expect_equal(sample2_row$raw_value, 70)
  })
  
  # Test error when no drug name provided
  expect_error(processDrugData(""), "Please select a drug.")
})

# Test for annotateDrugData function
test_that("annotateDrugData correctly annotates drug data", {
  # Create mock drug data
  drug_data <- data.frame(
    sampleid = c("sample1", "sample2", "sample3", "sample4"),
    zscore_value = c(0.1, 0.7, -0.3, 0.2),
    raw_value = c(10, 70, 30, 20),
    study = c("dataset1", "dataset1", "dataset1", "dataset1")
  )
  
  # Create mock sample annotations
  sample_annotations <- data.frame(
    SampleID = c("sample1", "sample2", "sample3", "sample4"),
    ProjectID = c("proj1", "proj2", "proj1", "proj2"),
    TumorType = c("Breast", "Lung", "Breast", "Colon"),
    ModelType = c("CellLine", "CellLine", "PDX", "PDO")
  )
  
  # Annotate drug data
  annotated_data <- annotateDrugData(drug_data, sample_annotations)
  
  # Verify results
  expect_true(is.data.frame(annotated_data))
  expect_true(all(c("TumorType", "ModelType") %in% colnames(annotated_data)))
  expect_equal(nrow(annotated_data), 4)
  
  # Check that annotations were correctly joined
  sample1_row <- annotated_data[annotated_data$sampleid == "sample1", ]
  expect_equal(sample1_row$TumorType, "Breast")
  expect_equal(sample1_row$ModelType, "CellLine")
  
  # Test handling of NULL input
  expect_null(annotateDrugData(NULL))
  
  # Test with NULL annotations but global variable available
  # This is complex to test since it uses .GlobalEnv, so we'll skip detailed testing
})

# Test for formatDrugTable function
test_that("formatDrugTable creates a formatted datatable", {
  # Create mock drug data
  drug_data <- data.frame(
    sampleid = c("sample1", "sample2", "sample3", "sample4"),
    zscore_value = c(0.1, 0.7, -0.3, 0.2),
    raw_value = c(10, 70, 30, 20),
    study = c("dataset1", "dataset1", "dataset1", "dataset1"),
    TumorType = c("Breast", "Lung", "Breast", "Colon"),
    ModelType = c("CellLine", "CellLine", "PDX", "PDO")
  )
  
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
  # This test requires mocking processDrugData and annotateDrugData
  
  # Setup the test environment with mock functions
  temp_env <- environment()
  
  # Mock processDrugData
  mock_processDrugData <- function(drug_name, data_type, tumor_type) {
    data.frame(
      sampleid = c("sample1", "sample2", "sample3", "sample4"),
      zscore_value = c(0.1, 0.7, -0.3, 0.2),
      raw_value = c(10, 70, 30, 20),
      study = c("dataset1", "dataset1", "dataset1", "dataset1")
    )
  }
  
  # Mock annotateDrugData
  mock_annotateDrugData <- function(drug_data, sample_annotations) {
    cbind(drug_data, 
          TumorType = c("Breast", "Lung", "Breast", "Colon"),
          ModelType = c("CellLine", "CellLine", "PDX", "PDO"))
  }
  
  # Assign mock functions in the environment
  assign("processDrugData", mock_processDrugData, envir = temp_env)
  assign("annotateDrugData", mock_annotateDrugData, envir = temp_env)
  
  # Call the function within this environment
  local({
    # Get drug data
    drug_data <- getDrugSensitivityData("Paclitaxel", 
                                       data_type = "all",
                                       tumor_type = "all")
    
    # Verify results
    expect_true(is.data.frame(drug_data))
    expect_true(all(c("sampleid", "zscore_value", "raw_value", "study", 
                       "TumorType", "ModelType") %in% colnames(drug_data)))
    expect_equal(nrow(drug_data), 4)
    
    # Test without annotations
    drug_data_no_anno <- getDrugSensitivityData("Paclitaxel", 
                                              include_annotations = FALSE)
    
    expect_true(is.data.frame(drug_data_no_anno))
    expect_false("TumorType" %in% colnames(drug_data_no_anno))
  })
})

# Test visualization functions
test_that("plotContinuousComparison creates a scatter plot", {
  # Create mock data
  data <- data.frame(
    value = c(0.1, 0.7, -0.3, 0.2, 0.5, -0.4),
    expression = c(1.5, 2.7, 0.8, 1.2, 2.0, 0.5)
  )
  
  # Create plot
  plot <- plotContinuousComparison(data, "expression", "value", "Drug IC50")
  
  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotContinuousGroups creates a boxplot with binned groups", {
  # Create mock data
  data <- data.frame(
    value = c(0.1, 0.7, -0.3, 0.2, 0.5, -0.4, 0.3, -0.2),
    expression = c(1.5, 2.7, 0.8, 1.2, 2.0, 0.5, 1.8, 0.9)
  )
  
  # Create plot
  plot <- plotContinuousGroups(data, "expression", "value", "Drug IC50", num_bins = 2)
  
  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotCategoryComparison creates a categorical boxplot", {
  # Create mock data
  data <- data.frame(
    value = c(0.1, 0.7, -0.3, 0.2, 0.5, -0.4, 0.3, -0.2),
    category = c("A", "B", "A", "B", "C", "A", "C", "B")
  )
  
  # Create plot
  plot <- plotCategoryComparison(data, "category", "value", "Drug IC50")
  
  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("createDrugComparisonPlot creates appropriate plot by variable type", {
  # Create mock data for continuous comparison
  data_cont <- data.frame(
    value = c(0.1, 0.7, -0.3, 0.2, 0.5, -0.4, 0.3, -0.2),
    expression = c(1.5, 2.7, 0.8, 1.2, 2.0, 0.5, 1.8, 0.9)
  )
  
  # Create mock data for categorical comparison
  data_cat <- data.frame(
    value = c(0.1, 0.7, -0.3, 0.2, 0.5, -0.4, 0.3, -0.2),
    category = c("A", "B", "A", "B", "C", "A", "C", "B")
  )
  
  # Test with continuous variable
  plot_cont <- createDrugComparisonPlot(data_cont, "expression", "value", "Drug IC50", show_groups_boxplot = FALSE)
  expect_true(inherits(plot_cont, "ggplot"))
  
  # Test with continuous variable and grouped boxplot
  plot_cont_group <- createDrugComparisonPlot(data_cont, "expression", "value", "Drug IC50", show_groups_boxplot = TRUE)
  expect_true(inherits(plot_cont_group, "grob") || inherits(plot_cont_group, "gtable"))
  
  # Test with categorical variable
  plot_cat <- createDrugComparisonPlot(data_cat, "category", "value", "Drug IC50")
  expect_true(inherits(plot_cat, "ggplot"))
  
  # Test with no data
  empty_data <- data.frame(value = numeric(0), var = character(0))
  plot_empty <- createDrugComparisonPlot(empty_data, "var", "value", "Drug IC50")
  expect_true(inherits(plot_empty, "ggplot"))
}) 