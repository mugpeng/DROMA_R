#!/usr/bin/env Rscript

# Test script for FuncDrugFeature.R functions
# Tests process_drug_data, annotate_drug_data, format_drug_table, get_drug_sensitivity_data

# Load required packages
library(testthat)
library(dplyr)
library(DT)
library(htmltools)

# Source the function file
source("Package_Function/FuncDrugFeature.R")

# Test for process_drug_data function
test_that("process_drug_data returns correct data structure", {
  drug_name <- "5-Fluorouracil"
  data_type <- "all"
  tumor_type <- "all"
  
  drug_data <- process_drug_data(drug_name, data_type, tumor_type)
  
  # Check result structure
  expect_true(is.data.frame(drug_data) || is.null(drug_data))
  
  if (!is.null(drug_data)) {
    # Check required columns
    expect_true(all(c("sampleid", "zscore_value", "raw_value", "study") %in% colnames(drug_data)))
    
    # Check data types
    expect_true(is.character(drug_data$sampleid))
    expect_true(is.numeric(drug_data$zscore_value))
    expect_true(is.numeric(drug_data$raw_value))
    expect_true(is.character(drug_data$study))
    
    # Check no NA values
    expect_false(any(is.na(drug_data)))
  }
})

test_that("process_drug_data handles empty drug name", {
  expect_error(process_drug_data(""), "Please select a drug")
})

# Test for annotate_drug_data function
test_that("annotate_drug_data adds annotations correctly", {
  drug_name <- "5-Fluorouracil"
  drug_data <- process_drug_data(drug_name)
  
  if (!is.null(drug_data)) {
    annotated_data <- annotate_drug_data(drug_data)
    
    # Check result is a data frame
    expect_true(is.data.frame(annotated_data))
    
    # Check it's not empty
    expect_gt(nrow(annotated_data), 0)
    
    # Check sample annotation columns were added
    expect_true("DataType" %in% colnames(annotated_data))
    
    # Check "ProjectID" was removed if it exists in sample_anno
    if ("ProjectID" %in% colnames(sample_anno)) {
      expect_false("ProjectID" %in% colnames(annotated_data))
    }
    
    # Check joined correctly
    expect_equal(nrow(drug_data), nrow(annotated_data))
  }
})

test_that("annotate_drug_data handles NULL input", {
  expect_null(annotate_drug_data(NULL))
})

# Test for format_drug_table function
test_that("format_drug_table creates a DT datatable", {
  drug_name <- "5-Fluorouracil"
  drug_data <- process_drug_data(drug_name)
  
  if (!is.null(drug_data)) {
    formatted_table <- format_drug_table(drug_data)
    
    # Check result is a datatable
    expect_s3_class(formatted_table, "datatables")
  }
})

test_that("format_drug_table handles NULL input", {
  expect_null(format_drug_table(NULL))
})

test_that("format_drug_table handles empty data frame", {
  empty_df <- data.frame()
  expect_null(format_drug_table(empty_df))
})

# Test for get_drug_sensitivity_data function
test_that("get_drug_sensitivity_data returns correct data", {
  drug_name <- "5-Fluorouracil"
  
  # With annotations
  with_anno <- get_drug_sensitivity_data(drug_name, include_annotations = TRUE)
  if (!is.null(with_anno)) {
    expect_true(is.data.frame(with_anno))
    expect_true("DataType" %in% colnames(with_anno))
  }
  
  # Without annotations
  without_anno <- get_drug_sensitivity_data(drug_name, include_annotations = FALSE)
  if (!is.null(without_anno)) {
    expect_true(is.data.frame(without_anno))
    expect_true(all(c("sampleid", "zscore_value", "raw_value", "study") %in% colnames(without_anno)))
  }
}) 