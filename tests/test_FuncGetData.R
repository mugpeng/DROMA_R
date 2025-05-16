#!/usr/bin/env Rscript

# Test script for FuncGetData.R functions
# Tests selectFeatures

# Load required packages
library(testthat)

# Source the function file
source("R/FuncGetData.R")

# Load necessary data files
load("data/search_vec.Rda")  # Contains fea_list required by selectFeatures
load("data/drug.Rda")
load("data/mRNA.Rda")
load("data/anno.Rda")

context("Data Selection Functions")

# Setup mock data
select_drugs <- "Paclitaxel"  # Example drug
select_omics_type <- "mRNA"    # Example omics type
select_omics <- "ABCB1"        # Example omics feature

# Test for selectFeatures function
test_that("selectFeatures returns correct data for drug", {
  # Call function
  myDrugs <- selectFeatures("drug", select_drugs,
                        data_type = "all",
                        tumor_type = "all")

  # Verify results
  expect_true(is.list(myDrugs))
  expect_true(length(myDrugs) > 0)

  # Test for a specific dataset
  if ("gdsc1" %in% names(myDrugs)) {
    drug_data <- myDrugs[["gdsc1"]]
    expect_true(is.numeric(drug_data))
    expect_true(length(drug_data) > 0)
  }
})

test_that("selectFeatures returns correct data for mRNA", {
  # Call function
  myOmics <- selectFeatures(select_omics_type, select_omics,
                         data_type = "all",
                         tumor_type = "all")

  # Verify results
  expect_true(is.list(myOmics))
  expect_true(length(myOmics) > 0)

  # Test for a specific dataset
  if ("gdsc" %in% names(myOmics)) {
    omics_data <- myOmics[["gdsc"]]
    expect_true(is.numeric(omics_data))
    expect_true(length(omics_data) > 0)
  }
})

test_that("selectFeatures handles invalid feature type", {
  expect_error(selectFeatures("invalid_type", "ABCB1"),
               "The select feature type doesn't exist")
})

test_that("selectFeatures handles invalid data_type", {
  expect_error(selectFeatures("mRNA", "ABCB1", data_type = "invalid_type"),
               "Invalid data_type")
})
