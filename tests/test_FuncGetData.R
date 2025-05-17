#!/usr/bin/env Rscript

# Test script for FuncGetData.R functions
# Tests selectFeatures

# Load required packages
library(testthat)
library(DROMA)
# Source the function file - use proper package testing approach
# When running tests in the package environment, we should use the package functions directly
# rather than sourcing the files

# Load necessary data files
setupDROMA()
loadDROMA("drug")
loadDROMA("mRNA")
context("Data Selection Functions")

# Setup mock data
select_drugs <- "Paclitaxel"  # Example drug
select_omics_type <- "mRNA"    # Example omics type
select_omics <- "ABCB1"        # Example omics feature

# Test for selectFeatures function
test_that("selectFeatures returns correct data for drug", {
  skip_if_not_installed("DROMA")
  skip_if(!exists("fea_list"))

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
  skip_if_not_installed("DROMA")
  skip_if(!exists("fea_list"))

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
  skip_if_not_installed("DROMA")

  expect_error(selectFeatures("invalid_type", "ABCB1"),
               "The select feature type doesn't exist")
})

test_that("selectFeatures handles invalid data_type", {
  skip_if_not_installed("DROMA")

  expect_error(selectFeatures("mRNA", "ABCB1", data_type = "invalid_type"),
               "Invalid data_type")
})
