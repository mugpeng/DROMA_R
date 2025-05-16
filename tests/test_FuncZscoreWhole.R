#!/usr/bin/env Rscript

# Test script for FuncZscoreWhole.R functions
# Tests zscoreNormalize, zscoreNormalizeDrug, applyZscoreNormalization, and resetToOriginal

# Load required packages
library(testthat)

# Source the function file
source("R/FuncZscoreWhole.R")
source("R/FuncGetData.R")  # For selectFeatures function if needed

# Load necessary data files
load("data/search_vec.Rda")  # Contains fea_list required by selectFeatures
load("data/mRNA.Rda")
load("data/drug.Rda")

context("Z-score Normalization Functions")

# Test for zscoreNormalize function
test_that("zscoreNormalize correctly normalizes omics data", {
  # Get a small subset of the mRNA data for testing
  if (exists("gdsc_mRNA")) {
    test_data <- gdsc_mRNA[1:10, 1:10]

    # Apply Z-score normalization
    normalized_data <- zscoreNormalize(test_data)

    # Verify results
    expect_equal(dim(normalized_data), dim(test_data))
    expect_true(all(abs(rowMeans(normalized_data, na.rm = TRUE)) < 1e-10))
    expect_true(all(abs(apply(normalized_data, 1, sd, na.rm = TRUE) - 1) < 1e-10 |
                    apply(test_data, 1, sd, na.rm = TRUE) < 1e-10))

    # Test handling of constant rows
    test_data_with_constant <- test_data
    test_data_with_constant[1,] <- rep(5, ncol(test_data))
    normalized_with_constant <- zscoreNormalize(test_data_with_constant)
    expect_true(all(is.na(normalized_with_constant[1,])))
  } else {
    skip("gdsc_mRNA data not available")
  }
})

# Test for zscoreNormalizeDrug function
test_that("zscoreNormalizeDrug correctly normalizes drug data", {
  # Get a small subset of the drug data for testing
  if (exists("gdsc1_drug")) {
    test_data <- gdsc1_drug[1:5, 1:10]

    # Apply Z-score normalization
    normalized_data <- zscoreNormalizeDrug(test_data)

    # Verify results
    expect_equal(dim(normalized_data), dim(test_data))
    expect_true(all(abs(rowMeans(normalized_data, na.rm = TRUE)) < 1e-10))
    expect_true(all(abs(apply(normalized_data, 1, sd, na.rm = TRUE) - 1) < 1e-10 |
                    apply(test_data, 1, sd, na.rm = TRUE) < 1e-10))
  } else {
    skip("gdsc1_drug data not available")
  }
})

# Test for applyZscoreNormalization function
test_that("applyZscoreNormalization doesn't fail", {
  # This function modifies global environment objects
  # We can only test that it doesn't fail, not its actual behavior
  expect_no_error(applyZscoreNormalization())
})

# Test for resetToOriginal function
test_that("resetToOriginal doesn't fail", {
  # This function modifies global environment objects
  # We can only test that it doesn't fail, not its actual behavior
  expect_no_error(resetToOriginal())
})
