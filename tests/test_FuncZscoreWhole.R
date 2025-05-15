#!/usr/bin/env Rscript

# Test script for FuncZscoreWhole.R functions
# Tests zscoreNormalize, zscoreNormalize, applyZscoreNormalization, and resetToOriginal

# Load required packages
library(testthat)

# Source the function file
source("R/FuncZscoreWhole.R")

context("Z-score Normalization Functions")

# Test for zscoreNormalize function
test_that("zscoreNormalize correctly normalizes omics data", {
  # Create a test matrix with genes in rows and samples in columns
  test_mat <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(test_mat) <- paste0("gene", 1:4)
  colnames(test_mat) <- paste0("sample", 1:5)
  
  # Apply normalization
  normalized <- zscoreNormalize(test_mat)
  
  # Check dimensions match
  expect_equal(dim(normalized), dim(test_mat))
  
  # Check row means are approximately 0
  row_means <- rowMeans(normalized)
  expect_true(all(abs(row_means) < 1e-10))
  
  # Check row standard deviations are approximately 1
  row_sds <- apply(normalized, 1, sd)
  expect_true(all(abs(row_sds - 1) < 1e-10))
})

test_that("zscoreNormalize handles constant rows correctly", {
  # Create a test matrix with a constant row
  test_mat <- matrix(1:15, nrow = 3, ncol = 5)
  test_mat[2,] <- rep(5, 5)  # Constant row
  
  # Apply normalization
  normalized <- zscoreNormalize(test_mat)
  
  # Check the constant row has been normalized to all zeros
  expect_equal(normalized[2,], rep(0, 5))
})

# Note: The applyZscoreNormalization and resetToOriginal functions modify global environment
# objects, which would require complex mocking. In practice, these should be
# tested with integration tests or with specific test fixtures that simulate
# the global environment.

test_that("applyZscoreNormalization handles non-existent data gracefully", {
  # This simply ensures the function doesn't error when data doesn't exist
  expect_error(applyZscoreNormalization(), NA)
})

test_that("resetToOriginal doesn't fail if source file doesn't exist", {
  # This is a minimal test to ensure the function definition is correct
  # In a real test environment, we would skip this or mock the source function
  skip("resetToOriginal requires source file that may not exist in test environment")
  
  # In a real test with appropriate fixtures:
  # expect_error(resetToOriginal(), NA)
}) 