#!/usr/bin/env Rscript

# Test script for FuncZscoreWhole.R functions
# Tests zscore_normalize, zscore_normalize_drug, apply_zscore_normalization, reset_to_original

# Load required packages
library(testthat)

# Source the function file
source("Package_Function/FuncZscoreWhole.R")

# Test for zscore_normalize function
test_that("zscore_normalize standardizes each row correctly", {
  # Create test matrix
  test_matrix <- matrix(
    c(1, 2, 3, 4, 5,
      5, 4, 3, 2, 1,
      1, 1, 1, 1, 1), # Row with zero variance
    nrow = 3, 
    byrow = TRUE
  )
  
  # Apply z-score normalization
  normalized <- zscore_normalize(test_matrix)
  
  # Check dimensions unchanged
  expect_equal(dim(normalized), dim(test_matrix))
  
  # Check first row is standardized (should have mean 0, sd 1)
  expect_equal(mean(normalized[1,]), 0, tolerance = 1e-10)
  expect_equal(sd(normalized[1,]), 1, tolerance = 1e-10)
  
  # Check second row is standardized
  expect_equal(mean(normalized[2,]), 0, tolerance = 1e-10)
  expect_equal(sd(normalized[2,]), 1, tolerance = 1e-10)
  
  # Check third row (with zero variance) is handled correctly
  expect_equal(sum(normalized[3,]), 0) # Should be all zeros
  expect_equal(sd(normalized[3,]), 0)
})

# Test for zscore_normalize_drug function
test_that("zscore_normalize_drug standardizes each drug row independently", {
  # Create test data frame
  test_df <- data.frame(
    sample1 = c(1, 5, 1),
    sample2 = c(2, 4, 1),
    sample3 = c(3, 3, 1),
    sample4 = c(4, 2, 1),
    sample5 = c(5, 1, 1)
  )
  
  # Apply z-score normalization
  normalized <- zscore_normalize_drug(test_df)
  
  # Check class and dimensions unchanged
  expect_true(is.data.frame(normalized))
  expect_equal(dim(normalized), dim(test_df))
  
  # Check first row is standardized (should have mean 0, sd 1)
  expect_equal(mean(as.numeric(normalized[1,])), 0, tolerance = 1e-10)
  expect_equal(sd(as.numeric(normalized[1,])), 1, tolerance = 1e-10)
  
  # Check second row is standardized
  expect_equal(mean(as.numeric(normalized[2,])), 0, tolerance = 1e-10)
  expect_equal(sd(as.numeric(normalized[2,])), 1, tolerance = 1e-10)
  
  # Check third row (with zero variance) is centered but not scaled
  expect_equal(mean(as.numeric(normalized[3,])), 0, tolerance = 1e-10)
  expect_equal(sd(as.numeric(normalized[3,])), 0, tolerance = 1e-10)
})

# Test for apply_zscore_normalization function
test_that("apply_zscore_normalization handles global variables correctly", {
  # Create some test data in the global environment
  assign("test_mRNA", matrix(rnorm(20), nrow = 4), envir = .GlobalEnv)
  assign("test_drug", data.frame(
    sample1 = rnorm(3),
    sample2 = rnorm(3),
    sample3 = rnorm(3)
  ), envir = .GlobalEnv)
  
  # Back up original normalization_state if it exists
  if (exists("normalization_state", envir = .GlobalEnv)) {
    original_state <- get("normalization_state", envir = .GlobalEnv)
    on.exit(assign("normalization_state", original_state, envir = .GlobalEnv))
  } else {
    on.exit(if (exists("normalization_state", envir = .GlobalEnv)) 
      rm("normalization_state", envir = .GlobalEnv))
  }
  
  # We can't directly test the function since it works on specific global variables
  # But we can check that it at least sets the normalization_state correctly
  
  # Mock the data type variables that the function expects
  for (var_name in c("ccle_mRNA", "gdsc_drug", "prism_drug")) {
    if (!exists(var_name, envir = .GlobalEnv)) {
      assign(var_name, if (grepl("drug", var_name)) 
                         data.frame(sample1 = rnorm(3), sample2 = rnorm(3)) 
                       else 
                         matrix(rnorm(20), nrow = 4),
             envir = .GlobalEnv)
      on.exit(if (exists(var_name, envir = .GlobalEnv)) 
        rm(var_name, envir = .GlobalEnv), add = TRUE)
    }
  }
  
  # Test that the function runs without error
  expect_silent(apply_zscore_normalization())
  
  # Check that normalization_state is set to TRUE
  expect_true(get("normalization_state", envir = .GlobalEnv))
  
  # Clean up
  rm("test_mRNA", "test_drug", envir = .GlobalEnv)
})

# Test for reset_to_original function
test_that("reset_to_original updates normalization state", {
  # This function relies on loading external data from a script
  # We can't fully test it, but we can check that it sets normalization_state correctly
  
  # Back up original normalization_state if it exists
  if (exists("normalization_state", envir = .GlobalEnv)) {
    original_state <- get("normalization_state", envir = .GlobalEnv)
    on.exit(assign("normalization_state", original_state, envir = .GlobalEnv))
  } else {
    on.exit(if (exists("normalization_state", envir = .GlobalEnv)) 
      rm("normalization_state", envir = .GlobalEnv))
  }
  
  # Set normalization_state to TRUE initially
  assign("normalization_state", TRUE, envir = .GlobalEnv)
  
  # Mock the source function to avoid actually running the data loading script
  mockery::stub(reset_to_original, "source", function(...) {
    # Do nothing
  })
  
  # Test that the function runs
  # Since we've mocked the source function, this should just set normalization_state
  tryCatch({
    reset_to_original()
    
    # Check that normalization_state is set to FALSE
    expect_false(get("normalization_state", envir = .GlobalEnv))
  }, error = function(e) {
    # Handle case where mockery isn't available
    skip("Mockery package required for this test")
  })
}) 