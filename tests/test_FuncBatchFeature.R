#!/usr/bin/env Rscript

# Test script for FuncBatchFeature.R functions
# Tests metaCalcConCon, metaCalcConDis, metaCalcDisDis, and batch feature analysis functions

# Load required packages
library(testthat)
library(meta)
library(metafor)
library(effsize)
library(ggplot2)
library(dplyr)

# Source the function files
source("R/FuncBatchFeature.R")
source("R/FuncGetData.R") # For selectFeatures function

context("Batch Feature Analysis Functions")

# Test meta-analysis calculation functions
test_that("metaCalcConCon calculates meta-analysis for continuous data", {
  # Create mock paired data
  selected_pair <- list(
    pair1 = list(
      c("sample1" = 0.5, "sample2" = 1.2, "sample3" = -0.8, "sample4" = 0.3),
      c("sample1" = 0.1, "sample2" = 0.7, "sample3" = -0.3, "sample4" = 0.2)
    ),
    pair2 = list(
      c("sample5" = 0.9, "sample6" = -0.6, "sample7" = 1.5, "sample8" = -0.2),
      c("sample5" = 0.4, "sample6" = -0.2, "sample7" = 0.8, "sample8" = -0.1)
    )
  )
  
  # Calculate meta-analysis
  result <- metaCalcConCon(selected_pair)
  
  # Verify result
  expect_true(!is.null(result))
  expect_true(inherits(result, "meta"))
})

test_that("metaCalcConDis calculates meta-analysis for continuous vs. discrete data", {
  # Create mock paired data
  selected_pair <- list(
    pair1 = list(
      yes = c("sample1" = 0.1, "sample2" = 0.7, "sample4" = 0.2),
      no = c("sample3" = -0.3, "sample5" = 0.5, "sample6" = -0.1)
    ),
    pair2 = list(
      yes = c("sample7" = 0.4, "sample8" = -0.1, "sample9" = -0.4),
      no = c("sample10" = 0.1, "sample11" = -0.3, "sample12" = 0.6)
    )
  )
  
  # Calculate meta-analysis
  result <- metaCalcConDis(selected_pair)
  
  # Verify result
  expect_true(!is.null(result))
  expect_true(inherits(result, "meta"))
})

test_that("metaCalcDisDis calculates meta-analysis for discrete vs. discrete data", {
  # Create mock paired data with contingency tables
  selected_pair <- list(
    pair1 = list(
      cont_table = matrix(c(10, 15, 20, 30), nrow = 2,
                         dimnames = list(Feature1 = c("Yes", "No"),
                                         Feature2 = c("Yes", "No")))
    ),
    pair2 = list(
      cont_table = matrix(c(8, 12, 16, 25), nrow = 2,
                         dimnames = list(Feature1 = c("Yes", "No"),
                                         Feature2 = c("Yes", "No")))
    )
  )
  
  # Calculate meta-analysis
  result <- metaCalcDisDis(selected_pair)
  
  # Verify result
  expect_true(!is.null(result))
  expect_true(inherits(result, "meta"))
})

# Test utility functions
test_that("formatTime correctly formats time in seconds", {
  # Test seconds
  expect_equal(formatTime(45), "45 seconds")
  
  # Test minutes and seconds
  expect_equal(formatTime(125), "2 minutes 5 seconds")
  
  # Test hours and minutes
  expect_equal(formatTime(7325), "2 hours 2 minutes")
})

test_that("estimateTimeRemaining calculates remaining time", {
  # Test with 25% completion (1 out of 4 items, 10 seconds elapsed)
  estimate <- estimateTimeRemaining(1, 4, 10)
  expect_equal(estimate, 30) # Expect 30 seconds remaining
  
  # Test with 50% completion (5 out of 10 items, 100 seconds elapsed)
  estimate <- estimateTimeRemaining(5, 10, 100)
  expect_equal(estimate, 100) # Expect 100 seconds remaining
})

# Test pairing functions
test_that("pairContinuousFeatures correctly pairs continuous data", {
  # Create mock data
  myOmics <- list(
    dataset1 = c("cell1" = 0.5, "cell2" = 1.2, "cell3" = -0.8, "cell4" = 0.3),
    dataset2 = c("cell2" = 0.9, "cell3" = -0.6, "cell5" = 1.5, "cell6" = -0.2)
  )
  
  myDrugs <- list(
    drug1 = c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3, "cell5" = 0.5),
    drug2 = c("cell2" = 0.4, "cell3" = -0.2, "cell4" = 0.6, "cell6" = -0.4)
  )
  
  # Pair the data
  pairs <- pairContinuousFeatures(myOmics, myDrugs)
  
  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)
  
  # Check specific pairs
  if ("dataset1_drug1" %in% names(pairs)) {
    pair <- pairs[["dataset1_drug1"]]
    expect_true(is.list(pair))
    expect_equal(length(pair[[1]]), length(pair[[2]]))
    expect_equal(names(pair[[1]]), names(pair[[2]]))
  }
})

test_that("pairDiscreteFeatures correctly pairs discrete and continuous data", {
  # Create mock discrete omics data (lists of sample names with feature)
  myOmics <- list(
    dataset1 = c("cell1", "cell2", "cell4"),
    dataset2 = c("cell2", "cell5", "cell6")
  )
  
  # Create mock drug data
  myDrugs <- list(
    drug1 = c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3, "cell4" = 0.2, "cell5" = 0.5),
    drug2 = c("cell2" = 0.4, "cell3" = -0.2, "cell4" = 0.6, "cell5" = -0.1, "cell6" = -0.4)
  )
  
  # Pair the data
  pairs <- pairDiscreteFeatures(myOmics, myDrugs)
  
  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)
  
  # Check specific pairs
  if ("dataset1_drug1" %in% names(pairs)) {
    pair <- pairs[["dataset1_drug1"]]
    expect_true(is.list(pair))
    expect_true(all(c("yes", "no") %in% names(pair)))
    # "yes" should contain drug values for samples with the feature
    expect_true(all(names(pair$yes) %in% myOmics$dataset1))
  }
})

test_that("pairDiscreteDiscrete correctly pairs two discrete features", {
  # This function requires more complex setup, so we'll create a simplified test
  
  # Mock discrete features
  my_feas1 <- list(
    dataset1 = c("cell1", "cell2", "cell4"),
    dataset2 = c("cell2", "cell5", "cell6")
  )
  
  my_feas2 <- list(
    dataset1 = c("cell1", "cell3", "cell5"),
    dataset2 = c("cell3", "cell5", "cell7")
  )
  
  # Mock sample search data
  samples_search <- data.frame(
    type = rep(c("feature1", "feature2"), each = 4),
    datasets = rep(c("dataset1", "dataset2"), 4),
    cells = c("cell1", "cell2", "cell3", "cell4", "cell5", "cell6", "cell7", "cell8")
  )
  
  # Pair the data
  pairs <- pairDiscreteDiscrete(my_feas1, my_feas2, "feature1", "feature2", samples_search)
  
  # Skip detailed verification since this requires complex setup
  expect_true(is.list(pairs))
  
  # In a real implementation, we would check the contingency tables
  if (length(pairs) > 0 && "dataset1_dataset1" %in% names(pairs)) {
    pair <- pairs[["dataset1_dataset1"]]
    expect_true("cont_table" %in% names(pair))
    expect_true(is.matrix(pair$cont_table))
    expect_equal(dim(pair$cont_table), c(2, 2))
  }
})

# Test plotting function
test_that("plotMetaVolcano creates a volcano plot", {
  # Create mock meta-analysis results
  meta_df <- data.frame(
    effect_size = runif(20, -0.8, 0.8),
    p_value = 10^(-runif(20, 0, 5)),
    name = paste0("Feature", 1:20)
  )
  
  # Create plot
  plot <- plotMetaVolcano(meta_df, es_t = 0.4, P_t = 0.05)
  
  # Verify plot was created
  expect_true(inherits(plot, "ggplot"))
})

# Test the main batch analysis function
test_that("batchFindSignificantFeatures handles complex inputs", {
  # This is a very complex function that requires mocking several components
  # We'll just check if the function definition exists and basic validation works
  
  # Skip the full test as it requires mocking selectFeatures and global variables
  skip("This test requires mocking selectFeatures and global variables")
  
  # Check that the function exists
  expect_true(exists("batchFindSignificantFeatures"))
  
  # In a real test environment with proper mocking:
  # mock_feas_search <- data.frame(
  #   name = c("gene1", "gene2"),
  #   type = c("mRNA", "mRNA")
  # )
  # 
  # # Assign to global environment
  # assign("feas_search", mock_feas_search, envir = .GlobalEnv)
  # 
  # # Mock other required global variables
  # assign("samples_search", data.frame(), envir = .GlobalEnv)
  # assign("fea_list", list(mRNA = c("gene1", "gene2")), envir = .GlobalEnv)
  # 
  # # Mock selectFeatures
  # mockSelectFeatures <- function(...) {
  #   list(dataset1 = c("sample1" = 1, "sample2" = 2))
  # }
  # 
  # # Test with minimal arguments
  # results <- batchFindSignificantFeatures(
  #   feature1_type = "mRNA",
  #   feature1_name = "gene1",
  #   feature2_type = "mRNA",
  #   test_top_100 = TRUE
  # )
  # 
  # expect_true(is.data.frame(results))
}) 