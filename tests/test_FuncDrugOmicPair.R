#!/usr/bin/env Rscript

# Test script for FuncDrugOmicPair.R functions
# Tests pairDrugOmic, analyzeContinuousDrugOmic, pairDiscreteDrugOmic, 
# analyzeDiscreteDrugOmic, and plotting functions

# Load required packages
library(testthat)
library(ggplot2)
library(metafor)
library(meta)
library(effsize)
library(ggpubr)

# Source the function files
source("R/FuncDrugOmicPair.R")
source("R/FuncGetData.R") # Needed for selectFeatures function

context("Drug-Omic Pairing Functions")

# Setup mock data for continuous tests
test_that("pairDrugOmic correctly pairs drug and omics data", {
  # Create mock omics and drug data
  myOmics <- list(
    dataset1 = c("cell1" = 0.5, "cell2" = 1.2, "cell3" = -0.8, "cell4" = 0.3),
    dataset2 = c("cell2" = 0.9, "cell3" = -0.6, "cell5" = 1.5, "cell6" = -0.2)
  )
  
  myDrugs <- list(
    drug1 = c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3, "cell5" = 0.5),
    drug2 = c("cell2" = 0.4, "cell3" = -0.2, "cell4" = 0.6, "cell6" = -0.4)
  )
  
  # Pair the data
  pairs <- pairDrugOmic(myOmics, myDrugs)
  
  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)
  
  # Check one specific pair
  expect_true("dataset1_drug1" %in% names(pairs))
  if ("dataset1_drug1" %in% names(pairs)) {
    pair <- pairs[["dataset1_drug1"]]
    expect_true(is.list(pair))
    expect_true(all(c("omic", "drug") %in% names(pair)))
    expect_equal(length(pair$omic), length(pair$drug))
    expect_equal(names(pair$omic), names(pair$drug))
  }
  
  # Test with merged = TRUE
  merged_pairs <- pairDrugOmic(myOmics, myDrugs, merged = TRUE)
  expect_true("merged_dataset" %in% names(merged_pairs))
})

test_that("analyzeContinuousDrugOmic correctly analyzes pairs", {
  # Create mock pairs
  myPairs <- list(
    dataset1_drug1 = list(
      omic = c("cell1" = 0.5, "cell2" = 1.2, "cell3" = -0.8),
      drug = c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3)
    ),
    dataset2_drug1 = list(
      omic = c("cell2" = 0.9, "cell3" = -0.6, "cell5" = 1.5),
      drug = c("cell2" = 0.4, "cell3" = -0.2, "cell5" = 0.8)
    )
  )
  
  # Add a merged dataset (should be skipped in meta-analysis)
  myPairs[["merged_dataset"]] <- list(
    omic = c("cell1" = 0.5, "cell2" = 1.2, "cell3" = -0.8, "cell5" = 1.5),
    drug = c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3, "cell5" = 0.8)
  )
  
  # Analyze the pairs
  results <- analyzeContinuousDrugOmic(myPairs)
  
  # Verify results
  expect_true(!is.null(results))
  expect_true(inherits(results, "meta"))
})

# Setup mock data for discrete tests
test_that("pairDiscreteDrugOmic correctly pairs discrete data", {
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
  pairs <- pairDiscreteDrugOmic(myOmics, myDrugs)
  
  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)
  
  # Check one specific pair
  expect_true("dataset1_drug1" %in% names(pairs))
  if ("dataset1_drug1" %in% names(pairs)) {
    pair <- pairs[["dataset1_drug1"]]
    expect_true(is.list(pair))
    expect_true(all(c("yes", "no") %in% names(pair)))
    # "yes" should contain drug values for samples with the feature
    expect_true(all(names(pair$yes) %in% myOmics$dataset1))
    # "no" should contain drug values for samples without the feature
    expect_true(all(!names(pair$no) %in% myOmics$dataset1))
  }
  
  # Test with merged = TRUE
  merged_pairs <- pairDiscreteDrugOmic(myOmics, myDrugs, merged = TRUE)
  expect_true("merged_dataset" %in% names(merged_pairs))
})

test_that("analyzeDiscreteDrugOmic correctly analyzes pairs", {
  # Create mock pairs
  myPairs <- list(
    dataset1_drug1 = list(
      yes = c("cell1" = 0.1, "cell2" = 0.7, "cell4" = 0.2),
      no = c("cell3" = -0.3, "cell5" = 0.5, "cell6" = -0.1)
    ),
    dataset2_drug1 = list(
      yes = c("cell2" = 0.4, "cell5" = -0.1, "cell6" = -0.4),
      no = c("cell1" = 0.1, "cell3" = -0.3, "cell4" = 0.6)
    )
  )
  
  # Add a merged dataset (should be skipped in meta-analysis)
  myPairs[["merged_dataset"]] <- list(
    yes = c("cell1" = 0.1, "cell2" = 0.7, "cell4" = 0.2, "cell5" = -0.1, "cell6" = -0.4),
    no = c("cell3" = -0.3)
  )
  
  # Analyze the pairs
  results <- analyzeDiscreteDrugOmic(myPairs)
  
  # Verify results
  expect_true(!is.null(results))
  expect_true(inherits(results, "meta"))
})

# Test plotting functions (minimal tests just to check they run)
test_that("createForestPlot creates a forest plot", {
  # This is a minimal test to check the function runs without error
  # For complete testing, we would need actual meta-analysis results
  
  skip("This test requires actual meta-analysis results")
  
  # In a real test with appropriate fixtures:
  # meta_obj <- analyzeContinuousDrugOmic(myPairs)
  # expect_error(createForestPlot(meta_obj), NA)
})

test_that("plotContinuousDrugOmic creates a scatter plot", {
  # Create mock data
  omic_values <- c("cell1" = 0.5, "cell2" = 1.2, "cell3" = -0.8, "cell4" = 0.3)
  drug_values <- c("cell1" = 0.1, "cell2" = 0.7, "cell3" = -0.3, "cell4" = 0.2)
  
  # Create plot
  plot <- plotContinuousDrugOmic(omic_values, drug_values, "Test Study")
  
  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotDiscreteDrugOmic creates a boxplot", {
  # Create mock data
  yes_values <- c("cell1" = 0.1, "cell2" = 0.7, "cell4" = 0.2)
  no_values <- c("cell3" = -0.3, "cell5" = 0.5, "cell6" = -0.1)
  
  # Create plot
  plot <- plotDiscreteDrugOmic(yes_values, no_values, "Test Study")
  
  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

# Test the main analysis function
test_that("analyzeDrugOmicPair handles both continuous and discrete data", {
  # This is a complex function that requires mocking selectFeatures
  # We'll mock the function by temporarily redefining it
  
  skip("This test requires mocking the selectFeatures function")
  
  # In a real test environment with proper mocking:
  # mockSelectFeatures <- function(...) {
  #   # Return appropriate mock data based on args
  # }
  # 
  # # Temporarily replace the selectFeatures function
  # original <- selectFeatures
  # assign("selectFeatures", mockSelectFeatures, envir = environment(analyzeDrugOmicPair))
  # 
  # # Test the function
  # result <- analyzeDrugOmicPair("mRNA", "GENE1", "DRUG1")
  # expect_true(is.list(result))
  # expect_true(!is.null(result$plot))
  # 
  # # Restore the original function
  # assign("selectFeatures", original, envir = environment(analyzeDrugOmicPair))
}) 