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
library(DROMA)

# Load necessary data files
setupDROMA()
loadDROMA("drug")
loadDROMA("mRNA")
loadDROMA("mut")

# When running tests in the package environment, we should use the package functions directly
# rather than sourcing the files

context("Batch Feature Analysis Functions")

# Setup real data for tests
# Define example features to use
select_drugs <- "Paclitaxel"  # Example drug
select_omics_type_continuous <- "mRNA"    # Example continuous omics type
select_omics_continuous <- "ABCB1"        # Example continuous omics feature
select_omics_type_discrete <- "mutation_gene"  # Example discrete omics type
select_omics_discrete <- "TP53"           # Example discrete omics feature

# Test meta-analysis calculation functions
test_that("metaCalcConCon calculates meta-analysis for continuous data", {
  skip_if_not_installed("DROMA")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  selected_pair <- pairContinuousFeatures(myOmics, myDrugs)

  # Calculate meta-analysis
  result <- metaCalcConCon(selected_pair)

  # Verify result
  expect_true(!is.null(result))
  expect_true(inherits(result, "meta"))
})

test_that("metaCalcConDis calculates meta-analysis for continuous vs. discrete data", {
  skip_if_not_installed("DROMA")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_discrete, select_omics_discrete,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  selected_pair <- pairDiscreteFeatures(myOmics, myDrugs)

  # Calculate meta-analysis
  result <- metaCalcConDis(selected_pair)

  # Verify result
  expect_true(!is.null(result))
  expect_true(inherits(result, "meta"))
})

test_that("metaCalcDisDis calculates meta-analysis for discrete vs. discrete data", {
  skip_if_not_installed("DROMA")

  # Get real data for two discrete features
  my_feas1 <- selectFeatures("mutation_gene", "TP53",
                          data_type = "all", tumor_type = "all")
  my_feas2 <- selectFeatures("mutation_gene", "ABCB1",
                          data_type = "all", tumor_type = "all")

  # Pair the data
  selected_pair <- pairDiscreteDiscrete(my_feas1, my_feas2,
                                      feature1_type = "mutation_gene",
                                      feature2_type = "mutation_gene",
                                      samples_search = samples_search)

  # Skip if no valid pairs were created
  skip_if(length(selected_pair) == 0, "No valid discrete-discrete pairs could be created with the test data")

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

# Test pairing functions
test_that("pairContinuousFeatures correctly pairs continuous data", {
  skip_if_not_installed("DROMA")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  pairs <- pairContinuousFeatures(myOmics, myDrugs)

  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)

  # Check specific pairs
  first_pair_name <- names(pairs)[1]
  pair <- pairs[[first_pair_name]]
  expect_true(is.list(pair))
  expect_equal(length(pair[[1]]), length(pair[[2]]))
  expect_equal(names(pair[[1]]), names(pair[[2]]))
})

test_that("pairDiscreteFeatures correctly pairs discrete and continuous data", {
  skip_if_not_installed("DROMA")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_discrete, select_omics_discrete,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  pairs <- pairDiscreteFeatures(myOmics, myDrugs)

  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)

  # Check specific pairs
  first_pair_name <- names(pairs)[1]
  pair <- pairs[[first_pair_name]]
  expect_true(is.list(pair))
  expect_true(all(c("yes", "no") %in% names(pair)))
})

test_that("pairDiscreteDiscrete correctly pairs two discrete features", {
  skip_if_not_installed("DROMA")

  # Get real data for two discrete features
  my_feas1 <- selectFeatures("mutation_gene", "TP53",
                          data_type = "all", tumor_type = "all")
  my_feas2 <- selectFeatures("mutation_gene", "KRAS",
                          data_type = "all", tumor_type = "all")

  # Pair the data
  pairs <- pairDiscreteDiscrete(my_feas1, my_feas2,
                              feature1_type = "mutation_gene",
                              feature2_type = "mutation_gene",
                              samples_search = samples_search)

  # Skip if no valid pairs were created
  skip_if(length(pairs) == 0, "No valid discrete-discrete pairs could be created with the test data")

  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)

  # Check specific pairs
  first_pair_name <- names(pairs)[1]
  pair <- pairs[[first_pair_name]]
  expect_true("cont_table" %in% names(pair))
  expect_true(is.matrix(pair$cont_table))
  expect_equal(dim(pair$cont_table), c(2, 2))
})

# Test plotting function
test_that("plotMetaVolcano creates a volcano plot", {
  skip_if_not_installed("DROMA")
  skip_if_not_installed("ggplot2")

  # Create a sample meta-analysis results dataframe
  meta_df <- data.frame(
    effect_size = runif(20, -0.8, 0.8),
    p_value = 10^(-runif(20, 0, 5)),
    name = paste0("Gene_", 1:20)
  )

  # Create plot
  plot <- plotMetaVolcano(meta_df, es_t = 0.4, P_t = 0.05)

  # Verify plot was created
  expect_true(inherits(plot, "ggplot"))
})

# Test the main batch analysis function
test_that("batchFindSignificantFeatures handles real data", {
  # skip("This test takes too long to run in a standard test suite")
  skip_if_not_installed("DROMA")

  # Run batch analysis with limited scope
  results <- batchFindSignificantFeatures(
    feature1_type = "mRNA",
    feature1_name = "ABCB1",
    feature2_type = "mRNA",
    data_type = "all",
    tumor_type = "all",
    cores = 1,
    test_top_100 = TRUE
  )

  # Verify results structure
  expect_true(is.data.frame(results))
  expect_true(all(c("p_value", "effect_size", "name") %in% colnames(results)))
})
