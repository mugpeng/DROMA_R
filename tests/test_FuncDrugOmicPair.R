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
library(patchwork)

# Source the function files
source("R/FuncDrugOmicPair.R")
source("R/FuncGetData.R") # Needed for selectFeatures function

context("Drug-Omic Pairing Functions")

# Setup real data for tests
# Define example features to use
select_drugs <- "Paclitaxel"  # Example drug
select_omics_type_continuous <- "mRNA"    # Example continuous omics type
select_omics_continuous <- "ABCB1"        # Example continuous omics feature
select_omics_type_discrete <- "mutation_gene"  # Example discrete omics type
select_omics_discrete <- "TP53"           # Example discrete omics feature

# Test for pairDrugOmic function
test_that("pairDrugOmic correctly pairs drug and omics data", {
  # Load example data
  load("data/drug.Rda")
  load("data/mRNA.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  pairs <- pairDrugOmic(myOmics, myDrugs)

  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)

  # Check that pairs have the expected structure
  first_pair_name <- names(pairs)[1]
  pair <- pairs[[first_pair_name]]
  expect_true(is.list(pair))
  expect_true(all(c("omic", "drug") %in% names(pair)))
  expect_equal(length(pair$omic), length(pair$drug))
  expect_equal(names(pair$omic), names(pair$drug))

  # Test with merged = TRUE
  merged_pairs <- pairDrugOmic(myOmics, myDrugs, merged = TRUE)
  expect_true("merged_dataset" %in% names(merged_pairs))
})

test_that("analyzeContinuousDrugOmic correctly analyzes pairs", {
  # Load example data
  load("data/drug.Rda")
  load("data/mRNA.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  myPairs <- pairDrugOmic(myOmics, myDrugs)

  # Analyze the pairs
  results <- analyzeContinuousDrugOmic(myPairs)

  # Verify results
  expect_true(!is.null(results))
  expect_true(inherits(results, "meta"))
})

# Test for pairDiscreteDrugOmic function
test_that("pairDiscreteDrugOmic correctly pairs discrete data", {
  # Load example data
  load("data/drug.Rda")
  load("data/mut.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_discrete, select_omics_discrete,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  pairs <- pairDiscreteDrugOmic(myOmics, myDrugs)

  # Verify structure and content
  expect_true(is.list(pairs))
  expect_true(length(pairs) > 0)

  # Check that pairs have the expected structure
  first_pair_name <- names(pairs)[1]
  pair <- pairs[[first_pair_name]]
  expect_true(is.list(pair))
  expect_true(all(c("yes", "no") %in% names(pair)))

  # Test with merged = TRUE
  merged_pairs <- pairDiscreteDrugOmic(myOmics, myDrugs, merged = TRUE)
  expect_true("merged_dataset" %in% names(merged_pairs))
})

test_that("analyzeDiscreteDrugOmic correctly analyzes pairs", {
  # Load example data
  load("data/drug.Rda")
  load("data/mut.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_discrete, select_omics_discrete,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  myPairs <- pairDiscreteDrugOmic(myOmics, myDrugs)

  # Analyze the pairs
  results <- analyzeDiscreteDrugOmic(myPairs)

  # Verify results
  expect_true(!is.null(results))
  expect_true(inherits(results, "meta"))
})

# Test plotting functions with real data
test_that("createForestPlot creates a forest plot", {
  # Load example data
  load("data/drug.Rda")
  load("data/mRNA.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data and analyze
  myPairs <- pairDrugOmic(myOmics, myDrugs)
  meta_obj <- analyzeContinuousDrugOmic(myPairs)

  # Skip if meta-analysis failed
  if (is.null(meta_obj)) {
    skip("Meta-analysis did not produce results with the test data")
  }

  # Create forest plot
  plot <- createForestPlot(meta_obj)

  # Verify plot was created (no easy way to test the plot itself)
  expect_true(TRUE)
})

test_that("plotContinuousDrugOmic creates a scatter plot", {
  # Load example data
  load("data/drug.Rda")
  load("data/mRNA.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_continuous, select_omics_continuous,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  myPairs <- pairDrugOmic(myOmics, myDrugs)

  # Get first pair for plotting
  first_pair_name <- names(myPairs)[1]
  pair <- myPairs[[first_pair_name]]

  # Create plot
  plot <- plotContinuousDrugOmic(pair$omic, pair$drug, first_pair_name)

  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

test_that("plotDiscreteDrugOmic creates a boxplot", {
  # Load example data
  load("data/drug.Rda")
  load("data/mut.Rda")
  load("data/anno.Rda")

  # Get real data using selectFeatures
  myOmics <- selectFeatures(select_omics_type_discrete, select_omics_discrete,
                         data_type = "all", tumor_type = "all")
  myDrugs <- selectFeatures("drug", select_drugs,
                         data_type = "all", tumor_type = "all")

  # Pair the data
  myPairs <- pairDiscreteDrugOmic(myOmics, myDrugs)

  # Get first pair for plotting
  first_pair_name <- names(myPairs)[1]
  pair <- myPairs[[first_pair_name]]

  # Create plot
  plot <- plotDiscreteDrugOmic(pair$yes, pair$no, first_pair_name)

  # Verify plot is created
  expect_true(inherits(plot, "ggplot"))
})

# Test the main analysis function
test_that("analyzeDrugOmicPair handles continuous data", {
  # Load example data
  load("data/drug.Rda")
  load("data/mRNA.Rda")
  load("data/anno.Rda")

  # Call the main analysis function
  result <- analyzeDrugOmicPair(
    select_omics_type = select_omics_type_continuous,
    select_omics = select_omics_continuous,
    select_drugs = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Verify results
  expect_true(is.list(result))
  expect_true(!is.null(result$plot))
  expect_true(!is.null(result$data))
})

test_that("analyzeDrugOmicPair handles discrete data", {
  # Load example data
  load("data/drug.Rda")
  load("data/mut.Rda")
  load("data/anno.Rda")

  # Call the main analysis function
  result <- analyzeDrugOmicPair(
    select_omics_type = select_omics_type_discrete,
    select_omics = select_omics_discrete,
    select_drugs = select_drugs,
    data_type = "all",
    tumor_type = "all"
  )

  # Verify results
  expect_true(is.list(result))
  expect_true(!is.null(result$plot))
  expect_true(!is.null(result$data))
})
