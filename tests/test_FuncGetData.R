#!/usr/bin/env Rscript

# Test script for FuncGetData.R functions
# Tests selFeatures, mergeDrugFeatures, and annoMergeFeatures

# Load required packages
library(testthat)

# Source the function file
source("Package_Function/FuncGetData.R")

# Test for selFeatures function
test_that("selFeatures returns correct data for drug", {
  select_drugs = "5-Fluorouracil"
  data_type = "all"
  tumor_type = "all"
  
  myDrugs <- selFeatures("drug", select_drugs, 
                         data_type = data_type, 
                         tumor_type = tumor_type)
  
  # Check that result is a list
  expect_true(is.list(myDrugs))
  
  # Check that it's not empty
  expect_gt(length(myDrugs), 0)
  
  # Check that each element is a numeric vector
  for(i in 1:length(myDrugs)) {
    expect_true(is.numeric(myDrugs[[i]]))
  }
})

test_that("selFeatures returns correct data for mRNA", {
  select_omics_type = "mRNA"
  select_omics = "ABCB1"
  data_type = "all"
  tumor_type = "all"
  
  myOmics <- selFeatures(select_omics_type, select_omics,
                         data_type = data_type, 
                         tumor_type = tumor_type)
  
  # Check that result is a list
  expect_true(is.list(myOmics))
  
  # Check that it's not empty
  expect_gt(length(myOmics), 0)
  
  # Check that each element is a numeric vector
  for(i in 1:length(myOmics)) {
    expect_true(is.numeric(myOmics[[i]]))
  }
})

test_that("selFeatures handles invalid feature type", {
  expect_error(selFeatures("invalid_type", "ABCB1"), 
               "The select feature type doesn't exsit")
})

test_that("selFeatures handles invalid data_type", {
  expect_error(selFeatures("mRNA", "ABCB1", data_type = "invalid_type"), 
               "Invalid data_type")
})
