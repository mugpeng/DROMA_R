#!/usr/bin/env Rscript

# Example script for DROMA package: Data Loading Scenarios
# This example demonstrates different ways to load data in the DROMA package

# Load required libraries
library(DROMA)

######################################
# Scenario 1: Basic setup and selective loading
######################################

# Setup DROMA environment - loads search vectors and sample annotations only
setupDROMA()
cat("After setupDROMA() - available objects in global environment:\n")
print(ls(pattern = "search|anno"))

# Selectively load only drug data
loadDROMA("drug")
cat("\nAfter loadDROMA('drug') - drug-related objects in global environment:\n")
print(ls(pattern = "drug"))

# Selectively load only mRNA data
loadDROMA("mRNA")
cat("\nAfter loadDROMA('mRNA') - mRNA-related objects in global environment:\n")
print(ls(pattern = "mRNA"))

######################################
# Scenario 2: Load mutation data with normalization disabled
######################################

# Normalization has no effect on mutation data as it's discrete
cat("\nLoading mutation data (normalization has no effect on discrete data):\n")
loadDROMA("mut", apply_normalization = FALSE)
cat("Mutation-related objects in global environment:\n")
print(ls(pattern = "mut"))

######################################
# Scenario 3: Load all data at once
######################################

# Clear environment first for demonstration (not necessary in actual use)
# Note: In a real script, you wouldn't typically do this
rm(list = ls())
library(DROMA)

# Setup environment
setupDROMA()

# Load all data types at once
cat("\nLoading all data types at once:\n")
loadDROMA("all")

# Check what objects were loaded
cat("\nAfter loadDROMA('all') - data objects in global environment:\n")
data_objects <- c(
  "Drug data:" = length(ls(pattern = "drug")),
  "mRNA data:" = length(ls(pattern = "mRNA")),
  "Mutation data:" = length(ls(pattern = "mut")),
  "Methylation data:" = length(ls(pattern = "meth")),
  "Protein data:" = length(ls(pattern = "protein")),
  "CNV data:" = length(ls(pattern = "cnv")),
  "Fusion data:" = length(ls(pattern = "fusion"))
)
print(data_objects)

######################################
# Scenario 4: Check available drugs after loading
######################################

cat("\nAvailable drugs (first 5):\n")
if(exists("drugs_search")) {
  print(head(drugs_search$name, 5))
} else {
  cat("drugs_search not found\n")
}

cat("\nAvailable omics types:\n")
if(exists("omics_search")) {
  print(unique(omics_search$type))
} else {
  cat("omics_search not found\n")
}
