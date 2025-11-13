#!/usr/bin/env Rscript

# Example script for DROMA.R package: Feature Selection with DromaSet Objects
# This example demonstrates how to retrieve omics and drug data using DromaSet and MultiDromaSet objects

# Load required libraries
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions
library(dplyr)

# Setup: Create DromaSet Objects ----

# Note: Replace with your actual database path
db_path <- "../Data/droma.sqlite"
# db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"

# Connect to DROMA database
connectDROMADatabase(db_path)

# Create individual DromaSet objects
cat("Creating DromaSet objects...\n")
gCSI <- createDromaSetFromDatabase("gCSI", db_path)
CCLE <- createDromaSetFromDatabase("CCLE", db_path)

# Create a MultiDromaSet for cross-project analysis
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE", "GDSC"),
  db_path = db_path
)

all_sets <- createMultiDromaSetFromAllProjects(db_path = db_path)

cat("DromaSet objects created successfully!\n")

# Basic Feature Selection Examples ----

# Example 1: Load drug data (Paclitaxel) from a single project
cat("\nExample 1: Loading drug data for Paclitaxel from gCSI\n")
paclitaxel_data <- loadTreatmentResponse(gCSI,
                                         select_drugs = "Paclitaxel",
                                         return_data = TRUE)

if (is.matrix(paclitaxel_data) && "Paclitaxel" %in% rownames(paclitaxel_data)) {
  drug_vector <- as.numeric(paclitaxel_data["Paclitaxel", ])
  names(drug_vector) <- colnames(paclitaxel_data)
  drug_vector <- drug_vector[!is.na(drug_vector)]
  cat(sprintf("Retrieved Paclitaxel data: %d samples\n", length(drug_vector)))
} else {
  cat("Paclitaxel not found in treatment response data\n")
}

# Example 2: Load mRNA data (ABCB1) from a single project
cat("\nExample 2: Loading mRNA data for ABCB1 from gCSI\n")
abcb1_data <- loadMolecularProfiles(CCLE,
                                    molecular_type = "mRNA",
                                    features = "ABCB1",
                                    return_data = TRUE)

if (is.matrix(abcb1_data) && "ABCB1" %in% rownames(abcb1_data)) {
  mrna_vector <- as.numeric(abcb1_data["ABCB1", ])
  names(mrna_vector) <- colnames(abcb1_data)
  mrna_vector <- mrna_vector[!is.na(mrna_vector)]
  cat(sprintf("Retrieved ABCB1 data: %d samples\n", length(mrna_vector)))
} else {
  cat("ABCB1 not found in mRNA data\n")
}

# Example 3: Load mutation data (TP53) from a single project
cat("\nExample 3: Loading mutation data for TP53 from gCSI\n")
tp53_data <- loadMolecularProfiles(gCSI,
                                   molecular_type = "mutation_gene",
                                   features = "TP53",
                                   return_data = TRUE)

if (is.data.frame(tp53_data) && "features" %in% colnames(tp53_data)) {
  tp53_samples <- tp53_data$samples[tp53_data$features == "TP53"]
  cat(sprintf("Retrieved TP53 mutation data: %d samples with mutations\n", length(tp53_samples)))
} else {
  cat("TP53 mutation data not found\n")
}

# Multi-Project Feature Selection ----

# Example 4: Load drug data across multiple projects
cat("\nExample 4: Loading Paclitaxel data across multiple projects\n")
multi_paclitaxel <- loadMultiProjectTreatmentResponseNormalized(multi_set,
                                                     drugs = "Paclitaxel")

cat("Retrieved Paclitaxel data from", length(multi_paclitaxel), "projects:\n")
for (project_name in names(multi_paclitaxel)) {
  if (is.matrix(multi_paclitaxel[[project_name]]) && "Paclitaxel" %in% rownames(multi_paclitaxel[[project_name]])) {
    drug_vector <- as.numeric(multi_paclitaxel[[project_name]]["Paclitaxel", ])
    drug_vector <- drug_vector[!is.na(drug_vector)]
    cat(sprintf("- %s: %d samples\n", project_name, length(drug_vector)))
  }
}

# Example 5: Load mRNA data across multiple projects
cat("\nExample 5: Loading ABCB1 mRNA data across multiple projects\n")
multi_abcb1 <- loadMultiProjectMolecularProfilesNormalized(multi_set,
                                                 molecular_type = "mRNA",
                                                 features = "ABCB1")

cat("Retrieved ABCB1 data from", length(multi_abcb1), "projects:\n")
for (project_name in names(multi_abcb1)) {
  if (is.matrix(multi_abcb1[[project_name]]) && "ABCB1" %in% rownames(multi_abcb1[[project_name]])) {
    mrna_vector <- as.numeric(multi_abcb1[[project_name]]["ABCB1", ])
    mrna_vector <- mrna_vector[!is.na(mrna_vector)]
    cat(sprintf("- %s: %d samples\n", project_name, length(mrna_vector)))
  }
}

# Filtering by Data Type and Tumor Type ----

# Example 6: Get Paclitaxel data only for cell lines
cat("\nExample 6: Loading Paclitaxel data only for cell lines\n")
paclitaxel_cellline <- loadTreatmentResponseNormalized(gCSI,
                                            drugs = "Paclitaxel",
                                            data_type = "CellLine",
                                            return_data = TRUE)

if (is.matrix(paclitaxel_cellline) && "Paclitaxel" %in% rownames(paclitaxel_cellline)) {
  drug_vector <- as.numeric(paclitaxel_cellline["Paclitaxel", ])
  drug_vector <- drug_vector[!is.na(drug_vector)]
  cat(sprintf("Retrieved Paclitaxel cell line data: %d samples\n", length(drug_vector)))
} else {
  cat("No Paclitaxel cell line data found\n")
}

# Example 7: Get ABCB1 data for a specific tumor type
cat("\nExample 7: Loading ABCB1 data for breast cancer\n")
abcb1_breast <- loadMolecularProfiles(gCSI,
                                     molecular_type = "mRNA",
                                     features = "ABCB1",
                                     tumor_type = "breast cancer",
                                     return_data = TRUE)

if (is.matrix(abcb1_breast) && "ABCB1" %in% rownames(abcb1_breast)) {
  mrna_vector <- as.numeric(abcb1_breast["ABCB1", ])
  mrna_vector <- mrna_vector[!is.na(mrna_vector)]
  cat(sprintf("Retrieved ABCB1 breast cancer data: %d samples\n", length(mrna_vector)))
} else {
  cat("No ABCB1 breast cancer data found\n")
}

# Using DROMA.R Functions with DromaSet Objects ----

# Example 8: Use processDrugData with DromaSet
cat("\nExample 8: Using processDrugData with DromaSet object\n")
drug_sensitivity_data <- processDrugData(
  gCSI,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

if (!is.null(drug_sensitivity_data)) {
  cat(sprintf("Processed drug sensitivity data: %d samples\n", nrow(drug_sensitivity_data)))
  cat("Columns:", paste(colnames(drug_sensitivity_data), collapse = ", "), "\n")
} else {
  cat("No drug sensitivity data processed\n")
}

# Example 9: Use getDrugSensitivityData with MultiDromaSet
cat("\nExample 9: Using getDrugSensitivityData with MultiDromaSet object\n")
multi_drug_data <- getDrugSensitivityData(
  multi_set,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  include_annotations = FALSE
)

if (!is.null(multi_drug_data)) {
  cat(sprintf("Retrieved drug sensitivity data from multiple projects: %d samples\n", nrow(multi_drug_data)))
  cat("Studies included:", paste(unique(multi_drug_data$study), collapse = ", "), "\n")
} else {
  cat("No multi-project drug sensitivity data retrieved\n")
}

# Error Handling Examples ----

# Example 10: Handle non-existent features gracefully
cat("\nExample 10: Handling non-existent features\n")
tryCatch({
  nonexistent_data <- loadTreatmentResponseNormalized(gCSI,
                                            drugs = "NonExistentDrug",
                                            return_data = TRUE)

  if (is.matrix(nonexistent_data) && "NonExistentDrug" %in% rownames(nonexistent_data)) {
    cat("Found non-existent drug (unexpected)\n")
  } else {
    cat("Non-existent drug not found in data (as expected)\n")
  }
}, error = function(e) {
  cat("Error caught:", conditionMessage(e), "\n")
})

# Example 11: Handle invalid molecular types
cat("\nExample 11: Handling invalid molecular types\n")
tryCatch({
  invalid_data <- loadMolecularProfiles(gCSI,
                                        molecular_type = "invalid_type",
                                        features = "ABCB1",
                                        return_data = TRUE)
}, error = function(e) {
  cat("Error caught (as expected):", conditionMessage(e), "\n")
})

cat("\nFeature selection examples completed!\n")
cat("Key takeaways:\n")
cat("1. Use createDromaSetFromDatabase() to create single project objects\n")
cat("2. Use createMultiDromaSetFromDatabase() for multi-project analysis\n")
cat("3. Use loadTreatmentResponseNormalized() for drug data\n")
cat("4. Use loadMolecularProfiles() for omics data\n")
cat("5. Use DROMA.R functions like processDrugData() with DromaSet objects\n")
