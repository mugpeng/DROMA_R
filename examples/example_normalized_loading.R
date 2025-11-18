#!/usr/bin/env Rscript

# Example script for DROMA.R package: Normalized Data Loading with DromaSet Objects
# This example demonstrates how to load data with automatic z-score normalization

# Load required libraries
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions

######################################
# Setup: Create DromaSet Objects
######################################

# Note: Replace with your actual database path
# db_path <- "path/to/your/droma.sqlite"
db_path <- "sql_db/droma.sqlite"

# Connect to DROMA database
connectDROMADatabase(db_path)

# Create DromaSet objects
cat("Creating DromaSet objects...\n")
gCSI <- createDromaSetFromDatabase("gCSI", db_path)

# Create a MultiDromaSet for cross-project analysis
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE"),
  db_path = db_path
)

cat("DromaSet objects created successfully!\n\n")

######################################
# Example 1: Load molecular profiles with z-score normalization (default)
######################################

cat("Example 1: Loading mRNA data with z-score normalization (default)\n")

# Load ABCB1 mRNA data with z-score normalization (default behavior)
abcb1_normalized <- loadMolecularProfiles(
  CCLE,
  molecular_type = "mRNA",
  features = "ABCB1"
)

if (!is.null(abcb1_normalized)) {
  cat("Loaded ABCB1 mRNA data with z-score normalization\n")
  cat("Data dimensions:", dim(abcb1_normalized), "\n")
  cat("Is normalized:", isZscoreNormalized(abcb1_normalized), "\n")
  cat("Sample of normalized values:", head(abcb1_normalized[1, ], 3), "\n\n")
} else {
  cat("No ABCB1 mRNA data found\n\n")
}

######################################
# Example 2: Load molecular profiles without z-score normalization
######################################

cat("Example 2: Loading mRNA data without z-score normalization\n")

# Load ABCB1 mRNA data without z-score normalization
abcb1_raw <- loadMolecularProfiles(
  gCSI,
  molecular_type = "mRNA",
  features = "ABCB1",
  zscore = FALSE
)

if (!is.null(abcb1_raw)) {
  cat("Loaded ABCB1 mRNA data without z-score normalization\n")
  cat("Data dimensions:", dim(abcb1_raw), "\n")
  cat("Is normalized:", isZscoreNormalized(abcb1_raw), "\n")
  cat("Sample of raw values:", head(abcb1_raw[1, ], 3), "\n\n")
} else {
  cat("No ABCB1 mRNA data found\n\n")
}

######################################
# Example 3: Compare normalized vs raw data
######################################

if (!is.null(abcb1_normalized) && !is.null(abcb1_raw)) {
  cat("Example 3: Comparing normalized vs raw data\n")

  # Extract values for the same samples
  common_samples <- intersect(colnames(abcb1_normalized), colnames(abcb1_raw))

  if (length(common_samples) > 0) {
    normalized_values <- abcb1_normalized[1, common_samples]
    raw_values <- abcb1_raw[1, common_samples]

    cat("Raw data statistics:\n")
    cat("  Mean:", round(mean(raw_values, na.rm = TRUE), 3), "\n")
    cat("  SD:", round(sd(raw_values, na.rm = TRUE), 3), "\n")
    cat("  Range:", round(range(raw_values, na.rm = TRUE), 3), "\n")

    cat("Normalized data statistics:\n")
    cat("  Mean:", round(mean(normalized_values, na.rm = TRUE), 3), "\n")
    cat("  SD:", round(sd(normalized_values, na.rm = TRUE), 3), "\n")
    cat("  Range:", round(range(normalized_values, na.rm = TRUE), 3), "\n\n")
  }
}

######################################
# Example 4: Load treatment response data with normalization
######################################

cat("Example 4: Loading treatment response data with z-score normalization\n")

# Load Paclitaxel drug data with z-score normalization (default)
paclitaxel_normalized <- loadTreatmentResponse(
  gCSI,
  drugs = "Paclitaxel"
)

if (!is.null(paclitaxel_normalized)) {
  cat("Loaded Paclitaxel data with z-score normalization\n")
  cat("Data dimensions:", dim(paclitaxel_normalized), "\n")
  cat("Is normalized:", isZscoreNormalized(paclitaxel_normalized), "\n")
  cat("Sample of normalized values:", head(paclitaxel_normalized[1, ], 3), "\n\n")
} else {
  cat("No Paclitaxel data found\n\n")
}

######################################
# Example 5: Multi-project loading with normalization
######################################

cat("Example 5: Loading multi-project data with z-score normalization\n")

# Load ABCB1 mRNA data across multiple projects with normalization
multi_abcb1_normalized <- loadMultiProjectMolecularProfiles(
  multi_set,
  molecular_type = "mRNA",
  features = "ABCB1",
  overlap_only = FALSE
)

if (!is.null(multi_abcb1_normalized) && length(multi_abcb1_normalized) > 0) {
  cat("Loaded ABCB1 data from", length(multi_abcb1_normalized), "projects with normalization\n")

  for (project_name in names(multi_abcb1_normalized)) {
    project_data <- multi_abcb1_normalized[[project_name]]
    if (!is.null(project_data)) {
      cat("  ", project_name, "- Dimensions:", dim(project_data),
          "- Normalized:", isZscoreNormalized(project_data), "\n")
    }
  }
  cat("\n")
} else {
  cat("No multi-project ABCB1 data found\n\n")
}

######################################
# Example 6: Multi-project treatment response with normalization
######################################

cat("Example 6: Loading multi-project treatment response with z-score normalization\n")

# Load Paclitaxel data across multiple projects with normalization
multi_paclitaxel_normalized <- loadMultiProjectTreatmentResponse(
  multi_set,
  drugs = "Paclitaxel",
  overlap_only = FALSE
)

if (!is.null(multi_paclitaxel_normalized) && length(multi_paclitaxel_normalized) > 0) {
  cat("Loaded Paclitaxel data from", length(multi_paclitaxel_normalized), "projects with normalization\n")

  for (project_name in names(multi_paclitaxel_normalized)) {
    project_data <- multi_paclitaxel_normalized[[project_name]]
    if (!is.null(project_data)) {
      cat("  ", project_name, "- Dimensions:", dim(project_data),
          "- Normalized:", isZscoreNormalized(project_data), "\n")
    }
  }
  cat("\n")
} else {
  cat("No multi-project Paclitaxel data found\n\n")
}

######################################
# Example 7: Apply normalization to existing data
######################################

cat("Example 7: Applying z-score normalization to existing data\n")

# Load raw data first
raw_mrna_data <- loadMolecularProfiles(
  gCSI,
  molecular_type = "mRNA",
  features = c("ABCB1", "TP53", "BRCA1"),
  zscore = FALSE
)

if (!is.null(raw_mrna_data)) {
  cat("Loaded raw mRNA data for multiple genes\n")
  cat("Original data - Normalized:", isZscoreNormalized(raw_mrna_data), "\n")

  # Apply normalization to existing data
  manually_normalized <- applyZscoreNormalization(raw_mrna_data)

  cat("After manual normalization - Normalized:", isZscoreNormalized(manually_normalized), "\n")

  # Compare statistics
  cat("Raw data mean (first gene):", round(mean(raw_mrna_data[1, ], na.rm = TRUE), 3), "\n")
  cat("Normalized data mean (first gene):", round(mean(manually_normalized[1, ], na.rm = TRUE), 3), "\n\n")
} else {
  cat("No mRNA data found for manual normalization example\n\n")
}

######################################
# Example 8: Handle discrete data types (mutations)
######################################

cat("Example 8: Attempting normalization on discrete data (mutations)\n")

# Try to load mutation data with normalization (should give warning)
tp53_mutations <- loadMolecularProfiles(
  gCSI,
  molecular_type = "mutation_gene",
  features = "TP53",
  zscore = TRUE  # This should trigger a warning
)

if (!is.null(tp53_mutations)) {
  cat("Loaded TP53 mutation data\n")
  cat("Data type:", class(tp53_mutations), "\n")
  cat("Is normalized:", isZscoreNormalized(tp53_mutations), "\n\n")
} else {
  cat("No TP53 mutation data found\n\n")
}

######################################
# Example 9: Load all features with normalization
######################################

cat("Example 9: Loading all mRNA features with z-score normalization\n")

# Load all mRNA features (this might take time with large datasets)
# In practice, you might want to limit this or use specific features
all_mrna_normalized <- loadMolecularProfiles(
  gCSI,
  molecular_type = "mRNA",
  # features = NULL means load all features
  zscore = TRUE
)

if (!is.null(all_mrna_normalized)) {
  cat("Loaded all mRNA features with normalization\n")
  cat("Data dimensions:", dim(all_mrna_normalized), "\n")
  cat("Is normalized:", isZscoreNormalized(all_mrna_normalized), "\n")
  cat("Sample gene names:", head(rownames(all_mrna_normalized), 3), "\n\n")
} else {
  cat("No mRNA data found\n\n")
}

######################################
# Example 10: Using normalized data in analysis
######################################

cat("Example 10: Using normalized data in drug-omics analysis\n")

# The existing analyzeDrugOmicPair function can work with normalized data
# by using the normalized loading functions internally or by passing pre-normalized data

if (!is.null(multi_abcb1_normalized) && !is.null(multi_paclitaxel_normalized)) {
  cat("Normalized data is ready for use in analyzeDrugOmicPair or other analysis functions\n")
  cat("The normalization ensures data is on comparable scales across projects\n")

  # Example of how this would integrate with existing analysis
  cat("Example usage: analyzeDrugOmicPair(multi_set, 'mRNA', 'ABCB1', 'Paclitaxel')\n")
  cat("(The analysis functions can be updated to use normalized loading internally)\n\n")
} else {
  cat("Normalized data not available for analysis example\n\n")
}

cat("Normalized data loading examples completed!\n")
cat("Key takeaways:\n")
cat("1. Use loadMolecularProfiles() for normalized molecular data\n")
cat("2. Use loadTreatmentResponse() for normalized drug response data\n")
cat("3. Use loadMultiProject*Normalized() for multi-project normalized data\n")
cat("4. Set zscore=FALSE to disable normalization when needed\n")
cat("5. Use isZscoreNormalized() to check if data has been normalized\n")
cat("6. Use applyZscoreNormalization() to normalize existing data\n")
cat("7. Z-score normalization is automatically skipped for discrete data types\n")
cat("8. Normalized data ensures comparable scales across different projects\n")
