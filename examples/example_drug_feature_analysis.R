#!/usr/bin/env Rscript

# Example script for DROMA.R package: Drug Feature Analysis with DromaSet Objects
# This example demonstrates drug data processing and visualization functions using DromaSet objects

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DT)
library(htmltools)
library(patchwork)
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions

######################################
# Setup: Create DromaSet Objects
######################################

# Note: Replace with your actual database path
db_path <- "sql_db/droma.sqlite"

# Connect to DROMA database
connectDROMADatabase(db_path)

# Create DromaSet objects
cat("Creating DromaSet objects...\n")
gCSI <- createDromaSetFromDatabase("gCSI", db_path)
CCLE <- createDromaSetFromDatabase("CCLE", db_path)

# Create a MultiDromaSet for cross-project analysis
multi_set <- createMultiDromaSetFromDatabase(
  project_names = c("gCSI", "CCLE"),
  db_path = db_path
)

cat("DromaSet objects created successfully!\n\n")

######################################
# Example 1: Basic Drug Data Processing
######################################

# Get processed drug data for Paclitaxel using DromaSet
cat("Example 1: Processing data for Paclitaxel using DromaSet\n")
paclitaxel_data <- processDrugData(
  gCSI,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display summary information
if (!is.null(paclitaxel_data)) {
  cat("Retrieved Paclitaxel data for", nrow(paclitaxel_data), "samples\n")
  cat("Data columns:", paste(colnames(paclitaxel_data), collapse = ", "), "\n")
  cat("Number of studies:", length(unique(paclitaxel_data$study)), "\n")
  cat("Studies:", paste(unique(paclitaxel_data$study), collapse = ", "), "\n\n")
} else {
  cat("No Paclitaxel data found\n\n")
}

######################################
# Example 2: Annotating Drug Data
######################################

# Annotate the drug data with sample information
cat("Example 2: Annotating Paclitaxel data with sample information\n")
if (!is.null(paclitaxel_data)) {
  cat("\nDemonstrating direct database loading of sample annotations:\n")
  annotated_data_db <- annotateDrugData(paclitaxel_data, db_path = db_path)

  # Display summaries with annotations
  if (!is.null(annotated_data_db) && "TumorType" %in% colnames(annotated_data_db)) {
    tumor_counts <- table(annotated_data_db$TumorType)
    cat("Tumor Types in the dataset:\n")
    print(tumor_counts)
  }

  if (!is.null(annotated_data_db) && "ModelType" %in% colnames(annotated_data_db)) {
    model_counts <- table(annotated_data_db$ModelType)
    cat("\nModel Types in the dataset:\n")
    print(model_counts)
  }
} else {
  cat("No data available for annotation\n")
}

######################################
# Example 3: Multi-Project Drug Data Processing
######################################

# Get processed drug data for Paclitaxel using MultiDromaSet
cat("\nExample 3: Processing data for Paclitaxel using MultiDromaSet\n")
multi_paclitaxel_data <- processDrugData(
  multi_set,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

# Display summary information
if (!is.null(multi_paclitaxel_data)) {
  cat("Retrieved Paclitaxel data from multiple projects:", nrow(multi_paclitaxel_data), "samples\n")
  cat("Studies included:", paste(unique(multi_paclitaxel_data$study), collapse = ", "), "\n")
} else {
  cat("No multi-project Paclitaxel data found\n")
}

######################################
# Example 4: Formatters and Wrappers
######################################

# Format drug data as a datatable (for interactive use)
cat("\nExample 4: Using the drug data table formatter and wrappers\n")
cat("In an interactive R session, this would display a formatted datatable\n")

# Example of using formatDrugTable
cat("Demonstrating formatDrugTable function (outputs HTML table in interactive session):\n")
if (!is.null(drug_data_with_db)) {
  # Select a subset of columns for display
  display_data <- drug_data_with_db %>%
    select(sampleid, study, zscore_value, raw_value, TumorType, ModelType) %>%
    head(10)  # Just show a few rows for the example

  # In interactive session, this would render an HTML table
  # formatted_table <- formatDrugTable(display_data, caption = "Paclitaxel sensitivity data with annotations")
  cat("Table would show", nrow(display_data), "rows with columns:",
      paste(colnames(display_data), collapse = ", "), "\n")
}

# Using the wrapper with database path for annotations
cat("Using getDrugSensitivityData with database path for annotations:\n")
drug_data_with_db <- getDrugSensitivityData(
  gCSI,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  include_annotations = TRUE,
  db_path = db_path
)

if (!is.null(drug_data_with_db)) {
  cat("Retrieved", nrow(drug_data_with_db), "samples with annotations\n\n")
} else {
  cat("No drug sensitivity data retrieved\n\n")
}

######################################
# Example 5: Visualization Functions
######################################

cat("\nExample 5: Visualization Functions - Continuous and Categorical Comparisons\n")

# Check if we have the required annotations in the data
if (!is.null(drug_data_with_db)) {
  # Continuous comparison: Drug sensitivity vs Age
  if ("Age" %in% colnames(drug_data_with_db)) {
    cat("Creating continuous comparison: Paclitaxel sensitivity vs Age\n")

    # Filter out missing values
    age_data <- drug_data_with_db[!is.na(drug_data_with_db$Age), ]

    if (nrow(age_data) >= 10) {  # Need enough samples for meaningful visualization
      cat("Found", nrow(age_data), "samples with Age data\n")

      # In interactive session, this would display a scatter plot
      p_age_scatter <- plotContinuousComparison(age_data,
                                               cont_column = "Age",
                                               value_column = "zscore_value",
                                               value_label = "Paclitaxel Sensitivity")

      # In interactive session, this would display a boxplot with binned age groups
      p_age_boxplot <- plotContinuousGroups(age_data,
                                            cont_column = "Age",
                                            value_column = "zscore_value",
                                            value_label = "Paclitaxel Sensitivity",
                                            num_bins = 3)

      # In interactive session, this would display a combined plot
      p_age_combined <- createDrugComparisonPlot(age_data,
                                                comparison_var = "Age",
                                                value_column = "zscore_value",
                                                value_label = "Paclitaxel Sensitivity")

      cat("Visualizations would include:\n")
      cat("1. Scatter plot with Spearman correlation between Age and drug sensitivity\n")
      cat("2. Boxplot with Age grouped into bins to compare drug sensitivity\n")
      cat("3. Combined visualization with both plots\n")
    } else {
      cat("Not enough samples with Age data for visualization\n")
    }
  } else {
    cat("Age data not available in the annotations\n")
  }

  # Categorical comparison: Drug sensitivity vs TumorType
  if ("TumorType" %in% colnames(drug_data_with_db)) {
    cat("\nCreating categorical comparison: Paclitaxel sensitivity vs TumorType\n")

    # Count samples per tumor type
    tumor_counts <- table(drug_data_with_db$TumorType)
    valid_tumors <- names(tumor_counts)[tumor_counts >= 3]

    if (length(valid_tumors) >= 2) {
      # Filter data to include only valid tumor types
      tumor_data <- drug_data_with_db[drug_data_with_db$TumorType %in% valid_tumors, ]
      cat("Found", length(valid_tumors), "tumor types with sufficient samples:",
          paste(valid_tumors, collapse = ", "), "\n")

      # In interactive session, this would display a categorical boxplot
       p_tumor <- plotCategoryComparison(tumor_data,
                                       category_column = "TumorType",
                                       value_column = "zscore_value",
                                       value_label = "Paclitaxel Sensitivity")

      cat("Visualization would include a boxplot comparing drug sensitivity across tumor types\n")
      cat("with statistical test to determine significance of differences\n")
    } else {
      cat("Not enough tumor types with sufficient samples for comparison\n")
    }
  } else {
    cat("TumorType data not available in the annotations\n")
  }
} else {
  cat("No annotated drug data available for visualization\n")
}

cat("\nDrug feature analysis examples completed!\n")
cat("Key takeaways:\n")
cat("1. Use processDrugData() with DromaSet objects for drug data processing\n")
cat("2. Use getDrugSensitivityData() for comprehensive drug data retrieval\n")
cat("3. Use loadMolecularProfilesNormalized() to get omics data for visualization\n")
cat("4. Combine drug and omics data for correlation analysis\n")
cat("5. Use MultiDromaSet for cross-project comparisons\n")
