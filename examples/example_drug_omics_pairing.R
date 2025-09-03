#!/usr/bin/env Rscript

# Example script for DROMA.R package: Drug-Omics Pairing Analysis with DromaSet Objects
# This example demonstrates how to analyze associations between drug sensitivity and omics features using DromaSet objects

# Load required libraries
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions
library(ggplot2)
library(metafor)
library(meta)
library(effsize)
library(ggpubr)
library(patchwork)
library(grid)

# Setup: Create DromaSet Objects ----

# Note: Replace with your actual database path
# db_path <- "path/to/your/droma.sqlite"
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

# Check available data types
cat("Available treatment responses in gCSI:\n")
available_drugs <- availableTreatmentResponses(gCSI)
print(available_drugs)

cat("\nAvailable molecular profiles in gCSI:\n")
available_omics <- availableMolecularProfiles(gCSI)
print(available_omics)

# Example 1: Analyze association between Paclitaxel and ABCB1 gene expression (single project) ----

cat("\nExample 1: Analyzing association between Paclitaxel and ABCB1 gene expression (single project)...\n")
result_single <- analyzeDrugOmicPair(
  gCSI,
  select_omics_type = "cnv",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Example 2: Analyze association between Paclitaxel and ABCB1 gene expression (multi-project) ----

cat("\nExample 2: Analyzing association between Paclitaxel and ABCB1 gene expression (multi-project)...\n")
result_multi <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

result_multi <- analyzeDrugOmicPair(
  multi_set_all,
  select_omics_type = "mRNA",
  select_omics = "PSMB5",
  select_drugs = "Bortezomib",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

# Overlap
result_multi2 <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = TRUE
)

# Display results
if (!is.null(result_multi) && !is.null(result_multi$meta)) {
  cat("Meta-analysis results for multi-project analysis:\n")
  print(result_multi$meta)

  # Create forest plot
  cat("Creating forest plot for multi-project analysis...\n")
  createForestPlot(result_multi$meta)

  # Display correlation plots
  if (!is.null(result_multi$plot)) {
    print(result_multi$plot)
  }
} else {
  cat("No significant results found for multi-project analysis\n")
}

# Example 3: Analyze association between Paclitaxel and TP53 mutation status ----

cat("\nExample 3: Analyzing association between Paclitaxel and TP53 mutation status...\n")
result_mut <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mutation_gene",
  select_omics = "TP53",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

# Display results
if (!is.null(result_mut) && !is.null(result_mut$meta)) {
  cat("Meta-analysis results for mutation association:\n")
  print(result_mut$meta)

  # Create forest plot
  cat("Creating forest plot for mutation analysis...\n")
  createForestPlot(result_mut$meta)

  # Display boxplots
  if (!is.null(result_mut$plot)) {
    print(result_mut$plot)
  }
} else {
  cat("No significant results found for mutation analysis\n")
}

# Example 4: Compare drug sensitivity across different data types ----

cat("\nExample 4: Analyzing Paclitaxel sensitivity across different data types...\n")

# Analyze cell lines only
result_cellline <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "CellLine",
  tumor_type = "all",
  overlap_only = FALSE
)

if (!is.null(result_cellline) && !is.null(result_cellline$meta)) {
  cat("Results for cell lines only:\n")
  print(result_cellline$meta)
} else {
  cat("No significant results for cell lines\n")
}

# Analyze PDX models only (if available)
result_pdx <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "PDX",
  tumor_type = "all",
  overlap_only = FALSE
)

if (!is.null(result_pdx) && !is.null(result_pdx$meta)) {
  cat("Results for PDX models only:\n")
  print(result_pdx$meta)
} else {
  cat("No significant results for PDX models\n")
}

# Example 5: Get drug sensitivity data for visualization ----

cat("\nExample 5: Retrieving drug sensitivity data for Paclitaxel...\n")
drug_data <- getDrugSensitivityData(
  multi_set,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

# Display summary of retrieved data
if (!is.null(drug_data)) {
  cat("Number of samples with Paclitaxel sensitivity data:", nrow(drug_data), "\n")
  cat("Studies with Paclitaxel data:\n")
  print(table(drug_data$study))

  # Show data structure
  cat("\nData structure:\n")
  cat("Columns:", paste(colnames(drug_data), collapse = ", "), "\n")
  cat("First few rows:\n")
  print(head(drug_data, 3))
} else {
  cat("No drug sensitivity data retrieved\n")
}

# Example 6: Analyze drug-drug associations----

cat("\nExample 6: Analyzing association between two drugs (Paclitaxel vs Cisplatin)...\n")
drug_drug_result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "drug",
  select_omics = "Cisplatin",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

if (!is.null(drug_drug_result) && !is.null(drug_drug_result$meta)) {
  cat("Drug-drug association results:\n")
  print(drug_drug_result$meta)

  # Create forest plot
  cat("Creating forest plot for drug-drug association...\n")
  createForestPlot(drug_drug_result$meta)

  if (!is.null(drug_drug_result$plot)) {
    print(drug_drug_result$plot)
  }
} else {
  cat("No significant drug-drug association found\n")
}

# Example 7: Tumor type-specific analysis ----

cat("\nExample 7: Tumor type-specific analysis (breast cancer)...\n")
breast_result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "breast cancer",
  overlap_only = FALSE
)

if (!is.null(breast_result) && !is.null(breast_result$meta)) {
  cat("Results for breast cancer samples:\n")
  print(breast_result$meta)

  if (!is.null(breast_result$plot)) {
    print(breast_result$plot)
  }
} else {
  cat("No significant results for breast cancer samples\n")
}

# Example 8: Using different omics types ----

cat("\nExample 8: Analyzing different omics types with Paclitaxel...\n")

# CNV analysis
cnv_result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "cnv",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

if (!is.null(cnv_result) && !is.null(cnv_result$meta)) {
  cat("CNV analysis results:\n")
  print(cnv_result$meta)
} else {
  cat("No significant CNV associations found\n")
}

# Methylation analysis (if available)
meth_result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "meth",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)

if (!is.null(meth_result) && !is.null(meth_result$meta)) {
  cat("Methylation analysis results:\n")
  print(meth_result$meta)
} else {
  cat("No significant methylation associations found\n")
}

# Example 9: Use create_plot_with_common_axes ----
plot_with_axis <- create_plot_with_common_axes(meth_result$plot,
                                               x_title = "Molecular State(mRNA expression or Mutation status)",
                                               y_title = "drug sensitivity (higher indicates resistance)")

# Example 10: Use create_plot_with_common_axes ----
multi_set_all <- createMultiDromaSetFromAllProjects(db_path = db_path)
result_multi_all <- analyzeDrugOmicPair(
  multi_set_all,
  select_omics_type = "mRNA",
  select_omics = "PSMB5",
  select_drugs = "Bortezomib",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)
result_multi_all$plot

result_multi_all2 <- analyzeDrugOmicPair(
  multi_set_all,
  select_omics_type = "mutation_gene",
  select_omics = "PSMB5",
  select_drugs = "Bortezomib",
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE
)


# mRNA
plot_sel <- result_multi_all$plot
plot_sel2 <- as.list(plot_sel)[1:9]
# plot_sel2 <- lapply(plot_sel[1:9], function(x) x[[1]])
plot_sel2 <- wrap_plots(plot_sel2, ncol = 3)
create_plot_with_common_axes(plot_sel2,
                             x_title = "Molecular State(mRNA expression or Mutation status)",
                             y_title = "drug sensitivity (higher indicates resistance)")()

createForestPlot(result_multi_all$meta)

# mut
plot_mut <- as.list(result_multi_all2$plot)[1:2]
plot_mut2 <- wrap_plots(plot_mut, ncol = 2)

# Summary ----
cat("\nDrug-omics pairing analysis examples completed!\n")
cat("Key takeaways:\n")
cat("1. Use analyzeDrugOmicPair() with DromaSet objects for drug-omics analysis\n")
cat("2. MultiDromaSet enables meta-analysis across multiple projects\n")
cat("3. Use createForestPlot() to visualize meta-analysis results\n")
cat("4. Filter by data_type and tumor_type for targeted analyses\n")
cat("5. Analyze different omics types (mRNA, mutations, CNV, methylation)\n")
cat("6. Compare drug-drug associations using the same framework\n")
cat("7. Use getDrugSensitivityData() for detailed drug sensitivity information\n")
