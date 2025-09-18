# Example: Statistical Plots with MultiDromaSet ----

# This example demonstrates how to use the enhanced statistical plotting functions
# that work with MultiDromaSet objects from the DROMA.Set package

# Load required packages ----
library(DROMA.Set)
library(DROMA.R)
library(ggpubr)
library(gg.gap)
library(dplyr)
library(treemapify)

# Example 1: Using MultiDromaSet ----

# Connect to database and create individual DromaSet objects first
# Replace with your actual database path
# db_path <- "path/to/your/droma.sqlite"
# db_path <- "250513-DROMA_R/sql_db/droma.sqlite"
connectDROMADatabase(db_path)

# Create individual DromaSet objects
gCSI <- createDromaSetFromDatabase("gCSI", db_path)
CCLE <- createDromaSetFromDatabase("CCLE", db_path)
GDSC1 <- createDromaSetFromDatabase("GDSC1", db_path)

# Create MultiDromaSet using the correct constructor
multi_set <- MultiDromaSet(
  name = c("gCSI", "CCLE", "GDSC1"),
  DromaSets = list(gCSI = gCSI, CCLE = CCLE, GDSC1 = GDSC1)
)

# Generate all statistical plots
cat("Generating comprehensive statistical plots...\n")
all_plots <- generateStatisticalPlots(
  project_names = "all",
  plot_types = "all",
  use_gap_plots = TRUE
)

# Display individual plots
if (!is.null(all_plots$counts)) {
  print("Drug and sample count plots:")
  print(all_plots$counts$detailed)
  print(all_plots$counts$summary)
}

if (!is.null(all_plots$overlaps)) {
  print("Overlap plots:")
  print(all_plots$overlaps$samples)
  print(all_plots$overlaps$drugs)
}

if (!is.null(all_plots$molecular)) {
  print("Molecular characteristics plot:")
  print(all_plots$molecular)
}

if (!is.null(all_plots$drug_moa)) {
  print("Drug MOA distribution plot:")
  print(all_plots$drug_moa)
}

if (!is.null(all_plots$tumor_types)) {
  print("Tumor type distribution plot:")
  print(all_plots$tumor_types)
}

# Example 2: Generate specific plot types ----

# Generate only count and overlap plots
count_overlap_plots <- generateStatisticalPlots(
  data_source = multi_set,
  plot_types = c("counts", "overlaps")
)

# Example 3: Using individual data objects (backward compatibility) ----

# If you have individual data matrices (legacy approach)
# This shows backward compatibility with existing code

# Example data matrices (replace with your actual data)
# Load drug response data from DromaSet objects
gCSI_drug <- loadTreatmentResponse(gCSI, return_data = TRUE)
CCLE_drug <- loadTreatmentResponse(CCLE, return_data = TRUE)
GDSC1_drug <- loadTreatmentResponse(GDSC1, return_data = TRUE)

data_list <- list(
  gCSI_drug = gCSI_drug,
  CCLE_drug = CCLE_drug,
  GDSC1_drug = GDSC1_drug
)

# Generate plots from individual data objects
individual_plots <- generateStatisticalPlots(
  data_source = data_list,
  plot_types = c("counts", "overlaps")
)

# Example 4: Create dashboard layout ----

# Combine plots into a dashboard
if (requireNamespace("patchwork", quietly = TRUE)) {
  dashboard <- createStatisticalDashboard(
    plot_list = all_plots,
    layout = "grid"
  )
  print("Statistical dashboard created!")
  # print(dashboard)
}

# Example 5: Advanced usage with custom annotations ----

# Access metadata from MultiDromaSet slots
sample_annotations <- multi_set@sampleMetadata
drug_annotations <- multi_set@treatmentMetadata

# Generate plots with explicit annotations
custom_plots <- generateStatisticalPlots(
  data_source = multi_set,
  plot_types = "all",
  sample_annotations = sample_annotations,
  drug_annotations = drug_annotations
)

# Example 6: Extracting data for custom analysis ----

# Access the underlying data structures
count_data <- all_plots$counts$data
print("Count data structure:")
print(head(count_data))

if (!is.null(all_plots$overlaps)) {
  drug_overlaps <- all_plots$overlaps$drug_overlaps
  sample_overlaps <- all_plots$overlaps$sample_overlaps

  print("Number of unique drugs per project:")
  print(sapply(drug_overlaps, length))

  print("Number of unique samples per project:")
  print(sapply(sample_overlaps, length))
}

# Example 7: Working with MultiDromaSet properties ----

# Explore MultiDromaSet structure
print("MultiDromaSet project names:")
print(multi_set@name)

print("MultiDromaSet dataset types:")
print(multi_set@datasetType)

print("Sample metadata columns:")
print(colnames(multi_set@sampleMetadata))

print("Treatment metadata columns:")
print(colnames(multi_set@treatmentMetadata))

# Example 8: Error handling and troubleshooting ----

# The functions include built-in error handling
# If a plot type fails, it will return NULL with a warning

# Check which plots were successfully generated
successful_plots <- names(all_plots)[!sapply(all_plots, is.null)]
cat("Successfully generated plots:", paste(successful_plots, collapse = ", "), "\n")

failed_plots <- c("counts", "overlaps", "molecular", "drug_moa", "tumor_types")[
  !c("counts", "overlaps", "molecular", "drug_moa", "tumor_types") %in% successful_plots
]
if (length(failed_plots) > 0) {
  cat("Failed to generate plots:", paste(failed_plots, collapse = ", "), "\n")
}

# Example 9: Saving plots ----

# Save individual plots
if (!is.null(all_plots$counts)) {
  ggsave("drug_sample_counts.png", all_plots$counts$detailed,
         width = 12, height = 8, dpi = 300)
  ggsave("total_counts_summary.png", all_plots$counts$summary,
         width = 6, height = 6, dpi = 300)
}

if (!is.null(all_plots$molecular)) {
  ggsave("molecular_characteristics.png", all_plots$molecular,
         width = 10, height = 8, dpi = 300)
}

if (!is.null(all_plots$tumor_types)) {
  ggsave("tumor_type_distribution.png", all_plots$tumor_types,
         width = 12, height = 8, dpi = 300)
}

# Example 10: Integration with existing DROMA.R workflows ----

# These statistical plots can be used alongside existing DROMA.R analysis
# For example, after performing batch feature analysis:

batch_results <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mRNA"
)

volcano_plot <- plotMetaVolcano(batch_results, es_t = 0.3, P_t = 0.05)

# Combined analysis report
combined_report <- list(
  statistical_overview = all_plots,
  batch_analysis = batch_results,
  volcano_plot = volcano_plot
)

cat("Example completed! Check your working directory for saved plots.\n")

# Tips for usage:
# 1. Ensure all required packages are installed:
#    - UpSetR (for overlap plots)
#    - treemapify (for MOA plots)
#    - ggrepel (for tumor type plots)
#    - patchwork (for dashboard layout)
#
# 2. The functions are designed to handle missing data gracefully
#
# 3. Large datasets may take some time to process - be patient!
#
# 4. If certain plot types fail, check your data structure and annotations
#
# 5. For troubleshooting, run plot types individually to identify issues
#
# 6. MultiDromaSet objects automatically merge metadata from constituent DromaSet objects
