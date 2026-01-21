# GSVA Analysis Example - Batch Analysis with Hallmark Gene Sets
#
# This example demonstrates:
# 1. Loading predefined Hallmark gene sets
# 2. Calculating GSVA scores for all gene sets
# 3. Batch analyzing multiple pathways against a single feature (drug/omics)
# 4. Batch analyzing a single pathway across multiple features
# 5. Visualizing results with a dotplot

library(DROMA.R)
library(DROMA.Set)

# ============================================================================
# Step 1: Load DromaSet or MultiDromaSet
# ============================================================================
# Example with a single dataset
# gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Or with multiple datasets
# multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")

# For this example, we'll load datasets from the database
db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
# db_path <- "/home/data/denglab/bigData/DROMA/droma.sqlite"

con <- connectDROMADatabase(db_path)

# Option 1: Load only PDO datasets
project_anno <- DROMA.Set::listDROMAProjects()
pdo_names <- project_anno[project_anno$dataset_type %in% "PDO",]$project_name

pdo_sets <- createMultiDromaSetFromAllProjects(db_path = db_path,
                                               include_projects = pdo_names,
                                               con = con)

# Option 2: Load all datasets
all_sets <- createMultiDromaSetFromAllProjects(db_path = db_path,
                                               con = con)

# ============================================================================
# Step 2: Get available predefined gene sets
# ============================================================================
# Get all available Hallmark gene sets from the default GMT file
# Option 1: Get only gene set names (default behavior)
available_gene_sets <- getAvailableGeneSets(gmt_file = "../Data/h.all.v2025.1.Hs.symbols.gmt")
cat("Found", length(available_gene_sets), "predefined gene sets\n")
cat("First 10 gene sets:", paste(head(available_gene_sets, 10), collapse = ", "), "\n")

# Option 2: Get full gene sets with genes (return_full_list = TRUE)
full_gene_sets <- getAvailableGeneSets(gmt_file = "../Data/h.all.v2025.1.Hs.symbols.gmt",
                                       return_full_list = TRUE)
cat("\nFull gene sets loaded. Each gene set contains its genes.\n")
cat("Example: First gene set '", names(full_gene_sets)[1], "' contains ",
    length(full_gene_sets[[1]]), " genes\n")
cat("First 10 genes in '", names(full_gene_sets)[1], "': ",
    paste(head(full_gene_sets[[1]], 10), collapse = ", "), "\n")

# View gene counts for all gene sets
gene_counts <- sapply(full_gene_sets, length)
cat("\nGene set sizes (min/max/mean): ", min(gene_counts), "/", max(gene_counts), "/",
    round(mean(gene_counts)), " genes\n")

# ============================================================================
# Step 3: Calculate GSVA scores for all gene sets
# ============================================================================
# Calculate GSVA scores using the GSVA method
# This will compute enrichment scores for all gene sets
cat("\nCalculating GSVA scores for all gene sets...\n")

# Option A: Use full_gene_sets directly (recommended - already contains gene lists)
gsva_scores <- calculateGSVAScores(
  dromaset_object = pdo_sets,
  gene_sets = full_gene_sets,  # Use first 2 gene sets with their gene lists
  method = "gsva",  # Options: "gsva", "ssgsea", "zscore", "plage"
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  min_size = 5,
  max_size = Inf
)
# gsva_scores2 <- gsva_scores[!(names(gsva_scores) %in% "HKUPDO")]

# Option B: Use gene set names with gmt_file parameter
# gsva_scores <- calculateGSVAScores(
#   dromaset_object = all_sets,
#   gene_sets = available_gene_sets[1:2],  # Gene set names
#   method = "gsva",
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   min_size = 5,
#   max_size = Inf,
#   gmt_file = "../Data/h.all.v2025.1.Hs.symbols.gmt"  # Must specify gmt_file for name lookup
# )

cat("GSVA scores calculated for", nrow(gsva_scores[[1]]), "gene sets\n")

# ============================================================================
# Step 4: Batch analyze multiple pathways against a single feature
# ============================================================================
# Analyze associations between GSVA scores and a single feature (drug or omics)
# This function tests multiple pathways against the same feature

# Example 4A: Analyze associations with a drug
cat("\nBatch analyzing associations with Bortezomib (drug)...\n")
results_drug <- batchAnalyzePathwaysWithFeature(
  dromaset_object = pdo_sets,  # Use same dataset as GSVA calculation
  gsva_scores = gsva_scores,
  feature_name = "Bortezomib",
  feature_type = "drug",  # Can be "drug", "mRNA", "mutation_gene", etc.
  gene_set_names = NULL,  # NULL means analyze all gene sets in gsva_scores
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  zscore = TRUE,
  show_progress = TRUE
)

cat("\nAnalysis complete!\n")
cat("Total gene sets analyzed:", nrow(results_drug), "\n")
cat("Significant associations (p < 0.05):", sum(results_drug$p_value < 0.05), "\n")
cat("Significant associations (q < 0.05):", sum(results_drug$q_value < 0.05), "\n")

# Example 4B: Analyze associations with a gene expression (mRNA)
# cat("\nBatch analyzing associations with TP53 expression...\n")
# results_mRNA <- batchAnalyzePathwaysWithFeature(
#   dromaset_object = pdo_sets,
#   gsva_scores = gsva_scores,
#   feature_name = "TP53",
#   feature_type = "mRNA",
#   gene_set_names = NULL,
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   zscore = TRUE,
#   show_progress = TRUE
# )

# Example 4C: Analyze associations with mutation status
# cat("\nBatch analyzing associations with KRAS mutation...\n")
# results_mutation <- batchAnalyzePathwaysWithFeature(
#   dromaset_object = pdo_sets,
#   gsva_scores = gsva_scores,
#   feature_name = "KRAS",
#   feature_type = "mutation_gene",
#   gene_set_names = NULL,
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   zscore = TRUE,
#   show_progress = TRUE
# )

# View top results
cat("\nTop 10 gene sets by p-value:\n")
print(head(results, 10))

# ============================================================================
# Step 5: Visualize results with volcano plot
# ============================================================================
# Create a volcano plot showing gene set associations
# Point size represents number of datasets, y-axis is -log10(q-value)
# Remove "HALLMARK_" prefix from pathway names for cleaner labels
cat("\nCreating volcano plot visualization...\n")

# Remove "HALLMARK_" prefix from pathway names (optional, for cleaner labels)
results_drug$name <- gsub("^HALLMARK_", "", results_drug$name)

# Optionally filter to top N pathways by p-value before plotting
# top_n <- 20
# if (nrow(results_drug) > top_n) {
#   results_drug <- results_drug[order(results_drug$p_value), ][1:top_n, ]
# }

# Create volcano plot using plotMetaVolcano
# Note: plotMetaVolcano uses q_value by default (set use_p_value=TRUE to use p_value)
volcano_plot <- plotMetaVolcano(
  meta_df = results_drug,
  es_t = 0,  # Effect size threshold (0 = no threshold, or use 0.4 for stricter filtering)
  P_t = 0.05,  # P-value threshold (or q-value threshold if use_p_value=FALSE)
  n_datasets_t = NULL,  # Optional: minimum number of datasets threshold
  label = TRUE,  # Show labels for significant pathways
  top_label_each = 5,  # Number of top pathways in each direction to label
  custom_labels = NULL,  # Optional: specific pathway names to label
  label_size = 4,
  point_size = 2.5,  # Base point size (scaled by n_datasets)
  point_alpha = 0.6,
  title = "GSVA Gene Set Associations with Bortezomib",
  custom_colors = NULL,  # NULL = use defaults (Down=blue, NS=grey, Up=red)
  use_p_value = TRUE  # FALSE = use q_value (recommended), TRUE = use p_value
)

# Display the plot
print(volcano_plot)

# ============================================================================
# Optional: Save results
# ============================================================================
# Save results to CSV
# write.csv(results, file = "gsva_bortezomib_results.csv", row.names = FALSE)

# Save plot
# ggsave("gsva_bortezomib_dotplot.png", plot = dotplot, width = 12, height = 8, dpi = 300)

# ============================================================================
# Step 6: Batch analyze a single pathway across multiple features
# ============================================================================
# This function takes a single gene list (pathway) and tests it against
# multiple features (drugs or omics). Similar to batchFindSignificantFeatures
# but uses pathway-level GSVA scores instead of single gene expression.

# Example 6A: Define a custom pathway and test against all drugs
cat("\n\n=== Testing a custom pathway across multiple drugs ===\n")
my_genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51", "PALB2",
              "FANCA", "FANCD2", "FANCG", "RAD51C", "RAD51D", "BRIP1")

cat("Testing pathway 'DNA_Repair' with", length(my_genes), "genes\n")
cat("Genes:", paste(my_genes, collapse = ", "), "\n")

# Test against all drugs (or use feature_names to specify specific drugs)
# Note: You can use pdo_sets or all_sets depending on your analysis needs
pathway_results <- batchAnalyzePathwayAcrossFeatures(
  dromaset_object = pdo_sets,  # Use PDO datasets, or change to all_sets for all datasets
  gene_list = my_genes,
  pathway_name = "DNA_Repair",
  feature_type = "drug",
  feature_names = NULL,  # NULL = test all drugs, or specify: c("Olaparib", "Cisplatin")
  gsva_method = "gsva",  # Options: "gsva", "ssgsea", "zscore", "plage"
  data_type = "all",
  tumor_type = "all",
  overlap_only = FALSE,
  zscore = TRUE,
  min_size = 3,  # Minimum genes that must be present in expression data
  show_progress = TRUE,
  test_top_n = NULL,  # NULL = test all features, or specify number to limit
  cores = 1  # Number of CPU cores for parallel processing (1 = sequential, >1 = parallel)
)

cat("\nPathway analysis complete!\n")
cat("Total features tested:", nrow(pathway_results), "\n")
cat("Significant associations (q < 0.05):", sum(pathway_results$q_value < 0.05), "\n")

# View top results
cat("\nTop 10 features associated with DNA_Repair pathway:\n")
print(head(pathway_results, 10))

# Example 6B: Test pathway against specific drugs only (with parallel processing)
# pathway_results_specific <- batchAnalyzePathwayAcrossFeatures(
#   dromaset_object = pdo_sets,  # or all_sets
#   gene_list = my_genes,
#   pathway_name = "DNA_Repair",
#   feature_type = "drug",
#   feature_names = c("Olaparib", "Cisplatin", "Carboplatin", "Temozolomide"),
#   gsva_method = "ssgsea",
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   zscore = TRUE,
#   min_size = 5,
#   show_progress = TRUE,
#   cores = 8  # Use 8 cores for faster processing
# )

# Example 6C: Test pathway against gene expressions (mRNA)
# pathway_results_mRNA <- batchAnalyzePathwayAcrossFeatures(
#   dromaset_object = pdo_sets,  # or all_sets
#   gene_list = my_genes,
#   pathway_name = "DNA_Repair",
#   feature_type = "mRNA",
#   feature_names = NULL,  # Test all genes, or specify: c("CDK1", "CDK2", "CCNB1")
#   gsva_method = "gsva",
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   zscore = TRUE,
#   min_size = 5,
#   show_progress = TRUE,
#   test_top_n = 1000  # Limit to top 1000 genes for faster testing
# )

# Example 6D: Test pathway against mutation status (with parallel processing)
# pathway_results_mutation <- batchAnalyzePathwayAcrossFeatures(
#   dromaset_object = pdo_sets,  # or all_sets
#   gene_list = my_genes,
#   pathway_name = "DNA_Repair",
#   feature_type = "mutation_gene",
#   feature_names = NULL,  # Test all mutations, or specify: c("TP53", "BRCA1", "BRCA2")
#   gsva_method = "gsva",
#   data_type = "all",
#   tumor_type = "all",
#   overlap_only = FALSE,
#   zscore = TRUE,
#   min_size = 5,
#   show_progress = TRUE,
#   cores = 4  # Use parallel processing for faster analysis
# )

# ============================================================================
# Optional: Analyze specific gene set in detail
# ============================================================================
# If you want to see detailed analysis for a specific gene set
# top_gene_set <- results$name[1]
# detailed_result <- analyzeGSVAAssociation(
#   dromaset_object = all_sets,
#   gsva_scores = gsva_scores2,
#   select_gene_set = top_gene_set,
#   select_drugs = "Bortezomib",
#   merged_enabled = TRUE,
#   meta_enabled = TRUE
# )
# print(detailed_result$plot)
