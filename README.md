# DROMA.R: Drug Omics Association Analysis Extension for DROMA.Set

[![Website](https://img.shields.io/website?url=https%3A//droma01.github.io/DROMA/)](https://droma01.github.io/DROMA/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: mpl-2-0](https://img.shields.io/badge/MPL-2.0-yellow.svg)](https://opensource.org/licenses/mpl-2-0)

## Overview

**DROMA.R** is an R package that provides advanced analysis functions for drug-omics associations using DromaSet and MultiDromaSet objects from the **DROMA.Set** package. It supports meta-analysis of drug-omics associations across multiple datasets, comprehensive visualization tools, and batch processing of features. This package extends DROMA.Set with statistical analysis capabilities for biomarker discovery in precision medicine. **All data loading functions now apply z-score normalization by default** for improved analysis consistency.

**CTRDB Support**: Includes specialized functions for analyzing Clinical Trial Database (CTRDB) data, enabling patient-level drug response analysis and cross-drug signature stratification.

### Core Design Principles

1. **Modular Architecture**: Functions organized into data loading, pairing, meta-analysis, and visualization modules
2. **Z-score Normalization**: All continuous data normalized by default for cross-dataset comparability
3. **Statistical Rigor**: Appropriate methods for continuous (Spearman) and discrete (Wilcoxon + Cliff's Delta) data
4. **Flexible Analysis**: Supports both single dataset and multi-dataset meta-analysis workflows

It is a part of [DROMA project](https://github.com/mugpeng/DROMA). Visit the [official DROMA website](https://droma01.github.io/DROMA/) for comprehensive documentation and interactive examples.

## Features

- **🔗 DROMA.Set Integration**: Works seamlessly with DromaSet and MultiDromaSet objects
- **📊 Meta-analysis**: Advanced statistical analysis across multiple datasets
- **🏥 Clinical Trial Database (CTRDB)**: Specialized functions for clinical drug response analysis
- **🎨 Comprehensive Visualization**: Forest plots, volcano plots, comparison plots with consistent theming
- **⚡ Batch Processing**: Efficient analysis of multiple features simultaneously
- **🧮 Multiple Statistical Methods**: Spearman correlation, Wilcoxon tests, Cliff's Delta effect sizes
- **🚀 Performance Optimization**: Parallel processing support for large datasets
- **🔄 Z-score Normalization**: All data loading functions apply z-score normalization by default

## Installation

### Prerequisites

First, install the DROMA.Set package:

```r
# Install DROMA.Set (replace with actual installation method)
# devtools::install_github("mugpeng/DROMA_Set")
```

### Install DROMA.R

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install DROMA.R from GitHub
devtools::install_github("mugpeng/DROMA_R")
```

## Quick Start

### 1. Load Required Packages

```r
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions
```

### 2. Create DromaSet Objects

```r
# Connect to DROMA database
connectDROMADatabase("path/to/your/droma.sqlite")

# Create a single DromaSet for one project
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Create a MultiDromaSet for multiple projects
multi_set <- createMultiDromaSetFromDatabase(
    project_names = c("gCSI", "CCLE"),
    db_path = "path/to/droma.sqlite"
)
```

### 3. Load Data with Z-score Normalization

```r
# Load molecular profiles (automatically z-score normalized)
mrna_data <- loadMolecularProfilesNormalized(
  gCSI,
  molecular_type = "mRNA",
  features = "ABCB1"
)

# Load drug response data (automatically z-score normalized)
drug_data <- loadTreatmentResponseNormalized(
  gCSI,
  drugs = "Paclitaxel"
)
```

### 4. Analyze Drug-Omics Associations

```r
# Single project analysis
result <- analyzeDrugOmicPair(
  gCSI,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Multi-project meta-analysis
multi_result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  overlap_only = FALSE
)

# Display results
print(result$meta)
plot(result$plot)
```

### 5. Batch Feature Analysis

```r
# Find genes associated with drug response
batch_results <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mRNA",
  overlap_only = FALSE
)

# Create volcano plot
volcano_plot <- plotMetaVolcano(batch_results, es_t = 0.3, P_t = 0.05)
print(volcano_plot)
```

## Core Functions by Module

### 🔧 Data Loading Module (Z-score Normalized by Default)
- **`loadMolecularProfilesNormalized()`**: Load molecular profiles with automatic z-score normalization
- **`loadTreatmentResponseNormalized()`**: Load drug response data with automatic z-score normalization
- **`loadMultiProjectMolecularProfilesNormalized()`**: Load multi-project molecular profiles
- **`loadMultiProjectTreatmentResponseNormalized()`**: Load multi-project drug response data
- **`applyZscoreNormalization()`**: Apply z-score normalization to existing data
- **`isZscoreNormalized()`**: Check if data has been z-score normalized

### 🔗 Data Pairing Module
- **`pairDrugOmic()`**: Pair continuous drug and omics data
- **`pairDiscreteDrugOmic()`**: Pair discrete omics with drug data
- **`pairContinuousFeatures()`**: General function for pairing continuous features
- **`pairDiscreteFeatures()`**: General function for pairing discrete and continuous features

### 📊 Meta-Analysis Module
- **`metaCalcConCon()`**: Meta-analysis for continuous vs continuous data (Spearman correlation)
- **`metaCalcConDis()`**: Meta-analysis for continuous vs discrete data (Wilcoxon + Cliff's Delta)
- **`analyzeContinuousDrugOmic()`**: Wrapper for continuous drug-omics analysis
- **`analyzeDiscreteDrugOmic()`**: Wrapper for discrete drug-omics analysis

### 🎨 Visualization Module
- **`createForestPlot()`**: Forest plots for meta-analysis results
- **`plotMetaVolcano()`**: Volcano plots for batch analysis results
- **`plotContinuousDrugOmic()`**: Scatter plots with correlation statistics
- **`plotDiscreteDrugOmic()`**: Box plots for group comparisons

### 🚀 High-Level Analysis Functions
- **`analyzeDrugOmicPair()`**: Complete workflow for single drug-omics pair analysis
- **`batchFindSignificantFeatures()`**: Batch screening of multiple features
- **`processDrugData()`**: Process drug sensitivity data with normalization
- **`getDrugSensitivityData()`**: Combined processing and annotation
- **`analyzeStratifiedDrugOmic()`**: Stratified analysis by another drug's response
- **`createStatisticalDashboard()`**: Interactive dashboard for statistical results
- **`generateStatisticalPlots()`**: Generate comprehensive statistical overview plots

### 🏥 Clinical Trial Database (CTRDB) Module
- **`analyzeClinicalDrugResponse()`**: Analyze clinical drug response with omics data
- **`analyzeStratifiedCTRDB()`**: Stratified analysis across different drugs
- **`getPatientExpressionData()`**: Retrieve patient expression data from CTRDB
- **`analyzeClinicalMeta()`**: Meta-analysis across clinical datasets
- **`getClinicalSummary()`**: Summary statistics for clinical analysis
- **`getStratifiedCTRDBSummary()`**: Summary for stratified CTRDB analysis

### 🛠 Utility Functions
- **`bright_palette_26`**: Pre-defined color palette for visualizations
- **`formatTime()`**: Format processing time outputs
- **`estimateTimeRemaining()`**: Estimate batch processing time
- **`formatDrugTable()`**: Format drug sensitivity data for display
- **`annotateDrugData()`**: Add clinical annotations to drug data
- **`loadFeatureData()`**: Load feature data for batch analysis
- **`filterFeatureData()`**: Filter features by minimum sample size

## Examples

### Example 1: Loading Data with Z-score Normalization

```r
# Load required packages
library(DROMA.Set)
library(DROMA.R)

# Create DromaSet
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Load mRNA data with z-score normalization (default)
mrna_normalized <- loadMolecularProfilesNormalized(
  gCSI,
  molecular_type = "mRNA",
  features = "ABCB1"
)

# Load without normalization
mrna_raw <- loadMolecularProfilesNormalized(
  gCSI,
  molecular_type = "mRNA",
  features = "ABCB1",
  zscore = FALSE
)

# Check if data is normalized
cat("Normalized:", isZscoreNormalized(mrna_normalized))
cat("Raw:", isZscoreNormalized(mrna_raw))

# Load drug data with normalization
drug_normalized <- loadTreatmentResponseNormalized(
  gCSI,
  drugs = "Paclitaxel"
)
```

### Example 2: Multi-Project Normalized Loading

```r
# Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))

# Load normalized data across projects
multi_mrna <- loadMultiProjectMolecularProfilesNormalized(
  multi_set,
  molecular_type = "mRNA",
  features = "ABCB1",
  overlap_only = FALSE
)

# Load normalized drug data across projects
multi_drugs <- loadMultiProjectTreatmentResponseNormalized(
  multi_set,
  drugs = "Paclitaxel",
  overlap_only = FALSE
)

# Check normalization status for each project
for (project in names(multi_mrna)) {
  cat(project, "normalized:", isZscoreNormalized(multi_mrna[[project]]), "\n")
}
```

### Example 3: Single Project Analysis

```r
# Load required packages
library(DROMA.Set)
library(DROMA.R)

# Create DromaSet
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Analyze Paclitaxel vs ABCB1 expression
result <- analyzeDrugOmicPair(
  gCSI,
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel"
)

# View results
print(result$meta)
plot(result$plot)
```

### Example 4: Multi-Project Meta-Analysis

```r
# Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))

# Analyze across projects with overlapping samples
result <- analyzeDrugOmicPair(
  multi_set,
  select_omics_type = "mutation_gene",
  select_omics = "TP53",
  select_drugs = "Cisplatin",
  overlap_only = FALSE
)

# Create forest plot
createForestPlot(result$meta)
```

### Example 5: Batch Analysis

```r
# Find all genes associated with Paclitaxel response
batch_results <- batchFindSignificantFeatures(
  multi_set,
  feature1_type = "drug",
  feature1_name = "Paclitaxel",
  feature2_type = "mRNA",
  cores = 4  # Use parallel processing
)

# Sort by significance
sorted_results <- batch_results[order(batch_results$p_value), ]
print(head(sorted_results, 10))

# Create volcano plot
volcano_plot <- plotMetaVolcano(batch_results)
print(volcano_plot)
```

### Example 6: Stratified Analysis

```r
# Analyze drug response stratified by another drug's sensitivity
# This helps identify context-dependent biomarkers

stratified_result <- analyzeStratifiedDrugOmic(
  dromaset_object = multi_set,
  stratification_drug = "Cisplatin",     # Drug for stratification
  select_omics_type = "mRNA",            # Omics data type
  select_omics = "ERCC1",                # Target gene
  select_drugs = "Bortezomib",           # Drug to analyze
  stratify_by = "response_median",       # Stratification method
  tumor_type = "all"
)

# View stratified results
print(stratified_result$statistics)
print(stratified_result$comparison)
```

### Example 7: Drug Sensitivity Data Processing

```r
# Process drug data with annotations
drug_data <- getDrugSensitivityData(
  dromaset_object = gCSI,
  drug_name = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  db_path = "path/to/droma.sqlite"  # For loading annotations
)

# View processed data
head(drug_data)

# Create sensitivity rank plot
formatDrugTable(drug_data, drug_name = "Paclitaxel")
```

### Example 8: Clinical Trial Database (CTRDB) Analysis

```r
# Connect to CTRDB database
con <- connectCTRDatabase("path/to/ctrdb.sqlite")

# Analyze clinical drug response
result <- analyzeClinicalDrugResponse(
  select_omics = "EGFR",
  select_drugs = "Erlotinib",
  data_type = "all",
  tumor_type = "all",
  connection = con,
  meta_enabled = TRUE
)

# View individual patient plots
print(result$plot)

# View meta-analysis forest plot
print(result$forest_plot)

# Get summary statistics
summary <- getClinicalSummary(result)
print(summary)
```

### Example 9: Stratified CTRDB Analysis

```r
# Stratified analysis: Drug B signature applied to Drug A
result <- analyzeStratifiedCTRDB(
  drug_b_name = "Cisplatin",      # Signature generation
  drug_a_name = "Paclitaxel",     # Signature application
  select_omics = "EGFR",          # Feature to analyze
  connection = con,
  top_n_genes = 100,
  data_type = "all",
  tumor_type = "all"
)

# View signature genes
print(result$signature_genes)

# View forest plot with meta-analyzed correlations
print(result$correlation_results$forest_plot)

# View combined scatter plots
print(result$correlation_results$combined_scatter_plot)

# Get comprehensive summary
summary <- getStratifiedCTRDBSummary(result)
print(summary)
```

## Data Types Supported

### Molecular Profiles
- **mRNA**: Gene expression data
- **cnv**: Copy number variation data
- **mutation_gene**: Gene-level mutation data
- **mutation_site**: Site-specific mutation data
- **fusion**: Gene fusion data
- **meth**: DNA methylation data
- **proteinrppa**: Reverse-phase protein array data
- **proteinms**: Mass spectrometry proteomics data

### Treatment Response
- **drug**: Drug sensitivity/response data

## Statistical Methods

### Continuous vs Continuous Analysis
When analyzing relationships between two continuous variables (e.g., gene expression vs drug response):
- **Spearman Correlation**: Non-parametric rank correlation coefficient
- **Fisher's Z-transformation**: Converts correlation coefficients to normally distributed values for meta-analysis
- **Random-Effects Model**: Combines results across multiple studies accounting for heterogeneity

### Continuous vs Discrete Analysis
When comparing continuous values between discrete groups (e.g., mutated vs wild-type):
- **Wilcoxon Rank-Sum Test**: Non-parametric test for group differences
- **Cliff's Delta**: Effect size measurement for group differences (range: -1 to 1)
  - |d| < 0.147: Negligible effect
  - 0.147 ≤ |d| < 0.33: Small effect
  - 0.33 ≤ |d| < 0.474: Medium effect
  - |d| ≥ 0.474: Large effect

### Meta-Analysis Approach
1. **Within-study analysis**: Calculate effect sizes and p-values for each dataset
2. **Quality control**: Filter results based on sample size and data quality
3. **Cross-study combination**: Use random-effects meta-analysis to combine results
4. **Heterogeneity assessment**: Evaluate variability between studies

## Performance Tips

1. For large-scale batch analyses, use `cores > 1` to enable parallel processing
2. Consider setting `overlap_only = TRUE` with MultiDromaSet when you need comparable sample analyses across all datasets
3. Use the `db_path` parameter with `annotateDrugData()` or `getDrugSensitivityData()` to automatically load sample annotations from the database when needed
4. For visualization of large datasets, consider filtering to specific tumor types or data types

## Contributing

Contributions to DROMA.R are welcome! Please feel free to submit issues or pull requests on GitHub.

## License

This project is licensed under the MPL-2 License - see the LICENSE file for details.

## Citation

If you use DROMA.R in your research, please cite:

```
Li, S., Peng, Y., Chen, M. et al. Facilitating integrative and personalized oncology omics analysis with UCSCXenaShiny. Commun Biol 7, 1200 (2024). https://doi.org/10.1038/s42003-024-06891-2

```

## Contact

For questions and feedback, please contact Peng Yu Zhong at yc47680@um.edu.mo.

---

**DROMA.R** - Advanced drug-omics association analysis powered by DROMA.Set 🧬💊📊
