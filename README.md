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

### 4. Analyze Drug-Omics Associations

```r
# Single project analysis
result <- analyzeDrugOmicPair(
  gCSI,
  feature_type = "mRNA",
  select_features = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Multi-project meta-analysis
multi_result <- analyzeDrugOmicPair(
  multi_set,
  feature_type = "mRNA",
  select_features = "ABCB1",
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

### 🔧 Data Loading and Processing
- **`loadFeatureData()`**: Unified feature loader used by the analysis workflows
- **`processDrugData()`**: Process one drug across a `DromaSet` or `MultiDromaSet`
- **`getDrugSensitivityData()`**: Process and optionally annotate one drug
- **`getAllDrugSensitivityData()`**: Collect all drugs with optional summarization
- **`annotateDrugData()`**: Attach sample annotations to processed drug data
- **`formatDrugTable()`**: Present processed drug data with raw and z-score values

### 📊 Pairwise and Batch Analysis
- **`analyzeDrugOmicPair()`**: Single feature vs single drug analysis for `DromaSet` or `MultiDromaSet`
- **`batchFindSignificantFeatures()`**: Batch screen many features against one feature/drug
- **`getSignificantFeatures()`**: Filter batch results by p-value, q-value, effect size, and dataset count
- **`getIntersectSignificantFeatures()`**: Intersect significant feature sets across analyses

### 🏥 Clinical Trial Database (CTRDB) Module
- **`analyzeClinicalDrugResponse()`**: Clinical response vs omics analysis across CTRDB cohorts
- **`batchFindClinicalSigResponse()`**: Batch screen many omics features against one clinical drug
- **`getInVitroCandidateFeatures()`**: Candidate discovery pipeline combining in vitro evidence
- **`extractPathwayGeneSupport()`**: Summarize pathway support for candidate genes
- **`prioritizeIntegratedCandidates()`**: Combine in vitro, PDO/PDX, and clinical evidence

### 🧬 GSVA and Pathway Analysis
- **`loadGeneSetsFromGMT()`**, **`getAvailableGeneSets()`**, **`getGeneSetByName()`**
- **`calculateGSVAScores()`**: Compute GSVA/ssGSEA/zscore/plage pathway scores
- **`analyzeGSVAAssociation()`**: Pathway vs drug/feature association analysis
- **`batchAnalyzePathwaysWithFeature()`**: Screen pathways against one drug or feature
- **`batchFindClinicalPathwaySigResponse()`**: Clinical response screening at pathway level
- **`batchAnalyzePathwayAcrossFeatures()`**: Pathway-by-feature batch analysis

### 🎨 Visualization and Prioritization
- **`createForestPlot()`**, **`plotMetaVolcano()`**, **`createPlotWithCommonAxes()`**
- **`plotClinicalDrugResponseRoc()`**: ROC/AUC plots for CTRDB response classification
- **`plotDrugSensitivityRank()`**: Ranked drug sensitivity view with raw or z-score values
- **`plotGSVAHeatmap()`**, **`prepareGSVAHeatmapMatrix()`**
- **`plotDrugPathwayEffectHeatmap()`**, **`prepareDrugPathwayEffectSizeMatrix()`**
- **`buildPriorityTable()`**, **`preparePriorityVisualizationData()`**
- **`plotPriorityTopBar()`**, **`plotPriorityEvidenceStacked()`**, **`plotPriorityEvidenceHeatmap()`**
- **`plotClinicallySupportedCandidateSelection()`**, **`plotPriorityClinicalPreclinicalBubble()`**
- **`plotContinuousComparison()`**, **`plotContinuousGroups()`**, **`plotCategoryComparison()`**

### 🛠 Enrichment Utilities
- **`runGO()`**: GO enrichment analysis with optional plotting
- **`runKEGG()`**: KEGG enrichment analysis with optional plotting

> **Note**: `getPatientExpressionData()` has been moved to `DROMA.Set` package (CTRDB_SQLManager.R) for better separation of database operations and analysis functions.

### 🛠 Utility Functions
- **`formatTime()`**: Format processing time outputs
- **`estimateTimeRemaining()`**: Estimate batch processing time
- **`loadFeatureData()`**: Load feature data for batch analysis

## Examples

### Example 1: Loading Feature Data for Analysis

```r
# Load required packages
library(DROMA.Set)
library(DROMA.R)

# Create DromaSet
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")

# Load one continuous feature (z-score by default)
mrna_feature <- loadFeatureData(
  gCSI,
  feature_type = "mRNA",
  select_features = "ABCB1",
  is_continuous = TRUE
)

# Load the same feature without z-score normalization
mrna_raw <- loadFeatureData(
  gCSI,
  feature_type = "mRNA",
  select_features = "ABCB1",
  is_continuous = TRUE,
  zscore = FALSE
)

# Load one drug profile
drug_feature <- loadFeatureData(
  gCSI,
  feature_type = "drug",
  select_features = "Paclitaxel",
  is_continuous = TRUE
)
```

### Example 2: Multi-Project Drug Processing

```r
# Create MultiDromaSet
multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))

# Collect one drug across projects with annotations
multi_drug_data <- getDrugSensitivityData(
  multi_set,
  select_drugs = "Paclitaxel",
  overlap_only = FALSE
)

# Summarize all drugs across projects
all_drugs <- getAllDrugSensitivityData(
  multi_set,
  summarize = "median",
  summarize_by = "drug"
)
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
  feature_type = "mRNA",
  select_features = "ABCB1",
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
  feature_type = "mutation_gene",
  select_features = "TP53",
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

### Example 6: Candidate Prioritization and Visualization

```r
# Visualize an already prepared priority table
top_bar_plot <- plotPriorityTopBar(priority_df)
print(top_bar_plot)
```

### Example 7: Drug Sensitivity Data Processing

```r
# Process drug data with annotations
drug_data <- getDrugSensitivityData(
  dromaset_object = gCSI,
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all",
  db_path = "path/to/droma.sqlite"  # For loading annotations
)

# View processed data
head(drug_data)

# Create a formatted data table
formatDrugTable(drug_data)
```

### Example 8: Clinical Trial Database (CTRDB) Analysis

```r
# Connect to CTRDB database
con <- connectCTRDatabase("path/to/ctrdb.sqlite")

# Analyze clinical drug response across CTRDB cohorts
result <- analyzeClinicalDrugResponse(
  select_omics = "EGFR",
  select_drugs = "Erlotinib",
  data_type = "all",
  tumor_type = "all",
  connection = con,
  meta_enabled = TRUE,
  merged_enabled = TRUE,
  normalization = "combat",  # choose from "none", "zscore", "combat"
  roc_plot = TRUE
)

# View individual patient plots
print(result$plot)

# View merged comparison and pooled ROC
print(result$merged_plot)
print(result$roc$merged$roc_plot)

# View meta-analysis object when multiple cohorts are available
print(result$meta)
```

### Example 9: Clinically Supported Candidate Selection Plot

```r
# Build priority summary tables first, then visualize retained vs filtered candidates
candidate_selection_plot <- plotClinicallySupportedCandidateSelection(
  candidate_selection_df,
  n_retained = 10,
  n_filtered = 4,
  add = c("TP53", "EGFR"),
  subtitle = NULL
)
print(candidate_selection_plot)
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
