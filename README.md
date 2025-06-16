# DROMA.R: Drug Omics Association Analysis Extension for DROMA.Set

[![Version](https://img.shields.io/badge/version-0.9.0-blue.svg)](https://github.com/mugpeng/DROMA_R)

## Overview

**DROMA.R** is an R package that provides advanced analysis functions for drug-omics associations using DromaSet and MultiDromaSet objects from the **DROMA.Set** package. It supports meta-analysis of drug-omics associations across multiple datasets, comprehensive visualization tools, and batch processing of features. This package extends DROMA.Set with statistical analysis capabilities for biomarker discovery in precision medicine. **All data loading functions now apply z-score normalization by default** for improved analysis consistency.

## Features

- **🔗 DROMA.Set Integration**: Works seamlessly with DromaSet and MultiDromaSet objects
- **📊 Meta-analysis**: Advanced statistical analysis across multiple datasets
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

### 3. Analyze Drug-Omics Associations

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

# Multi-project analysis
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

### 4. Batch Feature Analysis

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

## Core Functions

### Analysis Functions
- **`analyzeDrugOmicPair()`**: Analyze association between drug response and omics feature
- **`batchFindSignificantFeatures()`**: Batch analysis of multiple features
- **`analyzeContinuousDrugOmic()`**: Meta-analysis for continuous data
- **`analyzeDiscreteDrugOmic()`**: Meta-analysis for discrete data

### Data Loading Functions (Z-score Normalized by Default)
- **`loadMolecularProfilesNormalized()`**: Load molecular profiles with z-score normalization (default: TRUE)
- **`loadTreatmentResponseNormalized()`**: Load treatment response data with z-score normalization (default: TRUE)
- **`loadMultiProjectMolecularProfilesNormalized()`**: Load multi-project molecular profiles with normalization
- **`loadMultiProjectTreatmentResponseNormalized()`**: Load multi-project treatment response with normalization
- **`applyZscoreNormalization()`**: Apply z-score normalization to existing data
- **`isZscoreNormalized()`**: Check if data has been z-score normalized

### Data Pairing Functions
- **`pairDrugOmic()`**: Pair continuous drug and omics data
- **`pairDiscreteDrugOmic()`**: Pair discrete omics with drug data
- **`pairContinuousFeatures()`**: Pair continuous feature data
- **`pairDiscreteFeatures()`**: Pair discrete with continuous features

### Visualization Functions
- **`createForestPlot()`**: Create forest plots for meta-analysis results
- **`plotMetaVolcano()`**: Create volcano plots for batch analysis results
- **`plotContinuousDrugOmic()`**: Scatter plots for continuous associations
- **`plotDiscreteDrugOmic()`**: Box plots for discrete associations

### Utility Functions
- **`bright_palette_26`**: Pre-defined palette of 26 distinct colors
- **`formatTime()`**: Format time durations
- **`estimateTimeRemaining()`**: Estimate remaining processing time

### Data Processing Functions
- **`processDrugData()`**: Process drug sensitivity data using DromaSet objects
- **`annotateDrugData()`**: Add sample annotations to drug sensitivity data (now with database loading support)
- **`getDrugSensitivityData()`**: Wrapper for processing and annotating drug data

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

## Changelog
### Version 0.4.3
Update DESCRIPTION and R functions for drug sensitivity analysis
- Updated DESCRIPTION file to reflect new version (0.4.3) and changed R dependency to require R version 4.0.0.
- Modified `processDrugData` function to standardize column names (`sampleid` to `SampleID`, `study` to `ProjectID`).
- Enhanced `annotateDrugData` function to merge drug data with annotations using updated column names.
- Reintroduced `formatDrugTable` function with improved formatting for drug sensitivity data.
- Updated example scripts to demonstrate new functionalities in drug sensitivity rank plotting.

### Version 0.4.2
Add new functions for statistical analysis and update examples including createStatisticalDashboard and generateStatisticalPlots.
Removed FuncStat.R as it is no longer needed.
Updated example scripts to include new multi-set functionality for drug feature analysis and drug-omics pairing analysis.

---

**DROMA.R** - Advanced drug-omics association analysis powered by DROMA.Set 🧬💊📊
