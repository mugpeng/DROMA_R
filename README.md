# DROMA: Drug Omics Association Map


[![Version](https://img.shields.io/badge/version-0.3.1-blue.svg)](https://github.com/mugpeng/DROMA_R)

## Overview

DROMA (Drug Omics Association) is an R package that provides tools for analyzing associations between drug responses and various omics data types. It supports meta-analysis of drug-omics associations across multiple datasets, visualization of results, and batch processing of features.

## Features

- Meta-analysis of drug-omics associations across multiple datasets
- Support for various omics data types (mutation, gene expression, methylation, copy number variation, protein expression, fusion data)
- Comprehensive visualization tools (forest plots, comparison plots, volcano plots)
- Batch processing for analyzing multiple features simultaneously
- Statistical analyses with multiple methodologies based on data types
- Performance optimization for large datasets

## Installation

### From GitHub

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install DROMA from GitHub
devtools::install_github("mugpeng/DROMA_R")
```

### From source

```r
# Clone the repository
git clone https://github.com/mugpeng/DROMA_R.git

# Install dependencies
install.packages(c("meta", "metafor", "effsize", "parallel", "snowfall", 
                  "ggplot2", "ggpubr", "dplyr", "DT", "htmltools", 
                  "patchwork", "ggrepel", "gridExtra", "grid"))

# Build and install the package
devtools::install("./DROMA_R")
```

## Data

The package includes sample datasets for testing and demonstration purposes. 
To access the full datasets, please download them from Zenodo:

**Full datasets:** [https://zenodo.org/records/15392760](https://zenodo.org/records/15392760)

After downloading, you can load the data using:

```r
# Example of loading external data files
load("path/to/downloaded/mRNA.Rda")
load("path/to/downloaded/cnv.Rda")
load("path/to/downloaded/meth.Rda")
load("path/to/downloaded/protein.Rda")
```

## Quick Start

```r
library(DROMA)

# Setup DROMA environment
setupDROMA()

# Load required data types
loadDROMA("drug")   # Load drug data
loadDROMA("mRNA")   # Load mRNA data
loadDROMA("mut")    # Load mutation data

# Analyze association between Paclitaxel and ABCB1 gene expression
result <- analyzeDrugOmicPair(
  select_omics_type = "mRNA",
  select_omics = "ABCB1",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# Display results
print(result$data)

# Visualize the results
plot(result$plot)
```

## Examples

### 1. Drug-Omics Association Analysis

Analyze the association between a drug and an omics feature:

```r
# Load necessary data
loadDROMA("drug")
loadDROMA("mut")

# Analyze association between Paclitaxel and TP53 mutation
result <- analyzeDrugOmicPair(
  select_omics_type = "mutation_gene",
  select_omics = "TP53",
  select_drugs = "Paclitaxel",
  data_type = "all",
  tumor_type = "all"
)

# View results
print(result$data)

# Create forest plot
createForestPlot(result$data)
```

### 2. Batch Feature Analysis

Analyze multiple omics features in relation to drug response:

```r
# Load necessary data
loadDROMA("drug")
loadDROMA("mRNA")

# Batch analysis of multiple genes
results <- batchFindSignificantFeatures(
  drug_name = "Paclitaxel",
  omics_type = "mRNA",
  feature_list = c("ABCB1", "ERBB2", "ESR1", "BRCA1"),
  data_type = "all",
  tumor_type = "all"
)

# View results table
print(results$result_table)

# Create volcano plot
plotMetaVolcano(results$result_table)
```



### Full example

check scripts under `examples/`



## Documentation

For detailed documentation and function references, use the built-in R help system:

```r
?DROMA
help(package = "DROMA")
```

## Contributing

Contributions to DROMA are welcome! Please feel free to submit issues or pull requests on GitHub.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use DROMA in your research, please cite:

```
Zhong, P. Y. (2023). DROMA: Drug Omics Association Map. 
https://github.com/mugpeng/DROMA_R
```

## Contact

For questions and feedback, please contact Peng Yu Zhong at yc47680@um.edu.mo. 
