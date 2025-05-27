# DROMA.R Package Examples

This directory contains example scripts demonstrating how to use DROMA.R with DROMA.Set objects for drug-omics association analysis.

## Getting Started

Before running these examples, make sure you have both packages installed:

```r
# Install DROMA.Set first (replace with actual installation method)
# devtools::install_github("mugpeng/DROMA_Set")

# Install DROMA.R
# devtools::install_github("mugpeng/DROMA_R")

# Load the packages
library(DROMA.Set)  # For data management
library(DROMA.R)    # For analysis functions
```

## Available Examples

### Data Management with DROMA.Set

- `example_data_loading.R`: Demonstrates creating DromaSet and MultiDromaSet objects from databases
- `example_feature_selection.R`: Shows how to load and extract specific molecular profiles and treatment response data

### Drug-Omics Analysis with DROMA.R

- `example_drug_omics_pairing.R`: Demonstrates analyzing associations between drug sensitivity and omics features using DromaSet objects
- `example_drug_feature_analysis.R`: Shows drug data processing and visualization functions
- `example_batch_feature_analysis.R`: Demonstrates batch analysis to find features associated with drug response

## Running Examples

You can run the examples directly from R:

```r
# Ensure your working directory is set to the package root
source("examples/example_data_loading.R")
```

Or from the command line:

```bash
Rscript examples/example_data_loading.R
```

## Required Data

These examples require access to a DROMA SQLite database. The database should contain:

- Sample annotations (`sample_anno` table)
- Drug annotations (`drug_anno` table)
- Project-specific data tables (e.g., `gCSI_mRNA`, `CCLE_drug`, etc.)

## Example Workflow

A typical workflow with DROMA.R and DROMA.Set involves:

1. **Creating DromaSet objects** with `createDromaSetFromDatabase()` or `createMultiDromaSetFromDatabase()`
2. **Loading specific data** using `loadMolecularProfiles()` and `loadTreatmentResponse()` methods
3. **Analyzing relationships** between drugs and omics data with `analyzeDrugOmicPair()`
4. **Batch processing** multiple features with `batchFindSignificantFeatures()`
5. **Visualizing results** with the provided plotting functions

## Benefits of New Approach

1. **Better data management**: DromaSet objects encapsulate data and metadata
2. **Cross-project analysis**: MultiDromaSet enables seamless meta-analysis
3. **Efficient data loading**: Load only what you need, when you need it
4. **Sample filtering**: Built-in support for data type and tumor type filtering
5. **Overlap detection**: Automatic handling of overlapping samples across projects 