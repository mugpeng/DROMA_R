# DROMA Package Examples

This directory contains example scripts demonstrating how to use various features of the DROMA package.

## Getting Started

Before running these examples, make sure you have the DROMA package installed. If you installed it from GitHub:

```r
# Install devtools if not already installed
# install.packages("devtools")

# Install DROMA package from GitHub
# devtools::install_github("your-github-username/DROMA")

# Load the package
library(DROMA)
```

## Available Examples

### Data Management

- `example_data_loading.R`: Demonstrates various ways to load data using `setupDROMA()` and `loadDROMA()` functions.
- `example_feature_selection.R`: Shows how to use the `selectFeatures()` function to retrieve specific omics and drug features.

### Drug-Omics Analysis

- `example_drug_omics_pairing.R`: Demonstrates analyzing associations between drug sensitivity and omics features.
- `example_drug_feature_analysis.R`: Shows drug data processing and visualization functions.
- `example_batch_feature_analysis.R`: Demonstrates batch analysis to find features associated with drug response.

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

These examples require access to the DROMA datasets, which should be installed with the package. The main datasets used include:

- `drug.Rda`: Drug sensitivity data
- `mRNA.Rda`: Gene expression data
- `mut.Rda`: Mutation data
- `search_vec.Rda`: Search vectors for features
- `anno.Rda`: Sample annotations

## Example Workflow

A typical workflow with DROMA involves:

1. Setting up the environment with `setupDROMA()`
2. Loading necessary data with `loadDROMA()`
3. Selecting features of interest with `selectFeatures()`
4. Analyzing relationships between drugs and omics data
5. Visualizing results with the provided plotting functions 