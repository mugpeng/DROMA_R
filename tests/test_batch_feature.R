# Test script for FuncBatchFeature.R functions
library(meta)
library(metafor)
library(effsize)
library(parallel)
library(snowfall)
library(ggplot2)
library(ggpubr)

# Test 1: Basic batch feature analysis for continuous vs continuous features
# Example: Compare a drug (continuous) with mRNA expression (continuous)
test_batch_con_con <- function() {
  # Parameters
  feature1_type <- "drug"
  feature1_name <- "fenofibrate"
  feature2_type <- "proteinrppa"
  data_type <- "all"
  tumor_type <- "all"
  
  # Run batch analysis with 2 cores
  results <- BatchFindSigFeaturesPlus(
    feature1_type = feature1_type,
    feature1_name = feature1_name,
    feature2_type = feature2_type,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = 2,
    test_top_100 = F
  )
  
  # Display summary of results
  cat("Total features analyzed:", nrow(results), "\n")
  cat("Significant features (p < 0.05):", sum(results$p_value < 0.05), "\n")
  
  # Create and save volcano plot
  p <- plotMetaVolcano(
    results,
    es_t = 0.2,
    P_t = 0.05,
    label = TRUE,
    top_label_each = 5,
    title = paste(feature1_name, "vs", feature2_type)
  )
  
  print(p)
  
  # Return results for further analysis
  return(results)
}

# Test 2: Continuous vs Discrete features
# Example: Compare a drug (continuous) with mutations (discrete)
test_batch_con_dis <- function() {
  # Parameters
  feature1_type <- "drug"
  feature1_name <- "Bortezomib"
  feature2_type <- "mutation_gene"
  data_type <- "all"
  tumor_type <- "all"
  
  # Run batch analysis
  results <- BatchFindSigFeaturesPlus(
    feature1_type = feature1_type,
    feature1_name = feature1_name,
    feature2_type = feature2_type,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = 1
  )
  
  # Display summary of results
  cat("Total features analyzed:", nrow(results), "\n")
  cat("Significant features (p < 0.05):", sum(results$p_value < 0.05), "\n")
  
  # Create volcano plot with p-value adjustment
  p <- plotMetaVolcano(
    results,
    es_t = 0.3,
    P_t = 0.05,
    label = TRUE,
    top_label_each = 8,
    p_adj_method = "BH",
    title = paste(feature1_name, "vs", feature2_type, "(BH adjusted)")
  )
  
  print(p)
  
  return(results)
}

# Test 3: Test with smaller datasets by filtering on tumor type
test_batch_tumor_specific <- function() {
  # Parameters
  feature1_type <- "drug"
  feature1_name <- "Paclitaxel"
  feature2_type <- "proteinms"
  data_type <- "all"
  tumor_type <- "breast cancer"  # Filter to a specific tumor type
  
  # Run batch analysis
  results <- BatchFindSigFeaturesPlus(
    feature1_type = feature1_type,
    feature1_name = feature1_name,
    feature2_type = feature2_type,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = 1
  )
  
  # Display summary of results
  if (!is.null(results)) {
    cat("Total features analyzed:", nrow(results), "\n")
    cat("Significant features (p < 0.05):", sum(results$p_value < 0.05), "\n")
    
    # Create volcano plot
    p <- plotMetaVolcano(
      results,
      es_t = 0.3,
      P_t = 0.05,
      label = TRUE,
      top_label_each = 5,
      title = paste(feature1_name, "vs", feature2_type, "in", tumor_type)
    )
    
    print(p)
  } else {
    cat("No valid results found for the selected feature combination.\n")
  }
  
  return(results)
}

# Test 4: Testing the progress callback functionality
test_with_progress <- function() {
  # Parameters
  feature1_type <- "drug"
  feature1_name <- "Gemcitabine"
  feature2_type <- "mRNA"
  data_type <- "all"
  tumor_type <- "all"
  
  # Define a custom progress callback function
  progress_callback <- function(current, total, elapsed_time) {
    percent_done <- round(current / total * 100)
    
    # Calculate estimated time remaining
    if (current > 0) {
      rate <- as.numeric(elapsed_time) / current
      remaining <- (total - current) * rate
      time_remaining <- format_time(remaining)
    } else {
      time_remaining <- "Calculating..."
    }
    
    cat(sprintf("\rProgress: %d/%d (%d%%) - Est. remaining: %s", 
                current, total, percent_done, time_remaining))
    
    # Flush the output for real-time updates
    flush.console()
  }
  
  # Run batch analysis with progress callback
  results <- BatchFindSigFeaturesPlus(
    feature1_type = feature1_type,
    feature1_name = feature1_name,
    feature2_type = feature2_type,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = 1,
    progress_callback = progress_callback
  )
  
  cat("\n\nAnalysis complete!\n")
  
  if (!is.null(results)) {
    # Create volcano plot with custom colors
    p <- plotMetaVolcano(
      results,
      es_t = 0.3,
      P_t = 0.05,
      label = TRUE,
      title = paste(feature1_name, "vs", feature2_type),
      custom_colors = c("Down" = "#3498db", "NS" = "#95a5a6", "Up" = "#e74c3c")
    )
    
    print(p)
  }
  
  return(results)
}

# Test 5: Testing the pair data functions directly
test_pair_functions <- function() {
  # Test pairDrugOmic_batch1 (continuous vs continuous)
  drug_name <- "5-Fluorouracil"
  gene_name <- "ABCB1"
  data_type <- "all"
  tumor_type <- "all"
  
  # Get the data
  myDrugs <- selFeatures("drug", drug_name, data_type, tumor_type)
  myOmics <- selFeatures("mRNA", gene_name, data_type, tumor_type)
  
  # Create pairs
  pairs <- pairDrugOmic_batch1(myOmics, myDrugs)
  
  # Perform meta-analysis on the pairs
  meta_result <- metaCalConCon(pairs)
  
  # Display results
  if (!is.null(meta_result)) {
    cat("Meta-analysis results for", drug_name, "vs", gene_name, ":\n")
    cat("Effect size (random):", meta_result$TE.random, "\n")
    cat("P-value (random):", meta_result$pval.random, "\n")
    cat("Number of studies:", length(meta_result$studlab), "\n")
    
    # Print forest plot
    meta::forest(meta_result, 
                 xlab = "Correlation coefficient", 
                 slab = "study")
  } else {
    cat("Insufficient data for meta-analysis.\n")
  }
}

# Run tests (uncomment as needed)
# results1 <- test_batch_con_con()
# results2 <- test_batch_con_dis()
# results3 <- test_batch_tumor_specific()
# results4 <- test_with_progress()
# test_pair_functions()

# Example usage
if (FALSE) {
  # Simple example workflow
  # 1. Run batch analysis
  results <- BatchFindSigFeaturesPlus(
    feature1_type = "drug",
    feature1_name = "Paclitaxel",
    feature2_type = "mRNA",
    data_type = "all",
    tumor_type = "all",
    cores = 2
  )
  
  # 2. Create and save volcano plot
  p <- plotMetaVolcano(results, 
                      es_t = 0.3, 
                      P_t = 0.05, 
                      label = TRUE,
                      title = "Paclitaxel vs mRNA")
  
  # Save plot to file
  # ggsave("paclitaxel_mRNA_volcano.pdf", p, width = 10, height = 8)
  
  # 3. Get top significant features
  top_hits <- results[results$p_value < 0.05 & abs(results$effect_size) > 0.3, ]
  top_hits <- top_hits[order(top_hits$p_value), ]
  head(top_hits, 20)
} 