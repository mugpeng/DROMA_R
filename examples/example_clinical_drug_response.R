#!/usr/bin/env Rscript

# Example script for DROMA.R package: Clinical Drug Response Analysis with CTRDB
# This example demonstrates how to analyze associations between clinical drug response and omics features using CTRDB data

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

# Setup: Connect to DROMA Database ----

# Note: Replace with your actual database path
# db_path <- "path/to/your/droma.sqlite"
db_path <- "Data/droma.sqlite"

# Connect to DROMA database
cat("Connecting to DROMA database...\n")
connectDROMADatabase(db_path)
cat("Database connected successfully!\n\n")

connectCTRDatabase("Data/ctrdb.sqlite")

# Check available clinical data
cat("Checking available clinical data in CTRDB...\n")
sample_anno <- getDROMAAnnotation("sample")

# Check if CTRDB data exists
if ("ProjectID" %in% colnames(sample_anno)) {
  ctrdb_samples <- sample_anno[sample_anno$ProjectID == "CTRDB" & !is.na(sample_anno$ProjectID), ]
  cat("Found", nrow(ctrdb_samples), "CTRDB samples\n")

  # Check available drugs in clinical data
  if ("CliUsedDrug" %in% colnames(ctrdb_samples)) {
    available_drugs <- unique(unlist(strsplit(ctrdb_samples$CliUsedDrug[!is.na(ctrdb_samples$CliUsedDrug)], "[;,|]")))
    available_drugs <- trimws(available_drugs[available_drugs != "" & available_drugs != "0"])
    cat("Available clinical drugs:\n")
    print(head(available_drugs, 10))  # Show first 10 drugs
  }

  # Check available datasets
  if ("PatientID" %in% colnames(ctrdb_samples)) {
    unique_datasets <- unique(ctrdb_samples$PatientID[!is.na(ctrdb_samples$PatientID)])
    cat("Number of unique datasets:", length(unique_datasets), "\n")
    cat("Example patient IDs:\n")
    print(head(unique_datasets, 5))
  }
} else {
  cat("No CTRDB data found in the database\n")
}

cat("\n")

# Example 1: Analyze clinical response to Erlotinib with EGFR expression ----

cat("Example 1: Analyzing clinical response to Erlotinib with EGFR expression...\n")
tryCatch({
  result_egfr_cis <- analyzeClinicalDrugResponse(
    select_omics = "TP53",
    select_drugs = "Cisplatin",
    data_type = "all",
    tumor_type = "all"
  )

  # Display results
  if (!is.null(result_egfr_erlotinib)) {
    cat("Analysis completed successfully!\n")

    # Show summary statistics
    if (!is.null(result_egfr_erlotinib$summary)) {
      cat("Summary of analyzed dataset:\n")
      print(result_egfr_erlotinib$summary)
    }

    # Display individual patient plots
    if (!is.null(result_egfr_erlotinib$plot)) {
      cat("Displaying individual patient plots...\n")
      print(result_egfr_erlotinib$plot)
    }

    # Display meta-analysis results
    if (!is.null(result_egfr_erlotinib$meta)) {
      cat("Meta-analysis results:\n")
      print(result_egfr_erlotinib$meta)

      # Create forest plot
      cat("Creating forest plot for meta-analysis...\n")
      result_egfr_erlotinib$forest_plot
    } else {
      cat("No meta-analysis performed (insufficient dataset or data)\n")
    }

    # Get detailed summary statistics
    detailed_summary <- getClinicalSummary(result_egfr_erlotinib)
    if (!is.null(detailed_summary)) {
      cat("Detailed summary statistics:\n")
      print(detailed_summary)
    }
  } else {
    cat("No results obtained for EGFR-Erlotinib analysis\n")
  }
}, error = function(e) {
  cat("Error in EGFR-Erlotinib analysis:", e$message, "\n")
})

cat("\n")

# Example 2: Analyze clinical response to Cisplatin with TP53 expression ----

cat("Example 2: Analyzing clinical response to Cisplatin with TP53 expression...\n")
tryCatch({
  result_tp53_cisplatin <- analyzeClinicalDrugResponse(
    select_omics = "TP53",
    select_drugs = "Cisplatin",
    data_type = "all",
    tumor_type = "all"
  )

  # Display results
  if (!is.null(result_tp53_cisplatin)) {
    cat("TP53-Cisplatin analysis completed!\n")

    # Show summary
    if (!is.null(result_tp53_cisplatin$summary)) {
      cat("Number of datasets analyzed:", nrow(result_tp53_cisplatin$summary), "\n")
      print(result_tp53_cisplatin$summary)
    }

    # Display plots
    if (!is.null(result_tp53_cisplatin$plot)) {
      print(result_tp53_cisplatin$plot)
    }

    # Meta-analysis results
    if (!is.null(result_tp53_cisplatin$meta)) {
      cat("TP53-Cisplatin meta-analysis results:\n")
      cat("Overall effect size:", round(result_tp53_cisplatin$meta$TE.random, 3), "\n")
      cat("P-value:", format(result_tp53_cisplatin$meta$pval.random, scientific = TRUE), "\n")

      # Forest plot
      result_tp53_cisplatin$forest_plot
    }
  } else {
    cat("No results for TP53-Cisplatin analysis\n")
  }
}, error = function(e) {
  cat("Error in TP53-Cisplatin analysis:", e$message, "\n")
})

cat("\n")

# Example 3: Analyze clinical response to Paclitaxel with BRCA1 expression ----

cat("Example 3: Analyzing clinical response to Paclitaxel with BRCA1 expression...\n")
tryCatch({
  result_brca1_paclitaxel <- analyzeClinicalDrugResponse(
    select_omics = "BRCA1",
    select_drugs = "Paclitaxel",
    data_type = "all",
    tumor_type = "all"
  )

  if (!is.null(result_brca1_paclitaxel)) {
    cat("BRCA1-Paclitaxel analysis completed!\n")

    # Display individual patient results
    if (!is.null(result_brca1_paclitaxel$data)) {
      cat("Analysis results for individual dataset:\n")
      for (patient_id in names(result_brca1_paclitaxel$data)) {
        patient_data <- result_brca1_paclitaxel$data[[patient_id]]
        cat("Patient", patient_id, ":\n")
        cat("  Response samples:", length(patient_data$response), "\n")
        cat("  Non-response samples:", length(patient_data$non_response), "\n")
        cat("  Mean expression (Response):", round(mean(patient_data$response), 3), "\n")
        cat("  Mean expression (Non-response):", round(mean(patient_data$non_response), 3), "\n")
      }
    }

    # Show plots
    if (!is.null(result_brca1_paclitaxel$plot)) {
      print(result_brca1_paclitaxel$plot)
    }
  }
}, error = function(e) {
  cat("Error in BRCA1-Paclitaxel analysis:", e$message, "\n")
})

cat("\n")

# Example 4: Analyze clinical response with different tumor types ----

cat("Example 4: Analyzing clinical response in specific tumor types...\n")

# Lung cancer specific analysis
tryCatch({
  result_lung <- analyzeClinicalDrugResponse(
    select_omics = "EGFR",
    select_drugs = "Erlotinib",
    data_type = "all",
    tumor_type = "lung cancer"
  )

  if (!is.null(result_lung) && !is.null(result_lung$summary)) {
    cat("Lung cancer specific analysis:\n")
    cat("dataset analyzed:", nrow(result_lung$summary), "\n")

    if (!is.null(result_lung$meta)) {
      cat("Effect size:", round(result_lung$meta$TE.random, 3), "\n")
      cat("P-value:", format(result_lung$meta$pval.random, scientific = TRUE), "\n")
    }
  } else {
    cat("No results for lung cancer specific analysis\n")
  }
}, error = function(e) {
  cat("Error in lung cancer analysis:", e$message, "\n")
})

cat("\n")

# Example 5: Compare multiple omics features for the same drug ----

cat("Example 5: Comparing multiple omics features for Cisplatin response...\n")

omics_features <- c("TP53", "BRCA1", "ERCC1", "XRCC1")
cisplatin_results <- list()

for (feature in omics_features) {
  tryCatch({
    cat("Analyzing", feature, "with Cisplatin...\n")
    result <- analyzeClinicalDrugResponse(
      select_omics = feature,
      select_drugs = "Cisplatin",
      data_type = "all",
      tumor_type = "all",
      meta_enabled = FALSE  # Skip meta-analysis for speed
    )

    if (!is.null(result) && !is.null(result$summary)) {
      cisplatin_results[[feature]] <- result
      cat("  dataset analyzed:", nrow(result$summary), "\n")
    } else {
      cat("  No results for", feature, "\n")
    }
  }, error = function(e) {
    cat("  Error analyzing", feature, ":", e$message, "\n")
  })
}

# Summary of multi-feature analysis
if (length(cisplatin_results) > 0) {
  cat("Summary of Cisplatin multi-feature analysis:\n")
  for (feature in names(cisplatin_results)) {
    result <- cisplatin_results[[feature]]
    if (!is.null(result$summary)) {
      total_datasets <- nrow(result$summary)
      cat("  ", feature, ": ", total_datasets, " datasets\n")
    }
  }
}

cat("\n")

# Example 6: Analyze clinical response with meta-analysis disabled ----

cat("Example 6: Single patient analysis (meta-analysis disabled)...\n")
tryCatch({
  result_single <- analyzeClinicalDrugResponse(
    select_omics = "EGFR",
    select_drugs = "Gefitinib",
    data_type = "all",
    tumor_type = "all",
    meta_enabled = FALSE
  )

  if (!is.null(result_single)) {
    cat("Single patient analysis completed!\n")

    # Show individual patient statistics
    if (!is.null(result_single$data)) {
      for (patient_id in names(result_single$data)) {
        patient_data <- result_single$data[[patient_id]]

        # Perform statistical test for this patient
        if (length(patient_data$response) >= 2 && length(patient_data$non_response) >= 2) {
          test_result <- wilcox.test(patient_data$non_response, patient_data$response)
          cat("Patient", patient_id, ":\n")
          cat("  Wilcoxon test p-value:", format(test_result$p.value, digits = 4), "\n")

          # Calculate effect size
          effect_size <- tryCatch({
            effsize::cliff.delta(patient_data$non_response, patient_data$response)$estimate
          }, error = function(e) NA)

          if (!is.na(effect_size)) {
            cat("  Effect size (Cliff's Delta):", round(effect_size, 3), "\n")
          }
        }
      }
    }

    # Display plots
    if (!is.null(result_single$plot)) {
      print(result_single$plot)
    }
  }
}, error = function(e) {
  cat("Error in single patient analysis:", e$message, "\n")
})

cat("\n")

# Example 7: Error handling and data validation ----

cat("Example 7: Demonstrating error handling...\n")

# Test with non-existent drug
tryCatch({
  result_error1 <- analyzeClinicalDrugResponse(
    select_omics = "EGFR",
    select_drugs = "NonExistentDrug",
    data_type = "all",
    tumor_type = "all"
  )
}, error = function(e) {
  cat("Expected error for non-existent drug:", e$message, "\n")
})

# Test with non-existent omics feature
tryCatch({
  result_error2 <- analyzeClinicalDrugResponse(
    select_omics = "NonExistentGene",
    select_drugs = "Cisplatin",
    data_type = "all",
    tumor_type = "all"
  )
}, error = function(e) {
  cat("Expected error for non-existent gene:", e$message, "\n")
})

# Test with empty parameters
tryCatch({
  result_error3 <- analyzeClinicalDrugResponse(
    select_omics = "",
    select_drugs = "",
    data_type = "all",
    tumor_type = "all"
  )
}, error = function(e) {
  cat("Expected error for empty parameters:", e$message, "\n")
})

cat("\n")

# Example 8: Advanced plotting and visualization ----

cat("Example 8: Advanced visualization techniques...\n")

# Try to create a comprehensive analysis with visualization
tryCatch({
  comprehensive_result <- analyzeClinicalDrugResponse(
    select_omics = "KRAS",
    select_drugs = "Cetuximab",
    data_type = "all",
    tumor_type = "all"
  )

  if (!is.null(comprehensive_result)) {
    cat("Comprehensive analysis completed for KRAS-Cetuximab!\n")

    # Create custom visualization combining individual plots
    if (!is.null(comprehensive_result$plot)) {
      # Add custom title and styling
      final_plot <- comprehensive_result$plot +
        plot_annotation(
          title = "Clinical Response Analysis: KRAS Expression vs Cetuximab Response",
          subtitle = "Individual patient analyses with statistical comparisons",
          caption = "Data source: CTRDB clinical database",
          theme = theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            plot.caption = element_text(size = 10, hjust = 1)
          )
        )

      print(final_plot)
    }

    # Create summary visualization
    if (!is.null(comprehensive_result$summary)) {
      summary_data <- comprehensive_result$summary

      # Create a summary plot showing sample sizes
      summary_plot <- ggplot(summary_data, aes(x = PatientID)) +
        geom_col(aes(y = N_Response), fill = "#FB8072FF", alpha = 0.7, position = "dodge") +
        geom_col(aes(y = -N_Non_Response), fill = "#BEBADAFF", alpha = 0.7, position = "dodge") +
        coord_flip() +
        labs(
          title = "Sample Distribution Across datasets",
          subtitle = "Response (positive) vs Non-Response (negative)",
          x = "Patient ID",
          y = "Number of Samples"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

      print(summary_plot)
    }
  }
}, error = function(e) {
  cat("Error in comprehensive analysis:", e$message, "\n")
})

cat("\n")

# Summary and Best Practices ----

cat("Clinical drug response analysis examples completed!\n\n")
cat("Key takeaways and best practices:\n")
cat("1. Use analyzeClinicalDrugResponse() for clinical drug-omics analysis\n")
cat("2. The function automatically handles CTRDB data filtering and patient grouping\n")
cat("3. Meta-analysis is performed across datasets when multiple valid datasets exist\n")
cat("4. Use getClinicalSummary() for detailed statistical summaries\n")
cat("5. Forest plots visualize meta-analysis results across datasets\n")
cat("6. Individual patient plots show response vs non-response comparisons\n")
cat("7. Error handling ensures robust analysis even with missing data\n")
cat("8. Filter by tumor_type for cancer-specific analyses\n")
cat("9. Set meta_enabled=FALSE for single patient or exploratory analyses\n")
cat("10. Always validate data availability before running large-scale analyses\n\n")

cat("Clinical analysis workflow:\n")
cat("  1. Connect to database with connectDROMADatabase()\n")
cat("  2. Check available clinical data with getDROMAAnnotation('sample')\n")
cat("  3. Run analysis with analyzeClinicalDrugResponse()\n")
cat("  4. Visualize results with built-in plots and forest plots\n")
cat("  5. Extract detailed statistics with getClinicalSummary()\n")
cat("  6. Interpret results in clinical context\n\n")

cat("Statistical methods used:\n")
cat("  - Wilcoxon rank-sum test for group comparisons\n")
cat("  - Cliff's Delta for effect size calculation\n")
cat("  - Random-effects meta-analysis across datasets\n")
cat("  - Forest plots for meta-analysis visualization\n")
cat("  - Boxplots with statistical annotations for individual datasets\n\n")

cat("Example completed successfully!\n")
