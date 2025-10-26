# Stratified Analysis Functions for ACC Data ----

# Core principle: Stratified analysis should follow the same pairing logic as
# the original analysis functions, filtering samples within each pairing operation
# rather than creating pre-filtered datasets
#
# NOTE: This implementation is designed for ACC (Activity Area Above the Curve) data
# where higher values indicate greater drug sensitivity

#' Load stratified drug response data
#'
#' @description Loads drug response data and stratifies samples into sensitive/resistant groups.
#' For ACC (Activity Area Above the Curve) data where higher values indicate greater sensitivity.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Name of the drug used for stratification
#' @param strata_quantile Quantile threshold for creating strata (default: 0.33)
#' @param min_samples Minimum samples required in each stratum (default: 5)
#' @return List containing stratification information for each project
#' @export
getStratificationInfo <- function(dromaset_object,
                                stratification_drug,
                                strata_quantile = 0.33,
                                min_samples = 5) {

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object")
  }

  # Load stratification drug response data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    strat_data <- loadTreatmentResponse(dromaset_object,
                                                 select_drugs = stratification_drug,
                                                 data_type = "all",
                                                 tumor_type = "all",
                                                 return_data = TRUE,
                                                 zscore = TRUE)

    if (!is.matrix(strat_data) || !stratification_drug %in% rownames(strat_data)) {
      stop("Stratification drug '", stratification_drug, "' not found in dataset")
    }

    # Extract drug response vector
    drug_response <- as.numeric(strat_data[stratification_drug, ])
    names(drug_response) <- colnames(strat_data)
    drug_response <- drug_response[!is.na(drug_response)]

    if (length(drug_response) < (2 * min_samples)) {
      stop("Insufficient samples (", length(drug_response), ") for stratification. Minimum required: ", 2 * min_samples)
    }

    # Calculate quantiles
    lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
    upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

    # Create strata - ACC data: higher values = more sensitive
    sensitive_samples <- names(drug_response[drug_response >= upper_quantile])
    resistant_samples <- names(drug_response[drug_response <= lower_quantile])

    # Validate sample sizes
    if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
      stop("Insufficient samples in strata. Sensitive: ", length(sensitive_samples),
           ", Resistant: ", length(resistant_samples), ". Minimum required: ", min_samples)
    }

    result <- list()
    result[[dromaset_object@name]] <- list(
      sensitive = sensitive_samples,
      resistant = resistant_samples,
      thresholds = c(lower = lower_quantile, upper = upper_quantile)
    )

    attr(result, "type") <- "single"

  } else {
    # MultiDromaSet - handle each project separately
    strat_data_list <- loadMultiProjectTreatmentResponse(
      dromaset_object,
      select_drugs = stratification_drug,
      overlap_only = FALSE,
      data_type = "all",
      tumor_type = "all",
      zscore = TRUE
    )

    # Check if any projects have the stratification drug
    projects_with_drug <- names(strat_data_list)[
      sapply(strat_data_list, function(mat) is.matrix(mat) && stratification_drug %in% rownames(mat))
    ]

    if (length(projects_with_drug) == 0) {
      stop("Stratification drug '", stratification_drug, "' not found in any project")
    }

    # Warn about projects without the drug
    all_projects <- names(dromaset_object@DromaSets)
    if (length(all_projects) > length(projects_with_drug)) {
      projects_without <- setdiff(all_projects, projects_with_drug)
      warning("Stratification drug not found in projects: ",
              paste(projects_without, collapse = ", "))
    }

    stratification_info <- list()

    for (projects in names(strat_data_list)) {
      strat_data <- strat_data_list[[projects]]

      if (!is.matrix(strat_data) || !stratification_drug %in% rownames(strat_data)) {
        # Skip projects without the drug
        next
      }

      # Extract drug response vector
      drug_response <- as.numeric(strat_data[stratification_drug, ])
      names(drug_response) <- colnames(strat_data)
      drug_response <- drug_response[!is.na(drug_response)]

      if (length(drug_response) < (2 * min_samples)) {
        warning("Project ", projects, ": Insufficient samples (", length(drug_response),
                ") for stratification. Minimum required: ", 2 * min_samples)
        next
      }

      # Calculate quantiles
      lower_quantile <- quantile(drug_response, strata_quantile, na.rm = TRUE)
      upper_quantile <- quantile(drug_response, 1 - strata_quantile, na.rm = TRUE)

      # Create strata - ACC data: higher values = more sensitive
      sensitive_samples <- names(drug_response[drug_response >= upper_quantile])
      resistant_samples <- names(drug_response[drug_response <= lower_quantile])

      # Skip if insufficient samples
      if (length(sensitive_samples) < min_samples || length(resistant_samples) < min_samples) {
        warning("Project ", projects, ": Insufficient samples in strata. Sensitive: ",
                length(sensitive_samples), ", Resistant: ", length(resistant_samples),
                ". Minimum required: ", min_samples)
        next
      }

      stratification_info[[projects]] <- list(
        sensitive = sensitive_samples,
        resistant = resistant_samples,
        thresholds = c(lower = lower_quantile, upper = upper_quantile)
      )
    }

    if (length(stratification_info) == 0) {
      stop("No projects had sufficient samples for stratification")
    }

    result <- stratification_info
    attr(result, "type") <- "multi"
  }

  return(result)
}

#' Analyze stratified drug-omic associations
#'
#' @description Performs drug-omic association analysis within strata defined by
#' response to a stratification drug. This approach follows the pairing principle
#' by filtering samples within each pairing operation rather than creating pre-filtered datasets.
#' Designed for ACC (Activity Area Above the Curve) data where higher values indicate greater sensitivity.
#' @param dromaset_object A DromaSet or MultiDromaSet object
#' @param stratification_drug Drug used to define strata (e.g., "cisplatin")
#' @param strata_quantile Quantile threshold for stratification (default: 0.33)
#' @param feature_type Type of omics data to analyze
#' @param select_features Name of the specific omics feature
#' @param select_drugs Name of the target drug for association analysis
#' @param min_samples Minimum samples per stratum (default: 5)
#' @param ... Additional parameters passed to the analysis functions
#' @return List containing stratified analysis results and between-stratum comparisons
#' @export
#' @examples
#' \dontrun{
#' # Stratified analysis of ERCC1 expression and bortezomib response
#' # based on cisplatin sensitivity
#' result <- analyzeStratifiedDrugOmic(
#'   dromaset_object = multi_set,
#'   stratification_drug = "cisplatin",
#'   strata_quantile = 0.33,
#'   feature_type = "mRNA",
#'   select_features = "ERCC1",
#'   select_drugs = "bortezomib",
#'   data_type = "CellLine",
#'   tumor_type = "all"
#' )
#' }
analyzeStratifiedDrugOmic <- function(dromaset_object,
                                    stratification_drug,
                                    strata_quantile = 0.33,
                                    feature_type,
                                    select_features,
                                    select_drugs,
                                    min_samples = 5,
                                    ...) {

  # Validate inputs
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object")
  }

  if (missing(feature_type) || missing(select_features) || missing(select_drugs)) {
    stop("feature_type, select_features, and select_drugs are required")
  }

  # Extract additional parameters
  extra_params <- list(...)
  data_type <- extra_params$data_type %||% "all"
  tumor_type <- extra_params$tumor_type %||% "all"
  overlap_only <- extra_params$overlap_only %||% FALSE
  merged_enabled <- extra_params$merged_enabled %||% TRUE
  meta_enabled <- extra_params$meta_enabled %||% TRUE

  # Step 1: Get stratification information (sensitive/resistant sample names for each project)
  cat("Getting stratification info...\n")
  stratification_info <- getStratificationInfo(
    dromaset_object = dromaset_object,
    stratification_drug = stratification_drug,
    strata_quantile = strata_quantile,
    min_samples = min_samples
  )

  # Step 2: Load original data (same as analyzeDrugOmicPair)
  cat("Loading drug data...\n")
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    drug_data <- loadTreatmentResponse(dromaset_object,
                                              select_drugs = select_drugs,
                                              data_type = data_type,
                                              tumor_type = tumor_type,
                                              return_data = TRUE,
                                              zscore = TRUE)

    # Convert matrix to list format
    if (is.matrix(drug_data) && select_drugs %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[select_drugs, ])
      names(drug_vector) <- colnames(drug_data)
      drug_data_list <- list()
      drug_data_list[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    } else {
      drug_data_list <- list()
    }
  } else {
    # MultiDromaSet
    drug_data_list <- loadMultiProjectTreatmentResponse(
      dromaset_object,
      select_drugs = select_drugs,
      overlap_only = overlap_only,
      data_type = data_type,
      tumor_type = tumor_type,
      zscore = TRUE
    )

    # Extract specific drug from each project
    drug_data_list <- lapply(drug_data_list, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })
    drug_data_list <- drug_data_list[!sapply(drug_data_list, is.null)]
  }

  # Load omics data
  cat("Loading omics data...\n")
  if (feature_type %in% c("drug", "drug_raw")) {
    # Another drug
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadTreatmentResponse(dromaset_object,
                                                   select_drugs = select_features,
                                                   data_type = data_type,
                                                   tumor_type = tumor_type,
                                                   return_data = TRUE,
                                                   zscore = TRUE)

      if (is.matrix(omics_data) && select_features %in% rownames(omics_data)) {
        omics_vector <- as.numeric(omics_data[select_features, ])
        names(omics_vector) <- colnames(omics_data)
        omics_data_list <- list()
        omics_data_list[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
      } else {
        omics_data_list <- list()
      }
    } else {
      omics_data_list <- loadMultiProjectTreatmentResponse(
        dromaset_object,
        select_drugs = select_features,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type,
        zscore = TRUE
      )

      omics_data_list <- lapply(omics_data_list, function(drug_matrix) {
        if (is.matrix(drug_matrix) && select_features %in% rownames(drug_matrix)) {
          drug_vector <- as.numeric(drug_matrix[select_features, ])
          names(drug_vector) <- colnames(drug_matrix)
          return(drug_vector[!is.na(drug_vector)])
        }
        return(NULL)
      })
      omics_data_list <- omics_data_list[!sapply(omics_data_list, is.null)]
    }
  } else {
    # Molecular profiles
    if (inherits(dromaset_object, "DromaSet")) {
      omics_data <- loadMolecularProfiles(dromaset_object,
                                                  feature_type = feature_type,
                                                  select_features = select_features,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  return_data = TRUE,
                                                  zscore = TRUE)

      if (is.matrix(omics_data) && select_features %in% rownames(omics_data)) {
        omics_data_list <- list()
        if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
          # Continuous
          omics_vector <- as.numeric(omics_data[select_features, ])
          names(omics_vector) <- colnames(omics_data)
          omics_data_list[[dromaset_object@name]] <- omics_vector[!is.na(omics_vector)]
        } else {
          # Discrete
          gene_row <- omics_data[select_features, ]
          present_samples <- names(gene_row)[gene_row != 0]
          omics_data_list[[dromaset_object@name]] <- present_samples
        }
      } else {
        omics_data_list <- list()
      }
    } else {
      omics_data_list <- loadMultiProjectMolecularProfiles(
        dromaset_object,
        feature_type = feature_type,
        select_features = select_features,
        overlap_only = overlap_only,
        data_type = data_type,
        tumor_type = tumor_type,
        zscore = TRUE
      )

      omics_data_list <- lapply(omics_data_list, function(omics_matrix) {
        if (is.matrix(omics_matrix) && select_features %in% rownames(omics_matrix)) {
          if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
            # Continuous
            omics_vector <- as.numeric(omics_matrix[select_features, ])
            names(omics_vector) <- colnames(omics_matrix)
            return(omics_vector[!is.na(omics_vector)])
          } else {
            # Discrete
            gene_row <- omics_matrix[select_features, ]
            present_samples <- names(gene_row)[gene_row != 0]
            return(present_samples)
          }
        }
        return(NULL)
      })
      omics_data_list <- omics_data_list[!sapply(omics_data_list, is.null)]
    }
  }

  # Check if we have data
  if (length(drug_data_list) == 0 || length(omics_data_list) == 0) {
    stop("No data found for the specified drug-omics pair")
  }

  # Step 3: Filter data for each stratum and analyze
  cat("Analyzing sensitive group...\n")
  sensitive_loaded_data <- list(drugs = drug_data_list, omics = omics_data_list)
  sensitive_data <- filterStratifiedData(sensitive_loaded_data, stratification_info, "sensitive")
  sensitive_result <- analyzeStratifiedData(
    omics_data = sensitive_data$omics,
    drug_data = sensitive_data$drugs,
    feature_type = feature_type,
    merged_enabled = merged_enabled,
    meta_enabled = meta_enabled
  )

  cat("Analyzing resistant group...\n")
  resistant_loaded_data <- list(drugs = drug_data_list, omics = omics_data_list)
  resistant_data <- filterStratifiedData(resistant_loaded_data, stratification_info, "resistant")
  resistant_result <- analyzeStratifiedData(
    omics_data = resistant_data$omics,
    drug_data = resistant_data$drugs,
    feature_type = feature_type,
    merged_enabled = merged_enabled,
    meta_enabled = meta_enabled
  )

  # Compare results between strata
  comparison_result <- compareStratifiedResults(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    stratification_info = stratification_info,
    feature_type = feature_type,
    select_features = select_features,
    select_drugs = select_drugs
  )

  # Create omics expression comparison plot
  cat("Creating omics expression comparison plot...\n")
  expression_plot <- createStratifiedOmicExpressionPlot(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    feature_type = feature_type,
    select_features = select_features,
    stratification_drug = stratification_drug
  )

  # Add expression plot to comparison results
  comparison_result$expression_plot <- expression_plot

  # Compile final results
  result <- list(
    sensitive = sensitive_result,
    resistant = resistant_result,
    comparison = comparison_result,
    stratification_info = stratification_info,
    analysis_parameters = list(
      stratification_drug = stratification_drug,
      strata_quantile = strata_quantile,
      feature_type = feature_type,
      select_features = select_features,
      select_drugs = select_drugs,
      min_samples = min_samples,
      extra_params = extra_params
    )
  )

  class(result) <- "StratifiedDrugOmicResult"

  cat("Stratified analysis complete!\n")
  return(result)
}

#' Compare results between strata
#'
#' @description Compares association results between sensitive and resistant groups
#' @param sensitive_result Analysis result from sensitive group
#' @param resistant_result Analysis result from resistant group
#' @param stratification_info Information about stratification
#' @param feature_type Type of omics data analyzed
#' @param select_features Name of the omics feature
#' @param select_drugs Name of the drug analyzed
#' @return Comparison results with statistical tests
compareStratifiedResults <- function(sensitive_result, resistant_result, stratification_info,
                                   feature_type, select_features, select_drugs) {

  comparison <- list()

  # Extract correlation coefficients if available
  if (!is.null(sensitive_result$meta) && !is.null(resistant_result$meta)) {
    if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
      # Continuous-continuous analysis
      sens_cor <- sensitive_result$meta$TE.random
      sens_se <- sensitive_result$meta$seTE.random
      res_cor <- resistant_result$meta$TE.random
      res_se <- resistant_result$meta$seTE.random

      # Test for difference in correlations
      z_diff <- (sens_cor - res_cor) / sqrt(sens_se^2 + res_se^2)
      p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)

      comparison$correlation_comparison <- data.frame(
        Group = c("Sensitive", "Resistant", "Difference"),
        Correlation = c(sens_cor, res_cor, sens_cor - res_cor),
        SE = c(sens_se, res_se, sqrt(sens_se^2 + res_se^2)),
        P_value = c(sensitive_result$meta$pval.random, resistant_result$meta$pval.random, p_diff)
      )
    } else {
      # Discrete analysis - compare effect sizes
      sens_es <- sensitive_result$meta$TE.random
      sens_se <- sensitive_result$meta$seTE.random
      res_es <- resistant_result$meta$TE.random
      res_se <- resistant_result$meta$seTE.random

      # Test for difference in effect sizes
      z_diff <- (sens_es - res_es) / sqrt(sens_se^2 + res_se^2)
      p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)

      comparison$effect_comparison <- data.frame(
        Group = c("Sensitive", "Resistant", "Difference"),
        Effect_Size = c(sens_es, res_es, sens_es - res_es),
        SE = c(sens_se, res_se, sqrt(sens_se^2 + res_se^2)),
        P_value = c(sensitive_result$meta$pval.random, resistant_result$meta$pval.random, p_diff)
      )
    }
  }

  # Create comparison plot
  comparison$comparison_plot <- createStratifiedComparisonPlot(
    sensitive_result = sensitive_result,
    resistant_result = resistant_result,
    select_features = select_features,
    select_drugs = select_drugs,
    feature_type = feature_type
  )

  return(comparison)
}

#' Create stratified comparison plot
#'
#' @description Creates a visualization comparing results between strata
#' @param sensitive_result Results from sensitive group
#' @param resistant_result Results from resistant group
#' @param select_features Name of omics feature
#' @param select_drugs Name of drug
#' @param feature_type Type of omics data
#' @param p_value_digits Number of decimal places for p-values (default: 3)
#' @param add_difference_row Whether to add a row showing difference between groups (default: TRUE)
#' @return ggplot object
createStratifiedComparisonPlot <- function(sensitive_result, resistant_result,
                                          select_features, select_drugs,
                                          feature_type,
                                          p_value_digits = 3,
                                          add_difference_row = TRUE) {

  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for plotting")
    return(NULL)
  }

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    warning("gridExtra package is required for plotting")
    return(NULL)
  }

  # Check if we have valid results
  if (is.null(sensitive_result) || is.null(resistant_result)) {
    warning("Missing results for one or both strata")
    return(NULL)
  }

  # Format p-value for display
  format_pvalue <- function(p, digits = 3) {
    if (is.na(p) || is.null(p)) return("NA")
    if (p < 0.001) return("<0.001")
    if (p < 0.01) return(sprintf("<%.2f", p))
    if (p < 0.05) return(sprintf("<%.3f", p))
    return(sprintf("%.3f", round(p, digits)))
  }

  # Create plot based on omics type
  if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
    # Continuous-continuous analysis - create forest plot
    plot_list <- list()

    # Extract correlations if available
    if (!is.null(sensitive_result$plot) && !is.null(resistant_result$plot)) {
      # Create a comparison forest plot with p-values
      if (!is.null(sensitive_result$meta) && !is.null(resistant_result$meta)) {
        # Extract values
        sens_cor <- sensitive_result$meta$TE.random
        sens_lower <- sensitive_result$meta$lower.random
        sens_upper <- sensitive_result$meta$upper.random
        sens_p <- sensitive_result$meta$pval.random

        res_cor <- resistant_result$meta$TE.random
        res_lower <- resistant_result$meta$lower.random
        res_upper <- resistant_result$meta$upper.random
        res_p <- resistant_result$meta$pval.random

        # Calculate difference and p-value if requested
        if (add_difference_row && !is.na(sens_cor) && !is.na(res_cor)) {
          sens_se <- sensitive_result$meta$seTE.random
          res_se <- resistant_result$meta$seTE.random
          z_diff <- (sens_cor - res_cor) / sqrt(sens_se^2 + res_se^2)
          p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)
        } else {
          p_diff <- NA
        }

        # Create plot data
        groups <- c("Sensitive", "Resistant")
        correlations <- c(sens_cor, res_cor)
        ci_lowers <- c(sens_lower, res_lower)
        ci_uppers <- c(sens_upper, res_upper)
        p_values <- c(sens_p, res_p)

        # Add difference row if requested
        if (add_difference_row && !is.na(p_diff)) {
          groups <- c(groups, "Difference")
          correlations <- c(correlations, sens_cor - res_cor)
          ci_lowers <- c(ci_lowers, NA)
          ci_uppers <- c(ci_uppers, NA)
          p_values <- c(p_values, p_diff)
        }

        # Create group labels with p-values using newline
        p_labels <- sapply(p_values, format_pvalue, digits = p_value_digits)
        group_labels <- paste0(groups, "\np=", p_labels)

        plot_data <- data.frame(
          Group = group_labels,
          Group_Name = groups,  # Keep original group names for reference
          Correlation = correlations,
          CI_lower = ci_lowers,
          CI_upper = ci_uppers,
          P_value = p_values,
          stringsAsFactors = FALSE
        )

        # Create forest plot with p-values in group labels
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Correlation)) +
          ggplot2::geom_point(size = 3) +
          ggplot2::geom_errorbar(
            data = plot_data[!is.na(plot_data$CI_lower), ],
            ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
            width = 0.2
          ) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          ggplot2::coord_flip(ylim = range(plot_data$Correlation, plot_data$CI_lower, plot_data$CI_upper, na.rm = TRUE) * c(1.1, 1.1)) +
          ggplot2::labs(
            title = paste("Stratified Analysis:", select_features, "vs", select_drugs),
            subtitle = "Correlation with 95% Confidence Intervals",
            y = "Correlation Coefficient",
            x = ""
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            axis.text.y = ggplot2::element_text(hjust = 1, face = "bold", lineheight = 0.8),
            panel.grid.major.y = ggplot2::element_blank()
          )

        # Note: Significance indicators removed as requested

        plot_list$forest_plot <- p
      }
    }

    # Combine individual plots if available
    if (!is.null(sensitive_result$plot) && !is.null(resistant_result$plot)) {
      # Arrange plots side by side
      if (requireNamespace("patchwork", quietly = TRUE)) {
        # Add p-values to individual plots if available
        if (!is.null(sensitive_result$meta$pval.random)) {
          sens_p_label <- format_pvalue(sensitive_result$meta$pval.random, p_value_digits)
          sensitive_plot <- sensitive_result$plot +
            ggplot2::labs(caption = paste0("p = ", sens_p_label))
        } else {
          sensitive_plot <- sensitive_result$plot
        }

        if (!is.null(resistant_result$meta$pval.random)) {
          res_p_label <- format_pvalue(resistant_result$meta$pval.random, p_value_digits)
          resistant_plot <- resistant_result$plot +
            ggplot2::labs(caption = paste0("p = ", res_p_label))
        } else {
          resistant_plot <- resistant_result$plot
        }

        combined_plot <- patchwork::wrap_plots(
          list(Sensitive = sensitive_plot,
               Resistant = resistant_plot),
          ncol = 2
        ) +
          patchwork::plot_annotation(
            title = paste("Stratified Analysis:", select_features, "vs", select_drugs),
            subtitle = "Individual stratum scatter plots"
          )

        plot_list$individual_plots <- combined_plot
      } else {
        # Fallback: return plots in a list
        plot_list$individual_plots <- list(
          Sensitive = sensitive_result$plot,
          Resistant = resistant_result$plot
        )
      }
    }

    # Return both main plot and individual plots in a named list
    result <- list()
    if (!is.null(plot_list$forest_plot)) {
      result$forest_plot <- plot_list$forest_plot
    }
    if (!is.null(plot_list$individual_plots)) {
      result$individual_plots <- plot_list$individual_plots
    }

    if (length(result) == 0) {
      return(NULL)
    } else {
      # Add class for method dispatch
      class(result) <- "StratifiedComparisonPlots"
      return(result)
    }

  } else {
    # Discrete analysis - create effect size comparison
    if (!is.null(sensitive_result$meta) && !is.null(resistant_result$meta)) {
      # Extract values
      sens_es <- sensitive_result$meta$TE.random
      sens_lower <- sensitive_result$meta$lower.random
      sens_upper <- sensitive_result$meta$upper.random
      sens_p <- sensitive_result$meta$pval.random

      res_es <- resistant_result$meta$TE.random
      res_lower <- resistant_result$meta$lower.random
      res_upper <- resistant_result$meta$upper.random
      res_p <- resistant_result$meta$pval.random

      # Calculate difference and p-value if requested
      if (add_difference_row && !is.na(sens_es) && !is.na(res_es)) {
        sens_se <- sensitive_result$meta$seTE.random
        res_se <- resistant_result$meta$seTE.random
        z_diff <- (sens_es - res_es) / sqrt(sens_se^2 + res_se^2)
        p_diff <- 2 * pnorm(abs(z_diff), lower.tail = FALSE)
      } else {
        p_diff <- NA
      }

      # Create plot data
      groups <- c("Sensitive", "Resistant")
      effect_sizes <- c(sens_es, res_es)
      ci_lowers <- c(sens_lower, res_lower)
      ci_uppers <- c(sens_upper, res_upper)
      p_values <- c(sens_p, res_p)

      # Add difference row if requested
      if (add_difference_row && !is.na(p_diff)) {
        groups <- c(groups, "Difference")
        effect_sizes <- c(effect_sizes, sens_es - res_es)
        ci_lowers <- c(ci_lowers, NA)
        ci_uppers <- c(ci_uppers, NA)
        p_values <- c(p_values, p_diff)
      }

      # Create group labels with p-values using newline
      p_labels <- sapply(p_values, format_pvalue, digits = p_value_digits)
      group_labels <- paste0(groups, "\np=", p_labels)

      plot_data <- data.frame(
        Group = group_labels,
        Group_Name = groups,  # Keep original group names for reference
        Effect_Size = effect_sizes,
        CI_lower = ci_lowers,
        CI_upper = ci_uppers,
        P_value = p_values,
        stringsAsFactors = FALSE
      )

      # Create forest plot for effect sizes with p-values
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Effect_Size)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(
          data = plot_data[!is.na(plot_data$CI_lower), ],
          ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
          width = 0.2
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::coord_flip(ylim = range(plot_data$Effect_Size, plot_data$CI_lower, plot_data$CI_upper, na.rm = TRUE) * c(1.1, 1.1)) +
        ggplot2::labs(
          title = paste("Stratified Analysis:", select_features, "vs", select_drugs),
          subtitle = "Effect Sizes with 95% Confidence Intervals",
          y = "Effect Size (Mean Difference)",
          x = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          plot.subtitle = ggplot2::element_text(hjust = 0.5),
          axis.text.y = ggplot2::element_text(hjust = 1, face = "bold", lineheight = 0.8),
          panel.grid.major.y = ggplot2::element_blank()
        )

      # Note: Significance indicators removed as requested

      return(p)
    }

    # If no meta-analysis, try to combine individual plots
    if (!is.null(sensitive_result$plot) && !is.null(resistant_result$plot)) {
      if (requireNamespace("patchwork", quietly = TRUE)) {
        # Add p-values to individual plots if available
        if (!is.null(sensitive_result$meta$pval.random)) {
          sens_p_label <- format_pvalue(sensitive_result$meta$pval.random, p_value_digits)
          sensitive_plot <- sensitive_result$plot +
            ggplot2::labs(caption = paste0("p = ", sens_p_label))
        } else {
          sensitive_plot <- sensitive_result$plot
        }

        if (!is.null(resistant_result$meta$pval.random)) {
          res_p_label <- format_pvalue(resistant_result$meta$pval.random, p_value_digits)
          resistant_plot <- resistant_result$plot +
            ggplot2::labs(caption = paste0("p = ", res_p_label))
        } else {
          resistant_plot <- resistant_result$plot
        }

        return(patchwork::wrap_plots(
          list(Sensitive = sensitive_plot,
               Resistant = resistant_plot),
          ncol = 2
        ) +
          patchwork::plot_annotation(
            title = paste("Stratified Analysis:", select_features, "vs", select_drugs)
          ))
      } else {
        return(list(
          Sensitive = sensitive_result$plot,
          Resistant = resistant_result$plot
        ))
      }
    }
  }

  return(NULL)
}

#' Filter drug data for stratified analysis
#'
#' @description Filters loaded drug data to keep only samples from specified stratum.
#' Only filters datasets that are present in stratification_info.
#' @param loaded_data List containing loaded omics and drug data
#' @param stratification_info Stratification information from getStratificationInfo
#' @param stratum Which stratum ("sensitive" or "resistant")
#' @return List with filtered drug data (omics data is returned unchanged)
#' @export
filterStratifiedData <- function(loaded_data, stratification_info, stratum) {

  # Validate stratum
  stratum <- match.arg(stratum, c("sensitive", "resistant"))

  # Only filter drug data, keep omics data as is
  filtered_drugs <- list()

  # Filter only for projects that have stratification info
  for (projects in names(stratification_info)) {
    if (projects %in% names(loaded_data$drugs)) {
      # Get stratum samples for this project
      stratum_samples <- stratification_info[[projects]][[stratum]]

      # Filter drug data
      drug_data <- loaded_data$drugs[[projects]]
      if (!is.null(drug_data)) {
        filtered_drug_data <- drug_data[names(drug_data) %in% stratum_samples]
        if (length(filtered_drug_data) >= 3) {
          filtered_drugs[[projects]] <- filtered_drug_data
        }
      }
    }
  }

  # Return omics data unchanged and filtered drug data
  return(list(omics = loaded_data$omics, drugs = filtered_drugs))
}

#' Analyze stratified data
#'
#' @description Performs drug-omics analysis on pre-filtered stratified data
#' @param omics_data Filtered omics data
#' @param drug_data Filtered drug data
#' @param feature_type Type of omics data
#' @param merged_enabled Whether to create merged dataset
#' @param meta_enabled Whether to perform meta-analysis
#' @return Analysis results
#' @export
analyzeStratifiedData <- function(omics_data, drug_data, feature_type,
                                 merged_enabled = TRUE, meta_enabled = TRUE) {

  # Check if we have data
  if (length(drug_data) == 0 || length(omics_data) == 0) {
    return(list())
  }

  # Perform analysis based on omics type
  if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
    # Continuous omics data
    pairs <- pairContinuousFeatures(omics_data, drug_data, merged = merged_enabled)

    if (!is.null(pairs)) {
      result <- list(
        plot = plotMultipleCorrelations(
          lapply(pairs, function(p) {
            if (!is.null(p$feature1) && !is.null(p$feature2)) {
              list(feature1 = p$feature1, feature2 = p$feature2)
            } else {
              NULL
            }
          })[!sapply(lapply(pairs, function(p) {
            if (!is.null(p$feature1) && !is.null(p$feature2)) {
              list(feature1 = p$feature1, feature2 = p$feature2)
            } else {
              NULL
            }
          }), is.null)],
          x_label = "Omic Feature",
          y_label = "Drug Response"
        ),
        data = pairs
      )

      # Perform meta-analysis if multiple studies
      if (length(pairs) > 1 && meta_enabled) {
        # Remove merged dataset for meta-analysis
        if ("merged_dataset" %in% names(pairs)) {
          meta_pairs <- pairs[names(pairs) != "merged_dataset"]
        } else {
          meta_pairs <- pairs
        }

        if (length(meta_pairs) > 0) {
          result$meta <- metaCalcConCon(meta_pairs)
        }
      }
    } else {
      result <- list()
    }
  } else {
    # Discrete omics data
    pairs <- pairDiscreteFeatures(omics_data, drug_data, merged = merged_enabled)

    if (!is.null(pairs)) {
      result <- list(
        plot = plotMultipleGroupComparisons(pairs, y_label = "Drug Response"),
        data = pairs
      )

      # Perform meta-analysis if multiple studies
      if (length(pairs) > 1 && meta_enabled) {
        # Remove merged dataset for meta-analysis
        if ("merged_dataset" %in% names(pairs)) {
          meta_pairs <- pairs[names(pairs) != "merged_dataset"]
        } else {
          meta_pairs <- pairs
        }

        if (length(meta_pairs) > 0) {
          result$meta <- metaCalcConDis(meta_pairs)
        }
      }
    } else {
      result <- list()
    }
  }

  return(result)
}

# Print method for StratifiedComparisonPlots
#' @method print StratifiedComparisonPlots
print.StratifiedComparisonPlots <- function(x, ...) {
  cat("Stratified Comparison Plots\n")
  cat("========================\n")
  if (!is.null(x$forest_plot)) {
    cat("- Forest plot: Available (x$forest_plot)\n")
  }
  if (!is.null(x$individual_plots)) {
    cat("- Individual plots: Available (x$individual_plots)\n")
  }
  cat("\nUse plot() to display the forest plot, or access plots directly.\n")
  invisible(x)
}

# Plot method for StratifiedComparisonPlots
#' @method plot StratifiedComparisonPlots
plot.StratifiedComparisonPlots <- function(x, which = "forest", ...) {
  if (which == "forest" && !is.null(x$forest_plot)) {
    return(x$forest_plot)
  } else if (which == "individual" && !is.null(x$individual_plots)) {
    return(x$individual_plots)
  } else if (is.null(x$forest_plot) && is.null(x$individual_plots)) {
    warning("No plots available")
    return(NULL)
  } else if (which == "forest" && is.null(x$forest_plot)) {
    warning("Forest plot not available. Use which = 'individual' for individual plots.")
    return(x$individual_plots)
  } else if (which == "individual" && is.null(x$individual_plots)) {
    warning("Individual plots not available. Use which = 'forest' for forest plot.")
    return(x$forest_plot)
  }
}

#' Create stratified omics expression comparison plot
#'
#' @description Creates a boxplot comparing omics expression between sensitive and resistant groups
#' defined by response to a stratification drug. Uses existing data from sensitive_result and resistant_result.
#' Mimics the style from FuncPairVisualization.
#' @param sensitive_result Analysis result from sensitive group
#' @param resistant_result Analysis result from resistant group
#' @param feature_type Type of omics data to visualize
#' @param select_features Name of the specific omics feature
#' @param stratification_drug Name of the drug used for stratification
#' @return A ggplot2 object with boxplot comparing omics expression between strata
createStratifiedOmicExpressionPlot <- function(sensitive_result,
                                               resistant_result,
                                               feature_type,
                                               select_features,
                                               stratification_drug) {

  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for plotting")
    return(NULL)
  }

  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    warning("ggpubr package is required for plotting")
    return(NULL)
  }

  # Extract omics data from sensitive_result
  sensitive_omics <- NULL
  if (!is.null(sensitive_result$data)) {
    # Get omics data from the paired data
    if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
      # Continuous data
      if ("merged_dataset" %in% names(sensitive_result$data)) {
        sensitive_omics <- sensitive_result$data$merged_dataset$feature1
      } else if (length(sensitive_result$data) > 0) {
        # Take the first project with data
        first_project <- names(sensitive_result$data)[1]
        sensitive_omics <- sensitive_result$data[[first_project]]$feature1
      }
    } else {
      # Discrete data - convert to presence/absence
      # For discrete data, we need to reconstruct which samples have the feature
      all_sensitive_samples <- c()
      for (projects in names(sensitive_result$data)) {
        if (projects != "merged_dataset") {
          yes_samples <- sensitive_result$data[[projects]]$yes
          all_sensitive_samples <- c(all_sensitive_samples,
                                    paste0(projects, "_", yes_samples))
        }
      }
      # Create binary vector (1 = present, 0 = absent)
      sensitive_omics <- rep(0, length(all_sensitive_samples))
      names(sensitive_omics) <- all_sensitive_samples
      sensitive_omics[all_sensitive_samples %in% unlist(lapply(sensitive_result$data[names(sensitive_result$data) != "merged_dataset"], function(x) x$yes))] <- 1
    }
  }

  # Extract omics data from resistant_result
  resistant_omics <- NULL
  if (!is.null(resistant_result$data)) {
    # Get omics data from the paired data
    if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms", "drug")) {
      # Continuous data
      if ("merged_dataset" %in% names(resistant_result$data)) {
        resistant_omics <- resistant_result$data$merged_dataset$feature1
      } else if (length(resistant_result$data) > 0) {
        # Take the first project with data
        first_project <- names(resistant_result$data)[1]
        resistant_omics <- resistant_result$data[[first_project]]$feature1
      }
    } else {
      # Discrete data - convert to presence/absence
      all_resistant_samples <- c()
      for (projects in names(resistant_result$data)) {
        if (projects != "merged_dataset") {
          yes_samples <- resistant_result$data[[projects]]$yes
          all_resistant_samples <- c(all_resistant_samples,
                                   paste0(projects, "_", yes_samples))
        }
      }
      # Create binary vector (1 = present, 0 = absent)
      resistant_omics <- rep(0, length(all_resistant_samples))
      names(resistant_omics) <- all_resistant_samples
      resistant_omics[all_resistant_samples %in% unlist(lapply(resistant_result$data[names(resistant_result$data) != "merged_dataset"], function(x) x$yes))] <- 1
    }
  }

  # Check if we have omics data
  if (is.null(sensitive_omics) || is.null(resistant_omics) ||
      length(sensitive_omics) == 0 || length(resistant_omics) == 0) {
    warning("No omics data found in sensitive_result or resistant_result")
    return(NULL)
  }

  # Prepare data for plotting
  plot_data <- data.frame(
    expression = c(resistant_omics, sensitive_omics),
    stratum = c(rep("Resistant", length(resistant_omics)),
                rep("Sensitive", length(sensitive_omics))),
    stringsAsFactors = FALSE
  )

  # Create plot title
  plot_title <- paste0(select_features, " Expression by ", stratification_drug, " Response")

  # Create y-axis label based on omics type
  if (feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")) {
    y_label <- "Expression Level"
  } else {
    y_label <- "Presence"
  }

  # Create boxplot mimicking FuncPairVisualization style
  p <- ggpubr::ggboxplot(data = plot_data,
                        x = "stratum",
                        y = "expression",
                        fill = "stratum",
                        palette = c("#BEBADAFF", "#FB8072FF"),
                        add = "jitter",
                        add.params = list(alpha = 0.15)) +
    ggpubr::stat_compare_means(size = 6,
                              label.x = 0.8,
                              label.y = (max(plot_data$expression) - max(plot_data$expression)/8),
                              label = "p.format") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      title = ggplot2::element_text(size = 15, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      legend.position = "none"
    ) +
    ggplot2::coord_cartesian(ylim = c(NA, max(plot_data$expression) + max(plot_data$expression)/20)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::ylab(y_label) +
    ggplot2::xlab("Response Group")

  return(p)
}

# Helper function for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x
