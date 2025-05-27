# Data Processing Functions ----

#' Process Drug Sensitivity Data using DromaSet objects
#'
#' @description Creates a combined dataframe with raw and normalized drug sensitivity values from DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE).
#'   TRUE: Use only sample types present in all projects (recommended for meta-analysis)
#'   FALSE: Use all available samples from each project (may increase power but introduce bias)
#' @return A dataframe with combined raw and normalized drug sensitivity values
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' drug_data <- processDrugData(gCSI, "Paclitaxel")
#'
#' # Using MultiDromaSet with overlapping samples (recommended for meta-analysis)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))
#' drug_data <- processDrugData(multi_set, "Paclitaxel", overlap_only = FALSE)
#'
#' # Using MultiDromaSet with all samples (maximum sample size)
#' drug_data_all <- processDrugData(multi_set, "Paclitaxel", overlap_only = FALSE)
#' }
processDrugData <- function(dromaset_object, drug_name, data_type = "all", tumor_type = "all", overlap_only = FALSE) {
  if (drug_name == "") {
    stop("Please select a drug.")
  }

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }

  # Load drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet - Load normalized data
    drug_data <- loadTreatmentResponseNormalized(dromaset_object,
                                      drugs = drug_name,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE)

    # Convert matrix to list format for compatibility
    if (is.matrix(drug_data) && drug_name %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[drug_name, ])
      names(drug_vector) <- colnames(drug_data)
      drug_data_list <- list()
      drug_data_list[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    } else {
      stop("Drug not found in treatment response data")
    }

    # Load raw data (non-normalized)
    raw_drug_data <- loadTreatmentResponseNormalized(dromaset_object,
                                          drugs = drug_name,
                                          data_type = data_type,
                                          tumor_type = tumor_type,
                                          zscore = FALSE,
                                          return_data = TRUE)

    # Convert raw matrix to list format for compatibility
    if (is.matrix(raw_drug_data) && drug_name %in% rownames(raw_drug_data)) {
      raw_drug_vector <- as.numeric(raw_drug_data[drug_name, ])
      names(raw_drug_vector) <- colnames(raw_drug_data)
      raw_data_list <- list()
      raw_data_list[[dromaset_object@name]] <- raw_drug_vector[!is.na(raw_drug_vector)]
    } else {
      stop("Drug not found in raw treatment response data")
    }

  } else {
    # MultiDromaSet - Load normalized data
    drug_data <- loadMultiProjectTreatmentResponseNormalized(dromaset_object,
                                                  drugs = drug_name,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type)

    # Extract specific drug from each project
    drug_data_list <- lapply(drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && drug_name %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[drug_name, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    drug_data_list <- drug_data_list[!sapply(drug_data_list, is.null)]

    # Load raw data (non-normalized)
    raw_drug_data <- loadMultiProjectTreatmentResponseNormalized(dromaset_object,
                                                      drugs = drug_name,
                                                      overlap_only = overlap_only,
                                                      data_type = data_type,
                                                      tumor_type = tumor_type,
                                                      zscore = FALSE)

    # Extract specific drug from each project for raw data
    raw_data_list <- lapply(raw_drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && drug_name %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[drug_name, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    raw_data_list <- raw_data_list[!sapply(raw_data_list, is.null)]
  }

  # Get all study names (using union of keys from both lists)
  all_studies <- union(names(drug_data_list), names(raw_data_list))

  # Create an empty list to hold combined data for each study
  combined_data <- list()

  # For each study, process and combine the data
  for (study in all_studies) {
    if (study %in% names(drug_data_list) && study %in% names(raw_data_list)) {
      # Get values for this study
      drug_values <- drug_data_list[[study]]
      raw_values <- raw_data_list[[study]]

      # Get common sample IDs
      common_samples <- intersect(names(drug_values), names(raw_values))

      if (length(common_samples) > 0) {
        # Create dataframe with both values
        study_df <- data.frame(
          sampleid = common_samples,
          zscore_value = as.numeric(drug_values[common_samples]),
          raw_value = as.numeric(raw_values[common_samples]),
          study = study,
          stringsAsFactors = FALSE
        )

        combined_data[[study]] <- study_df
      }
    }
  }

  # Combine all studies into a single dataframe
  if (length(combined_data) > 0) {
    result_df <- bind_rows(combined_data)
    result_df <- na.omit(result_df)

    return(result_df)
  } else {
    return(NULL)
  }
}

#' Annotate Drug Sensitivity Data
#'
#' @description Adds sample annotations to drug sensitivity data
#' @details This function merges drug sensitivity data with sample annotations.
#'   It can obtain sample annotations from three sources, in order of preference:
#'   1. From the provided `sample_annotations` parameter
#'   2. From the global environment variable `sample_anno`
#'   3. From the SQLite database specified by `db_path`
#'
#' @param drug_data Dataframe containing drug sensitivity data from processDrugData
#' @param sample_annotations Dataframe containing sample annotations (default: uses global sample_anno)
#' @param db_path Optional path to SQLite database for loading sample annotations if not provided and not in global environment
#' @return A dataframe with drug sensitivity data and sample annotations
#' @export
#' @examples
#' \dontrun{
#' # Using existing sample_anno in global environment
#' drug_data <- processDrugData(gCSI, "Paclitaxel")
#' annotated_data <- annotateDrugData(drug_data)
#'
#' # Using provided sample annotations
#' my_annotations <- data.frame(SampleID = c("sample1", "sample2"),
#'                             TumorType = c("BRCA", "LUAD"))
#' annotated_data <- annotateDrugData(drug_data, sample_annotations = my_annotations)
#'
#' # Loading directly from database
#' annotated_data <- annotateDrugData(drug_data, db_path = "path/to/droma.sqlite")
#' }
annotateDrugData <- function(drug_data, sample_annotations = NULL, db_path = NULL) {
  if (is.null(drug_data)) {
    return(NULL)
  }

  # Use provided sample annotations or global sample_anno
  annotations <- sample_annotations
  if (is.null(annotations)) {
    if (!is.null(db_path) && file.exists(db_path)) {
      # Try to load sample annotations from the database
      tryCatch({
        # Load required package
        if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) {
          warning("DBI and/or RSQLite packages not available. Cannot load sample annotations from database.")
          return(drug_data)
        }

        # Connect to database
        con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
        on.exit(DBI::dbDisconnect(con), add = TRUE)

        # Query sample annotations
        annotations <- DBI::dbGetQuery(con, "SELECT * FROM sample_anno")
      }, error = function(e) {
        warning(paste("Failed to load sample annotations from database:", e$message,
                      "\nReturning non-annotated data."))
        return(drug_data)
      })
    } else {
      warning("sample_anno object not found and no valid db_path provided. Returning non-annotated data.")
      return(drug_data)
    }
  }

  # Merge drug data with annotations
  merged_df <- left_join(drug_data,
                         annotations %>% select(-ProjectID), # Remove ProjectID from the join
                         by = c("sampleid" = "SampleID"))

  # Clean up
  if ("ProjectRawName" %in% colnames(merged_df)) {
    merged_df$ProjectRawName <- NULL
  }

  merged_df <- unique(merged_df)
  return(merged_df)
}

#' Format Drug Data Table
#'
#' @description Creates a formatted datatable for drug sensitivity data
#' @param drug_data Dataframe containing drug sensitivity data
#' @param caption Caption text for the table
#' @return A formatted DT::datatable object
#' @export
formatDrugTable <- function(drug_data, caption = "Drug sensitivity data - showing both raw and Z-score normalized values") {
  if (is.null(drug_data) || nrow(drug_data) == 0) {
    return(NULL)
  }

  DT::datatable(
    drug_data,
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left; color: black; font-size: 14px;',
      htmltools::strong(caption)
    ),
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv')
    ),
    extensions = 'Buttons',
    rownames = FALSE,
    filter = 'top'
  ) %>%
  # Format numeric columns to 3 decimal places
  DT::formatRound(columns = c('zscore_value', 'raw_value'), digits = 3)
}

# Visualization Functions ----
#' Get Drug Sensitivity Data using DromaSet objects
#'
#' @description Wrapper function that processes and returns drug sensitivity data from DromaSet or MultiDromaSet objects
#' @details This function combines the functionality of `processDrugData()` and `annotateDrugData()`.
#'   When `include_annotations = TRUE`, it will add sample annotations to the drug data.
#'   Sample annotations can be provided directly via the `sample_annotations` parameter,
#'   or loaded from the global `sample_anno` variable, or retrieved from the SQLite
#'   database specified by `db_path`.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE).
#'   TRUE: Use only sample types present in all projects (recommended for meta-analysis)
#'   FALSE: Use all available samples from each project (may increase power but introduce bias)
#' @param include_annotations Logical indicating whether to include sample annotations
#' @param sample_annotations Optional dataframe containing sample annotations
#' @param db_path Optional path to SQLite database for loading sample annotations if not provided and not in global environment
#' @return A dataframe with drug sensitivity data
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel")
#'
#' # Using MultiDromaSet with overlapping samples (recommended)
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"))
#' drug_data <- getDrugSensitivityData(multi_set, "Paclitaxel",
#'                                    include_annotations = TRUE)
#'
#' # Using MultiDromaSet with all available samples
#' drug_data_all <- getDrugSensitivityData(multi_set, "Paclitaxel",
#'                                        overlap_only = FALSE,
#'                                        include_annotations = TRUE)
#'
#' # Using database path to load sample annotations
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel",
#'                                    include_annotations = TRUE,
#'                                    db_path = "path/to/droma.sqlite")
#'
#' # Using custom sample annotations
#' my_annotations <- data.frame(SampleID = c("sample1", "sample2"),
#'                             TumorType = c("BRCA", "LUAD"))
#' drug_data <- getDrugSensitivityData(gCSI, "Paclitaxel",
#'                                    include_annotations = TRUE,
#'                                    sample_annotations = my_annotations)
#' }
getDrugSensitivityData <- function(dromaset_object,
                                   drug_name,
                                   data_type = "all",
                                   tumor_type = "all",
                                   overlap_only = FALSE,
                                   include_annotations = TRUE,
                                   sample_annotations = NULL,
                                   db_path = NULL) {
  # Process drug data
  drug_data <- processDrugData(dromaset_object, drug_name, data_type, tumor_type, overlap_only)

  # Add annotations if requested
  if (include_annotations) {
    drug_data <- annotateDrugData(drug_data, sample_annotations, db_path)
  }

  return(drug_data)
}

#' Plot continuous variable comparison
#'
#' @description Creates a scatter plot with correlation statistics for a continuous variable
#' @param data Data frame containing the variables to plot
#' @param cont_column Name of the continuous column to use for x-axis
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @return A ggplot2 object with the comparison plot
#' @export
plotContinuousComparison <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity") {
  # Create scatter plot with correlation information
  p <- ggscatter(data, x = cont_column, y = value_column, alpha = 0.2) +
    stat_cor(size = 6, method = "spearman") +
    stat_smooth(formula = y ~ x, method = "lm") +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12)
    ) +
    ggtitle(paste(value_label, "vs", cont_column))


  return(p)
}

#' Plot continuous variable as groups
#'
#' @description Creates boxplots for binned groups of a continuous variable
#' @param data Data frame containing the variables to plot
#' @param cont_column Name of the continuous column to bin into groups
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param num_bins Number of bins to create from the continuous variable (default: 4)
#' @return A ggplot2 object with the grouped boxplot
#' @export
plotContinuousGroups <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", num_bins = 4) {
  # Create bins for the continuous variable
  cont_values <- data[[cont_column]]

  # Create bins
  cont_bins <- cut(cont_values,
                   breaks = num_bins,
                   include.lowest = TRUE,
                   labels = FALSE)

  # Create labels for groups
  cont_range <- range(cont_values, na.rm = TRUE)
  bin_width <- diff(cont_range) / num_bins
  group_labels <- sapply(1:num_bins, function(i) {
    lower <- cont_range[1] + (i-1) * bin_width
    upper <- cont_range[1] + i * bin_width
    paste0(round(lower), "-", round(upper))
  })

  # Add group information to data
  group_column <- paste0(cont_column, "_group")
  label_column <- paste0(cont_column, "_group_label")

  data[[group_column]] <- cont_bins
  data[[label_column]] <- group_labels[data[[group_column]]]

  # Create boxplot with statistical test
  p <- ggboxplot(data, x = label_column, y = value_column,
                 fill = label_column, palette = "jco",
                 add = "jitter", add.params = list(alpha = 0.15)) +
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(data[[value_column]]) - max(data[[value_column]])/8),
                       label = "p.format") +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggtitle(paste(value_label, "by", cont_column, "Group"))


  return(p)
}

#' Plot category comparison
#'
#' @description Creates boxplots comparing values across categories
#' @param data Data frame containing the variables to plot
#' @param category_column Name of the categorical column to use for grouping
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @return A ggplot2 object with the category comparison plot
#' @export
plotCategoryComparison <- function(data, category_column, value_column = "value", value_label = "Drug Sensitivity") {
  # Count observations per category and filter out categories with too few samples
  category_counts <- table(data[[category_column]])
  valid_categories <- names(category_counts)[category_counts >= 3]

  if (length(valid_categories) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Not enough samples per category for comparison") +
             theme_void())
  }

  data_filtered <- data[data[[category_column]] %in% valid_categories, ]

  # Create improved boxplot with consistent styling
  p <- ggboxplot(data_filtered, x = category_column, y = value_column,
                 fill = category_column,
                 palette = bright_palette_26,
                 add = "jitter",
                 add.params = list(alpha = 0.15)) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggtitle(paste(value_label, "by", category_column))


  # Add statistical comparison with appropriate method
  if (length(valid_categories) >= 2) {
    max_y <- max(data_filtered[[value_column]], na.rm = TRUE)
    label_x_pos <- length(valid_categories) / 2.5
    if (length(valid_categories) == 2) {
      # For two groups, use wilcoxon with clear label positioning
      p <- p + stat_compare_means(size = 6,
                                  label.x = label_x_pos,
                                  label.y = (max_y - max_y/8),
                                  label = "p.format")
    } else {
      # For more than two groups, use Kruskal-Wallis
      # Add global p-value at top
      p <- p + stat_compare_means(method = "kruskal.test",
                                  size = 6,
                                  label.x = label_x_pos,
                                  label.y = max_y + (max_y * 0.15),
                                  label = "p.format")
    }
  }

  return(p)
}

#' Create drug comparison plot
#'
#' @description Creates appropriate comparison plots based on variable type
#' @param data Data frame containing the variables to plot
#' @param comparison_var Name of the variable to compare against drug sensitivity
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param num_bins Number of bins to create from continuous variables (default: 4)
#' @param show_groups_boxplot Logical, whether to show grouped boxplot for continuous variables
#' @return A ggplot2 object or grid with comparison plots
#' @export
createDrugComparisonPlot <- function(data, comparison_var, value_column = "value", value_label = "Drug Sensitivity",
                                     num_bins = 4, show_groups_boxplot = TRUE) {
  # Handle missing values in the comparison variable
  data <- data[!is.na(data[[comparison_var]]), ]

  if (nrow(data) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No data available for this comparison") +
             theme_void())
  }

  # Check if the comparison variable is numeric/continuous
  if (is.numeric(data[[comparison_var]])) {
    # For continuous variables
    p1 <- plotContinuousComparison(data, cont_column = comparison_var,
                                   value_column = value_column, value_label = value_label)

    # Also create a boxplot with bins if requested
    if (show_groups_boxplot) {
      # Create grouped boxplot
      p2 <- plotContinuousGroups(data, cont_column = comparison_var,
                                 value_column = value_column, value_label = value_label,
                                 num_bins = num_bins)

      # Return grid of both plots
      return(patchwork::wrap_plots(p1,p2, ncol = 2))
    }

    return(p1)
  } else {
    # For categorical variables
    return(plotCategoryComparison(data, category_column = comparison_var,
                                  value_column = value_column, value_label = value_label))
  }
}
