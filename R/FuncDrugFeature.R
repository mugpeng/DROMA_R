# Data Processing Functions ----

#' Process Drug Sensitivity Data using DromaSet objects
#'
#' @description Creates a combined dataframe with raw and normalized drug sensitivity values from DromaSet or MultiDromaSet objects
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_drugs Character string specifying the drug name
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
processDrugData <- function(dromaset_object, select_drugs, data_type = "all", tumor_type = "all", overlap_only = FALSE) {
  if (select_drugs == "") {
    stop("Please select a drug.")
  }

  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }

  # Load drug data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet - Load normalized data
    drug_data <- loadTreatmentResponse(dromaset_object,
                                      select_drugs = select_drugs,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = TRUE)

    # Convert matrix to list format for compatibility
    if (is.matrix(drug_data) && select_drugs %in% rownames(drug_data)) {
      drug_vector <- as.numeric(drug_data[select_drugs, ])
      names(drug_vector) <- colnames(drug_data)
      drug_data_list <- list()
      drug_data_list[[dromaset_object@name]] <- drug_vector[!is.na(drug_vector)]
    } else {
      stop("Drug not found in treatment response data")
    }

    # Load raw data (non-normalized)
    raw_drug_data <- loadTreatmentResponse(dromaset_object,
                                          select_drugs = select_drugs,
                                          data_type = data_type,
                                          tumor_type = tumor_type,
                                          zscore = FALSE,
                                          return_data = TRUE)

    # Convert raw matrix to list format for compatibility
    if (is.matrix(raw_drug_data) && select_drugs %in% rownames(raw_drug_data)) {
      raw_drug_vector <- as.numeric(raw_drug_data[select_drugs, ])
      names(raw_drug_vector) <- colnames(raw_drug_data)
      raw_data_list <- list()
      raw_data_list[[dromaset_object@name]] <- raw_drug_vector[!is.na(raw_drug_vector)]
    } else {
      stop("Drug not found in raw treatment response data")
    }

  } else {
    # MultiDromaSet - Load normalized data
    drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                  select_drugs = select_drugs,
                                                  overlap_only = overlap_only,
                                                  data_type = data_type,
                                                  tumor_type = tumor_type,
                                                  zscore = TRUE)

    # Extract specific drug from each project
    drug_data_list <- lapply(drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
        names(drug_vector) <- colnames(drug_matrix)
        return(drug_vector[!is.na(drug_vector)])
      }
      return(NULL)
    })

    # Remove NULL entries
    drug_data_list <- drug_data_list[!sapply(drug_data_list, is.null)]

    # Load raw data (non-normalized)
    raw_drug_data <- loadMultiProjectTreatmentResponse(dromaset_object,
                                                      select_drugs = select_drugs,
                                                      overlap_only = overlap_only,
                                                      data_type = data_type,
                                                      tumor_type = tumor_type,
                                                      zscore = FALSE)

    # Extract specific drug from each project for raw data
    raw_data_list <- lapply(raw_drug_data, function(drug_matrix) {
      if (is.matrix(drug_matrix) && select_drugs %in% rownames(drug_matrix)) {
        drug_vector <- as.numeric(drug_matrix[select_drugs, ])
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
          SampleID = common_samples,
          zscore_value = as.numeric(drug_values[common_samples]),
          raw_value = as.numeric(raw_values[common_samples]),
          ProjectID = study,
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
                         annotations,
                         by = c("SampleID", "ProjectID"))

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

#' Get Drug Sensitivity Data using DromaSet objects
#'
#' @description Wrapper function that processes and returns drug sensitivity data from DromaSet or MultiDromaSet objects
#' @details This function combines the functionality of `processDrugData()` and `annotateDrugData()`.
#'   When `include_annotations = TRUE`, it will add sample annotations to the drug data.
#'   Sample annotations can be provided directly via the `sample_annotations` parameter,
#'   or loaded from the global `sample_anno` variable, or retrieved from the SQLite
#'   database specified by `db_path`.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_drugs Character string specifying the drug name
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
                                   select_drugs,
                                   data_type = "all",
                                   tumor_type = "all",
                                   overlap_only = FALSE,
                                   include_annotations = TRUE,
                                   sample_annotations = NULL,
                                   db_path = NULL) {
  # Process drug data
  drug_data <- processDrugData(dromaset_object, select_drugs, data_type, tumor_type, overlap_only)

  # Add annotations if requested
  if (include_annotations) {
    drug_data <- annotateDrugData(drug_data, sample_annotations, db_path)
  }

  return(drug_data)
}

# Visualization Functions ----
#' Plot continuous variable comparison
#'
#' @description Creates a scatter plot with correlation statistics for a continuous variable
#' @param data Data frame containing the variables to plot
#' @param cont_column Name of the continuous column to use for x-axis
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with the comparison plot
#' @export
plotContinuousComparison <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", title = NULL) {
  # Create scatter plot with correlation information
  # Determine y-axis label based on value_column
  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  # Determine plot title
  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "vs", cont_column)
  }

  p <- ggscatter(data, x = cont_column, y = value_column, alpha = 0.2) +
    stat_cor(size = 6, method = "spearman") +
    stat_smooth(formula = y ~ x, method = "lm") +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    ggtitle(plot_title) +
    xlab(cont_column) +
    ylab(y_label)

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
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with the grouped boxplot
#' @export
plotContinuousGroups <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", num_bins = 4, title = NULL) {
  # Remove rows with NA in either column
  data <- data[!is.na(data[[cont_column]]) & !is.na(data[[value_column]]), ]

  # Determine y-axis label based on value_column
  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  # Determine plot title
  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "by", cont_column, "Group")
  }

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
  data[[label_column]] <- factor(group_labels[data[[group_column]]], levels = group_labels)

  # Calculate overall median for reference line
  overall_median <- median(data[[value_column]], na.rm = TRUE)

  # Create boxplot with statistical test
  p <- ggboxplot(data, x = label_column, y = value_column,
                 fill = label_column, palette = soft_palette_26,
                 add = "jitter", add.params = list(alpha = 0.15)) +
    geom_hline(yintercept = overall_median, linetype = "dashed",
               color = "gray40", size = 0.8, alpha = 0.7) +
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(data[[value_column]]) - max(data[[value_column]])/8),
                       label = "p.format") +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggtitle(plot_title) +
    xlab(paste(cont_column, "Group")) +
    ylab(y_label)

  return(p)
}

#' Plot category comparison
#'
#' @description Creates boxplots comparing values across categories
#' @param data Data frame containing the variables to plot
#' @param category_column Name of the categorical column to use for grouping
#' @param value_column Name of the value column to use for y-axis (default: "value")
#' @param value_label Label for the value variable (default: "Drug Sensitivity")
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot2 object with the category comparison plot
#' @export
plotCategoryComparison <- function(data, category_column, value_column = "value", value_label = "Drug Sensitivity", title = NULL) {
  # Remove rows with NA in either column
  data <- data[!is.na(data[[category_column]]) & !is.na(data[[value_column]]), ]

  # Determine y-axis label based on value_column
  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  # Determine plot title
  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "by", category_column)
  }

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

  # Calculate median for each category and sort by median value
  category_medians <- tapply(data_filtered[[value_column]], 
                             data_filtered[[category_column]], 
                             median, na.rm = TRUE)
  sorted_categories <- names(sort(category_medians))
  
  # Special handling for TumorType: put "non-cancer" at the end if it exists
  if (category_column == "TumorType" && "non-cancer" %in% sorted_categories) {
    # Reorder categories to put "non-cancer" last
    sorted_categories <- c(sorted_categories[sorted_categories != "non-cancer"], "non-cancer")
  }
  
  data_filtered[[category_column]] <- factor(data_filtered[[category_column]], levels = sorted_categories)

  # Calculate overall median for reference line
  overall_median <- median(data_filtered[[value_column]], na.rm = TRUE)

  # Create improved boxplot with consistent styling
  p <- ggboxplot(data_filtered, x = category_column, y = value_column,
                 fill = category_column,
                 palette = soft_palette_26,
                 add = "jitter",
                 add.params = list(alpha = 0.15)) +
    geom_hline(yintercept = overall_median, linetype = "dashed",
               color = "gray40", size = 0.8, alpha = 0.7) +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(NA, max(data_filtered[[value_column]], na.rm = TRUE) + max(data_filtered[[value_column]], na.rm = TRUE)/20)) +
    ggtitle(plot_title) +
    xlab(NULL) +
    ylab(y_label)

  # Add statistical comparison with appropriate method
  if (length(valid_categories) >= 2) {
    max_y <- max(data_filtered[[value_column]], na.rm = TRUE)
    label_x_pos <- length(valid_categories) / 2.5
    if (length(valid_categories) == 2) {
      # For two groups, use wilcoxon with clear label positioning
      p <- p + stat_compare_means(size = 6,
                                  label.x.npc = "left",
                                  label.y.npc = "top",
                                  label = "p.format")
    } else {
      # For more than two groups, use Kruskal-Wallis
      # Add global p-value at top
      p <- p + stat_compare_means(method = "kruskal.test",
                                  size = 6,
                                  label.x.npc = "left",
                                  label.y.npc = "top",
                                  label = "p.format")
    }
  }

  return(p)
}


#' Create Drug Sensitivity Rank Plot
#'
#' @description Creates a rank plot showing drug sensitivity values ordered from most sensitive (highest values) to least sensitive
#' @details This function creates a rank plot where samples are ordered by drug sensitivity values.
#'   For AUC data, higher values indicate lower sensitivity, so samples are ranked with highest values on the left.
#'   The function supports highlighting specific samples and coloring by categorical variables.
#'   For MultiDromaSet objects with zscore=TRUE and merge=TRUE, it can combine data from multiple projects
#'   when the drug appears in at least two projects. When merge=FALSE, separate plots are returned for each project.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param select_drugs Character string specifying the drug name
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDC", "PDO", "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or specific tumor type
#' @param highlight Specifying samples to highlight. Can be:
#'   - A numeric value (e.g., 10) to highlight the top N samples by rank
#'   - A character value from data_type (e.g., "CellLine") to highlight all samples of that type
#'   - A character value from tumor_type (e.g., "breast cancer") to highlight all samples of that tumor type
#'   - A character vector of specific sample IDs to highlight
#'   Note: If more than 20 samples are highlighted, only the top 20 will be labeled
#' @param color Character string specifying the variable to use for coloring points.
#'   Options: NULL (default, no coloring), "data_type", "tumor_type", or any column name in sample annotations
#' @param zscore Logical, whether to use z-score normalized values (default: FALSE)
#' @param merge Logical, only applicable when zscore=TRUE and using MultiDromaSet.
#'   If TRUE, merges data from multiple projects when drug appears in at least 2 projects (default: FALSE)
#'   If FALSE, returns separate plots for each project
#' @param point_size Numeric, size of points in the plot (default: 2)
#' @param highlight_alpha Numeric, alpha transparency for non-highlighted points (default: 0.6)
#' @param sample_annotations Optional dataframe containing sample annotations
#' @param db_path Optional path to SQLite database for loading sample annotations
#' @return For DromaSet: A ggplot2 object. For MultiDromaSet with merge=FALSE: A list of ggplot2 objects (one per project).
#'   For MultiDromaSet with merge=TRUE: A single ggplot2 object with combined data.
#' @export
#' @examples
#' \dontrun{
#' # Basic rank plot
#' rank_plot <- plotDrugSensitivityRank(gCSI, "Paclitaxel")
#'
#' # Highlight top 10 most sensitive samples
#' rank_plot <- plotDrugSensitivityRank(gCSI, "Paclitaxel", highlight = 10)
#'
#' # Highlight cell line samples and color by tumor type
#' rank_plot <- plotDrugSensitivityRank(
#'   gCSI, "Paclitaxel",
#'   highlight = "CellLine",
#'   color = "TumorType"
#' )
#'
#' # MultiDromaSet with separate plots for each project
#' rank_plots <- plotDrugSensitivityRank(
#'   multi_set, "Paclitaxel",
#'   merge = FALSE,
#'   color = "TumorType"
#' )
#'
#' # MultiDromaSet with merged data
#' rank_plot <- plotDrugSensitivityRank(
#'   multi_set, "Paclitaxel",
#'   zscore = TRUE,
#'   merge = TRUE,
#'   color = "ProjectID"
#' )
#' }
plotDrugSensitivityRank <- function(dromaset_object,
                                   select_drugs,
                                   data_type = "all",
                                   tumor_type = "all",
                                   overlap_only = FALSE,
                                   highlight = NULL,
                                   color = NULL,
                                   zscore = FALSE,
                                   merge = FALSE,
                                   point_size = 2,
                                   highlight_alpha = 0.6,
                                   sample_annotations = NULL,
                                   db_path = NULL) {

  # Validate inputs
  if (missing(select_drugs) || is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
  }

  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("dromaset_object must be a DromaSet or MultiDromaSet object")
  }

  # For MultiDromaSet, check if merge is appropriate
  # if (inherits(dromaset_object, "MultiDromaSet") && merge && !zscore) {
  #   warning("Without z-score normalization (zscore=FALSE), merging data from different studies may not be appropriate. Consider setting zscore=TRUE. Setting merge=FALSE.")
  #   merge <- FALSE
  # }

  # Get drug sensitivity data
  drug_data <- getDrugSensitivityData(
    dromaset_object = dromaset_object,
    select_drugs = select_drugs,
    data_type = data_type,
    tumor_type = tumor_type,
    overlap_only = overlap_only,
    include_annotations = TRUE,
    sample_annotations = sample_annotations,
    db_path = db_path
  )

  if (is.null(drug_data) || nrow(drug_data) == 0) {
    stop("No drug sensitivity data found for the specified parameters")
  }

  # Determine which value column to use
  value_column <- if (zscore) "zscore_value" else "raw_value"
  if (!value_column %in% colnames(drug_data)) {
    stop("Required value column '", value_column, "' not found in data")
  }

  # Remove rows with missing values in the value column
  drug_data <- drug_data[!is.na(drug_data[[value_column]]), ]

  if (nrow(drug_data) == 0) {
    stop("No valid drug sensitivity values found after removing missing data")
  }

  # Handle MultiDromaSet cases
  if (inherits(dromaset_object, "MultiDromaSet")) {
    if (merge && zscore) {
      # Check if drug appears in multiple studies
      study_counts <- table(drug_data$ProjectID)
      if (length(study_counts) < 2) {
        warning("Drug '", select_drugs, "' found in only ", length(study_counts),
                " project(s). merge=TRUE requires at least 2 projects. Using individual project data.")
        merge <- FALSE
      } else {
        # Create a combined dataset indicator for merged analysis
        # drug_data$ProjectID <- "Combined"
        return(createSingleRankPlot(drug_data, select_drugs, highlight, color, zscore,
                                  value_column, point_size, highlight_alpha, merge))
      }
    }

    if (!merge) {
      # Create separate plots for each project
      projects <- unique(drug_data$ProjectID)
      plot_list <- list()

      for (project in projects) {
        project_data <- drug_data[drug_data$ProjectID == project, ]
        if (nrow(project_data) > 0) {
          plot_list[[project]] <- createSingleRankPlot(
            project_data, select_drugs, highlight, color, zscore,
            value_column, point_size, highlight_alpha, merge, project
          )
        }
      }

      return(plot_list)
    }
  }

  # For DromaSet or single project analysis
  return(createSingleRankPlot(drug_data, select_drugs, highlight, color, zscore,
                            value_column, point_size, highlight_alpha, merge))
}

#' Create a single rank plot (internal function)
#' @param drug_data Data frame with drug sensitivity data
#' @param select_drugs Name of the drug
#' @param highlight Highlighting specification
#' @param color Coloring specification
#' @param zscore Whether using zscore values
#' @param value_column Which value column to use
#' @param point_size Size of points
#' @param highlight_alpha Alpha for non-highlighted points
#' @param merge Whether this is a merged plot
#' @param projects Optional project name for title
#' @return ggplot2 object
createSingleRankPlot <- function(drug_data, select_drugs, highlight, color, zscore,
                               value_column, point_size, highlight_alpha, merge,
                               projects = NULL) {

  # Rank samples by drug sensitivity (highest values = most sensitive = rank 1)
  drug_data <- drug_data[order(drug_data[[value_column]], decreasing = TRUE), ]
  drug_data$rank <- 1:nrow(drug_data)

  # Determine highlighting
  highlight_samples <- character(0)
  is_numeric_highlight <- FALSE
  if (!is.null(highlight)) {
    if (length(highlight) == 1 && is.numeric(highlight)) {
      # If highlight is a number, highlight top N unique samples by rank
      is_numeric_highlight <- TRUE
      unique_samples_ordered <- unique(drug_data$SampleID)
      n_highlight <- min(as.integer(highlight), length(unique_samples_ordered))
      highlight_samples <- unique_samples_ordered[1:n_highlight]
    } else if (length(highlight) == 1) {
      # Check if highlight value matches data_type values
      if ("DataType" %in% colnames(drug_data) && highlight %in% drug_data$DataType) {
        highlight_samples <- drug_data$SampleID[drug_data$DataType == highlight]
      }
      # Check if highlight value matches tumor_type values
      else if ("TumorType" %in% colnames(drug_data) && highlight %in% drug_data$TumorType) {
        highlight_samples <- drug_data$SampleID[drug_data$TumorType == highlight]
      }
      # Check other annotation columns
      else {
        for (col in colnames(drug_data)) {
          if (col != "SampleID" && highlight %in% drug_data[[col]]) {
            highlight_samples <- drug_data$SampleID[drug_data[[col]] == highlight]
            break
          }
        }
      }
      # If no matches found, treat as specific sample IDs
      if (length(highlight_samples) == 0) {
        highlight_samples <- highlight[highlight %in% drug_data$SampleID]
      }
    } else {
      # Multiple values - treat as specific sample IDs
      highlight_samples <- highlight[highlight %in% drug_data$SampleID]
    }
  }

  # Add highlighting indicator
  drug_data$highlighted <- drug_data$SampleID %in% highlight_samples

  # Determine coloring variable
  color_data <- NULL
  color_title <- ""
  if (!is.null(color)) {
    if (color == "data_type" && "DataType" %in% colnames(drug_data)) {
      color_data <- drug_data$DataType
      color_title <- "Data Type"
    } else if (color == "tumor_type" && "TumorType" %in% colnames(drug_data)) {
      color_data <- drug_data$TumorType
      color_title <- "Tumor Type"
    } else if (color %in% colnames(drug_data)) {
      color_data <- drug_data[[color]]
      color_title <- gsub("_", " ", tools::toTitleCase(color))
    } else {
      warning("Color variable '", color, "' not found in data. Proceeding without coloring.")
    }
  }

  # Add points with conditional aesthetics
  if (!is.null(color_data)) {
    drug_data$color_var <- color_data
    p <- ggplot(drug_data, aes(x = rank, y = .data[[value_column]]))
    p <- p + geom_point(aes(color = color_var),
                       size = point_size,
                       alpha = ifelse(drug_data$highlighted, 1, highlight_alpha))
  } else {
    p <- ggplot(drug_data, aes(x = rank, y = .data[[value_column]]))
    p <- p + geom_point(size = point_size,
                       alpha = ifelse(drug_data$highlighted, 1, highlight_alpha),
                       color = "steelblue")
  }

  # Add highlighting circles and labels for highlighted samples
  if (length(highlight_samples) > 0) {
    highlighted_data <- drug_data[drug_data$highlighted, ]
    
    # For numeric highlight, work with unique samples; otherwise use all highlighted data
    if (is_numeric_highlight) {
      n_unique_highlighted <- length(unique(highlighted_data$SampleID))
      
      # Limit labeling to top 20 if too many highlighted samples
      if (n_unique_highlighted > 20) {
        # Sort by rank and take top 20 unique samples
        highlighted_data <- highlighted_data[order(highlighted_data$rank), ]
        unique_samples_to_label <- unique(highlighted_data$SampleID)[1:20]
        highlighted_data_to_label <- highlighted_data[highlighted_data$SampleID %in% unique_samples_to_label, ]
        highlighted_data_to_label <- highlighted_data_to_label[!duplicated(highlighted_data_to_label$SampleID), ]
        
        # Warn user about limitation
        message("Note: ", n_unique_highlighted, " unique samples were highlighted, but only the top 20 will be labeled to avoid overcrowding.")
      } else {
        # For numeric highlight, only label first occurrence of each sample
        highlighted_data_to_label <- highlighted_data[!duplicated(highlighted_data$SampleID), ]
      }
      highlighted_data_for_circle <- highlighted_data[!duplicated(highlighted_data$SampleID), ]
    } else {
      # For non-numeric highlight, keep all data (including duplicates in merged mode)
      if (nrow(highlighted_data) > 20) {
        # Sort by rank and take top 20
        highlighted_data <- highlighted_data[order(highlighted_data$rank), ]
        highlighted_data_to_label <- highlighted_data[1:20, ]
        
        # Warn user about limitation
        message("Note: ", nrow(highlighted_data), " samples were highlighted, but only the top 20 will be labeled to avoid overcrowding.")
      } else {
        highlighted_data_to_label <- highlighted_data
      }
      highlighted_data_for_circle <- highlighted_data
    }
    
    # Create labels with project info if merged
    if (merge && "ProjectID" %in% colnames(highlighted_data_to_label)) {
      highlighted_data_to_label$label <- paste0(highlighted_data_to_label$SampleID, " (", highlighted_data_to_label$ProjectID, ")")
    } else {
      highlighted_data_to_label$label <- highlighted_data_to_label$SampleID
    }

    # Add larger circles around highlighted points
    p <- p + geom_point(data = highlighted_data_for_circle,
                       aes(x = rank, y = .data[[value_column]]),
                       size = point_size + 0.5,
                       shape = 1,
                       color = "black",
                       stroke = 1.2)

    # Add labels for highlighted samples (limited to top 20)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = highlighted_data_to_label,
        aes(x = rank, y = .data[[value_column]], label = label),
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        max.overlaps = Inf,
        min.segment.length = 0
      )
    } else {
      warning("ggrepel package not available. Labels will not be added to highlighted points.")
    }
  }

  # Customize plot appearance
  value_label <- if (zscore) "Z-score Drug Sensitivity(Area Above Curve)" else "Drug Sensitivity(Area Above Curve)"
  plot_title <- paste("Drug Sensitivity Rank Plot:", select_drugs)
  if (!is.null(projects)) {
    plot_title <- paste(plot_title, "-", projects)
  } else if (merge) {
    plot_title <- paste(plot_title, "(Merged)")
  }

  # Create labs - only include color if it's used
  lab_list <- list(
    title = plot_title,
    x = "Rank (1 = Most Sensitive)",
    y = value_label
  )
  
  if (!is.null(color_data) && color_title != "") {
    lab_list$color <- color_title
  }
  
  p <- p +
    do.call(labs, lab_list) +
    theme_bw() +
    theme(
      title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  # Add color scale if coloring is used - using colors from theme_utils.R
  if (!is.null(color_data)) {
    if (is.factor(color_data) || is.character(color_data)) {
      # For categorical variables - use soft_palette_26 from theme_utils.R
      n_colors <- length(unique(color_data))
      if (n_colors <= length(soft_palette_26)) {
        p <- p + scale_color_manual(values = soft_palette_26[1:n_colors])
      } else {
        # Use colorRampPalette to extend the palette if needed
        extended_colors <- colorRampPalette(soft_palette_26)(n_colors)
        p <- p + scale_color_manual(values = extended_colors)
      }
    } else {
      # For continuous variables - use soft colors from soft_palette_26
      p <- p + scale_color_gradient(low = "#80B1D3", high = "#FB8072")  # Soft blue to soft coral
    }
  }

  # Add summary statistics as subtitle
  n_samples <- nrow(drug_data)
  n_highlighted <- sum(drug_data$highlighted)
  median_value <- median(drug_data[[value_column]], na.rm = TRUE)

  subtitle_text <- paste0(
    "N = ", n_samples, " samples",
    if (n_highlighted > 0) paste0(", ", n_highlighted, " highlighted"),
    ", Median = ", round(median_value, 3)
  )

  p <- p + labs(subtitle = subtitle_text)

  return(p)
}