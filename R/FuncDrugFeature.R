bright_palette_26 <- c(
  "#FF5733", # Bright red-orange
  "#33C4FF", # Bright sky blue
  "#FF33E9", # Bright pink
  "#33FF57", # Bright green
  "#FFD133", # Bright yellow
  "#8B33FF", # Bright purple
  "#FF8B33", # Bright orange
  "#33FFC4", # Bright teal
  "#FF3355", # Bright red
  "#33FFED", # Bright cyan
  "#C433FF", # Bright violet
  "#BBFF33", # Bright lime
  "#FF33BB", # Bright magenta
  "#33FFB2", # Bright mint
  "#3357FF", # Bright blue
  "#FF9933", # Bright amber
  "#33FF78", # Bright seafoam
  "#ED33FF", # Bright fuchsia
  "#57FF33", # Bright chartreuse
  "#FF33C4", # Bright hot pink
  "#33A5FF", # Bright azure
  "#FFB233", # Bright gold
  "#3378FF", # Bright royal blue
  "#FF5733", # Bright vermilion
  "#33FF9E", # Bright turquoise
  "#D6FF33"  # Bright yellow-green
)

# process function ----
#' Process Drug Sensitivity Data
#'
#' Creates a combined dataframe with raw and normalized drug sensitivity values
#'
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @return A dataframe with combined raw and normalized drug sensitivity values
#' @export
process_drug_data <- function(drug_name, data_type = "all", tumor_type = "all") {
  if (drug_name == "") {
    stop("Please select a drug.")
  }
  
  # Get z-score normalized data
  normalized_data_list <- selFeatures(
    select_feas_type = "drug", 
    select_feas = drug_name,
    data_type = data_type, 
    tumor_type = tumor_type
  )
  
  # Get raw data
  raw_data_list <- selFeatures(
    select_feas_type = "drug_raw", 
    select_feas = drug_name,
    data_type = data_type, 
    tumor_type = tumor_type
  )
  
  # Get all study names (using union of keys from both lists)
  all_studies <- union(names(normalized_data_list), names(raw_data_list))
  
  # Create an empty list to hold combined data for each study
  combined_data <- list()
  
  # For each study, process and combine the data
  for (study in all_studies) {
    if (study %in% names(normalized_data_list) && study %in% names(raw_data_list)) {
      # Get values for this study
      normalized_values <- normalized_data_list[[study]]
      raw_values <- raw_data_list[[study]]
      
      # Get common sample IDs
      common_samples <- intersect(names(normalized_values), names(raw_values))
      
      if (length(common_samples) > 0) {
        # Create dataframe with both values
        study_df <- data.frame(
          sampleid = common_samples,
          zscore_value = as.numeric(normalized_values[common_samples]),
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
#' Adds sample annotations to drug sensitivity data
#'
#' @param drug_data Dataframe containing drug sensitivity data from process_drug_data
#' @param sample_annotations Dataframe containing sample annotations (default: uses global sample_anno)
#' @return A dataframe with drug sensitivity data and sample annotations
#' @export
annotate_drug_data <- function(drug_data, sample_annotations = NULL) {
  if (is.null(drug_data)) {
    return(NULL)
  }
  
  # Use provided sample annotations or global sample_anno
  annotations <- sample_annotations
  if (is.null(annotations)) {
    if (exists("sample_anno", envir = .GlobalEnv)) {
      annotations <- base::get("sample_anno", envir = .GlobalEnv)
    } else {
      warning("sample_anno object not found. Returning non-annotated data.")
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
#' Creates a formatted datatable for drug sensitivity data
#'
#' @param drug_data Dataframe containing drug sensitivity data
#' @param caption Caption text for the table
#' @return A formatted DT::datatable object
#' @export
format_drug_table <- function(drug_data, caption = "Drug sensitivity data - showing both raw and Z-score normalized values") {
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

#' Get Drug Sensitivity Data
#'
#' Wrapper function that processes and returns drug sensitivity data
#'
#' @param drug_name Character string specifying the drug name
#' @param data_type Filter by data type ("all", "CellLine", "PDC", "PDO", "PDX")
#' @param tumor_type Filter by tumor type (use "all" for all tumor types)
#' @param include_annotations Logical indicating whether to include sample annotations
#' @param sample_annotations Optional dataframe containing sample annotations
#' @return A dataframe with drug sensitivity data
#' @export
get_drug_sensitivity_data <- function(drug_name, 
                                     data_type = "all", 
                                     tumor_type = "all",
                                     include_annotations = TRUE,
                                     sample_annotations = NULL) {
  # Process drug data
  drug_data <- process_drug_data(drug_name, data_type, tumor_type)
  
  # Add annotations if requested
  if (include_annotations) {
    drug_data <- annotate_drug_data(drug_data, sample_annotations)
  }
  
  return(drug_data)
}


# plot ----
# Create a comparison plot for continuous variables (like Age)
plot_continuous_comparison <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity") {
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

# Create boxplots for continuous variable groups
plot_continuous_groups <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", num_bins = 4) {
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

# Create a categorical comparison plot
plot_category_comparison <- function(data, category_column, value_column = "value", value_label = "Drug Sensitivity") {
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

# Create a drug feature comparison plot
create_drug_comparison_plot <- function(data, comparison_var, value_column = "value", value_label = "Drug Sensitivity", 
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
    p1 <- plot_continuous_comparison(data, cont_column = comparison_var, 
                                     value_column = value_column, value_label = value_label)
    
    # Also create a boxplot with bins if requested
    if (show_groups_boxplot) {
      # Create grouped boxplot
      p2 <- plot_continuous_groups(data, cont_column = comparison_var, 
                                   value_column = value_column, value_label = value_label,
                                   num_bins = num_bins)
      
      # Return grid of both plots
      return(grid.arrange(p1, p2, ncol = 2))
    }
    
    return(p1)
  } else {
    # For categorical variables
    return(plot_category_comparison(data, category_column = comparison_var, 
                                    value_column = value_column, value_label = value_label))
  }
}
