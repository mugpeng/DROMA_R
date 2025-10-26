# Visualization Functions Module ----

#' Plot correlation between two continuous variables
#'
#' @description Creates a scatter plot with correlation statistics
#' @param x_values Vector of x-axis values
#' @param y_values Vector of y-axis values
#' @param x_label Label for x-axis
#' @param y_label Label for y-axis
#' @param title Plot title
#' @param method Correlation method ("spearman", "pearson")
#' @return A ggplot2 object with scatter plot and correlation statistics
#' @export
plotCorrelation <- function(x_values, y_values,
                           x_label = "Feature 1",
                           y_label = "Feature 2",
                           title = "Correlation Plot",
                           method = "spearman") {
  # Combine data into dataframe
  cor_df <- data.frame(
    x = x_values,
    y = y_values
  )

  # Create scatter plot with correlation statistics
  p <- ggscatter(cor_df, x = "x", y = "y", alpha = 0.2) +
    stat_cor(size = 6, method = method) +
    stat_smooth(formula = y ~ x, method = "lm") +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    ggtitle(title) +
    xlab(x_label) +
    ylab(y_label)

  return(p)
}

#' Plot continuous values by discrete groups
#'
#' @description Creates a boxplot comparing continuous values between discrete groups
#' @param no_values Values for group without feature present
#' @param yes_values Values for group with feature present
#' @param group_labels Labels for the two groups
#' @param title Plot title
#' @param y_label Label for y-axis
#' @return A ggplot2 object with boxplot and statistical test
#' @export
plotGroupComparison <- function(no_values, yes_values,
                                group_labels = c("Without", "With"),
                                title = "Group Comparison",
                                y_label = "Value") {
  # Combine data into dataframe
  box_df <- data.frame(
    values = c(no_values, yes_values),
    group = rep(group_labels, times = c(length(no_values), length(yes_values)))
  )

  # Create boxplot with statistical test
  p <- ggboxplot(data = box_df, x = "group", y = "values",
                 fill = "group", palette = c("#BEBADAFF", "#FB8072FF"),
                 add = "jitter", add.params = list(alpha = 0.15)) +
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(box_df$values) - max(box_df$values)/8),
                       label = "p.format") +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(NA, max(box_df$values) + max(box_df$values)/20)) +
    ggtitle(title) +
    ylab(y_label) + 
    xlab("")

  return(p)
}

#' Plot multiple correlation plots in a grid
#'
#' @description Creates and combines multiple correlation plots
#' @param pairs_list List of paired datasets with feature1 and feature2 values
#' @param x_label Common x-axis label
#' @param y_label Common y-axis label
#' @param ncol Number of columns in the grid
#' @return A combined plot with all correlations or NULL if no valid pairs
#' @export
plotMultipleCorrelations <- function(pairs_list,
                                    x_label = "Feature 1",
                                    y_label = "Feature 2",
                                    ncol = 3) {
  # Initialize list to store plots
  p_list <- list()

  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    # Skip merged dataset for individual plots
    if (names(pairs_list)[i] == "merged_dataset") next

    tryCatch({
      feat1_vals <- pairs_list[[i]]$feature1
      feat2_vals <- pairs_list[[i]]$feature2

      # Ensure adequate data for plotting
      if (length(feat1_vals) < 3 || length(feat2_vals) < 3) next

      # Create plot and add to list (without axis labels to avoid redundancy)
      p_list[[i]] <- plotCorrelation(feat1_vals, feat2_vals,
                                     title = names(pairs_list)[i],
                                     x_label = "",
                                     y_label = "")
    }, error = function(e) {
      warning("Error creating plot for pair ", names(pairs_list)[i], ": ", e$message)
    })
  }

  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]

  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (length(p_list) <= 3) {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = length(p_list))
      } else {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = ncol)
      }

      # Add overall title if we have multiple plots
      if (length(p_list) > 1) {
        combined_plot <- combined_plot +
          patchwork::plot_annotation(
            title = paste("Multiple Correlations:", x_label, "vs", y_label),
            theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
          )
      }

      return(combined_plot)
    } else {
      warning("patchwork package not available. Returning list of individual plots.")
      return(p_list)
    }
  } else {
    return(NULL)
  }
}

#' Plot multiple group comparisons in a grid
#'
#' @description Creates and combines multiple group comparison plots
#' @param pairs_list List of paired datasets with yes/no groups
#' @param x_label Label for the grouping variable (discrete feature)
#' @param y_label Common y-axis label
#' @param group_labels Labels for the two groups
#' @param ncol Number of columns in the grid
#' @return A combined plot with all comparisons or NULL if no valid pairs
#' @export
plotMultipleGroupComparisons <- function(pairs_list,
                                        x_label = "Feature",
                                        y_label = "Value",
                                        group_labels = c("Without", "With"),
                                        ncol = 3) {
  # Initialize list to store plots
  p_list <- list()

  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    # Skip merged dataset for individual plots
    if (names(pairs_list)[i] == "merged_dataset") next

    tryCatch({
      yes_vals <- pairs_list[[i]]$yes
      no_vals <- pairs_list[[i]]$no

      # Ensure adequate data for plotting
      if (length(yes_vals) < 3 || length(no_vals) < 3) next

      # Create plot and add to list (without axis labels to avoid redundancy)
      p_list[[i]] <- plotGroupComparison(no_vals, yes_vals,
                                         group_labels = group_labels,
                                         title = names(pairs_list)[i],
                                         y_label = "")
    }, error = function(e) {
      warning("Error creating plot for pair ", names(pairs_list)[i], ": ", e$message)
    })
  }

  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]

    # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      if (length(p_list) <= 3) {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = length(p_list))
      } else {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = ncol)
      }

      # Add overall title if we have multiple plots
      if (length(p_list) > 1) {
        combined_plot <- combined_plot +
          patchwork::plot_annotation(
            title = paste("Multiple Group Comparisons:", x_label, "vs", y_label),
            theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
          )
      }

      return(combined_plot)
    } else {
      warning("patchwork package not available. Returning list of individual plots.")
      return(p_list)
    }
  } else {
    return(NULL)
  }
}

#' Create plot with common axis labels
#'
#' @description Creates a plot with common axis labels for multiple subplots from plotMultipleCorrelations or plotMultipleGroupComparisons
#' @param p A patchwork object from plotMultipleCorrelations or plotMultipleGroupComparisons
#' @param x_title Common x-axis title (optional, only for correlation plots)
#' @param y_title Common y-axis title
#' @return A patchwork/ggplot object with common axis labels
#' @export
createPlotWithCommonAxes <- function(p, x_title = NULL, y_title = NULL) {
  # Check if plot object is valid
  if (is.null(p)) {
    warning("Input plot is NULL. Returning NULL.")
    return(NULL)
  }
  
  # If p is a list (when patchwork is not available), return it as is
  if (is.list(p) && !inherits(p, "ggplot") && !inherits(p, "patchwork")) {
    warning("Input is a list of plots, not a combined plot. Returning original list.")
    return(p)
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package not available. Returning original plot.")
    return(p)
  }
  
  # Build caption text for axis labels
  caption_text <- ""
  if (!is.null(x_title) && !is.null(y_title)) {
    caption_text <- paste0("X-axis: ", x_title, "   Y-axis: ", y_title)
  } else if (!is.null(y_title)) {
    caption_text <- paste0("Y-axis: ", y_title)
  } else if (!is.null(x_title)) {
    caption_text <- paste0("X-axis: ", x_title)
  }
  
  # Add caption as axis labels using plot_annotation
  if (caption_text != "") {
    p <- p + patchwork::plot_annotation(
      caption = caption_text,
      theme = theme(
        plot.caption = element_text(size = 14, hjust = 0.5, face = "bold")
      )
    )
  }
  
  return(p)
}