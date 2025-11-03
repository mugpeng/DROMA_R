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
    stat_cor(size = 6, method = method, label.x.npc = "left", label.y.npc = "top") +
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
    stat_compare_means(size = 6,
                       label.x.npc = "left",
                       label.y.npc = "top",
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

#' Create a volcano plot from meta-analysis results
#'
#' @param meta_df Data frame containing meta-analysis results with columns:
#'                effect_size, q_value (or p_value), n_datasets, and name
#' @param es_t Effect size threshold to consider significant
#' @param P_t Q-value threshold to consider significant
#' @param n_datasets_t Minimum number of datasets threshold (NULL for no threshold)
#' @param label Whether to add labels to top points (TRUE/FALSE)
#' @param top_label_each Number of top points in each direction to label
#' @param label_size Size of text labels
#' @param point_size Size of points
#' @param point_alpha Alpha transparency of points
#' @param title Plot title (NULL for no title)
#' @param custom_colors Custom color vector for Up, NS, Down (NULL for defaults)
#' @param use_p_value Logical, whether to use p_value column instead of q_value (default: FALSE)
#' @return ggplot object with volcano plot
#' @export
plotMetaVolcano <- function(meta_df,
                            es_t = .4,
                            P_t = .01,
                            n_datasets_t = NULL,
                            label = TRUE,
                            top_label_each = 5,
                            label_size = 5,
                            point_size = 2.5,
                            point_alpha = 0.6,
                            title = NULL,
                            custom_colors = NULL,
                            use_p_value = FALSE) {

  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  required_cols <- c("effect_size", "name", "n_datasets")
  p_val_col <- if(use_p_value) "p_value" else "q_value"
  if(!all(c(required_cols, p_val_col) %in% colnames(meta_df))) {
    stop("meta_df must contain columns: effect_size, ", p_val_col, ", n_datasets, and name")
  }

  # Default colors
  if(is.null(custom_colors)) {
    custom_colors <- c("Down" = "#44bce4", "NS" = "grey", "Up" = "#fc7474")
  }

  # Determine which p-value column to use
  p_val_col <- if(use_p_value) "p_value" else "q_value"

  # Group the points based on thresholds
  if(is.null(n_datasets_t)) {
    # No n_datasets threshold
    meta_df$group <- dplyr::case_when(
      meta_df$effect_size > es_t & meta_df[[p_val_col]] < P_t ~ "Up",
      meta_df$effect_size < -es_t & meta_df[[p_val_col]] < P_t ~ "Down",
      TRUE ~ "NS"
    )
  } else {
    # With n_datasets threshold
    meta_df$group <- dplyr::case_when(
      meta_df$effect_size > es_t & meta_df[[p_val_col]] < P_t & meta_df$n_datasets >= n_datasets_t ~ "Up",
      meta_df$effect_size < -es_t & meta_df[[p_val_col]] < P_t & meta_df$n_datasets >= n_datasets_t ~ "Down",
      TRUE ~ "NS"
    )
  }

  # Count significant findings
  sig_counts <- table(meta_df$group)
  sig_text <- paste0(
    "Up: ", sum(meta_df$group == "Up"), ", ",
    "Down: ", sum(meta_df$group == "Down"), ", ",
    "Total: ", nrow(meta_df)
  )

  # Get min and max for legend breaks
  min_datasets <- min(meta_df$n_datasets)
  max_datasets <- max(meta_df$n_datasets)
  
  # Basic volcano plot
  # Using scale_size makes the radius proportional to n_datasets
  # This means n_datasets=30 will have 3x the radius of n_datasets=10
  p <- ggplot(data = meta_df,
              aes(x = effect_size,
                  y = -log10(.data[[p_val_col]]))) +
    geom_point(alpha = point_alpha,
               aes(color = group, size = n_datasets)) +
    scale_size(range = c(point_size * 0.5, point_size * 2), 
               name = "N Datasets",
               breaks = function(x) {
                 # Generate nice breaks for the legend
                 if(max_datasets - min_datasets <= 5) {
                   return(seq(min_datasets, max_datasets, by = 1))
                 } else {
                   pretty_breaks <- pretty(c(min_datasets, max_datasets), n = 4)
                   return(pretty_breaks[pretty_breaks >= min_datasets & 
                                       pretty_breaks <= max_datasets])
                 }
               }) +
    theme_bw() +
    theme(
      title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, colour = "black"),
      axis.text = element_text(size = 15, color = "black"),
      legend.title = element_text(size = 13, colour = "black", face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      legend.box = "vertical",
      text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black")
    ) +
    ylab(if(use_p_value) "-log10(Pvalue)" else "-log10(Qvalue)") +
    xlab("Effect Size") +
    scale_color_manual(values = custom_colors, guide = "none") +
    geom_vline(xintercept = c(-es_t, es_t), lty = 4, col = "black", lwd = 0.5) +
    geom_hline(yintercept = -log10(P_t), lty = 4, col = "black", lwd = 0.5) +
    annotate("text", x = 0,
             y = max(-log10(meta_df[[p_val_col]]), na.rm = TRUE) * 0.9,
             label = sig_text, hjust = 0.5, size = 5, fontface = "bold")

  # Add title if provided
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }

  # Add labels if requested
  if(label) {
    meta_df2 <- meta_df[meta_df$group != "NS",]

    # Skip labeling if there are no significant points
    if(nrow(meta_df2) > 0) {
      # Get top points to label, ensuring exact counts
      ordered_indices <- order(meta_df2$effect_size)
      n_low <- min(top_label_each, sum(meta_df2$group == "Down"))
      n_high <- min(top_label_each, sum(meta_df2$group == "Up"))
      low_indices <- head(ordered_indices, n_low)
      high_indices <- tail(ordered_indices, n_high)
      forlabel_names <- unique(c(meta_df2$name[low_indices], meta_df2$name[high_indices]))
      forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]

      p <- p +
        geom_point(aes(size = n_datasets), shape = 1, data = forlabel_df, show.legend = FALSE) +
        ggrepel::geom_text_repel(
          data = forlabel_df,
          aes(label = name),
          size = label_size,
          color = "black",
          box.padding = 0.5,
          point.padding = 0.3,
          force = 5,
          max.overlaps = 20
        )
    }
  }

  p
}

#' Create a forest plot for meta-analysis results
#'
#' @description Creates a standardized forest plot for visualizing meta-analysis results
#' @param meta_obj Meta-analysis object from metagen() function
#' @param xlab Label for x-axis
#' @param show_common Logical, whether to show common effect model
#' @return A forest plot visualization
#' @export
createForestPlot <- function(meta_obj,
                             xlab = "Effect Size (95% CI)",
                             show_common = FALSE) {
  # Validate input
  if (!inherits(meta_obj, "meta") && !inherits(meta_obj, "metagen")) {
    stop("Input must be a meta-analysis object from the 'meta' package")
  }

  # Format p-value text for random effects model
  p_val <- meta_obj$pval.random
  p_text <- if(p_val < 0.001) {
    paste("Random-Effects Model (p =", format(p_val, scientific = TRUE, digits = 3), ")")
  } else {
    paste("Random-Effects Model (p =", round(p_val, 3), ")")
  }

  # Create forest plot
  meta::forest(meta_obj,
               xlab = xlab,
               slab = "study",
               print.pval.common = show_common,
               boxsize = 0.2,
               lineheight = "auto",
               print.pval.Q = FALSE,
               print.I2 = FALSE,
               print.tau2 = FALSE,
               common = show_common,
               text.random = p_text
  )
}

#' Create plot with common axis labels
#'
#' @description Creates a plot with common axis labels for multiple subplots from plotMultipleCorrelations or plotMultipleGroupComparisons
#' @param p A patchwork object from plotMultipleCorrelations or plotMultipleGroupComparisons
#' @param x_title Common x-axis title (optional, only for correlation plots)
#' @param y_title Common y-axis title
#' @return A patchwork/ggplot object with common axis labels
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
