# Visualization Functions Module ----

# Correlation / Pair Comparison ----

.wrapVisualizationLabels <- function(x, width = 24) {
  vapply(x, function(label) {
    label <- as.character(label)
    if (is.na(label) || label == "") {
      return(label)
    }
    if (grepl("_", label, fixed = TRUE)) {
      parts <- strsplit(label, "_", fixed = TRUE)[[1]]
      wrapped <- paste(strwrap(paste(parts, collapse = " "), width = width), collapse = "\n")
      wrapped <- gsub(" ", "_", wrapped, fixed = TRUE)
      if (grepl("\n", wrapped, fixed = TRUE)) {
        wrapped <- paste0(wrapped, "\n")
      }
      return(wrapped)
    }
    paste(strwrap(label, width = width), collapse = "\n")
  }, character(1))
}

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
#' @param subtitle Optional subtitle (e.g. comma-separated dataset IDs for merged panels).
#' @param y_label Label for y-axis
#' @return A ggplot2 object with boxplot and statistical test
plotGroupComparison <- function(no_values, yes_values,
                                group_labels = c("Without", "With"),
                                title = "Group Comparison",
                                y_label = "Value",
                                subtitle = NULL) {
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
      plot.subtitle = element_text(size = 10, color = "gray35"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 15, hjust = 0.85),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(NA, max(box_df$values) + max(box_df$values)/20)) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ylab(y_label) +
    xlab("")

  return(p)
}


#' ROC / AUC plot for clinical drug response (expression vs response groups)
#'
#' @description
#' ggplot2 ROC curve comparable to the clinical validation panel style: diagonal
#' reference line, smooth curve, and annotated AUC.
#'
#' @param no_values Non-response expression values (same order as \code{plotGroupComparison}).
#' @param yes_values Response expression values.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param positive_direction One of \code{"auto"}, \code{"high"}, \code{"low"}; passed to the
#'   internal ROC builder in the clinical module (\code{.computeClinicalResponseRoc}).
#' @return A \code{ggplot} object, or \code{NULL} if ROC cannot be computed.
#'         Attribute \code{roc_auc} holds the numeric AUC when the plot is non-NULL.
#' @export
plotClinicalDrugResponseRoc <- function(no_values,
                                       yes_values,
                                       title = "Clinical readout",
                                       subtitle = "Response classification performance",
                                       positive_direction = c("auto", "high", "low")) {
  positive_direction <- match.arg(positive_direction)
  prep <- .computeClinicalResponseRoc(no_values, yes_values,
                                     positive_direction = positive_direction)
  if (is.null(prep$df)) {
    return(NULL)
  }
  auc_val <- prep$auc
  df <- prep$df
  # Panel D–style annotation position (lower-right of diagonal)
  txt_x <- 0.55
  txt_y <- 0.12
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$FPR, .data$TPR)) +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                        color = "#B0A99F", linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_line(color = "#4B74F2", linewidth = 1.2) +
    ggplot2::annotate("text", x = txt_x, y = txt_y,
                     label = sprintf("AUC = %.2f", auc_val),
                     hjust = 0, size = 3.8) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "False positive rate",
      y = "True positive rate"
    ) +
    ggplot2::coord_equal() +
    ggplot2::lims(x = c(0, 1), y = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#6F685E"),
      panel.grid.minor = ggplot2::element_blank()
    )
  attr(p, "roc_auc") <- auc_val
  p
}


#' Plot multiple correlation plots in a grid
#'
#' @description Creates and combines multiple correlation plots
#' @param pairs_list List of paired datasets with feature1 and feature2 values
#' @param x_label Common x-axis label
#' @param y_label Common y-axis label
#' @param ncol Number of columns in the grid
#' @param data_type_anno Optional character string to add as annotation in plot titles (e.g., "Cell lines"). If provided, it will be appended to titles as "(annotation)"
#' @return A combined plot with all correlations or NULL if no valid pairs
plotMultipleCorrelations <- function(pairs_list,
                                    x_label = "Feature 1",
                                    y_label = "Feature 2",
                                    ncol = 3,
                                    data_type_anno = NULL) {
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
        # Create title with optional data type annotation
        title_suffix <- if (!is.null(data_type_anno) && nchar(data_type_anno) > 0) {
          paste0(" (", data_type_anno, ")")
        } else {
          ""
        }
        combined_plot <- combined_plot +
          patchwork::plot_annotation(
            title = paste0("Multiple Correlations: ", x_label, " vs ", y_label, title_suffix),
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
#' @param data_type_anno Optional character string to add as annotation in plot titles (e.g., "Cell lines"). If provided, it will be appended to titles as "(annotation)"
#' @return A combined plot with all comparisons or NULL if no valid pairs
plotMultipleGroupComparisons <- function(pairs_list,
                                        x_label = "Feature",
                                        y_label = "Value",
                                        group_labels = c("Without", "With"),
                                        ncol = 3,
                                        data_type_anno = NULL) {
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
      # Determine ncol: use 2 columns when there are 4 plots, otherwise follow existing logic
      if (length(p_list) <= 3) {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = length(p_list))
      } else if (length(p_list) == 4) {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = 2)
      } else {
        combined_plot <- patchwork::wrap_plots(p_list, ncol = ncol)
      }

      # Add overall title if we have multiple plots
      if (length(p_list) > 1) {
        # Create title with optional data type annotation
        title_suffix <- if (!is.null(data_type_anno) && nchar(data_type_anno) > 0) {
          paste0(" (", data_type_anno, ")")
        } else {
          ""
        }
        combined_plot <- combined_plot +
          patchwork::plot_annotation(
            title = paste0("Multiple Group Comparisons: ", x_label, " vs ", y_label, title_suffix),
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
#' @param size_col Column name to use for point size (e.g., "n_datasets"). NULL for fixed size.
#' @param label Whether to add labels to top points (TRUE/FALSE)
#' @param top_label_each Number of top points in each direction to label
#' @param custom_labels Character vector of gene names to label (overrides top_label_each if provided)
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
                            size_col = "n_datasets",
                            label = TRUE,
                            top_label_each = 5,
                            custom_labels = NULL,
                            label_size = 5,
                            point_size = 2.5,
                            point_alpha = 0.6,
                            title = NULL,
                            custom_colors = NULL,
                            use_p_value = FALSE) {

  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  required_cols <- c("effect_size", "name")
  p_val_col <- if(use_p_value) "p_value" else "q_value"
  required_cols <- c(required_cols, p_val_col)
  if(!is.null(n_datasets_t)) required_cols <- c(required_cols, "n_datasets")
  if(!is.null(size_col)) required_cols <- c(required_cols, size_col)
  if(!all(required_cols %in% colnames(meta_df))) {
    stop("meta_df must contain columns: ", paste(required_cols, collapse = ", "))
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

  # Build aesthetic mapping and point layer based on size_col
  # Use abs() for size when size_col has negative values
  if (is.null(size_col)) {
    aes_base <- aes(x = effect_size, y = -log10(.data[[p_val_col]]), color = group)
    use_abs_size <- FALSE
  } else {
    sz_vals <- meta_df[[size_col]]
    use_abs_size <- any(sz_vals < 0, na.rm = TRUE)
    aes_base <- aes(x = effect_size, y = -log10(.data[[p_val_col]]), color = group,
                    size = if (use_abs_size) abs(.data[[size_col]]) else .data[[size_col]])
    # For scale_size we need the effective values (abs when negative exists)
    sz_vals <- if (use_abs_size) abs(sz_vals) else sz_vals
  }

  # Basic volcano plot
  p <- ggplot(data = meta_df, aes_base)
  if (is.null(size_col)) {
    p <- p + geom_point(alpha = point_alpha, size = point_size)
  } else {
    min_sz <- min(sz_vals, na.rm = TRUE)
    max_sz <- max(sz_vals, na.rm = TRUE)
    range_sz <- max_sz - min_sz
    # Legend breaks: different logic for [0,1] fraction vs integer-like vs large range
    if (max_sz <= 1 || (range_sz < 1 && max_sz <= 1)) {
      # Fraction/proportion range (e.g. I2, correlation): use pretty decimal breaks
      size_breaks <- pretty(c(min_sz, max_sz), n = 5)
      size_breaks <- size_breaks[size_breaks >= min_sz & size_breaks <= max_sz]
      if (length(size_breaks) < 2) {
        size_breaks <- unique(c(min_sz, max_sz))
      }
    } else if (range_sz <= 5 && all(sz_vals == floor(sz_vals), na.rm = TRUE)) {
      # Small integer range (e.g. n_datasets): use integer seq
      size_breaks <- seq(min_sz, max_sz, by = 1)
    } else {
      # Large or mixed range: use pretty()
      size_breaks <- pretty(c(min_sz, max_sz), n = 4)
      size_breaks <- size_breaks[size_breaks >= min_sz & size_breaks <= max_sz]
    }
    p <- p + geom_point(alpha = point_alpha) +
      scale_size(range = c(point_size * 0.5, point_size * 2), name = size_col, breaks = size_breaks)
  }
  p <- p + theme_bw() +
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
      # Determine which genes to label
      if(!is.null(custom_labels)) {
        # Use custom labels if provided
        forlabel_df <- meta_df2[meta_df2$name %in% custom_labels,]
      } else {
        # Get top points to label, ensuring exact counts
        ordered_indices <- order(meta_df2$effect_size)
        n_low <- min(top_label_each, sum(meta_df2$group == "Down"))
        n_high <- min(top_label_each, sum(meta_df2$group == "Up"))
        low_indices <- head(ordered_indices, n_low)
        high_indices <- tail(ordered_indices, n_high)
        forlabel_names <- unique(c(meta_df2$name[low_indices], meta_df2$name[high_indices]))
        forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]
      }

      p <- p +
        (if (is.null(size_col))
          geom_point(shape = 1, size = point_size, data = forlabel_df, show.legend = FALSE)
        else
          geom_point(aes(size = if (use_abs_size) abs(.data[[size_col]]) else .data[[size_col]]), shape = 1, data = forlabel_df, show.legend = FALSE)) +
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

# Layout / Shared Utilities ----

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

# Priority Visualization ----

#' Plot top priority candidates as grouped horizontal bars
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @param show_priority_only Logical, whether to show only the priority bars
#'   without direction annotations.
#' @return A ggplot2 object.
#' @export
plotPriorityTopBar <- function(priority_df,
                               title = "Top prioritized candidates",
                               show_priority_only = FALSE) {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "PriorityScore", "candidate_class")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.logical(show_priority_only) || length(show_priority_only) != 1) {
    stop("show_priority_only must be TRUE or FALSE")
  }

  priority_df$display_name_wrapped <- .wrapVisualizationLabels(as.character(priority_df$display_name), width = 24)
  priority_df$display_name_wrapped <- factor(
    priority_df$display_name_wrapped,
    levels = rev(priority_df$display_name_wrapped)
  )
  max_score <- max(priority_df$PriorityScore, na.rm = TRUE)
  if (!show_priority_only) {
    extra_cols <- c("direction_consensus", "direction_pattern_short")
    if (!all(extra_cols %in% colnames(priority_df))) {
      stop("priority_df must contain columns: ", paste(extra_cols, collapse = ", "),
           " when show_priority_only = FALSE")
    }
    priority_df$direction_consensus_label <- dplyr::case_when(
      priority_df$direction_consensus == "Up" ~ "Positive",
      priority_df$direction_consensus == "Down" ~ "Negative",
      priority_df$direction_consensus == "Mixed" ~ "Mixed",
      TRUE ~ "NA"
    )
    priority_df$direction_consensus_label <- factor(
      priority_df$direction_consensus_label,
      levels = c("Positive", "Negative", "Mixed", "NA")
    )
    priority_df$annotation_x <- priority_df$PriorityScore + 0.03 * max_score
  }

  p <- ggplot2::ggplot(
    priority_df,
    ggplot2::aes(x = PriorityScore, y = display_name_wrapped, fill = candidate_class)
  ) +
    ggplot2::geom_col(width = 0.72, alpha = 0.92) +
    ggplot2::scale_fill_manual(values = c("Gene" = "#E98371", "Pathway" = "#66A9C9")) +
    ggplot2::labs(
      title = title,
      x = "Priority score",
      y = NULL,
      fill = "Candidate class"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!show_priority_only) {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(shape = direction_consensus_label),
        x = 0.02,
        size = 2.8,
        stroke = 0.7,
        fill = "white",
        color = "#2F2A24"
      ) +
      ggplot2::scale_shape_manual(
        values = c("Positive" = 24, "Negative" = 25, "Mixed" = 23, "NA" = 21),
        name = "Biomarker type"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = annotation_x, label = direction_pattern_short),
        hjust = 0,
        size = 3.1,
        color = "#4B433A"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1.22 * max_score), clip = "off") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5.5, 70, 5.5, 5.5)
      )
  }

  p
}

# Priority Evidence Summary ----

#' Plot decomposed priority evidence as stacked bars
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPriorityEvidenceStacked <- function(priority_df,
                                        title = "Priority score decomposition") {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "S_cellsets", "S_pdopdx", "S_clinical", "S_direction_adj")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }

  stacked_df <- rbind(
    data.frame(display_name = priority_df$display_name, component = "Cell sets", value = 0.25 * priority_df$S_cellsets),
    data.frame(display_name = priority_df$display_name, component = "PDO/PDX", value = 0.30 * priority_df$S_pdopdx),
    data.frame(display_name = priority_df$display_name, component = "Clinical", value = 0.35 * priority_df$S_clinical),
    data.frame(display_name = priority_df$display_name, component = "Direction", value = 0.10 * priority_df$S_direction_adj)
  )
  stacked_df$component <- factor(
    stacked_df$component,
    levels = c("Cell sets", "PDO/PDX", "Clinical", "Direction")
  )
  stacked_df$display_name_wrapped <- .wrapVisualizationLabels(as.character(stacked_df$display_name), width = 24)
  stacked_df$display_name_wrapped <- factor(
    stacked_df$display_name_wrapped,
    levels = rev(unique(stacked_df$display_name_wrapped))
  )

  ggplot2::ggplot(
    stacked_df,
    ggplot2::aes(x = value, y = display_name_wrapped, fill = component)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_manual(
      values = c(
        "Cell sets" = "#D8B365",
        "PDO/PDX" = "#80B1D3",
        "Clinical" = "#FB8072",
        "Direction" = "#B3B3B3"
      )
    ) +
    ggplot2::labs(
      title = title,
      x = "Weighted contribution to priority score",
      y = NULL,
      fill = "Evidence component"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot multi-source evidence heatmap for prioritized candidates
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPriorityEvidenceHeatmap <- function(priority_df,
                                        title = "Evidence heatmap across sources") {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "S_cellsets", "S_pdopdx", "S_clinical", "dir_cell", "dir_pdopdx", "dir_clin")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }

  heatmap_df <- rbind(
    data.frame(display_name = priority_df$display_name, source = "Cell sets", score = priority_df$S_cellsets, direction = priority_df$dir_cell),
    data.frame(display_name = priority_df$display_name, source = "PDO/PDX", score = priority_df$S_pdopdx, direction = priority_df$dir_pdopdx),
    data.frame(display_name = priority_df$display_name, source = "Clinical", score = priority_df$S_clinical, direction = priority_df$dir_clin)
  )
  heatmap_df$display_name_wrapped <- .wrapVisualizationLabels(as.character(heatmap_df$display_name), width = 24)
  heatmap_df$display_name_wrapped <- factor(
    heatmap_df$display_name_wrapped,
    levels = rev(unique(heatmap_df$display_name_wrapped))
  )
  heatmap_df$source <- factor(heatmap_df$source, levels = c("Cell sets", "PDO/PDX", "Clinical"))
  heatmap_df$direction_symbol <- dplyr::case_when(
    heatmap_df$direction == "Up" ~ "U",
    heatmap_df$direction == "Down" ~ "D",
    TRUE ~ "."
  )

  ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = source, y = display_name_wrapped, fill = score)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = direction_symbol), size = 4.8, color = "#1F1F1F") +
    ggplot2::scale_fill_gradient(low = "#F6F1E8", high = "#C65D53", limits = c(0, 1), name = "Source score") +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid = ggplot2::element_blank()
    )
}

#' Plot a GSVA heatmap
#'
#' @param gsva_heatmap A numeric matrix prepared by \code{prepareGSVAHeatmapMatrix()} with
#'   pathways in rows and datasets in columns.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param low_color Low-value color.
#' @param mid_color Midpoint color.
#' @param high_color High-value color.
#' @param midpoint Midpoint for the diverging color scale.
#' @param border_color Tile border color.
#' @param cluster_rows Logical, whether to cluster pathway order.
#' @param cluster_cols Logical, whether to cluster dataset order.
#' @param show_x_text Logical, whether to show dataset labels.
#' @return A ggplot2 object.
#' @export
plotGSVAHeatmap <- function(gsva_heatmap,
                            title = "GSVA pathway activity",
                            subtitle = "Representative pathway score matrix",
                            low_color = "#6AAED6",
                            mid_color = "#FBFAF7",
                            high_color = "#E98484",
                            midpoint = 0,
                            border_color = "#F6F2EA",
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            show_x_text = FALSE) {
  if (!is.matrix(gsva_heatmap) || !is.numeric(gsva_heatmap)) {
    stop("gsva_heatmap must be a numeric matrix")
  }
  if (is.null(rownames(gsva_heatmap)) || is.null(colnames(gsva_heatmap))) {
    stop("gsva_heatmap must have rownames and colnames")
  }

  row_order <- rownames(gsva_heatmap)
  col_order <- colnames(gsva_heatmap)
  cluster_mat <- gsva_heatmap

  if (anyNA(cluster_mat)) {
    row_means <- apply(cluster_mat, 1, function(x) mean(x, na.rm = TRUE))
    row_means[!is.finite(row_means)] <- 0
    for (i in seq_len(nrow(cluster_mat))) {
      na_idx <- is.na(cluster_mat[i, ])
      if (any(na_idx)) {
        cluster_mat[i, na_idx] <- row_means[i]
      }
    }
  }

  if (cluster_rows && nrow(gsva_heatmap) > 1) {
    row_dist <- stats::dist(cluster_mat)
    if (length(row_dist) > 0 && all(is.finite(row_dist))) {
      row_order <- rownames(gsva_heatmap)[stats::hclust(row_dist)$order]
    }
  }

  if (cluster_cols && ncol(gsva_heatmap) > 1) {
    col_dist <- stats::dist(t(cluster_mat))
    if (length(col_dist) > 0 && all(is.finite(col_dist))) {
      col_order <- colnames(gsva_heatmap)[stats::hclust(col_dist)$order]
    }
  }

  heatmap_df <- as.data.frame(as.table(gsva_heatmap), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("Pathway", "Dataset", "Score")
  heatmap_df$Pathway <- factor(heatmap_df$Pathway, levels = rev(row_order))
  heatmap_df$Dataset <- factor(heatmap_df$Dataset, levels = col_order)

  ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Dataset, y = Pathway, fill = Score)) +
    ggplot2::geom_tile(color = border_color, linewidth = 0.45) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = midpoint,
      na.value = "#E7E1D7",
      name = "Score"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = if (show_x_text) {
        ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, color = "#2F2A24")
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(color = "#2F2A24"),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13, color = "#2F2A24"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#6F685E"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(color = "#2F2A24")
    )
}

#' Plot a drug-by-pathway effect size heatmap
#'
#' @param effect_size_matrix Numeric matrix with drugs in rows and pathways in columns.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param low_color Low-value color.
#' @param mid_color Midpoint color.
#' @param high_color High-value color.
#' @param midpoint Midpoint for the diverging color scale.
#' @param border_color Tile border color.
#' @param zscore Optional z-score transform before plotting: \code{NULL} (none),
#'   \code{"row"} (per drug), \code{"col"} (per pathway), or \code{"both"}
#'   (row z-score then column z-score). NAs are temporarily imputed for scaling
#'   then restored for display.
#' @param cluster_rows Logical, whether to cluster drugs.
#' @param cluster_cols Logical, whether to cluster pathways.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @return A ggplot2 object.
#' @export
plotDrugPathwayEffectHeatmap <- function(effect_size_matrix,
                                         title = "Drug-pathway effect size heatmap",
                                         subtitle = "Effect sizes from GSVA-feature association analysis",
                                         low_color = "#6AAED6",
                                         mid_color = "#FBFAF7",
                                         high_color = "#E98484",
                                         midpoint = 0,
                                         border_color = "#F6F2EA",
                                         zscore = NULL,
                                         cluster_rows = FALSE,
                                         cluster_cols = FALSE,
                                         x_label = NULL,
                                         y_label = NULL) {
  if (!is.matrix(effect_size_matrix) || !is.numeric(effect_size_matrix)) {
    stop("effect_size_matrix must be a numeric matrix")
  }
  if (is.null(rownames(effect_size_matrix)) || is.null(colnames(effect_size_matrix))) {
    stop("effect_size_matrix must have rownames and colnames")
  }

  plot_mat <- effect_size_matrix
  na_mask <- is.na(plot_mat)

  if (!is.null(zscore)) {
    zscore <- match.arg(zscore, c("row", "col", "both"))
    if (anyNA(plot_mat)) {
      fill_value <- mean(plot_mat, na.rm = TRUE)
      if (!is.finite(fill_value)) {
        fill_value <- 0
      }
      plot_mat[na_mask] <- fill_value
    }
    if (zscore %in% c("row", "both")) {
      plot_mat <- t(scale(t(plot_mat)))
    }
    if (zscore %in% c("col", "both")) {
      plot_mat <- scale(plot_mat)
    }
    plot_mat[is.nan(plot_mat)] <- 0
    if (any(na_mask)) {
      plot_mat[na_mask] <- NA
    }
  }

  fill_label <- if (is.null(zscore)) "Effect size" else "Z-score"
  row_order <- rownames(plot_mat)
  col_order <- colnames(plot_mat)
  cluster_mat <- plot_mat

  if (anyNA(cluster_mat)) {
    fill_value <- mean(cluster_mat, na.rm = TRUE)
    if (!is.finite(fill_value)) {
      fill_value <- 0
    }
    cluster_mat[is.na(cluster_mat)] <- fill_value
  }

  if (cluster_rows && nrow(cluster_mat) > 1) {
    row_order <- rownames(cluster_mat)[stats::hclust(stats::dist(cluster_mat))$order]
  }
  if (cluster_cols && ncol(cluster_mat) > 1) {
    col_order <- colnames(cluster_mat)[stats::hclust(stats::dist(t(cluster_mat)))$order]
  }

  heatmap_df <- as.data.frame(as.table(plot_mat), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("Drug", "Pathway", "EffectSize")
  heatmap_df$Drug <- factor(heatmap_df$Drug, levels = rev(row_order))
  heatmap_df$Pathway <- factor(heatmap_df$Pathway, levels = col_order)

  ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Pathway, y = Drug, fill = EffectSize)) +
    ggplot2::geom_tile(color = border_color, linewidth = 0.45) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = midpoint,
      na.value = "#E7E1D7",
      name = fill_label
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 28, hjust = 1, vjust = 1, color = "#2F2A24"),
      axis.text.y = ggplot2::element_text(color = "#2F2A24"),
      axis.title = ggplot2::element_text(color = "#2F2A24"),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13, color = "#2F2A24"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#6F685E"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(color = "#2F2A24")
    )
}

.wrapCandidateSelectionLabels <- function(x, width = 18) {
  vapply(x, function(label) {
    paste(strwrap(as.character(label), width = width), collapse = "\n")
  }, character(1))
}

#' Plot clinically supported candidate selection
#'
#' @param candidate_df Data frame with candidate-level preclinical and clinical support annotations.
#' @param n_retained Number of retained candidates to display.
#' @param n_filtered Number of filtered candidates to display.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param label_width Width used for wrapping candidate labels.
#' @return A patchwork plot object.
#' @export
plotClinicallySupportedCandidateSelection <- function(candidate_df,
                                                      n_retained = 6,
                                                      n_filtered = 4,
                                                      add = NULL,
                                                      title = "Clinically supported candidate selection",
                                                      subtitle = "Preclinical hits retained after clinical validation and concordant direction",
                                                      label_width = 18) {
  candidate_df <- as.data.frame(candidate_df)
  required_cols <- c(
    "name", "effect_size_clinical", "cell_supported", "pdopdx_supported",
    "clinical_supported", "direction_concordant", "retained"
  )
  if (!all(required_cols %in% colnames(candidate_df))) {
    stop("candidate_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(n_retained) || length(n_retained) != 1 || n_retained < 0) {
    stop("n_retained must be a single non-negative number")
  }
  if (!is.numeric(n_filtered) || length(n_filtered) != 1 || n_filtered < 0) {
    stop("n_filtered must be a single non-negative number")
  }
  if (!is.null(add) && !is.character(add)) {
    stop("add must be NULL or a character vector of candidate names")
  }

  candidate_df$clinical_supported[is.na(candidate_df$clinical_supported)] <- FALSE
  candidate_df$direction_concordant[is.na(candidate_df$direction_concordant)] <- FALSE
  candidate_df$retained[is.na(candidate_df$retained)] <- FALSE
  if (!"clinical_abs" %in% colnames(candidate_df)) {
    candidate_df$clinical_abs <- abs(ifelse(is.na(candidate_df$effect_size_clinical), 0, candidate_df$effect_size_clinical))
  }
  if (!"preclinical_abs" %in% colnames(candidate_df)) {
    if ("effect_size_preclinical" %in% colnames(candidate_df)) {
      candidate_df$preclinical_abs <- abs(ifelse(is.na(candidate_df$effect_size_preclinical), 0, candidate_df$effect_size_preclinical))
    } else {
      candidate_df$preclinical_abs <- 0
    }
  }

  retained_df <- candidate_df[candidate_df$retained, , drop = FALSE]
  retained_df <- retained_df[order(-retained_df$clinical_abs, -retained_df$preclinical_abs), , drop = FALSE]
  filtered_df <- candidate_df[!candidate_df$retained, , drop = FALSE]
  filtered_df <- filtered_df[order(-filtered_df$clinical_abs, -filtered_df$preclinical_abs), , drop = FALSE]

  add_candidates <- character()
  if (!is.null(add)) {
    add_candidates <- unique(add[!is.na(add) & nzchar(add)])
    add_candidates <- add_candidates[add_candidates %in% candidate_df$name]
  }

  plot_candidates <- unique(c(
    head(retained_df$name, n_retained),
    head(filtered_df$name, n_filtered),
    add_candidates
  ))
  if (length(plot_candidates) == 0) {
    plot_candidates <- head(candidate_df$name, min(nrow(candidate_df), n_retained + n_filtered))
  }

  plot_df <- candidate_df[candidate_df$name %in% plot_candidates, , drop = FALSE]
  plot_df <- plot_df[order(
    -plot_df$retained,
    -plot_df$clinical_supported,
    -plot_df$clinical_abs,
    -plot_df$preclinical_abs
  ), , drop = FALSE]

  plot_df$name_wrapped <- .wrapCandidateSelectionLabels(plot_df$name, width = label_width)
  plot_df$name_wrapped <- factor(plot_df$name_wrapped, levels = rev(plot_df$name_wrapped))
  plot_df$retained_label <- ifelse(plot_df$retained, "clinical biomarker", "preclinical biomarker")

  evidence_df <- rbind(
    data.frame(name_wrapped = plot_df$name_wrapped, Evidence = "Cell/PDC", Status = plot_df$cell_supported),
    data.frame(name_wrapped = plot_df$name_wrapped, Evidence = "PDO/PDX", Status = plot_df$pdopdx_supported),
    data.frame(name_wrapped = plot_df$name_wrapped, Evidence = "Clinical", Status = plot_df$clinical_supported),
    data.frame(name_wrapped = plot_df$name_wrapped, Evidence = "Concordant", Status = plot_df$direction_concordant)
  )
  evidence_df$Evidence <- factor(evidence_df$Evidence, levels = c("Cell/PDC", "PDO/PDX", "Clinical", "Concordant"))
  evidence_df$Marker <- ifelse(evidence_df$Status, "•", "")

  evidence_plot <- ggplot2::ggplot(evidence_df, ggplot2::aes(Evidence, name_wrapped, fill = Status)) +
    ggplot2::geom_tile(color = "#F6F2EA", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = Marker), size = 6.2, color = "#2F2A24") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#F6D3CD", `FALSE` = "#F4EFE7")) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 24, hjust = 1, vjust = 1, color = "#2F2A24"),
      axis.text.y = ggplot2::element_text(color = "#2F2A24"),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  effect_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(effect_size_clinical, name_wrapped, fill = retained_label)) +
    ggplot2::geom_vline(xintercept = 0, color = "#B0A99F", linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = effect_size_clinical, yend = name_wrapped),
      color = "#C9C1B6", linewidth = 0.95
    ) +
    ggplot2::geom_point(shape = 21, size = 3.6, color = "#6F685E", stroke = 0.55) +
    ggplot2::scale_fill_manual(values = c("clinical biomarker" = "#FB8A7A", "preclinical biomarker" = "#D9D3C8")) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "effect size",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15, color = "#2F2A24"),
      plot.subtitle = ggplot2::element_text(size = 12, color = "#6F685E"),
      axis.text.x = ggplot2::element_text(color = "#2F2A24"),
      axis.title.x = ggplot2::element_text(color = "#2F2A24"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(color = "#2F2A24", size = 12.5),
      legend.position = "top",
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )

  patchwork::wrap_plots(evidence_plot, effect_plot, nrow = 1, widths = c(1.1, 2.3))
}

#' Plot clinical versus preclinical evidence as a bubble chart
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPriorityClinicalPreclinicalBubble <- function(priority_df,
                                                  title = "Clinical vs preclinical support") {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "S_preclinical", "S_clinical", "candidate_class", "n_selected_sources", "direction_consensus", "PriorityScore")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  priority_df$display_name_wrapped <- .wrapVisualizationLabels(as.character(priority_df$display_name), width = 24)
  priority_df$candidate_class <- factor(priority_df$candidate_class, levels = c("Gene", "Pathway"))
  priority_df$direction_consensus_label <- dplyr::case_when(
    priority_df$direction_consensus == "Up" ~ "Positive",
    priority_df$direction_consensus == "Down" ~ "Negative",
    priority_df$direction_consensus == "Mixed" ~ "Mixed",
    TRUE ~ "NA"
  )
  priority_df$direction_consensus_label <- factor(
    priority_df$direction_consensus_label,
    levels = c("Positive", "Negative", "Mixed", "NA")
  )

  ggplot2::ggplot(
    priority_df,
    ggplot2::aes(
      x = S_preclinical,
      y = S_clinical,
      size = PriorityScore,
      shape = direction_consensus_label,
      color = candidate_class,
      alpha = n_selected_sources
    )
  ) +
    ggplot2::geom_point(stroke = 0.8, fill = "white") +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = display_name_wrapped),
      size = 3.5,
      color = "black",
      box.padding = 0.9,
      point.padding = 0.6,
      force = 8,
      force_pull = 0.5,
      min.segment.length = 0,
      max.overlaps = 80
    ) +
    ggplot2::scale_color_manual(values = c("Gene" = "#E98371", "Pathway" = "#66A9C9")) +
    ggplot2::scale_size_continuous(range = c(3.5, 10), name = "Priority score") +
    ggplot2::scale_shape_manual(values = c("Positive" = 24, "Negative" = 25, "Mixed" = 23, "NA" = 21),
                                name = "Biomarker type") +
    ggplot2::scale_alpha_continuous(range = c(0.65, 1), guide = "none") +
    ggplot2::labs(
      title = title,
      x = "Mean preclinical evidence score",
      y = "Clinical evidence score",
      color = "candidate class"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 11, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Drug Feature Visualization ----

#' Plot continuous comparison
#'
#' @description Creates a scatter plot showing the association between a
#'   continuous variable and a response value.
#' @param data Data frame containing the variables to plot.
#' @param cont_column Name of the continuous column to use for x-axis.
#' @param value_column Name of the value column to use for y-axis.
#' @param value_label Label for the value variable.
#' @param title Plot title.
#' @return A ggplot2 object with the comparison plot.
#' @export
plotContinuousComparison <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", title = NULL) {
  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "vs", cont_column)
  }

  ggscatter(data, x = cont_column, y = value_column, alpha = 0.2) +
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
}

#' Plot continuous groups
#'
#' @description Creates boxplots for binned groups of a continuous variable.
#' @param data Data frame containing the variables to plot.
#' @param cont_column Name of the continuous column to bin into groups.
#' @param value_column Name of the value column to use for y-axis.
#' @param value_label Label for the value variable.
#' @param num_bins Number of bins to create from the continuous variable.
#' @param title Plot title.
#' @return A ggplot2 object with the grouped boxplot.
#' @export
plotContinuousGroups <- function(data, cont_column, value_column = "value", value_label = "Drug Sensitivity", num_bins = 4, title = NULL) {
  data <- data[!is.na(data[[cont_column]]) & !is.na(data[[value_column]]), ]

  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "by", cont_column, "Group")
  }

  cont_values <- data[[cont_column]]
  cont_bins <- cut(cont_values, breaks = num_bins, include.lowest = TRUE, labels = FALSE)
  cont_range <- range(cont_values, na.rm = TRUE)
  bin_width <- diff(cont_range) / num_bins
  group_labels <- sapply(1:num_bins, function(i) {
    lower <- cont_range[1] + (i - 1) * bin_width
    upper <- cont_range[1] + i * bin_width
    paste0(round(lower), "-", round(upper))
  })

  group_column <- paste0(cont_column, "_group")
  label_column <- paste0(cont_column, "_group_label")
  data[[group_column]] <- cont_bins
  data[[label_column]] <- factor(group_labels[data[[group_column]]], levels = group_labels)
  overall_median <- median(data[[value_column]], na.rm = TRUE)

  ggboxplot(data, x = label_column, y = value_column,
            fill = label_column, palette = soft_palette_26,
            add = "jitter", add.params = list(alpha = 0.15)) +
    geom_hline(yintercept = overall_median, linetype = "dashed",
               color = "gray40", size = 0.8, alpha = 0.7) +
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(data[[value_column]]) - max(data[[value_column]]) / 8),
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
}

#' Plot category comparison
#'
#' @description Creates boxplots comparing values across categories.
#' @param data Data frame containing the variables to plot.
#' @param category_column Name of the categorical column to use for grouping.
#' @param value_column Name of the value column to use for y-axis.
#' @param value_label Label for the value variable.
#' @param title Plot title.
#' @return A ggplot2 object with the category comparison plot.
#' @export
plotCategoryComparison <- function(data, category_column, value_column = "value", value_label = "Drug Sensitivity", title = NULL) {
  data <- data[!is.na(data[[category_column]]) & !is.na(data[[value_column]]), ]

  y_label <- if (value_column == "zscore_value") {
    "Z-score Drug Sensitivity(Area Above Curve)"
  } else {
    "Drug Sensitivity(Area Above Curve)"
  }

  plot_title <- if (!is.null(title)) {
    title
  } else {
    paste(value_label, "by", category_column)
  }

  category_counts <- table(data[[category_column]])
  valid_categories <- names(category_counts)[category_counts >= 3]

  if (length(valid_categories) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Not enough samples per category for comparison") +
             theme_void())
  }

  data_filtered <- data[data[[category_column]] %in% valid_categories, ]
  category_medians <- tapply(data_filtered[[value_column]],
                             data_filtered[[category_column]],
                             median, na.rm = TRUE)
  sorted_categories <- names(sort(category_medians))

  if (category_column == "TumorType" && "non-cancer" %in% sorted_categories) {
    sorted_categories <- c(sorted_categories[sorted_categories != "non-cancer"], "non-cancer")
  }

  data_filtered[[category_column]] <- factor(data_filtered[[category_column]], levels = sorted_categories)
  overall_median <- median(data_filtered[[value_column]], na.rm = TRUE)

  p <- ggboxplot(data_filtered, x = category_column, y = value_column,
                 fill = category_column, palette = soft_palette_26,
                 add = "jitter", add.params = list(alpha = 0.15)) +
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
    coord_cartesian(ylim = c(NA, max(data_filtered[[value_column]], na.rm = TRUE) + max(data_filtered[[value_column]], na.rm = TRUE) / 20)) +
    ggtitle(plot_title) +
    xlab(NULL) +
    ylab(y_label)

  if (length(valid_categories) >= 2) {
    if (length(valid_categories) == 2) {
      p <- p + stat_compare_means(size = 6,
                                  label.x.npc = "left",
                                  label.y.npc = "top",
                                  label = "p.format")
    } else {
      p <- p + stat_compare_means(method = "kruskal.test",
                                  size = 6,
                                  label.x.npc = "left",
                                  label.y.npc = "top",
                                  label = "p.format")
    }
  }

  p
}

#' Create Drug Sensitivity Rank Plot
#'
#' @description Creates a rank plot showing drug sensitivity values ordered
#'   from most sensitive to least sensitive.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object.
#' @param select_drugs Character string specifying the drug name.
#' @param data_type Filter by data type.
#' @param tumor_type Filter by tumor type.
#' @param overlap_only Whether to use only overlapping samples.
#' @param highlight Highlighting specification.
#' @param color Variable to use for coloring points.
#' @param zscore Logical, whether to use z-score normalized values.
#' @param merge Logical, whether to merge multiple projects for plotting.
#' @param point_size Numeric, size of points in the plot.
#' @param highlight_alpha Numeric, alpha transparency for non-highlighted points.
#' @param sample_annotations Optional dataframe containing sample annotations.
#' @param db_path Optional path to SQLite database for loading sample annotations.
#' @return A ggplot2 object or list of ggplot2 objects.
#' @export
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
  if (missing(select_drugs) || is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
  }
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("dromaset_object must be a DromaSet or MultiDromaSet object")
  }

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

  value_column <- if (zscore) "zscore_value" else "raw_value"
  if (!value_column %in% colnames(drug_data)) {
    stop("Required value column '", value_column, "' not found in data")
  }

  drug_data <- drug_data[!is.na(drug_data[[value_column]]), ]
  if (nrow(drug_data) == 0) {
    stop("No valid drug sensitivity values found after removing missing data")
  }

  if (inherits(dromaset_object, "MultiDromaSet")) {
    if (merge && zscore) {
      study_counts <- table(drug_data$ProjectID)
      if (length(study_counts) < 2) {
        warning("Drug '", select_drugs, "' found in only ", length(study_counts),
                " project(s). merge=TRUE requires at least 2 projects. Using individual project data.")
        merge <- FALSE
      } else {
        return(createSingleRankPlot(drug_data, select_drugs, highlight, color, zscore,
                                    value_column, point_size, highlight_alpha, merge))
      }
    }

    if (!merge) {
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

  createSingleRankPlot(drug_data, select_drugs, highlight, color, zscore,
                       value_column, point_size, highlight_alpha, merge)
}

#' @keywords internal
createSingleRankPlot <- function(drug_data, select_drugs, highlight, color, zscore,
                               value_column, point_size, highlight_alpha, merge,
                               projects = NULL) {
  drug_data <- drug_data[order(drug_data[[value_column]], decreasing = TRUE), ]
  drug_data$rank <- seq_len(nrow(drug_data))

  highlight_samples <- character(0)
  is_numeric_highlight <- FALSE
  if (!is.null(highlight)) {
    if (length(highlight) == 1 && is.numeric(highlight)) {
      is_numeric_highlight <- TRUE
      unique_samples_ordered <- unique(drug_data$SampleID)
      n_highlight <- min(as.integer(highlight), length(unique_samples_ordered))
      highlight_samples <- unique_samples_ordered[seq_len(n_highlight)]
    } else if (length(highlight) == 1) {
      if ("DataType" %in% colnames(drug_data) && highlight %in% drug_data$DataType) {
        highlight_samples <- drug_data$SampleID[drug_data$DataType == highlight]
      } else if ("TumorType" %in% colnames(drug_data) && highlight %in% drug_data$TumorType) {
        highlight_samples <- drug_data$SampleID[drug_data$TumorType == highlight]
      } else {
        for (col in colnames(drug_data)) {
          if (col != "SampleID" && highlight %in% drug_data[[col]]) {
            highlight_samples <- drug_data$SampleID[drug_data[[col]] == highlight]
            break
          }
        }
      }
      if (length(highlight_samples) == 0) {
        highlight_samples <- highlight[highlight %in% drug_data$SampleID]
      }
    } else {
      highlight_samples <- highlight[highlight %in% drug_data$SampleID]
    }
  }

  drug_data$highlighted <- drug_data$SampleID %in% highlight_samples

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

  if (!is.null(color_data)) {
    drug_data$color_var <- color_data
    p <- ggplot(drug_data, aes(x = rank, y = .data[[value_column]])) +
      geom_point(aes(color = color_var), size = point_size,
                 alpha = ifelse(drug_data$highlighted, 1, highlight_alpha))
  } else {
    p <- ggplot(drug_data, aes(x = rank, y = .data[[value_column]])) +
      geom_point(size = point_size,
                 alpha = ifelse(drug_data$highlighted, 1, highlight_alpha),
                 color = "steelblue")
  }

  if (length(highlight_samples) > 0) {
    highlighted_data <- drug_data[drug_data$highlighted, ]
    if (is_numeric_highlight) {
      n_unique_highlighted <- length(unique(highlighted_data$SampleID))
      if (n_unique_highlighted > 20) {
        highlighted_data <- highlighted_data[order(highlighted_data$rank), ]
        unique_samples_to_label <- unique(highlighted_data$SampleID)[1:20]
        highlighted_data_to_label <- highlighted_data[highlighted_data$SampleID %in% unique_samples_to_label, ]
        highlighted_data_to_label <- highlighted_data_to_label[!duplicated(highlighted_data_to_label$SampleID), ]
        message("Note: ", n_unique_highlighted, " unique samples were highlighted, but only the top 20 will be labeled to avoid overcrowding.")
      } else {
        highlighted_data_to_label <- highlighted_data[!duplicated(highlighted_data$SampleID), ]
      }
      highlighted_data_for_circle <- highlighted_data[!duplicated(highlighted_data$SampleID), ]
    } else {
      if (nrow(highlighted_data) > 20) {
        highlighted_data <- highlighted_data[order(highlighted_data$rank), ]
        highlighted_data_to_label <- highlighted_data[1:20, ]
        message("Note: ", nrow(highlighted_data), " samples were highlighted, but only the top 20 will be labeled to avoid overcrowding.")
      } else {
        highlighted_data_to_label <- highlighted_data
      }
      highlighted_data_for_circle <- highlighted_data
    }

    if (merge && "ProjectID" %in% colnames(highlighted_data_to_label)) {
      highlighted_data_to_label$label <- paste0(highlighted_data_to_label$SampleID, " (", highlighted_data_to_label$ProjectID, ")")
    } else {
      highlighted_data_to_label$label <- highlighted_data_to_label$SampleID
    }

    p <- p + geom_point(data = highlighted_data_for_circle,
                        aes(x = rank, y = .data[[value_column]]),
                        size = point_size + 0.5,
                        shape = 1,
                        color = "black",
                        stroke = 1.2)

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

  value_label <- if (zscore) "Z-score Drug Sensitivity(Area Above Curve)" else "Drug Sensitivity(Area Above Curve)"
  plot_title <- paste("Drug Sensitivity Rank Plot:", select_drugs)
  if (!is.null(projects)) {
    plot_title <- paste(plot_title, "-", projects)
  } else if (merge) {
    plot_title <- paste(plot_title, "(Merged)")
  }

  lab_list <- list(title = plot_title, x = "Rank (1 = Most Sensitive)", y = value_label)
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

  if (!is.null(color_data)) {
    if (is.factor(color_data) || is.character(color_data)) {
      n_colors <- length(unique(color_data))
      if (n_colors <= length(soft_palette_26)) {
        p <- p + scale_color_manual(values = soft_palette_26[1:n_colors])
      } else {
        extended_colors <- colorRampPalette(soft_palette_26)(n_colors)
        p <- p + scale_color_manual(values = extended_colors)
      }
    } else {
      p <- p + scale_color_gradient(low = "#80B1D3", high = "#FB8072")
    }
  }

  n_samples <- nrow(drug_data)
  n_highlighted <- sum(drug_data$highlighted)
  median_value <- median(drug_data[[value_column]], na.rm = TRUE)
p + labs(subtitle = paste0("n = ", n_samples, "; highlighted = ", n_highlighted, "; median = ", round(median_value, 3)))
}

# Enrichment Visualization ----

#' Create a dotplot for enrichment results
#'
#' @description Creates a dotplot visualization for enrichResult objects.
#' @param enrich_result An enrichResult object from clusterProfiler.
#' @param use_padj Logical, whether to use adjusted p-value.
#' @param showCategory Number of top pathways to display.
#' @param title Plot title.
#' @param x_label Label for x-axis.
#' @param size_by Variable for point size.
#' @param color_low Low-significance color.
#' @param color_high High-significance color.
#' @param ... Additional parameters passed to ggplot2 functions.
#' @return A ggplot2 object with dotplot visualization.
#' @export
plotEnrichDotplot <- function(enrich_result,
                              use_padj = TRUE,
                              showCategory = 10,
                              title = NULL,
                              x_label = NULL,
                              size_by = c("Count", "GeneRatio"),
                              color_low = "#BEBADAFF",
                              color_high = "#FB8072FF",
                              ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it using: install.packages('ggplot2')")
  }
  if (missing(enrich_result) || is.null(enrich_result)) {
    stop("enrich_result is required and must be an enrichResult object")
  }
  if (!inherits(enrich_result, "enrichResult")) {
    stop("enrich_result must be an enrichResult object from clusterProfiler")
  }
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  size_by <- match.arg(size_by)
  result_df <- enrich_result@result
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("No enrichment results found in enrich_result object")
  }

  pval_col <- if (use_padj && "p.adjust" %in% colnames(result_df)) {
    "p.adjust"
  } else if ("pvalue" %in% colnames(result_df)) {
    "pvalue"
  } else {
    stop("Neither p.adjust nor pvalue column found in enrichment results")
  }

  if (is.numeric(showCategory) && showCategory > 0) {
    result_df <- result_df[order(result_df[[pval_col]]), , drop = FALSE]
    result_df <- head(result_df, min(showCategory, nrow(result_df)))
  } else if (is.character(showCategory)) {
    if (!"Description" %in% colnames(result_df)) {
      stop("Description column not found. Cannot filter by pathway names.")
    }
    result_df <- result_df[result_df$Description %in% showCategory, , drop = FALSE]
    if (nrow(result_df) == 0) {
      stop("None of the specified pathways found in enrichment results")
    }
  } else {
    stop("showCategory must be a positive integer or a character vector of pathway names")
  }

  plot_df <- result_df
  plot_df$x_value <- -log10(plot_df[[pval_col]])
  if (size_by == "Count" && "Count" %in% colnames(plot_df)) {
    plot_df$size_value <- plot_df$Count
    size_label <- "Gene Count"
  } else if (size_by == "GeneRatio" && "GeneRatio" %in% colnames(plot_df)) {
    plot_df$size_value <- sapply(plot_df$GeneRatio, function(x) {
      parts <- strsplit(as.character(x), "/")[[1]]
      if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA
    })
    size_label <- "Gene Ratio"
  } else {
    if ("Count" %in% colnames(plot_df)) {
      plot_df$size_value <- plot_df$Count
      size_label <- "Gene Count"
    } else {
      plot_df$size_value <- 1
      size_label <- "Size"
      warning("Neither Count nor GeneRatio found. Using constant size.")
    }
  }

  if (!"Description" %in% colnames(plot_df)) {
    stop("Description column not found in enrichment results")
  }
  plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description[order(plot_df[[pval_col]])]))

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$x_value,
                                    y = .data$Description,
                                    size = .data$size_value,
                                    color = .data[[pval_col]])) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_gradient(low = color_high, high = color_low,
                                  name = if (use_padj) "p.adjust" else "pvalue",
                                  trans = "log10") +
    ggplot2::scale_size_continuous(name = size_label, range = c(3, 8)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right"
    ) +
    ggplot2::xlab(if (!is.null(x_label)) x_label else paste0("-log10(", if (use_padj) "p.adjust" else "pvalue", ")")) +
    ggplot2::ylab("Pathway")

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  } else {
    enrich_type <- if ("ONTOLOGY" %in% colnames(result_df)) {
      paste0("GO ", result_df$ONTOLOGY[1], " Enrichment")
    } else {
      "KEGG Pathway Enrichment"
    }
    p <- p + ggplot2::ggtitle(enrich_type)
  }

  p
}

#' Create a barplot for enrichment results
#'
#' @description Creates a barplot visualization for enrichResult objects.
#' @param enrich_result An enrichResult object from clusterProfiler.
#' @param use_padj Logical, whether to use adjusted p-value.
#' @param showCategory Number of top pathways to display.
#' @param title Plot title.
#' @param x_label Label for x-axis.
#' @param color_low Low-significance color.
#' @param color_high High-significance color.
#' @param ... Additional parameters passed to ggplot2 functions.
#' @return A ggplot2 object with barplot visualization.
#' @export
plotEnrichBarplot <- function(enrich_result,
                              use_padj = TRUE,
                              showCategory = 10,
                              title = NULL,
                              x_label = NULL,
                              color_low = "#BEBADAFF",
                              color_high = "#FB8072FF",
                              ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it using: install.packages('ggplot2')")
  }
  if (missing(enrich_result) || is.null(enrich_result)) {
    stop("enrich_result is required and must be an enrichResult object")
  }
  if (!inherits(enrich_result, "enrichResult")) {
    stop("enrich_result must be an enrichResult object from clusterProfiler")
  }
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  result_df <- enrich_result@result
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("No enrichment results found in enrich_result object")
  }

  pval_col <- if (use_padj && "p.adjust" %in% colnames(result_df)) {
    "p.adjust"
  } else if ("pvalue" %in% colnames(result_df)) {
    "pvalue"
  } else {
    stop("Neither p.adjust nor pvalue column found in enrichment results")
  }

  if (is.numeric(showCategory) && showCategory > 0) {
    result_df <- result_df[order(result_df[[pval_col]]), , drop = FALSE]
    result_df <- head(result_df, min(showCategory, nrow(result_df)))
  } else if (is.character(showCategory)) {
    if (!"Description" %in% colnames(result_df)) {
      stop("Description column not found. Cannot filter by pathway names.")
    }
    result_df <- result_df[result_df$Description %in% showCategory, , drop = FALSE]
    if (nrow(result_df) == 0) {
      stop("None of the specified pathways found in enrichment results")
    }
  } else {
    stop("showCategory must be a positive integer or a character vector of pathway names")
  }

  plot_df <- result_df
  plot_df$x_value <- -log10(plot_df[[pval_col]])
  if (!"Description" %in% colnames(plot_df)) {
    stop("Description column not found in enrichment results")
  }
  plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description[order(plot_df[[pval_col]])]))

  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = .data$x_value,
                                    y = .data$Description,
                                    fill = .data[[pval_col]])) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.8) +
    ggplot2::scale_fill_gradient(low = color_high, high = color_low,
                                 name = if (use_padj) "p.adjust" else "pvalue",
                                 trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right"
    ) +
    ggplot2::xlab(if (!is.null(x_label)) x_label else paste0("-log10(", if (use_padj) "p.adjust" else "pvalue", ")")) +
    ggplot2::ylab("Pathway")

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  } else {
    enrich_type <- if ("ONTOLOGY" %in% colnames(result_df)) {
      paste0("GO ", result_df$ONTOLOGY[1], " Enrichment")
    } else {
      "KEGG Pathway Enrichment"
    }
    p <- p + ggplot2::ggtitle(enrich_type)
  }

  p
}
