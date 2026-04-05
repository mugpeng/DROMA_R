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
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 15, hjust = 0.85),
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

#' Prepare merged priority table for visualization
#'
#' @description Merges gene and pathway priority tables, annotates candidate type,
#'   adds shortened labels, and derives summary columns used by plotting helpers.
#' @param gene_priority Data frame produced by \code{buildPriorityTable()} for genes.
#' @param pathway_priority Data frame produced by \code{buildPriorityTable()} for pathways.
#' @param top_n Optional integer. If provided, keep only the top N candidates by
#'   \code{PriorityScore} after merging.
#' @param pathway_prefix Prefix to strip from pathway names for display.
#' @return A merged and annotated data frame ready for plotting.
#' @export
preparePriorityVisualizationData <- function(gene_priority,
                                             pathway_priority,
                                             top_n = NULL,
                                             pathway_prefix = "^HALLMARK_") {
  gene_priority <- as.data.frame(gene_priority)
  pathway_priority <- as.data.frame(pathway_priority)

  if (!"name" %in% colnames(gene_priority) || !"PriorityScore" %in% colnames(gene_priority)) {
    stop("gene_priority must contain at least 'name' and 'PriorityScore' columns")
  }
  if (!"name" %in% colnames(pathway_priority) || !"PriorityScore" %in% colnames(pathway_priority)) {
    stop("pathway_priority must contain at least 'name' and 'PriorityScore' columns")
  }

  gene_priority$candidate_class <- "Gene"
  pathway_priority$candidate_class <- "Pathway"

  merged_priority <- rbind(pathway_priority, gene_priority)
  merged_priority <- merged_priority[order(-merged_priority$PriorityScore, merged_priority$name), , drop = FALSE]
  rownames(merged_priority) <- NULL

  if (!is.null(top_n)) {
    top_n <- min(as.integer(top_n), nrow(merged_priority))
    merged_priority <- utils::head(merged_priority, top_n)
  }

  merged_priority$display_name <- ifelse(
    merged_priority$candidate_class == "Pathway",
    gsub(pathway_prefix, "", merged_priority$name),
    merged_priority$name
  )
  merged_priority$display_name <- gsub("_", " ", merged_priority$display_name)
  merged_priority$display_name <- ifelse(
    nchar(merged_priority$display_name) > 34,
    paste0(substr(merged_priority$display_name, 1, 31), "..."),
    merged_priority$display_name
  )

  merged_priority$rank <- seq_len(nrow(merged_priority))
  merged_priority$n_selected_sources <- if ("n_selected_sources" %in% colnames(merged_priority)) {
    merged_priority$n_selected_sources
  } else {
    rowSums(cbind(
      if ("selected_cellsets" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_cellsets), FALSE, merged_priority$selected_cellsets) else FALSE,
      if ("selected_pdopdx" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_pdopdx), FALSE, merged_priority$selected_pdopdx) else FALSE,
      if ("selected_clinical" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_clinical), FALSE, merged_priority$selected_clinical) else FALSE
    ))
  }
  merged_priority$S_preclinical <- rowMeans(cbind(
    if ("S_cellsets" %in% colnames(merged_priority)) merged_priority$S_cellsets else NA_real_,
    if ("S_pdopdx" %in% colnames(merged_priority)) merged_priority$S_pdopdx else NA_real_
  ), na.rm = TRUE)
  merged_priority$S_preclinical[!is.finite(merged_priority$S_preclinical)] <- 0

  merged_priority$direction_pattern <- apply(
    cbind(
      if ("dir_cell" %in% colnames(merged_priority)) merged_priority$dir_cell else NA_character_,
      if ("dir_pdopdx" %in% colnames(merged_priority)) merged_priority$dir_pdopdx else NA_character_,
      if ("dir_clin" %in% colnames(merged_priority)) merged_priority$dir_clin else NA_character_
    ),
    1,
    function(x) paste(ifelse(is.na(x) | x == "", "NA", x), collapse = " / ")
  )

  merged_priority$display_name <- factor(
    merged_priority$display_name,
    levels = rev(merged_priority$display_name)
  )

  merged_priority
}

#' Plot top priority candidates as grouped horizontal bars
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPriorityTopBar <- function(priority_df,
                               title = "Top prioritized candidates") {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "PriorityScore", "candidate_class", "n_selected_sources")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }

  ggplot2::ggplot(
    priority_df,
    ggplot2::aes(x = PriorityScore, y = display_name, fill = candidate_class)
  ) +
    ggplot2::geom_col(width = 0.72, alpha = 0.92) +
    ggplot2::geom_point(
      ggplot2::aes(shape = factor(n_selected_sources)),
      x = 0.02,
      size = 2.6,
      stroke = 0.6,
      fill = "white",
      color = "#2F2A24"
    ) +
    ggplot2::scale_fill_manual(values = c("Gene" = "#E98371", "Pathway" = "#66A9C9")) +
    ggplot2::scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24), name = "Selected sources") +
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
}

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

  ggplot2::ggplot(
    stacked_df,
    ggplot2::aes(x = value, y = display_name, fill = component)
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
  heatmap_df$source <- factor(heatmap_df$source, levels = c("Cell sets", "PDO/PDX", "Clinical"))
  heatmap_df$direction_symbol <- dplyr::case_when(
    heatmap_df$direction == "Up" ~ "U",
    heatmap_df$direction == "Down" ~ "D",
    TRUE ~ "."
  )

  ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = source, y = display_name, fill = score)
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

#' Plot clinical versus preclinical evidence as a bubble chart
#'
#' @param priority_df Data frame prepared by \code{preparePriorityVisualizationData()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPriorityClinicalPreclinicalBubble <- function(priority_df,
                                                  title = "Clinical vs preclinical support") {
  priority_df <- as.data.frame(priority_df)
  required_cols <- c("display_name", "S_preclinical", "S_clinical", "PriorityScore", "candidate_class", "n_selected_sources")
  if (!all(required_cols %in% colnames(priority_df))) {
    stop("priority_df must contain columns: ", paste(required_cols, collapse = ", "))
  }

  ggplot2::ggplot(
    priority_df,
    ggplot2::aes(
      x = S_preclinical,
      y = S_clinical,
      size = PriorityScore,
      fill = candidate_class,
      alpha = n_selected_sources
    )
  ) +
    ggplot2::geom_point(shape = 21, color = "#2F2A24", stroke = 0.4) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = as.character(display_name)),
      size = 3.5,
      color = "black",
      box.padding = 0.45,
      point.padding = 0.25,
      max.overlaps = 50
    ) +
    ggplot2::scale_fill_manual(values = c("Gene" = "#E98371", "Pathway" = "#66A9C9")) +
    ggplot2::scale_size_continuous(range = c(3, 10), name = "Priority score") +
    ggplot2::scale_alpha_continuous(range = c(0.65, 1), guide = "none") +
    ggplot2::labs(
      title = title,
      x = "Mean preclinical evidence score",
      y = "Clinical evidence score",
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
}
