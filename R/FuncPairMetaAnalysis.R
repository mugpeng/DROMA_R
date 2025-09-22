# Meta-analysis Functions Module ----

#' Meta-analysis for continuous vs continuous features
#'
#' @description Performs meta-analysis on correlation between two continuous features
#' @param selected_pair List of paired data with two continuous variables
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcConCon <- function(selected_pair) {
  if (length(selected_pair) < 1) return(NULL)

  # Test pairs one by one
  cal_list <- lapply(1:length(selected_pair), function(y) {
    fea1_sel <- selected_pair[[y]][[1]]
    fea2_sel <- selected_pair[[y]][[2]]

    # Check for minimum length
    if (length(fea1_sel) < 3 || length(fea2_sel) < 3) return(NULL)

    options(warn = -1)
    cor_re <- tryCatch(
      cor.test(fea1_sel, fea2_sel, method = "spearman"),
      error = function(x) { return(NULL) }
    )

    if (is.null(cor_re)) return(NULL)

    data.frame(
      p = cor_re$p.value,
      effect = cor_re$estimate,
      N = length(fea2_sel)
    )
  })

  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)
  
  # Handle extreme correlation values to prevent infinite z-scores
  # Clamp correlation coefficients to avoid numerical issues in Fisher's z transform
  cal_re$effect <- pmax(pmin(cal_re$effect, 0.999), -0.999)
  
  cal_re$se <- sqrt((1 - cal_re$effect^2) / (cal_re$N - 2))
  cal_re$z <- 0.5 * log((1 + cal_re$effect) / (1 - cal_re$effect))  # Fisher's z
  cal_re$se_z <- 1 / sqrt(cal_re$N - 3)

  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = z, seTE = se_z, data = cal_re, sm = "Z",
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )}),
    error = function(x) { return(NULL) }
  )
  cal_meta_re
}

#' Meta-analysis for continuous vs discrete features
#'
#' @description Performs meta-analysis comparing continuous values between discrete groups
#' @param selected_pair List of paired data with continuous and discrete variables
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcConDis <- function(selected_pair) {
  options(warn = -1)
  if (length(selected_pair) < 1) return(NULL)

  cal_list <- lapply(1:length(selected_pair), function(y) {
    yes_drugs <- selected_pair[[y]][[1]]
    no_drugs <- selected_pair[[y]][[2]]

    # Check for minimum length
    if (length(yes_drugs) < 3 || length(no_drugs) < 3) return(NULL)

    wilcox_re <- tryCatch(
      wilcox.test(no_drugs, yes_drugs),
      error = function(x) { return(NULL) }
    )
    if (is.null(wilcox_re)) return(NULL)

    cliff_delta <- tryCatch(
      cliff.delta(no_drugs, yes_drugs),
      error = function(x) { return(NULL) }
    )
    if (is.null(cliff_delta)) return(NULL)

    data.frame(
      p = wilcox_re$p.value,
      effect = cliff_delta$estimate,
      N = length(yes_drugs) + length(no_drugs),
      n1 = length(yes_drugs),
      n2 = length(no_drugs)
    )
  })

  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)
  # Calculate standard error for Cliff's Delta
  cal_re$se <- sqrt((1 - cal_re$effect^2) * (cal_re$n1 + cal_re$n2 + 1) /
                      (12 * cal_re$n1 * cal_re$n2))

  cal_meta_re <- tryCatch(
    suppressWarnings({
      meta_result <- metagen(TE = effect,
                            seTE = se,
                            data = cal_re,
                            control = list(maxiter = 2000,
                                           stepadj = 0.1,
                                           threshold = 0.000001),
                            sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
      )
    }),
    error = function(x) { return(NULL) }
  )
  cal_meta_re
}

#' Meta-analysis for discrete vs discrete features
#'
#' @description Performs meta-analysis on two discrete features using odds ratio
#' @param selected_pair List of paired data with contingency tables
#' @return Meta-analysis result object or NULL if insufficient data
#' @export
metaCalcDisDis <- function(selected_pair) {
  # Check if we have enough pairs for meta-analysis
  if (length(selected_pair) < 1) return(NULL)

  # Calculate statistics for each pair
  cal_list <- lapply(1:length(selected_pair), function(y) {
    cont_table <- selected_pair[[y]]$cont_table

    # Skip if any cell has too few observations (e.g., < 3)
    if (any(cont_table < 3)) return(NULL)

    # Calculate odds ratio and its standard error
    tryCatch({
      # Extract values from contingency table
      a <- cont_table[1,1] # yes-yes
      b <- cont_table[1,2] # yes-no
      c <- cont_table[2,1] # no-yes
      d <- cont_table[2,2] # no-no

      # Calculate log odds ratio and its standard error
      log_or <- log((a * d)/(b * c))
      se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)

      # Calculate Fisher's exact test p-value
      fisher_test <- fisher.test(cont_table)

      data.frame(
        log_or = log_or,
        se = se_log_or,
        p = fisher_test$p.value,
        N = sum(cont_table)
      )
    }, error = function(x) NULL)
  })

  # Remove NULL results and combine
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)

  # Perform meta-analysis using random effects model
  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = log_or,
              seTE = se,
              data = cal_re,
              sm = "OR", # Specify odds ratio as summary measure
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )
    }),
    error = function(x) NULL
  )

  cal_meta_re
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

#' Create a volcano plot from meta-analysis results
#'
#' @param meta_df Data frame containing meta-analysis results with columns:
#'                effect_size, p_value, and name
#' @param es_t Effect size threshold to consider significant
#' @param P_t P-value threshold to consider significant
#' @param label Whether to add labels to top points (TRUE/FALSE)
#' @param top_label_each Number of top points in each direction to label
#' @param label_size Size of text labels
#' @param point_size Size of points
#' @param point_alpha Alpha transparency of points
#' @param title Plot title (NULL for no title)
#' @param p_adj_method Method for p-value adjustment ("none", "BH", "bonferroni")
#' @param custom_colors Custom color vector for Up, NS, Down (NULL for defaults)
#' @return ggplot object with volcano plot
#' @export
plotMetaVolcano <- function(meta_df,
                            es_t = .4,
                            P_t = .001,
                            label = TRUE,
                            top_label_each = 5,
                            label_size = 5,
                            point_size = 2.5,
                            point_alpha = 0.6,
                            title = NULL,
                            p_adj_method = "none",
                            custom_colors = NULL) {

  # Input validation
  if(!is.data.frame(meta_df)) stop("meta_df must be a data frame")
  if(!all(c("effect_size", "p_value", "name") %in% colnames(meta_df))) {
    stop("meta_df must contain columns: effect_size, p_value, and name")
  }

  # Handle p-value adjustment if requested
  if(p_adj_method != "none") {
    meta_df$p_value <- p.adjust(meta_df$p_value, method = p_adj_method)
  }

  # Default colors
  if(is.null(custom_colors)) {
    custom_colors <- c("Down" = "#44bce4", "NS" = "grey", "Up" = "#fc7474")
  }

  # Group the points based on thresholds
  meta_df$group <- dplyr::case_when(
    meta_df$effect_size > es_t & meta_df$p_value < P_t ~ "Up",
    meta_df$effect_size < -es_t & meta_df$p_value < P_t ~ "Down",
    TRUE ~ "NS"
  )

  # Count significant findings
  sig_counts <- table(meta_df$group)
  sig_text <- paste0(
    "Up: ", sum(meta_df$group == "Up"), ", ",
    "Down: ", sum(meta_df$group == "Down"), ", ",
    "Total: ", nrow(meta_df)
  )

  # Basic volcano plot
  p <- ggplot(data = meta_df,
              aes(x = effect_size,
                  y = -log10(p_value))) +
    geom_point(size = point_size, alpha = point_alpha,
               aes(color = group)) +
    theme_bw() +
    theme(
      legend.position = "none",
      title = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 15, colour = "black"),
      axis.text = element_text(size = 15, color = "black"),
      legend.title = element_text(size = 15, colour = "black"),
      legend.text = element_text(size = 15),
      text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black")
    ) +
    ylab("-log10(Pvalue)") +
    xlab("Effect Size") +
    scale_color_manual(values = custom_colors) +
    geom_vline(xintercept = c(-es_t, es_t), lty = 4, col = "black", lwd = 0.5) +
    geom_hline(yintercept = -log10(P_t), lty = 4, col = "black", lwd = 0.5) +
    annotate("text", x = min(meta_df$effect_size, na.rm = TRUE) * 0.8,
             y = max(-log10(meta_df$p_value), na.rm = TRUE) * 0.9,
             label = sig_text, hjust = 0, size = 5)

  # Add title if provided
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }

  # Add labels if requested
  if(label) {
    meta_df2 <- meta_df[meta_df$group != "NS",]

    # Skip labeling if there are no significant points
    if(nrow(meta_df2) > 0) {
      # Get top points to label
      low_indices <- head(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      high_indices <- tail(order(meta_df2$effect_size), min(top_label_each, nrow(meta_df2)))
      forlabel_names <- c(meta_df2$name[low_indices], meta_df2$name[high_indices])
      forlabel_df <- meta_df2[meta_df2$name %in% forlabel_names,]

      p <- p +
        geom_point(size = point_size + 0.5, shape = 1, data = forlabel_df) +
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