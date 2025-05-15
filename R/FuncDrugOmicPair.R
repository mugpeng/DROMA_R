# Con ----
# Function to pair omics and drug data
pairDrugOmic <- function(myOmics, myDrugs, merged = FALSE){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      omic_sel <- na.omit(omic_sel); drug_sel <- na.omit(drug_sel)
      if(length(na.omit(omic_sel)) == 0 | length(na.omit(drug_sel)) == 0){ return(NULL) }
      intersected_cells <- intersect(names(omic_sel), names(drug_sel))
      omic_sel <- omic_sel[match(intersected_cells, names(omic_sel))]
      drug_sel <- drug_sel[match(intersected_cells, names(drug_sel))]
      if(length(na.omit(omic_sel)) < 3 | length(na.omit(drug_sel)) < 3){ return(NULL) }
      list("omic" = omic_sel,
           "drug" = drug_sel)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  if(length(pair_list2) < 1) {stop("Please try to another drug-omic pair. This pair do not have result.")}
  # If merged is TRUE and z-score normalization is enabled, create a merged dataset
  if(merged & length(pair_list2) > 1) {
    # Create a merged dataset by combining all pairs using the more efficient approach
    combined_list <- list(
      # Combine all omic vectors
      unlist(lapply(pair_list2, function(x) x$omic)),
      
      # Combine all drug vectors
      unlist(lapply(pair_list2, function(x) x$drug))
    )
    pair_list2[["merged_dataset"]] <- list(
      "omic" = combined_list[[1]],
      "drug" = combined_list[[2]]
    )
  }
  pair_list2
}

# Perform meta-analysis on continuous drug-omic pairs
analyze_continuous_drugomic <- function(myPairs) {
  # Initialize list to store correlation results
  test_list <- list()
  valid_indices <- c()
  
  # Analyze each pair
  for (x in seq_along(myPairs)) {
    # Skip merged dataset for meta-analysis
    if (names(myPairs)[x] == "merged_dataset") next
    
    # Try to perform correlation test
    tryCatch({
      omic_sel <- myPairs[[x]]$omic
      drug_sel <- myPairs[[x]]$drug
      
      # Check for valid data
      if (length(omic_sel) < 3 || length(drug_sel) < 3) next
      
      # Perform correlation test
      cor_re <- cor.test(omic_sel, drug_sel, method = "spearman")
      test_list[[x]] <- data.frame(
        p = cor_re$p.value,
        effect = cor_re$estimate,
        N = length(omic_sel)
      )
      # Track valid indices for proper study name mapping
      valid_indices <- c(valid_indices, x)
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Return NULL if no correlations could be calculated
  if (length(test_list) < 1) return(NULL)
  
  # Prepare data for meta-analysis
  meta_df <- do.call(rbind, test_list)
  # Use valid_indices to correctly map study names
  meta_df$study <- names(myPairs)[valid_indices]
  meta_df$se <- sqrt((1 - meta_df$effect^2) / (meta_df$N - 2))
  meta_df$z <- 0.5 * log((1 + meta_df$effect) / (1 - meta_df$effect))  # Fisher's z 
  meta_df$se_z <- 1 / sqrt(meta_df$N - 3)         
  
  # Perform meta-analysis
  tryCatch({
    # Only perform meta-analysis if we have at least 2 studies
    if (nrow(meta_df) >= 2) {
      cal_meta_re <- metagen(TE = z, 
                             seTE = se_z, 
                             data = meta_df, 
                             sm = "Z",
                             studlab = study)
      return(cal_meta_re)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

# Discrete ----
# Function to pair discrete omics and drug data
pairDrugOmic2 <- function(myOmics, myDrugs, merged = FALSE){
  pair_list2 <- lapply(1:length(myOmics), function(x){
    omic_sel <- myOmics[[x]]
    pair_list <- lapply(1:length(myDrugs), function(y){
      drug_sel <- myDrugs[[y]]
      yes_drugs <- na.omit(drug_sel[names(drug_sel) %in% omic_sel])
      no_drugs <- na.omit(drug_sel[!names(drug_sel) %in% omic_sel])
      if(length(yes_drugs) < 3 | length(no_drugs) < 3){ 
        return(NULL)
      }
      list(yes = yes_drugs,
           no = no_drugs)
    })
    names(pair_list) <- paste0(names(myOmics)[x], "_",
                               names(myDrugs))
    pair_list
  })
  pair_list2 <- unlist(pair_list2, recursive = F)
  pair_list2 <- pair_list2[!sapply(pair_list2, is.null)]
  if(length(pair_list2) < 1){stop("Please try to another drug-omic pair. This pair do not have result.")}
  # If merged is TRUE and z-score normalization is enabled, create a merged dataset
  if(merged & length(pair_list2) > 1) {
    # Create a merged dataset by combining all pairs using the more efficient approach
    combined_list <- list(
      # Combine all yes vectors
      unlist(lapply(pair_list2, function(x) x$yes)),
      
      # Combine all no vectors
      unlist(lapply(pair_list2, function(x) x$no))
    )
    
    # Only add merged dataset if we have enough data points
    pair_list2[["merged_dataset"]] <- list(
      "yes" = combined_list[[1]],
      "no" = combined_list[[2]]
    )
  }
  pair_list2
}

# Perform meta-analysis on discrete drug-omic pairs
analyze_discrete_drugomic <- function(myPairs) {
  # Initialize list to store test results
  test_list <- list()
  valid_indices <- c()

  # Analyze each pair
  for (x in seq_along(myPairs)) {
    # Skip merged dataset for meta-analysis
    if (names(myPairs)[x] == "merged_dataset") next
    
    # Try to perform Wilcoxon test and effect size calculation
    tryCatch({
      yes_drugs <- myPairs[[x]]$yes
      no_drugs <- myPairs[[x]]$no
      
      # Check for valid data
      if (length(yes_drugs) < 3 || length(no_drugs) < 3) next
      
      # Perform statistical test
      wilcox_test <- wilcox.test(no_drugs, yes_drugs)
      cliff_delta <- cliff.delta(no_drugs, yes_drugs)
      
      # Store results
      test_list[[x]] <- data.frame(
        p = wilcox_test$p.value,
        effect = cliff_delta$estimate,
        N = length(yes_drugs) + length(no_drugs),
        n1 = length(yes_drugs),
        n2 = length(no_drugs)
      )
      valid_indices <- c(valid_indices, x)
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Return NULL if no tests could be performed
  if (length(test_list) < 1) return(NULL)
  
  # Prepare data for meta-analysis
  meta_df <- do.call(rbind, test_list)
  meta_df$study <- names(myPairs)[valid_indices]
  
  # Calculate standard error for Cliff's Delta
  meta_df$se <- sqrt((1 - meta_df$effect^2) * (meta_df$n1 + meta_df$n2 + 1) / 
                       (12 * meta_df$n1 * meta_df$n2))
  
  # Perform meta-analysis
  tryCatch({
    # Only perform meta-analysis if we have at least 2 studies
    if (nrow(meta_df) >= 2) {
      meta_result <- metagen(TE = effect, 
                             seTE = se, 
                             data = meta_df,
                             sm = "CMD",  # Custom Mean Difference (using Cliff's Delta)
                             studlab = study)
      return(meta_result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

# Plot ----
# Create a standardized forest plot for meta-analysis results
create_forest_plot <- function(meta_obj, 
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

# Plot continuous drug-omic correlation for a single study
plot_continuous_drugomic <- function(omic_values, drug_values, study_name) {
  # Combine data into dataframe
  cor_df <- data.frame(
    genes = omic_values,
    drugs = drug_values
  )
  
  # Create scatter plot with correlation statistics
  ggscatter(cor_df, x = "genes", y = "drugs", alpha = 0.2) +
    stat_cor(size = 6, method = "spearman") + 
    stat_smooth(formula = y ~ x, method = "lm") + 
    theme_bw() +
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12)
    ) + 
    ggtitle(study_name)
}

# Create plots for all continuous drug-omic pairs
plot_all_continuous_drugomic <- function(pairs_list) {
  # Initialize list to store plots
  p_list <- list()
  
  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    
    # Try to create the plot, continue if error
    tryCatch({
      omic_sel <- pairs_list[[i]]$omic
      drug_sel <- pairs_list[[i]]$drug
      
      # Ensure adequate data for plotting
      if (length(omic_sel) < 3 || length(drug_sel) < 3) next
      
      # Create plot and add to list
      p_list[[i]] <- plot_continuous_drugomic(omic_sel, drug_sel, names(pairs_list)[i])
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]
  
  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    return(wrap_plots(p_list, ncol = 3))
  } else {
    return(NULL)
  }
}

# Plot discrete drug-omic comparison for a single study
plot_discrete_drugomic <- function(yes_values, no_values, study_name) {
  # Combine data into dataframe
  box_df <- data.frame(
    drugs = c(no_values, yes_values),
    events = rep(c("no", "yes"), times = c(length(no_values), length(yes_values)))
  )
  
  # Create boxplot with statistical test
  ggboxplot(data = box_df, x = "events", y = "drugs",
            fill = "events", palette = c("#BEBADAFF", "#FB8072FF"),
            add = "jitter", add.params = list(alpha = 0.15)) + 
    stat_compare_means(size = 6, label.x = 0.8,
                       label.y = (max(box_df$drugs) - max(box_df$drugs)/8),
                       label = "p.format") + 
    theme_bw() + 
    theme(
      axis.title = element_blank(),
      title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "none"
    ) + 
    coord_cartesian(ylim = c(NA, max(box_df$drugs) + max(box_df$drugs)/20)) + 
    ggtitle(study_name)
}

# Create plots for all discrete drug-omic pairs
plot_all_discrete_drugomic <- function(pairs_list) {
  # Initialize list to store plots
  p_list <- list()
  
  # Create plot for each pair
  for (i in seq_along(pairs_list)) {
    
    # Try to create the plot, continue if error
    tryCatch({
      yes_drugs <- pairs_list[[i]]$yes
      no_drugs <- pairs_list[[i]]$no
      
      # Ensure adequate data for plotting
      if (length(yes_drugs) < 3 || length(no_drugs) < 3) next
      
      # Create plot and add to list
      p_list[[i]] <- plot_discrete_drugomic(yes_drugs, no_drugs, names(pairs_list)[i])
    }, error = function(e) {
      # Continue to next pair on error
    })
  }
  
  # Remove NULL entries from list
  p_list <- p_list[!sapply(p_list, is.null)]
  
  # Combine plots using patchwork if plots exist
  if (length(p_list) > 0) {
    return(wrap_plots(p_list, ncol = 3))
  } else {
    return(NULL)
  }
}

create_plot_with_common_axes <- function(p, x_title = "Common X-Axis Title", 
                                         y_title = "Common Y-Axis Title") {
  # Create a function that will generate the plot when called
  function() {
    # Convert patchwork to a grob
    p_grob <- patchworkGrob(p)
    
    # Create a new plotting area
    grid.newpage()
    
    # Draw the patchwork
    grid.draw(p_grob)
    
    # Add common x-axis title
    grid.text(x_title,
              x = 0.5, y = 0.02,
              gp = gpar(fontsize = 18, fontface = "bold"))
    
    # Add common y-axis title (rotated)
    grid.text(y_title,
              x = 0.01, y = 0.5,
              rot = 90,
              gp = gpar(fontsize = 18, fontface = "bold"))
    
    # Return the grob for potential further use
    invisible(p_grob)
  }
}

# All ----
# Function to handle both continuous and discrete omics data in 
# one function
oneDrugOmicPair <- function(select_omics_type, select_omics,
                            select_drugs,
                            data_type = "all", tumor_type = "all",
                            merged_enabled = TRUE,
                            meta_enabled = TRUE){
  # Get drug data
  myDrugs <- selFeatures("drug", select_drugs, 
                        data_type = data_type, 
                        tumor_type = tumor_type)
  myOmics <- selFeatures(select_omics_type, select_omics,
                         data_type = data_type, 
                         tumor_type = tumor_type)
  # Initialize result list
  result <- list()
  
  # Handle continuous omics data
  if(select_omics_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")){
    # Pair data
    myPairs <- pairDrugOmic(myOmics, myDrugs, merged = merged_enabled)
    
    # Create plots
    result$plot <- plot_all_continuous_drugomic(myPairs)
    
    # Perform meta-analysis
    if(meta_enabled){
      meta_result <- analyze_continuous_drugomic(myPairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    # Store data
    result$data <- myPairs
  } else {
    # Pair data
    myPairs <- pairDrugOmic2(myOmics, myDrugs, merged = merged_enabled)
    
    # Create plots
    result$plot <- plot_all_discrete_drugomic(myPairs)
    
    # Perform meta-analysis
    if(meta_enabled){
      meta_result <- analyze_discrete_drugomic(myPairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    # Store data
    result$data <- myPairs
  }
  
  # Return results if there's a plot
  if (is.null(result$plot)) {
    return(list())
  } else {
    return(result)
  }
}

