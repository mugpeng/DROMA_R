# Enhanced Statistical Plotting Functions for MultiDromaSet ----

# Note: This file requires the following packages for full functionality:
# - ggplot2 (for plotting)
# - dplyr (for case_when function)
# - UpSetR (for overlap plots)
# - treemapify (for MOA treemap plots)
# - ggrepel (for tumor type bubble plots)
# - patchwork or gridExtra (for dashboard layout)

# Import case_when from dplyr if available
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required for this function. Please install with install.packages('dplyr')")
} else {
  # Import case_when from dplyr
  case_when <- dplyr::case_when
}

#' Generate comprehensive statistical plots from DROMA database projects
#'
#' @description Creates statistical overview plots showing drug/sample counts, overlaps,
#' molecular characteristics, and tumor type distributions from projects in the DROMA database
#' @param project_names Character vector of project names to include. Use "all" to include all projects from database
#' @param plot_types Character vector of plot types to generate: "counts", "overlaps",
#' "molecular", "drug_moa", "tumor_types", or "all"
#' @param data_type Character, filter by data type: "all" (default), "CellLine", "PDO", "PDC", or "PDX"
#' @param connection Optional database connection object. If NULL, uses global connection
#' @param use_gap_plots Logical, whether to use gap plots for large count differences
#' @return List of ggplot objects with statistical plots
#' @export
#' @examples
#' \dontrun{
#' # Connect to database
#' con <- connectDROMADatabase("path/to/droma.sqlite")
#'
#' # Generate plots for all projects
#' stat_plots <- generateStatisticalPlots("all", plot_types = "all")
#'
#' # Generate plots for specific projects
#' stat_plots <- generateStatisticalPlots(c("gCSI", "CCLE"), plot_types = "counts")
#'
#' # Generate plots for cell line projects only
#' stat_plots <- generateStatisticalPlots("all", data_type = "CellLine")
#' }
generateStatisticalPlots <- function(project_names = "all",
                                   plot_types = "all",
                                   data_type = "all",
                                   connection = NULL,
                                   use_gap_plots = TRUE) {

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Standardize plot types
  if ("all" %in% plot_types) {
    plot_types <- c("counts", "overlaps", "molecular", "drug_moa", "tumor_types")
  }

  # Get project information from database
  all_projects <- DROMA.Set::listDROMAProjects(connection = connection)
  if (nrow(all_projects) == 0) {
    stop("No projects found in database")
  }

  # Filter projects by data_type if specified
  if (data_type != "all") {
    all_projects <- all_projects[all_projects$dataset_type == data_type, ]
    if (nrow(all_projects) == 0) {
      stop("No projects found with data_type: ", data_type)
    }
  }

  # Determine which projects to use
  if ("all" %in% project_names) {
    selected_projects <- all_projects$project_name
  } else {
    # Validate provided project names
    invalid_projects <- setdiff(project_names, all_projects$project_name)
    if (length(invalid_projects) > 0) {
      warning("Invalid project names: ", paste(invalid_projects, collapse = ", "))
    }
    selected_projects <- intersect(project_names, all_projects$project_name)
    if (length(selected_projects) == 0) {
      stop("No valid projects found")
    }
  }

  # Filter project info to selected projects
  project_info <- all_projects[all_projects$project_name %in% selected_projects, ]

  # Get annotations from database
  sample_annotations <- NULL
  drug_annotations <- NULL

  sample_annotations <- tryCatch({
    DROMA.Set::getDROMAAnnotation("sample", connection = connection)
  }, error = function(e) {
    warning("Could not load sample annotations from database: ", e$message)
    NULL
  })

  drug_annotations <- tryCatch({
    DROMA.Set::getDROMAAnnotation("drug", connection = connection)
  }, error = function(e) {
    warning("Could not load drug annotations from database: ", e$message)
    NULL
  })

  # Initialize results list
  plot_results <- list()

  # Generate requested plots
  if ("counts" %in% plot_types) {
    plot_results$counts <- generateCountPlots(project_info, use_gap_plots)
  }

  if ("overlaps" %in% plot_types) {
    plot_results$overlaps <- generateOverlapPlots(selected_projects, connection)
  }

  if ("molecular" %in% plot_types) {
    plot_results$molecular <- generateCharacteristicsPlot(project_info)
  }

  if ("drug_moa" %in% plot_types && !is.null(drug_annotations)) {
    plot_results$drug_moa <- generateDrugMOAPlot(drug_annotations)
  }

  if ("tumor_types" %in% plot_types && !is.null(sample_annotations)) {
    plot_results$tumor_types <- generateTumorTypePlot(sample_annotations)
  }

  return(plot_results)
}

#' Generate count plots for drugs and samples
#'
#' @description Creates bar plots showing drug and sample counts across projects using database info
#' @param project_info Data frame with project information from database
#' @param use_gap_plots Logical, whether to use gap plots
#' @return List of ggplot objects for count visualizations
generateCountPlots <- function(project_info, use_gap_plots = TRUE) {

  # Prepare data for plotting using database information
  all_stat <- data.frame(
    counts = numeric(),
    source = character(),
    type = character(),
    group = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(project_info)) {
    project_row <- project_info[i, ]
    project_name <- project_row$project_name
    dataset_type <- if(is.na(project_row$dataset_type)) "Other" else project_row$dataset_type
    sample_count <- if(is.na(project_row$sample_count)) 0 else project_row$sample_count
    drug_count <- if(is.na(project_row$drug_count)) 0 else project_row$drug_count

    # Add dimensions to data frame
    all_stat <- rbind(all_stat, data.frame(
      counts = c(drug_count, sample_count),
      source = rep(project_name, 2),
      type = c("Drugs", "Samples"),
      group = rep(dataset_type, 2),
      stringsAsFactors = FALSE
    ))
  }

  # Calculate totals and create ordering
  source_totals <- aggregate(counts ~ source + group, data = all_stat, FUN = sum)

  # Ensure all expected group levels are present, ordered by typical hierarchy
  expected_groups <- c("CellLine", "PDC", "PDO", "PDX", "Other")
  present_groups <- intersect(expected_groups, unique(all_stat$group))
  all_stat$group <- factor(all_stat$group, levels = present_groups)

  # Create rank for ordering within groups
  source_order <- source_totals[order(source_totals$group, source_totals$counts),]
  source_order$rank <- ave(source_order$counts, source_order$group,
                          FUN = function(x) rank(x, ties.method = "first"))

  # Join rank back to original data
  all_stat <- merge(all_stat, source_order[,c("source", "rank")], by = "source")
  all_stat <- all_stat[order(all_stat$group, all_stat$rank, all_stat$type),]
  all_stat$source <- factor(all_stat$source, levels = unique(all_stat$source))

  # Add alternating background
  all_stat$group_id <- as.numeric(all_stat$group)
  all_stat$group_id <- all_stat$group_id %% 2

  # Create main count plot without labels (for gg.gap)
  p_count_no_labels <- ggplot(all_stat, aes(x = source, y = counts, fill = type)) +
    # Add alternating background
    geom_rect(data = subset(all_stat, group_id == 1),
              aes(xmin = as.numeric(source) - 0.5,
                  xmax = as.numeric(source) + 0.5,
                  ymin = -Inf, ymax = Inf),
              fill = "gray90", alpha = 0.3,
              inherit.aes = FALSE) +
    geom_col(position = "dodge") +
    geom_text(aes(label = counts), size = 4.5,
              position = position_dodge(0.9),
              vjust = -0.5) +
    scale_fill_manual(values = c("Drugs" = "#F7766D", "Samples" = "#31BFC4")) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      title = element_text(size = 17, color = "black"),
      axis.title.y = element_text(size = 17, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
      legend.title = element_text(size = 17, colour = "black"),
      legend.text = element_text(size = 17),
      legend.position = "none"
    ) +
    ylab("Count")
  p_count_main <- p_count_no_labels
  # Apply gg.gap if requested and available
  if(use_gap_plots){
    p_count_no_labels_gap <- gg.gap(
      plot = p_count_no_labels,
      ylim = c(0, max(all_stat$counts) * 1.05),
      segments = list(c(2000, 54000)),
      tick_width = c(250, 1000),
      rel_heights = c(0.7, 0.05, 0.25)
    )
    p_count_main <- p_count_no_labels_gap
  }
  # Create summary plot
  total_samples <- sum(project_info$sample_count, na.rm = TRUE)
  total_drugs <- sum(project_info$drug_count, na.rm = TRUE)

  summary_data <- data.frame(
    counts = c(total_samples, total_drugs),
    type = c("Samples", "Drugs")
  )

  # Create summary plot with legend (for extracting legend)
  p_count_summary_with_legend <- ggplot(summary_data, aes(x = type, y = counts, fill = type)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = counts), size = 5,
              position = position_dodge(0.9),
              vjust = -0.8) +
    scale_fill_manual(values = c("Drugs" = "#F7766D", "Samples" = "#31BFC4")) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      title = element_text(size = 17, color = "black"),
      axis.title.y = element_text(size = 17, color = "black"),
      axis.text = element_text(size = 17, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 17)
    ) +
    ylab("Count") +
    ggtitle("")

  # Extract legend
  legend <- NULL
  tryCatch({
    if (requireNamespace("cowplot", quietly = TRUE)) {
      legend <- cowplot::get_legend(p_count_summary_with_legend)
    }
  }, error = function(e) {
    warning("cowplot package not available for legend extraction")
  })

  # Create summary plot without legend
  p_count_summary_no_legend <- p_count_summary_with_legend +
    theme(legend.position = "none")
  p_count_summary <- p_count_summary_no_legend
  # Apply gg.gap to summary plot if available and appropriate
  if (use_gap_plots) {
    p_count_summary_gapped <- gg.gap(
      plot = p_count_summary_no_legend,
      ylim = c(0, 75000),
      segments = list(c(15000, 50000)),
      tick_width = c(4000, 10000),
      rel_heights = c(0.7, 0.05, 0.25)
    )
    p_count_summary <- p_count_summary_gapped
  }

  # Create combined summary plot with legend (mimicking p_count_drugandsample_sum2)
  p_count_summary_combined <- NULL
  if (!is.null(legend)) {
    tryCatch({
      if (requireNamespace("cowplot", quietly = TRUE)) {
        p_count_summary_combined <- cowplot::plot_grid(
          legend,
          p_count_summary_gapped,
          ncol = 1,
          rel_heights = c(0.2, 0.9)
        )
      }
    }, error = function(e) {
      warning("Could not create combined summary plot: ", e$message)
      p_count_summary_combined <- p_count_summary_gapped
    })
  } else {
    p_count_summary_combined <- p_count_summary_gapped
  }

  return(list(
    detailed = p_count_main,                    # Main plot with gap (p_count_drugandsample)
    summary = p_count_summary_combined,         # Summary with legend (p_count_drugandsample_sum2)
    data = all_stat                            # Raw data
  ))
}

#' Generate overlap plots for samples and drugs from database
#'
#' @description Creates UpSet plots showing overlaps between projects using database data
#' @param selected_projects Character vector of project names
#' @param connection Database connection object
#' @return List of UpSet plot objects
generateOverlapPlots <- function(selected_projects, connection) {

  # Create lists of drugs and samples for all datasets
  drug_list <- list()
  sample_list <- list()

  for (project_name in selected_projects) {
    # Get drugs for this project
    tryCatch({
      drugs <- DROMA.Set::listDROMAFeatures(project_name, "drug", connection = connection)
      if (length(drugs) > 0) {
        drug_list[[project_name]] <- drugs
      }
    }, error = function(e) {
      warning("Could not get drugs for project ", project_name, ": ", e$message)
    })

    # Get samples for this project
    tryCatch({
      samples <- DROMA.Set::listDROMASamples(project_name, connection = connection)
      if (length(samples) > 0) {
        sample_list[[project_name]] <- samples
      }
    }, error = function(e) {
      warning("Could not get samples for project ", project_name, ": ", e$message)
    })
  }

  # Create UpSet plots (requires UpSetR package)
  tryCatch({
    if (requireNamespace("UpSetR", quietly = TRUE)) {
      p_overlap_sample <- NULL
      p_overlap_drug <- NULL

      if (length(sample_list) > 1) {
        p_overlap_sample <- UpSetR::upset(
          UpSetR::fromList(sample_list),
          mainbar.y.label = "Sample Counts",
          text.scale = 2,
          nsets = length(sample_list)
        )
      }

      if (length(drug_list) > 1) {
        p_overlap_drug <- UpSetR::upset(
          UpSetR::fromList(drug_list),
          mainbar.y.label = "Drug Counts",
          text.scale = 2,
          nsets = length(drug_list)
        )
      }

      return(list(
        samples = p_overlap_sample,
        drugs = p_overlap_drug,
        sample_overlaps = sample_list,
        drug_overlaps = drug_list
      ))
    } else {
      warning("UpSetR package not available. Install with: install.packages('UpSetR')")
      return(NULL)
    }
  }, error = function(e) {
    warning("Could not create overlap plots: ", e$message)
    return(NULL)
  })
}

#' Generate characteristics plot using database information
#'
#' @description Creates a dot plot showing which data types (including drug and molecular) are available in each project
#' @param project_info Data frame with project information from database
#' @return ggplot object showing data availability
generateCharacteristicsPlot <- function(project_info) {

  # Extract data types for each project from data_types column
  plot_data <- data.frame(
    Project = character(),
    Feature = character(),
    Available = logical(),
    Dataset_Type = character(),
    stringsAsFactors = FALSE
  )

  # All possible features including drug and molecular data types
  all_features <- c("drug", "mRNA", "cnv", "mutation_gene", "mutation_site",
                   "fusion", "proteinrppa", "proteinms", "meth")

  for (i in 1:nrow(project_info)) {
    project_name <- project_info$project_name[i]
    data_types_str <- project_info$data_types[i]
    dataset_type <- if(is.na(project_info$dataset_type[i])) "Other" else project_info$dataset_type[i]

    if (!is.na(data_types_str)) {
      project_data_types <- unlist(strsplit(data_types_str, ","))
      project_data_types <- trimws(project_data_types)
    } else {
      project_data_types <- character(0)
    }

    for (feature in all_features) {
      plot_data <- rbind(plot_data, data.frame(
        Project = project_name,
        Feature = feature,
        Available = feature %in% project_data_types,
        Dataset_Type = dataset_type,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(plot_data) == 0) {
    warning("No data type information available")
    return(NULL)
  }

  # Create nice feature names
  feature_names <- c(
    "drug" = "Drug Response",
    "mRNA" = "Gene Expression",
    "cnv" = "Copy Number",
    "mutation_gene" = "Gene Mutation",
    "mutation_site" = "Site Mutation",
    "fusion" = "Gene Fusion",
    "meth" = "DNA Methylation",
    "proteinrppa" = "Protein RPPA",
    "proteinms" = "Protein MS"
  )

  plot_data$Feature_Display <- feature_names[plot_data$Feature]
  plot_data$Feature_Display[is.na(plot_data$Feature_Display)] <- plot_data$Feature[is.na(plot_data$Feature_Display)]

  # Calculate availability counts for ordering
  # Order features by availability (from most to least available)
  feature_availability <- aggregate(Available ~ Feature_Display,
                                   data = plot_data,
                                   FUN = sum)
  feature_order <- feature_availability$Feature_Display[order(-feature_availability$Available)]

  # Order projects by total availability (from most to least available)
  project_availability <- aggregate(Available ~ Project,
                                   data = plot_data,
                                   FUN = sum)
  project_order <- project_availability$Project[order(-project_availability$Available)]

  # Apply ordering to factors
  plot_data$Feature_Display <- factor(plot_data$Feature_Display, levels = feature_order)
  plot_data$Project <- factor(plot_data$Project, levels = rev(project_order)) # Reverse for y-axis

  # Define colors for dataset types
  dataset_colors <- c(
    "CellLine" = "#E31A1C",    # Red
    "PDC" = "#1F78B4",         # Blue
    "PDO" = "#33A02C",         # Green
    "PDX" = "#FF7F00",         # Orange
    "Unavailable" = "#CCCCCC"        # Gray
  )

  # Create the plot with conditional coloring
  # Available points use dataset type colors, unavailable points use grey80
  plot_data$Point_Color <- ifelse(plot_data$Available,
                                 as.character(plot_data$Dataset_Type),
                                 "Unavailable")

  # Ensure factor levels for Dataset_Type match our color scheme
  # plot_data$Dataset_Type <- factor(plot_data$Dataset_Type,
  #                                 levels = names(dataset_colors))
  # Create the plot
  p_characteristics <- ggplot(plot_data, aes(x = Feature_Display, y = Project)) +
    geom_point(aes(size = Available, color = Point_Color), alpha = 0.8) +
    scale_size_manual(values = c(2, 6), guide = "none") +
    scale_color_manual(values = dataset_colors,
                      name = "Dataset Type",
                      breaks = names(dataset_colors),
                      labels = names(dataset_colors)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 17, color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(x = "", y = "") +
    coord_fixed(ratio = 1)

  return(p_characteristics)
}

#' Generate drug MOA treemap plot from database annotations
#'
#' @description Creates a treemap visualization of drug mechanisms of action using database annotations
#' @param drug_annotations Data frame from getDROMAAnnotation("drug") with drug metadata
#' @return ggplot object of drug MOA treemap
generateDrugMOAPlot <- function(drug_annotations) {

  # Check if drug_annotations is provided and not empty
  if (is.null(drug_annotations) || nrow(drug_annotations) == 0) {
    warning("No drug annotations provided or drug annotations data frame is empty")
    return(NULL)
  }

  # Find the MOA column (check multiple possible names)
  moa_columns <- c("MOA Targets", "MOA", "moa_targets", "moa", "Targets", "targets")
  moa_col <- NULL
  for (col in moa_columns) {
    if (col %in% colnames(drug_annotations)) {
      moa_col <- col
      break
    }
  }

  if (is.null(moa_col)) {
    warning("No valid MOA column found in drug annotations. Expected columns: ",
            paste(moa_columns, collapse = ", "))
    return(NULL)
  }

  # Extract MOA data and process
  moa_data <- drug_annotations[[moa_col]]

  # Remove NA, empty, and "0" values
  moa_data <- moa_data[!is.na(moa_data) & moa_data != "" & moa_data != "0"]

  if (length(moa_data) == 0) {
    warning("No valid MOA data found after filtering")
    return(NULL)
  }

  # Count frequency of each MOA
  moa_counts <- table(moa_data)

  # Split compound MOAs and count individual targets
  all_targets <- character()
  for (moa in names(moa_counts)) {
    # Split by common separators
    targets <- unlist(strsplit(moa, "[;,|]"))
    targets <- trimws(targets)
    targets <- targets[targets != "" & targets != "0"]

    # Add each target the number of times this MOA appears
    all_targets <- c(all_targets, rep(targets, moa_counts[moa]))
  }

  if (length(all_targets) == 0) {
    warning("No valid targets found after processing MOA data")
    return(NULL)
  }

  # Count target frequencies
  target_counts <- table(all_targets)
  target_df <- data.frame(
    Target = names(target_counts),
    Frequency = as.numeric(target_counts),
    stringsAsFactors = FALSE
  )

  # Filter for top targets (more dynamic filtering)
  min_freq <- max(1, ceiling(length(all_targets) * 0.01))  # At least 1% frequency or minimum 1
  target_df <- target_df[target_df$Frequency >= min_freq, ]

  if (nrow(target_df) == 0) {
    warning("No targets meet the minimum frequency threshold")
    return(NULL)
  }

  # Sort by frequency and take top targets if still too many
  target_df <- target_df[order(-target_df$Frequency), ]
  if (nrow(target_df) > 20) {
    target_df <- target_df[1:20, ]
  }

  # Create treemap plot
  tryCatch({
    if (requireNamespace("treemapify", quietly = TRUE)) {
      p_moa <- ggplot(target_df, aes(area = Frequency, fill = Frequency, label = Target)) +
        treemapify::geom_treemap() +
        treemapify::geom_treemap_text(
          aes(label = paste0(Target, "\n(", Frequency, ")")),
          colour = "white", place = "centre",
          reflow = TRUE, grow = FALSE, size = 12) +
        scale_fill_gradient(low = "#56B1F7", high = "#1F3552") +
        theme_void() +
        theme(legend.position = "none",
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12)) +
        labs(fill = "Frequency", title = "Drug Mechanism of Action (MOA) Distribution")

      return(p_moa)
    } else {
      warning("treemapify package not available. Install with: install.packages('treemapify')")

      # Fallback to bar plot
      p_moa_bar <- ggplot(target_df, aes(x = reorder(Target, Frequency), y = Frequency)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        theme_bw() +
        labs(x = "Target", y = "Frequency", title = "Drug MOA Distribution") +
        theme(axis.text = element_text(size = 12))

      return(p_moa_bar)
    }
  }, error = function(e) {
    warning("Could not create MOA plot: ", e$message)
    return(NULL)
  })
}

#' Generate tumor type bubble plot from database annotations
#'
#' @description Creates a bubble plot visualization of tumor types using database sample annotations
#' @param sample_annotations Data frame from getDROMAAnnotation("sample") with sample metadata
#' @return ggplot object of tumor type bubble plot
generateTumorTypePlot <- function(sample_annotations) {

  # Check if sample_annotations is provided and not empty
  if (is.null(sample_annotations) || nrow(sample_annotations) == 0) {
    warning("No sample annotations provided or sample annotations data frame is empty")
    return(NULL)
  }

  # Find the tumor type column (check multiple possible names)
  tumor_type_columns <- c("TumorType", "Tumor_Type", "tumor_type", "TumorTypes",
                         "Cancer_Type", "cancer_type")
  tumor_col <- NULL
  for (col in tumor_type_columns) {
    if (col %in% colnames(sample_annotations)) {
      tumor_col <- col
      break
    }
  }

  if (is.null(tumor_col)) {
    warning("No valid tumor type column found in sample annotations. Expected columns: ",
            paste(tumor_type_columns, collapse = ", "))
    return(NULL)
  }

  # Extract and process tumor type data
  tumor_data_raw <- sample_annotations[[tumor_col]]
  tumor_data_raw <- tumor_data_raw[!is.na(tumor_data_raw) & tumor_data_raw != "" & tumor_data_raw != "Unknown"]

  if (length(tumor_data_raw) == 0) {
    warning("No valid tumor type data found after filtering")
    return(NULL)
  }

  # Create frequency table
  tumor_freq <- table(tumor_data_raw)
  tumor_data <- data.frame(
    TumorType = names(tumor_freq),
    Frequency = as.numeric(tumor_freq),
    stringsAsFactors = FALSE
  )

  # Clean tumor type names for better display
  tumor_data$TumorType <- gsub("_", " ", tumor_data$TumorType)
  tumor_data$TumorType <- tolower(tumor_data$TumorType)

  # Group tumor types by system
  tumor_data$System <- case_when(
    grepl("lung|aerodigestive|nasopharyngeal", tumor_data$TumorType) ~ "Respiratory",
    grepl("gastrointestinal|stomach|liver|pancreatic", tumor_data$TumorType) ~ "Digestive",
    grepl("breast|ovarian|cervical|endometrial|uterine|vulvar", tumor_data$TumorType) ~ "Reproductive - Female",
    grepl("prostate|testicular", tumor_data$TumorType) ~ "Reproductive - Male",
    grepl("haematopoietic|lymphoid", tumor_data$TumorType) ~ "Blood & Lymphatic",
    grepl("nervous", tumor_data$TumorType) ~ "Nervous System",
    grepl("skin", tumor_data$TumorType) ~ "Integumentary",
    grepl("kidney|bladder", tumor_data$TumorType) ~ "Urinary",
    grepl("sarcoma", tumor_data$TumorType) ~ "Connective Tissue",
    TRUE ~ "Other"
  )

  # Shorten names for better display
  tumor_data$TumorType <- gsub(" cancer", "", tumor_data$TumorType)

  # Filter out very small frequencies if there are too many tumor types
  if (nrow(tumor_data) > 20) {
    tumor_data <- tumor_data[tumor_data$Frequency > 1, ]
    if (nrow(tumor_data) > 20) {
      tumor_data <- tumor_data[order(-tumor_data$Frequency)[1:20], ]
    }
  }

  # Create bubble plot
  tryCatch({
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      # Create the bubble chart
      p_tumor_bubble <- ggplot(tumor_data, aes(x = System, y = Frequency, size = Frequency, color = System)) +
        geom_point(alpha = 0.7) +
        ggrepel::geom_text_repel(
          aes(label = paste0(TumorType, " (", Frequency, ")")),
          color = "black",
          size = 3.2,
          box.padding = 0.7,
          segment.color = "grey50",
          max.overlaps = Inf
        ) +
        scale_size(range = c(3, 15)) +
        scale_color_brewer(palette = "Set3") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size = 17, color = "black"),
          legend.position = "none",
          plot.title = element_text(size = 17, color = "black", hjust = 0.5)
        ) +
        labs(
          title = "Tumor Type Distribution by System",
          x = "",
          y = "Number of Samples"
        )

      return(p_tumor_bubble)
    } else {
      warning("ggrepel package not available. Install with: install.packages('ggrepel')")

      # Fallback to regular bar plot
      p_tumor_bar <- ggplot(tumor_data, aes(x = reorder(TumorType, Frequency),
                                         y = Frequency, fill = System)) +
        geom_col() +
        coord_flip() +
        scale_fill_brewer(palette = "Set3") +
        theme_bw() +
        labs(x = "Tumor Type", y = "Frequency", fill = "System",
             title = "Tumor Type Distribution") +
        theme(axis.text = element_text(size = 10))

      return(p_tumor_bar)
    }
  }, error = function(e) {
    warning("Could not create tumor type plot: ", e$message)
    return(NULL)
  })
}

#' Create comprehensive statistical dashboard
#'
#' @description Combines multiple statistical plots into a dashboard layout
#' @param plot_list List of plots from generateStatisticalPlots
#' @param layout Character specifying layout: "grid", "tabs", or "sequential"
#' @return Combined plot object or list of plots
#' @export
createStatisticalDashboard <- function(plot_list, layout = "grid") {

  # Filter out NULL plots
  valid_plots <- plot_list[!sapply(plot_list, is.null)]

  if (length(valid_plots) == 0) {
    warning("No valid plots to display")
    return(NULL)
  }

  tryCatch({
    if (layout == "grid" && requireNamespace("patchwork", quietly = TRUE)) {
      # Use patchwork for grid layout
      combined_plot <- Reduce(`+`, valid_plots)
      return(combined_plot)
    } else if (layout == "grid" && requireNamespace("gridExtra", quietly = TRUE)) {
      # Use gridExtra as fallback
      return(gridExtra::grid.arrange(grobs = valid_plots))
    } else {
      # Return individual plots
      return(valid_plots)
    }
  }, error = function(e) {
    warning("Could not create dashboard layout: ", e$message)
    return(valid_plots)
  })
}
