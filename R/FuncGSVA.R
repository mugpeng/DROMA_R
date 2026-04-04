# GSVA Analysis Functions ----

#' Load gene sets from GMT file
#'
#' @description Reads a GMT (Gene Matrix Transposed) file and returns a list of gene sets
#' @param gmt_file Path to the GMT file
#' @return A named list where each element is a character vector of gene symbols for a gene set
#' @export
#' @examples
#' \dontrun{
#' # Load Hallmark gene sets
#' gmt_path <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "DROMA")
#' gene_sets <- loadGeneSetsFromGMT(gmt_path)
#' }
loadGeneSetsFromGMT <- function(gmt_file) {
  if (!file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file)
  }
  
  gene_sets <- list()
  
  con <- file(gmt_file, "r")
  on.exit(close(con))
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      gene_set_name <- parts[1]
      # Skip description (parts[2]) and get genes (parts[3:end])
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]  # Remove empty strings
      if (length(genes) > 0) {
        gene_sets[[gene_set_name]] <- genes
      }
    }
  }
  
  return(gene_sets)
}

#' Get available predefined gene sets
#'
#' @description Returns a list of available predefined gene set names from the package's GMT file,
#'              or optionally returns the full gene sets with their genes
#' @param gmt_file Optional path to GMT file. If NULL, uses the default Hallmark gene sets file
#' @param return_full_list Logical, if TRUE returns a named list where each element is a character vector
#'                        of gene symbols for that gene set. If FALSE (default), returns only gene set names.
#' @return If return_full_list is FALSE: A character vector of gene set names.
#'         If return_full_list is TRUE: A named list where each element is a character vector of gene symbols.
#' @export
#' @examples
#' \dontrun{
#' # Get all available Hallmark gene sets (names only)
#' available_sets <- getAvailableGeneSets()
#' print(available_sets)
#'
#' # Get full gene sets with genes
#' full_gene_sets <- getAvailableGeneSets(return_full_list = TRUE)
#' # View genes in first gene set
#' print(full_gene_sets[[1]])
#' # View names and gene counts
#' sapply(full_gene_sets, length)
#' }
getAvailableGeneSets <- function(gmt_file = NULL, return_full_list = FALSE) {
  if (is.null(gmt_file)) {
    # Try multiple possible paths for the default GMT file
    possible_paths <- c(
      # Relative path from package root
      file.path("Data", "h.all.v2025.1.Hs.symbols.gmt"),
      # Absolute path if Data is in workspace root
      file.path(getwd(), "Data", "h.all.v2025.1.Hs.symbols.gmt"),
      # Path relative to package installation
      file.path(system.file(package = "DROMA"), "..", "..", "Data", "h.all.v2025.1.Hs.symbols.gmt")
    )
    
    default_path <- NULL
    for (path in possible_paths) {
      normalized_path <- normalizePath(path, mustWork = FALSE)
      if (file.exists(normalized_path)) {
        default_path <- normalized_path
        break
      }
    }
    
    if (is.null(default_path)) {
      stop("Default GMT file not found. Please specify gmt_file parameter. ",
           "Tried paths: ", paste(possible_paths, collapse = ", "))
    }
    gmt_file <- default_path
  }
  
  gene_sets <- loadGeneSetsFromGMT(gmt_file)
  
  if (return_full_list) {
    return(gene_sets)
  } else {
    return(names(gene_sets))
  }
}

#' Get gene set by name from predefined gene sets
#'
#' @description Retrieves a specific gene set from the predefined GMT file
#' @param gene_set_name Name of the gene set to retrieve
#' @param gmt_file Optional path to GMT file. If NULL, uses the default Hallmark gene sets file
#' @return A character vector of gene symbols in the gene set
#' @export
#' @examples
#' \dontrun{
#' # Get HALLMARK_APOPTOSIS gene set
#' apoptosis_genes <- getGeneSetByName("HALLMARK_APOPTOSIS")
#' }
getGeneSetByName <- function(gene_set_name, gmt_file = NULL) {
  if (is.null(gmt_file)) {
    # Try multiple possible paths for the default GMT file
    possible_paths <- c(
      # Relative path from package root
      file.path("Data", "h.all.v2025.1.Hs.symbols.gmt"),
      # Absolute path if Data is in workspace root
      file.path(getwd(), "Data", "h.all.v2025.1.Hs.symbols.gmt"),
      # Path relative to package installation
      file.path(system.file(package = "DROMA"), "..", "..", "Data", "h.all.v2025.1.Hs.symbols.gmt")
    )
    
    default_path <- NULL
    for (path in possible_paths) {
      normalized_path <- normalizePath(path, mustWork = FALSE)
      if (file.exists(normalized_path)) {
        default_path <- normalized_path
        break
      }
    }
    
    if (is.null(default_path)) {
      stop("Default GMT file not found. Please specify gmt_file parameter. ",
           "Tried paths: ", paste(possible_paths, collapse = ", "))
    }
    gmt_file <- default_path
  }
  
  gene_sets <- loadGeneSetsFromGMT(gmt_file)
  
  if (!gene_set_name %in% names(gene_sets)) {
    stop("Gene set '", gene_set_name, "' not found in GMT file. Available gene sets: ", 
         paste(head(names(gene_sets), 10), collapse = ", "), "...")
  }
  
  return(gene_sets[[gene_set_name]])
}

.prepareGSVAGeneSets <- function(gene_sets, gmt_file = NULL) {
  all_available_gene_sets <- tryCatch({
    getAvailableGeneSets(gmt_file = gmt_file, return_full_list = TRUE)
  }, error = function(e) {
    NULL
  })
  available_gene_set_names <- if (!is.null(all_available_gene_sets)) names(all_available_gene_sets) else character(0)

  if (is.character(gene_sets) && length(gene_sets) == 1) {
    if (gene_sets %in% available_gene_set_names) {
      gene_sets <- all_available_gene_sets[gene_sets]
    } else {
      warning("'", gene_sets, "' is not a predefined gene set name. Treating as a single gene.")
      gene_sets <- list(custom_gene_set = gene_sets)
    }
  } else if (is.character(gene_sets) && length(gene_sets) > 1) {
    matching_names <- gene_sets %in% available_gene_set_names

    if (all(matching_names)) {
      gene_sets <- all_available_gene_sets[gene_sets]
    } else if (any(matching_names)) {
      unrecognized <- gene_sets[!matching_names]
      warning("The following are not predefined gene set names and will be ignored: ",
              paste(unrecognized, collapse = ", "))
      gene_sets <- all_available_gene_sets[gene_sets[matching_names]]
    } else {
      message("No predefined gene set names matched. Treating input as a custom gene list for a single gene set.")
      gene_sets <- list(custom_gene_set = gene_sets)
    }
  } else if (is.list(gene_sets)) {
    if (is.null(names(gene_sets))) {
      names(gene_sets) <- paste0("gene_set_", seq_along(gene_sets))
    }
  } else {
    stop("gene_sets must be a character vector, a named list, or a predefined gene set name")
  }

  if (length(gene_sets) == 0) {
    stop("No valid gene sets provided")
  }

  gene_sets
}

.calculateMatrixGSVA <- function(expr_data, gene_sets, method, min_size, max_size, ..., verbose_value = FALSE) {
  expr_data <- as.matrix(expr_data)

  if (any(is.na(expr_data))) {
    expr_data[is.na(expr_data)] <- 0
  }

  available_genes <- rownames(expr_data)
  filtered_gene_sets <- lapply(gene_sets, function(gs) {
    intersect(gs, available_genes)
  })
  gene_set_sizes <- sapply(filtered_gene_sets, length)
  valid_gene_sets <- filtered_gene_sets[gene_set_sizes >= min_size & gene_set_sizes <= max_size]

  if (length(valid_gene_sets) == 0) {
    return(NULL)
  }

  dots <- list(...)

  if (method == "gsva") {
    kcdf_value <- if ("kcdf" %in% names(dots)) dots$kcdf else "Gaussian"
    tau_value <- if ("tau" %in% names(dots)) dots$tau else 1
    maxDiff_value <- if ("maxDiff" %in% names(dots)) dots$maxDiff else TRUE
    absRanking_value <- if ("absRanking" %in% names(dots)) dots$absRanking else FALSE

    gsva_param <- GSVA::gsvaParam(
      expr_data,
      valid_gene_sets,
      kcdf = kcdf_value,
      tau = tau_value,
      maxDiff = maxDiff_value,
      absRanking = absRanking_value,
      minSize = min_size,
      maxSize = max_size
    )
  } else if (method == "ssgsea") {
    alpha_value <- if ("alpha" %in% names(dots)) dots$alpha else 0.25
    normalize_value <- if ("normalize" %in% names(dots)) dots$normalize else TRUE

    gsva_param <- GSVA::ssgseaParam(
      expr_data,
      valid_gene_sets,
      alpha = alpha_value,
      normalize = normalize_value,
      minSize = min_size,
      maxSize = max_size
    )
  } else if (method == "zscore") {
    gsva_param <- GSVA::zscoreParam(
      expr_data,
      valid_gene_sets,
      minSize = min_size,
      maxSize = max_size
    )
  } else if (method == "plage") {
    gsva_param <- GSVA::plageParam(
      expr_data,
      valid_gene_sets,
      minSize = min_size,
      maxSize = max_size
    )
  } else {
    stop("Unknown method: ", method)
  }

  as.matrix(GSVA::gsva(gsva_param, verbose = verbose_value))
}

#' Calculate GSVA scores for gene sets
#'
#' @description Main function for calculating gene set enrichment scores using GSVA methods
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param gene_sets Either a named list of gene sets (each element is a character vector of gene symbols), 
#'                  or a character vector of gene symbols for a single gene set, 
#'                  or a character string naming a predefined gene set
#' @param method GSVA method to use: "gsva" (default), "ssgsea", "zscore", or "plage"
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param min_size Minimum number of genes in a gene set that must be present in the expression data (default: 5)
#' @param max_size Maximum number of genes in a gene set (default: Inf)
#' @param gmt_file Optional path to GMT file if using predefined gene set names
#' @param ... Additional parameters passed to GSVA::gsva()
#' @return A list containing GSVA scores (gene sets x samples matrix) for each dataset
#' @export
#' @examples
#' \dontrun{
#' # Using DromaSet with custom gene list
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' custom_genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' scores <- calculateGSVAScores(gCSI, custom_genes, method = "ssgsea")
#'
#' # Using predefined gene set name
#' scores <- calculateGSVAScores(gCSI, "HALLMARK_APOPTOSIS", method = "gsva")
#'
#' # Using multiple gene sets as a list
#' gene_sets_list <- list(
#'   apoptosis = c("TP53", "BAX", "BCL2"),
#'   cell_cycle = c("CDK1", "CDK2", "CCNB1")
#' )
#' scores <- calculateGSVAScores(gCSI, gene_sets_list, method = "plage")
#'
#' # Using MultiDromaSet
#' multi_set <- createMultiDromaSetFromDatabase(c("gCSI", "CCLE"), "path/to/droma.sqlite")
#' scores <- calculateGSVAScores(multi_set, "HALLMARK_APOPTOSIS")
#' }
calculateGSVAScores <- function(dromaset_object, gene_sets,
                                method = c("gsva", "ssgsea", "zscore", "plage"),
                                data_type = "all", tumor_type = "all",
                                overlap_only = FALSE,
                                min_size = 5, max_size = Inf,
                                gmt_file = NULL, ...) {
  
  # Check if GSVA package is available
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required. Please install it using: BiocManager::install('GSVA')")
  }
  
  # Validate input object
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  # Validate method
  method <- match.arg(method)
  
  gene_sets <- .prepareGSVAGeneSets(gene_sets, gmt_file = gmt_file)
  
  # Load mRNA expression data
  if (inherits(dromaset_object, "DromaSet")) {
    # Single DromaSet
    expr_data <- loadMolecularProfiles(dromaset_object,
                                      feature_type = "mRNA",
                                      select_features = NULL,
                                      data_type = data_type,
                                      tumor_type = tumor_type,
                                      return_data = TRUE,
                                      zscore = FALSE)  # GSVA will handle normalization
    
    if (is.null(expr_data) || !is.matrix(expr_data) || nrow(expr_data) == 0 || ncol(expr_data) == 0) {
      stop("No mRNA expression data found for the specified filters")
    }
    
    # Check for NA values and replace with 0
    if (any(is.na(expr_data))) {
      na_count <- sum(is.na(expr_data))
      total_count <- length(expr_data)
      warning("Found ", na_count, " NA values (", round(100 * na_count / total_count, 2), 
              "%) in expression data. Replacing with 0.")
      expr_data[is.na(expr_data)] <- 0
    }
    
    verbose_value <- if ("verbose" %in% names(list(...))) list(...)$verbose else FALSE
    gsva_result <- .calculateMatrixGSVA(
      expr_data = expr_data,
      gene_sets = gene_sets,
      method = method,
      min_size = min_size,
      max_size = max_size,
      verbose_value = verbose_value,
      ...
    )

    if (is.null(gsva_result)) {
      stop("No valid gene sets found after filtering. Check min_size and max_size parameters.")
    }
    
    result <- list()
    result[[dromaset_object@name]] <- gsva_result
    
  } else {
    # MultiDromaSet
    expr_data_list <- loadMultiProjectMolecularProfiles(dromaset_object,
                                                       feature_type = "mRNA",
                                                       select_features = NULL,
                                                       overlap_only = overlap_only,
                                                       data_type = data_type,
                                                       tumor_type = tumor_type,
                                                       zscore = FALSE)
    
    if (is.null(expr_data_list) || length(expr_data_list) == 0) {
      stop("No mRNA expression data found for the specified filters")
    }
    
    result <- list()
    
    for (dataset_name in names(expr_data_list)) {
      expr_data <- expr_data_list[[dataset_name]]
      
      if (is.null(expr_data) || !is.matrix(expr_data) || nrow(expr_data) == 0 || ncol(expr_data) == 0) {
        warning("Skipping dataset ", dataset_name, ": no expression data")
        next
      }
      
      # Check for NA values and replace with 0
      if (any(is.na(expr_data))) {
        na_count <- sum(is.na(expr_data))
        total_count <- length(expr_data)
        warning("Found ", na_count, " NA values (", round(100 * na_count / total_count, 2), 
                "%) in dataset '", dataset_name, "'. Replacing with 0.")
        expr_data[is.na(expr_data)] <- 0
      }
      
      verbose_value <- if ("verbose" %in% names(list(...))) list(...)$verbose else FALSE
      gsva_result <- .calculateMatrixGSVA(
        expr_data = expr_data,
        gene_sets = gene_sets,
        method = method,
        min_size = min_size,
        max_size = max_size,
        verbose_value = verbose_value,
        ...
      )

      if (is.null(gsva_result)) {
        warning("Skipping dataset ", dataset_name, ": no valid gene sets after filtering")
        next
      }
      
      result[[dataset_name]] <- gsva_result
    }
    
    if (length(result) == 0) {
      stop("No valid GSVA scores calculated for any dataset")
    }
  }
  
  return(result)
}

#' Analyze GSVA score associations with drugs or omics features
#'
#' @description Analyzes associations between GSVA scores and drug sensitivity or omics features,
#'              similar to analyzeDrugOmicPair but for gene set enrichment scores
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param gsva_scores GSVA scores from calculateGSVAScores(). Can be a matrix (single dataset) 
#'                    or a list of matrices (multiple datasets)
#' @param select_gene_set Name of the gene set to analyze (must be a rowname in gsva_scores)
#' @param select_drugs Name of the drug to analyze (for drug-GSVA association)
#' @param select_omics Optional: Name of the omics feature to analyze (for omics-GSVA association).
#'                     If provided, select_drugs will be ignored.
#' @param feature_type Type of omics data if select_omics is provided (e.g., "mRNA", "mutation_gene")
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param merged_enabled Logical, whether to create a merged dataset from all studies
#' @param meta_enabled Logical, whether to perform meta-analysis
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @param data_type_anno Optional character string to add as annotation in plot titles
#' @param plot_ncol Number of columns for multi-plot layout (default: 3)
#' @return A list containing plot (individual study plots), merged_plot (if merged_enabled=TRUE), 
#'         meta-analysis results, and data
#' @export
#' @examples
#' \dontrun{
#' # Calculate GSVA scores first
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' gsva_scores <- calculateGSVAScores(gCSI, "HALLMARK_APOPTOSIS", method = "ssgsea")
#'
#' # Analyze association with drug sensitivity
#' result <- analyzeGSVAAssociation(gCSI, gsva_scores, "HALLMARK_APOPTOSIS", "Paclitaxel")
#'
#' # Analyze association with omics feature
#' result <- analyzeGSVAAssociation(gCSI, gsva_scores, "HALLMARK_APOPTOSIS", 
#'                                  select_omics = "TP53", feature_type = "mutation_gene")
#' }
analyzeGSVAAssociation <- function(dromaset_object, gsva_scores, select_gene_set,
                                  select_drugs = NULL, select_omics = NULL,
                                  feature_type = NULL,
                                  data_type = "all", tumor_type = "all",
                                  overlap_only = FALSE,
                                  merged_enabled = TRUE,
                                  meta_enabled = TRUE,
                                  zscore = TRUE,
                                  data_type_anno = NULL,
                                  plot_ncol = 3) {
  
  # Validate inputs
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  if (is.null(select_drugs) && is.null(select_omics)) {
    stop("Either select_drugs or select_omics must be provided")
  }
  
  if (!is.null(select_omics) && is.null(feature_type)) {
    stop("feature_type must be provided when select_omics is specified")
  }
  
  # Process GSVA scores - convert to list format if needed
  if (is.matrix(gsva_scores)) {
    # Single dataset - convert to list
    dataset_name <- if (inherits(dromaset_object, "DromaSet")) {
      dromaset_object@name
    } else {
      names(gsva_scores)[1]  # Use first dataset name if available
    }
    gsva_scores <- list(gsva_scores)
    names(gsva_scores) <- dataset_name
  }
  
  # Check if gene set exists in scores
  if (!select_gene_set %in% rownames(gsva_scores[[1]])) {
    stop("Gene set '", select_gene_set, "' not found in GSVA scores. Available gene sets: ",
         paste(head(rownames(gsva_scores[[1]]), 10), collapse = ", "), "...")
  }
  
  # Extract GSVA scores for the selected gene set
  gsva_data <- lapply(gsva_scores, function(score_mat) {
    if (select_gene_set %in% rownames(score_mat)) {
      score_vec <- as.numeric(score_mat[select_gene_set, ])
      names(score_vec) <- colnames(score_mat)
      return(score_vec[!is.na(score_vec)])
    }
    return(NULL)
  })
  
  # Remove NULL entries
  gsva_data <- gsva_data[!sapply(gsva_data, is.null)]
  
  if (length(gsva_data) == 0) {
    stop("No valid GSVA scores found for gene set '", select_gene_set, "'")
  }
  
  # Load comparison data (drug or omics)
  if (!is.null(select_omics)) {
    # Analyze GSVA vs omics
    is_continuous_omics <- feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")
    
    comparison_data <- loadFeatureData(dromaset_object, feature_type, select_omics,
                                      data_type = data_type, tumor_type = tumor_type,
                                      overlap_only = overlap_only, 
                                      is_continuous = is_continuous_omics,
                                      return_all_samples = !is_continuous_omics,
                                      zscore = zscore)
    
    comparison_data <- filterFeatureData(comparison_data, min_samples = 3, 
                                       is_discrete_with_all = !is_continuous_omics)
    
    if (is.null(comparison_data) || length(comparison_data) == 0) {
      stop("No sufficient comparison data found for the specified omics feature")
    }
    
    # Pair GSVA scores with omics data
    # For continuous: GSVA (feature1/x) vs omics (feature2/y)
    # For discrete: omics groups (discrete) vs GSVA scores (continuous)
    if (is_continuous_omics) {
      myPairs <- pairContinuousFeatures(gsva_data, comparison_data, merged = merged_enabled)
    } else {
      # For discrete omics, pair omics (discrete) with GSVA (continuous)
      myPairs <- pairDiscreteFeatures(comparison_data, gsva_data, merged = merged_enabled)
    }
    
  } else {
    # Analyze GSVA vs drug
    drug_data <- loadFeatureData(dromaset_object, "drug", select_drugs,
                                data_type = data_type, tumor_type = tumor_type,
                                overlap_only = overlap_only, is_continuous = TRUE,
                                zscore = zscore)
    
    drug_data <- filterFeatureData(drug_data, min_samples = 3, is_discrete_with_all = FALSE)
    
    if (is.null(drug_data) || length(drug_data) == 0) {
      stop("No sufficient drug data found for the specified drug")
    }
    
    # Pair GSVA scores (feature1/x) with drug data (feature2/y)
    myPairs <- pairContinuousFeatures(gsva_data, drug_data, merged = merged_enabled)
  }
  
  # Check if we have sufficient data after pairing
  if (is.null(myPairs) || length(myPairs) == 0) {
    stop("No sufficient paired data found. Each dataset needs at least 3 samples.")
  }
  
  # Initialize result list
  result <- list()
  
  # Determine if we're dealing with continuous or discrete comparison
  is_continuous_comparison <- is.null(select_omics) || 
    feature_type %in% c("mRNA", "meth", "proteinrppa", "cnv", "proteinms")
  
  # Create title suffix
  title_suffix <- if (!is.null(data_type_anno) && nchar(data_type_anno) > 0) {
    paste0(" (", data_type_anno, ")")
  } else {
    ""
  }
  
  if (is_continuous_comparison) {
    # Handle continuous comparison (GSVA vs drug or continuous omics)
    individual_pairs <- myPairs
    merged_pair <- NULL
    
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }
    
    # Create plots for individual studies
    if (length(individual_pairs) == 1) {
      single_pair <- individual_pairs[[1]]
      # GSVA is feature1 (x-axis), comparison is feature2 (y-axis)
      x_label <- "GSVA Score"
      y_label <- if (!is.null(select_omics)) {
        paste(feature_type, "expression")
      } else {
        "drug sensitivity(Area Above Curve)"
      }
      
      title <- if (!is.null(select_omics)) {
        paste0(select_gene_set, " vs ", select_omics, " (", feature_type, ")", title_suffix)
      } else {
        paste0(select_gene_set, " vs ", select_drugs, title_suffix)
      }
      
      result$plot <- plotCorrelation(single_pair$feature1,
                                     single_pair$feature2,
                                     x_label = x_label,
                                     y_label = y_label,
                                     title = title,
                                     method = "spearman")
    } else if (length(individual_pairs) > 1) {
      # GSVA is feature1 (x-axis), comparison is feature2 (y-axis)
      x_label <- select_gene_set
      y_label <- if (!is.null(select_omics)) {
        paste0(select_omics, " (", feature_type, ")")
      } else {
        "Drug Response"
      }
      
      multi_plot <- plotMultipleCorrelations(individual_pairs,
                                            x_label = x_label,
                                            y_label = y_label,
                                            data_type_anno = data_type_anno)
      
      x_title <- "GSVA Score"
      y_title <- if (!is.null(select_omics)) {
        paste(feature_type, "expression")
      } else {
        "drug sensitivity(Area Above Curve)"
      }
      
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                            x_title = x_title,
                                            y_title = y_title)
    }
    
    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      # GSVA is feature1 (x-axis), comparison is feature2 (y-axis)
      x_label <- "GSVA Score"
      y_label <- if (!is.null(select_omics)) {
        paste(feature_type, "expression")
      } else {
        "drug sensitivity(Area Above Curve)"
      }
      
      title <- if (!is.null(select_omics)) {
        paste0(select_gene_set, " vs ", select_omics, " (", feature_type, ")", title_suffix)
      } else {
        paste0(select_gene_set, " vs ", select_drugs, title_suffix)
      }
      
      result$merged_plot <- plotCorrelation(merged_pair$feature1,
                                           merged_pair$feature2,
                                           x_label = x_label,
                                           y_label = y_label,
                                           title = title,
                                           method = "spearman")
    }
    
    # Perform meta-analysis
    if (meta_enabled && length(individual_pairs) > 1) {
      meta_result <- metaCalcConCon(individual_pairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    result$data <- myPairs
    
  } else {
    # Handle discrete comparison (GSVA vs discrete omics like mutations)
    individual_pairs <- myPairs
    merged_pair <- NULL
    
    if (merged_enabled && "merged_dataset" %in% names(myPairs)) {
      merged_pair <- myPairs[["merged_dataset"]]
      individual_pairs <- myPairs[names(myPairs) != "merged_dataset"]
    }
    
    # Create plots for individual studies
    if (length(individual_pairs) == 1) {
      single_pair <- individual_pairs[[1]]
      result$plot <- plotGroupComparison(single_pair$no, single_pair$yes,
                                        group_labels = c(paste("Without", select_omics), 
                                                         paste("With", select_omics)),
                                        title = paste0(select_gene_set, " vs ", select_omics, 
                                                      " (", feature_type, ")", title_suffix),
                                        y_label = "GSVA Score")
    } else if (length(individual_pairs) > 1) {
      ncol_value <- if (length(individual_pairs) == 4) 2 else plot_ncol
      multi_plot <- plotMultipleGroupComparisons(individual_pairs,
                                                group_labels = c(paste("Without", select_omics), 
                                                               paste("With", select_omics)),
                                                x_label = paste0(select_omics, " (", feature_type, ")"),
                                                y_label = "GSVA Score",
                                                ncol = ncol_value,
                                                data_type_anno = data_type_anno)
      result$plot <- createPlotWithCommonAxes(multi_plot,
                                            y_title = "GSVA Score")
    }
    
    # Create plot for merged dataset if available
    if (!is.null(merged_pair) && merged_enabled) {
      result$merged_plot <- plotGroupComparison(merged_pair$no, merged_pair$yes,
                                                group_labels = c(paste("Without", select_omics), 
                                                               paste("With", select_omics)),
                                                title = paste0(select_gene_set, " vs ", select_omics, 
                                                              " (", feature_type, ")", title_suffix),
                                                y_label = "GSVA Score")
    }
    
    # Perform meta-analysis
    if (meta_enabled && length(individual_pairs) > 1) {
      meta_result <- metaCalcConDis(individual_pairs)
      if (!is.null(meta_result)) {
        result$meta <- meta_result
      }
    }
    
    result$data <- myPairs
  }
  
  # Return results
  if (is.null(result$plot) && is.null(result$merged_plot)) {
    return(list())
  } else {
    return(result)
  }
}

#' Batch analyze GSVA scores for multiple pathways against a single feature
#'
#' @description Performs batch analysis of associations between GSVA scores for multiple gene sets
#'              and a single feature (drug or omics), returning a data frame with statistics for each gene set.
#'              Use this when you have pre-computed GSVA scores for multiple pathways and want to test
#'              each pathway against the same feature.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param gsva_scores GSVA scores from calculateGSVAScores(). Can be a matrix (single dataset) 
#'                    or a list of matrices (multiple datasets)
#' @param feature_name Name of the feature to analyze (drug name or omics feature name)
#' @param feature_type Type of feature: "drug" (default), or omics types like "mRNA", "mutation_gene", etc.
#' @param gene_set_names Optional: character vector of gene set names to analyze. 
#'                       If NULL, analyzes all gene sets in gsva_scores
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param zscore Logical, whether to apply z-score normalization (default: TRUE)
#' @param show_progress Logical, whether to show progress messages (default: TRUE)
#' @return A data frame with columns: name (gene set name), p_value, effect_size, n_datasets, q_value
#' @export
#' @examples
#' \dontrun{
#' # Calculate GSVA scores for all Hallmark gene sets
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' all_gene_sets <- getAvailableGeneSets()
#' gsva_scores <- calculateGSVAScores(gCSI, all_gene_sets, method = "gsva")
#'
#' # Batch analyze associations with a drug
#' results <- batchAnalyzePathwaysWithFeature(gCSI, gsva_scores, "Bortezomib", feature_type = "drug")
#'
#' # Batch analyze associations with a gene expression
#' results <- batchAnalyzePathwaysWithFeature(gCSI, gsva_scores, "TP53", feature_type = "mRNA")
#'
#' # Batch analyze associations with mutation status
#' results <- batchAnalyzePathwaysWithFeature(gCSI, gsva_scores, "KRAS", feature_type = "mutation_gene")
#' }
batchAnalyzePathwaysWithFeature <- function(dromaset_object, gsva_scores, feature_name,
                                            feature_type = "drug",
                                            gene_set_names = NULL,
                                            data_type = "all", tumor_type = "all",
                                            overlap_only = FALSE,
                                            zscore = TRUE,
                                            show_progress = TRUE) {
  
  # Validate inputs
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  # Validate feature_type
  valid_feature_types <- c("mRNA", "cnv", "meth", "proteinrppa", "proteinms",
                           "drug", "mutation_gene", "mutation_site", "fusion")
  if (!feature_type %in% valid_feature_types) {
    stop(paste0("Invalid feature_type. Please choose from: ",
                paste(valid_feature_types, collapse = ", ")))
  }
  
  # Determine if feature is continuous or discrete
  continuous_types <- c("drug", "cnv", "proteinrppa", "proteinms", "meth", "mRNA")
  is_continuous <- feature_type %in% continuous_types
  
  # Process GSVA scores - convert to list format if needed
  if (is.matrix(gsva_scores)) {
    dataset_name <- if (inherits(dromaset_object, "DromaSet")) {
      dromaset_object@name
    } else {
      "dataset1"
    }
    gsva_scores <- list(gsva_scores)
    names(gsva_scores) <- dataset_name
  }
  
  # Get gene set names to analyze
  if (is.null(gene_set_names)) {
    gene_set_names <- rownames(gsva_scores[[1]])
  }
  
  if (length(gene_set_names) == 0) {
    stop("No gene sets found in GSVA scores")
  }
  
  # Load feature data once
  feature_data <- loadFeatureData(dromaset_object, feature_type, feature_name,
                                  data_type = data_type, tumor_type = tumor_type,
                                  overlap_only = overlap_only, 
                                  is_continuous = is_continuous,
                                  return_all_samples = !is_continuous,
                                  zscore = zscore)
  
  feature_data <- filterFeatureData(feature_data, min_samples = 3, 
                                    is_discrete_with_all = !is_continuous)
  
  if (is.null(feature_data) || length(feature_data) == 0) {
    stop("No sufficient data found for the specified feature: ", feature_name)
  }
  
  # Initialize results
  results_list <- list()
  total_sets <- length(gene_set_names)
  
  # Process each gene set
  for (i in seq_along(gene_set_names)) {
    gene_set_name <- gene_set_names[i]
    
    if (show_progress && (i %% 10 == 0 || i == 1 || i == total_sets)) {
      message(sprintf("Processing gene set %d/%d: %s", i, total_sets, gene_set_name))
    }
    
    # Extract GSVA scores for this gene set
    gsva_data <- lapply(gsva_scores, function(score_mat) {
      if (gene_set_name %in% rownames(score_mat)) {
        score_vec <- as.numeric(score_mat[gene_set_name, ])
        names(score_vec) <- colnames(score_mat)
        return(score_vec[!is.na(score_vec)])
      }
      return(NULL)
    })
    
    # Remove NULL entries
    gsva_data <- gsva_data[!sapply(gsva_data, is.null)]
    
    if (length(gsva_data) == 0) {
      next
    }
    
    # Pair GSVA scores with feature data based on feature type
    tryCatch({
      if (is_continuous) {
        # Continuous feature: use correlation
        myPairs <- pairContinuousFeatures(gsva_data, feature_data, merged = FALSE)
        
        if (is.null(myPairs) || length(myPairs) == 0) {
          next
        }
        
        # Perform meta-analysis or single correlation
        if (length(myPairs) > 1) {
          meta_result <- metaCalcConCon(myPairs)
          
          if (!is.null(meta_result)) {
            p_val <- meta_result[["pval.random"]]
            eff_size <- meta_result[["TE.random"]]
            n_datasets <- length(meta_result[["studlab"]])
            
            # Replace NA values
            if (is.na(p_val)) p_val <- 1
            if (is.na(eff_size)) eff_size <- 0
            
            results_list[[gene_set_name]] <- data.frame(
              name = gene_set_name,
              p_value = p_val,
              effect_size = eff_size,
              n_datasets = n_datasets,
              stringsAsFactors = FALSE
            )
          }
        } else {
          # Single dataset - calculate correlation directly
          single_pair <- myPairs[[1]]
          if (length(single_pair$feature1) >= 3 && length(single_pair$feature2) >= 3) {
            cor_test <- cor.test(single_pair$feature1, single_pair$feature2, method = "spearman")
            
            results_list[[gene_set_name]] <- data.frame(
              name = gene_set_name,
              p_value = cor_test$p.value,
              effect_size = cor_test$estimate,
              n_datasets = 1,
              stringsAsFactors = FALSE
            )
          }
        }
      } else {
        # Discrete feature: use group comparison
        # For discrete: feature groups (discrete) vs GSVA scores (continuous)
        myPairs <- pairDiscreteFeatures(feature_data, gsva_data, merged = FALSE)
        
        if (is.null(myPairs) || length(myPairs) == 0) {
          next
        }
        
        # Perform meta-analysis or single comparison
        if (length(myPairs) > 1) {
          meta_result <- metaCalcConDis(myPairs)
          
          if (!is.null(meta_result)) {
            p_val <- meta_result[["pval.random"]]
            eff_size <- meta_result[["TE.random"]]
            n_datasets <- length(meta_result[["studlab"]])
            
            # Replace NA values
            if (is.na(p_val)) p_val <- 1
            if (is.na(eff_size)) eff_size <- 0
            
            results_list[[gene_set_name]] <- data.frame(
              name = gene_set_name,
              p_value = p_val,
              effect_size = eff_size,
              n_datasets = n_datasets,
              stringsAsFactors = FALSE
            )
          }
        } else {
          # Single dataset - calculate group comparison directly
          single_pair <- myPairs[[1]]
          if (length(single_pair$no) >= 3 && length(single_pair$yes) >= 3) {
            wilcox_test <- wilcox.test(single_pair$no, single_pair$yes)
            # Effect size: difference in medians / pooled SD (approximate)
            pooled_sd <- sd(c(single_pair$no, single_pair$yes))
            effect_size <- (median(single_pair$yes) - median(single_pair$no)) / pooled_sd
            
            results_list[[gene_set_name]] <- data.frame(
              name = gene_set_name,
              p_value = wilcox_test$p.value,
              effect_size = effect_size,
              n_datasets = 1,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }, error = function(e) {
      # Silently skip on error
      NULL
    })
  }
  
  # Combine results
  if (length(results_list) == 0) {
    return(data.frame(
      name = character(0),
      p_value = numeric(0),
      effect_size = numeric(0),
      n_datasets = integer(0),
      q_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # Add q_value (FDR correction)
  results_df$q_value <- p.adjust(results_df$p_value, method = "BH")
  
  # Sort by p-value
  results_df <- results_df[order(results_df$p_value), ]
  
  return(results_df)
}

#' Batch analysis of clinical drug response for multiple pathways
#'
#' @description Calculates patient-level GSVA scores from CTRDB expression data
#'   and tests pathway associations with clinical drug response using the same
#'   response-group and meta-analysis logic as \code{batchFindClinicalSigResponse()}.
#' @param select_pathways Character vector of pathway names to evaluate.
#' @param select_drugs Character string specifying the drug name.
#' @param gene_sets Optional named list of gene sets, or predefined gene set
#'   names accepted by \code{calculateGSVAScores()}.
#' @param method GSVA method to use: "gsva" (default), "ssgsea", "zscore", or "plage".
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX").
#' @param tumor_type Filter by tumor type ("all" or specific tumor types).
#' @param min_size Minimum number of genes required per pathway in each patient.
#' @param max_size Maximum number of genes allowed per pathway.
#' @param gmt_file Optional GMT file path for predefined pathway names.
#' @param cores Number of CPU cores (kept for API consistency; not used).
#' @param zscore Logical, whether to z-score pathway scores within each patient.
#' @param show_progress Logical, whether to print progress messages.
#' @param ... Additional parameters passed to GSVA.
#' @return Data frame with \code{p_value}, \code{effect_size}, \code{n_datasets},
#'   \code{name}, and \code{q_value}.
#' @export
batchFindClinicalPathwaySigResponse <- function(select_pathways,
                                                select_drugs,
                                                gene_sets = select_pathways,
                                                method = c("gsva", "ssgsea", "zscore", "plage"),
                                                data_type = "all",
                                                tumor_type = "all",
                                                min_size = 5,
                                                max_size = Inf,
                                                gmt_file = NULL,
                                                cores = 1,
                                                zscore = TRUE,
                                                show_progress = TRUE,
                                                ...) {
  if (is.null(select_pathways) || length(select_pathways) == 0) {
    stop("select_pathways must be a non-empty character vector")
  }
  if (is.null(select_drugs) || select_drugs == "") {
    stop("select_drugs must be specified")
  }
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("GSVA package is required. Please install it using: BiocManager::install('GSVA')")
  }

  method <- match.arg(method)
  gene_sets <- .prepareGSVAGeneSets(gene_sets, gmt_file = gmt_file)
  missing_pathways <- setdiff(select_pathways, names(gene_sets))
  if (length(missing_pathways) > 0) {
    warning("The following pathways were not found in gene_sets and will be ignored: ",
            paste(missing_pathways, collapse = ", "))
  }

  gene_sets <- gene_sets[intersect(select_pathways, names(gene_sets))]
  select_pathways <- names(gene_sets)
  if (length(select_pathways) == 0) {
    return(data.frame(
      p_value = numeric(0),
      effect_size = numeric(0),
      n_datasets = integer(0),
      name = character(0),
      q_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  if (!exists("ctrdb_connection", envir = .GlobalEnv)) {
    stop("No CTRDB database connection found. Connect first with connectCTRDatabase()")
  }
  connection <- get("ctrdb_connection", envir = .GlobalEnv)

  sample_anno <- tryCatch({
    getDROMAAnnotation("sample")
  }, error = function(e) {
    stop("Failed to load sample annotations from DROMA: ", e$message)
  })

  if (is.null(sample_anno) || nrow(sample_anno) == 0) {
    stop("No sample annotations found in database")
  }
  if (!"CliUsedDrug" %in% colnames(sample_anno)) {
    stop("CliUsedDrug column not found in sample annotations")
  }

  drug_samples <- sample_anno[!is.na(sample_anno$CliUsedDrug) &
                               grepl(select_drugs, sample_anno$CliUsedDrug, ignore.case = TRUE), ]
  if (nrow(drug_samples) == 0) {
    stop("No samples found with drug: ", select_drugs)
  }

  if ("ProjectID" %in% colnames(drug_samples)) {
    ctrdb_samples <- drug_samples[drug_samples$ProjectID == "CTRDB", ]
  } else {
    stop("ProjectID column not found in sample annotations")
  }
  if (nrow(ctrdb_samples) == 0) {
    stop("No CTRDB samples found with drug: ", select_drugs)
  }

  if (!is.null(data_type) && data_type != "all") {
    if ("DataType" %in% colnames(ctrdb_samples)) {
      ctrdb_samples <- ctrdb_samples[ctrdb_samples$DataType == data_type, ]
      if (nrow(ctrdb_samples) == 0) {
        stop("No samples found with data type: ", data_type)
      }
    } else {
      warning("DataType column not found. Skipping data type filter.")
    }
  }

  if (!is.null(tumor_type) && tumor_type != "all") {
    if ("TumorType" %in% colnames(ctrdb_samples)) {
      tumor_pattern <- paste0("^", tumor_type, "$", collapse = "|")
      ctrdb_samples <- ctrdb_samples[grepl(tumor_pattern, ctrdb_samples$TumorType, ignore.case = TRUE), ]
      if (nrow(ctrdb_samples) == 0) {
        stop("No samples found with tumor type: ", tumor_type)
      }
    } else {
      warning("TumorType column not found. Skipping tumor type filter.")
    }
  }

  if (!"PatientID" %in% colnames(ctrdb_samples)) {
    stop("PatientID column not found in sample annotations")
  }

  unique_datasets <- unique(ctrdb_samples$PatientID)
  unique_datasets <- unique_datasets[!is.na(unique_datasets)]
  if (length(unique_datasets) == 0) {
    stop("No valid PatientIDs found")
  }

  message("Found ", length(unique_datasets), " datasets with drug ", select_drugs)

  patient_context <- list()
  total_datasets <- length(unique_datasets)

  for (i in seq_along(unique_datasets)) {
    patient_id <- unique_datasets[i]

    if (show_progress && (i %% 10 == 0 || i == 1 || i == total_datasets)) {
      message(sprintf("Processing clinical GSVA dataset %d/%d: %s", i, total_datasets, patient_id))
    }

    tryCatch({
      patient_samples <- ctrdb_samples[ctrdb_samples$PatientID == patient_id, ]

      if (!"SampleID" %in% colnames(patient_samples) ||
          !"Response" %in% colnames(patient_samples)) {
        warning("Missing SampleID or Response columns for patient: ", patient_id)
        next
      }

      expr_data <- getPatientExpressionData(patient_id, connection, auto_log = TRUE)
      if (is.null(expr_data) || nrow(expr_data) == 0) {
        warning("No expression data found for patient: ", patient_id)
        next
      }

      common_samples <- intersect(patient_samples$SampleID, colnames(expr_data))
      if (length(common_samples) < 3) {
        warning("Insufficient samples for patient: ", patient_id, " (", length(common_samples), " samples)")
        next
      }

      patient_meta <- patient_samples[patient_samples$SampleID %in% common_samples, ]
      expr_subset <- expr_data[, common_samples, drop = FALSE]
      pathway_scores <- .calculateMatrixGSVA(
        expr_data = expr_subset,
        gene_sets = gene_sets,
        method = method,
        min_size = min_size,
        max_size = max_size,
        verbose_value = FALSE,
        ...
      )

      if (is.null(pathway_scores) || nrow(pathway_scores) == 0) {
        next
      }

      patient_context[[patient_id]] <- list(
        pathway_scores = pathway_scores,
        patient_meta = patient_meta,
        common_samples = common_samples
      )
    }, error = function(e) {
      warning("Error processing patient ", patient_id, ": ", e$message)
    })
  }

  if (length(patient_context) == 0) {
    return(data.frame(
      p_value = numeric(0),
      effect_size = numeric(0),
      n_datasets = integer(0),
      name = character(0),
      q_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  message("Successfully processed ", length(patient_context), " datasets")

  cal_re_list <- lapply(select_pathways, function(pathway_name) {
    patient_data_list <- list()

    for (patient_id in names(patient_context)) {
      ctx <- patient_context[[patient_id]]
      if (!pathway_name %in% rownames(ctx$pathway_scores)) {
        next
      }

      pathway_result <- tryCatch({
        pathway_values <- as.numeric(ctx$pathway_scores[pathway_name, ctx$common_samples])
        names(pathway_values) <- ctx$common_samples

        if (zscore) {
          mean_value <- mean(pathway_values, na.rm = TRUE)
          sd_value <- stats::sd(pathway_values, na.rm = TRUE)
          if (sd_value > 0) {
            pathway_values <- (pathway_values - mean_value) / sd_value
          }
        }

        response_groups <- split(
          pathway_values,
          ctx$patient_meta$Response[match(names(pathway_values), ctx$patient_meta$SampleID)]
        )

        if (!"Response" %in% names(response_groups) || !"Non_response" %in% names(response_groups)) {
          return(NULL)
        }

        response_values <- response_groups$Response
        non_response_values <- response_groups$Non_response

        if (length(response_values) < 2 || length(non_response_values) < 2) {
          return(NULL)
        }

        list(
          response = response_values,
          non_response = non_response_values,
          metadata = ctx$patient_meta
        )
      }, error = function(e) NULL)

      if (!is.null(pathway_result)) {
        patient_data_list[[patient_id]] <- pathway_result
      }
    }

    if (length(patient_data_list) == 0) {
      return(NULL)
    }

    if (length(patient_data_list) > 1) {
      meta <- analyzeClinicalMeta(patient_data_list)
      if (is.null(meta)) {
        return(NULL)
      }
      p_val <- meta[["pval.random"]]
      eff_size <- meta[["TE.random"]]
      n_datasets <- length(meta[["studlab"]])
    } else {
      dat <- patient_data_list[[1]]
      resp <- dat$response
      non_resp <- dat$non_response

      if (length(resp) < 2 || length(non_resp) < 2) {
        return(NULL)
      }

      wilcox_re <- tryCatch(
        suppressWarnings(stats::wilcox.test(resp, non_resp)),
        error = function(e) NULL
      )
      cliff_re <- tryCatch(
        suppressWarnings(effsize::cliff.delta(resp, non_resp)),
        error = function(e) NULL
      )

      if (is.null(wilcox_re) || is.null(cliff_re)) {
        return(NULL)
      }

      p_val <- wilcox_re$p.value
      eff_size <- cliff_re$estimate
      n_datasets <- 1L
    }

    if (is.na(p_val)) {
      p_val <- 1
    }
    if (is.na(eff_size)) {
      eff_size <- 0
    }

    data.frame(
      p_value = p_val,
      effect_size = eff_size,
      n_datasets = n_datasets,
      name = pathway_name,
      stringsAsFactors = FALSE
    )
  })

  valid <- !sapply(cal_re_list, is.null)
  cal_re_list <- cal_re_list[valid]

  if (length(cal_re_list) == 0) {
    return(data.frame(
      p_value = numeric(0),
      effect_size = numeric(0),
      n_datasets = integer(0),
      name = character(0),
      q_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  cal_re_df <- do.call(rbind, cal_re_list)
  rownames(cal_re_df) <- NULL
  cal_re_df$q_value <- p.adjust(cal_re_df$p_value, method = "BH")
  cal_re_df[order(cal_re_df$p_value), , drop = FALSE]
}

#' Batch analyze a single pathway (gene list) across multiple features
#'
#' @description Given a single gene list (pathway), calculates GSVA score and tests associations
#'              with multiple features (drugs or omics). Similar to batchFindSignificantFeatures
#'              but uses pathway-level GSVA scores instead of single gene expression.
#' @param dromaset_object Either a DromaSet or MultiDromaSet object
#' @param gene_list Character vector of gene symbols defining the pathway
#' @param pathway_name Name for this pathway (used in output)
#' @param feature_type Type of features to test against: "drug" (default), "mRNA", "cnv", 
#'                     "mutation_gene", "mutation_site", "fusion", "meth", "proteinrppa", "proteinms"
#' @param feature_names Optional: character vector of specific feature names to test.
#'                      If NULL, tests all available features of the specified type
#' @param gsva_method GSVA method to use: "gsva" (default), "ssgsea", "zscore", or "plage"
#' @param data_type Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")
#' @param tumor_type Filter by tumor type ("all" or specific tumor types)
#' @param overlap_only For MultiDromaSet, whether to use only overlapping samples (default: FALSE)
#' @param zscore Logical, whether to apply z-score normalization to feature data (default: TRUE)
#' @param min_size Minimum number of genes in gene_list that must be present in expression data (default: 5)
#' @param show_progress Logical, whether to show progress messages (default: TRUE)
#' @param test_top_n Integer, number of top features to test (default: NULL for all features)
#' @param cores Number of CPU cores to use for parallel processing (default: 1 for sequential).
#'               Uses future + furrr with task chunking for efficiency when cores > 1
#' @return A data frame with columns: name (feature name), p_value, effect_size, n_datasets, q_value
#' @export
#' @examples
#' \dontrun{
#' # Define a custom gene list
#' my_genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2", "RAD51", "PALB2")
#'
#' # Test pathway association with all drugs
#' gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
#' results <- batchAnalyzePathwayAcrossFeatures(gCSI, my_genes, "DNA_Repair", 
#'                                              feature_type = "drug")
#'
#' # Test pathway association with specific drugs
#' results <- batchAnalyzePathwayAcrossFeatures(gCSI, my_genes, "DNA_Repair",
#'                                              feature_type = "drug",
#'                                              feature_names = c("Olaparib", "Cisplatin"))
#'
#' # Test pathway association with all genes (mRNA)
#' results <- batchAnalyzePathwayAcrossFeatures(gCSI, my_genes, "DNA_Repair",
#'                                              feature_type = "mRNA", test_top_n = 1000)
#'
#' # Test pathway association with mutation status
#' results <- batchAnalyzePathwayAcrossFeatures(gCSI, my_genes, "DNA_Repair",
#'                                              feature_type = "mutation_gene")
#'
#' # With parallel processing
#' results <- batchAnalyzePathwayAcrossFeatures(gCSI, my_genes, "DNA_Repair",
#'                                              feature_type = "drug", cores = 8)
#' }
batchAnalyzePathwayAcrossFeatures <- function(dromaset_object,
                                              gene_list,
                                              pathway_name,
                                              feature_type = "drug",
                                              feature_names = NULL,
                                              gsva_method = c("gsva", "ssgsea", "zscore", "plage"),
                                              data_type = "all",
                                              tumor_type = "all",
                                              overlap_only = FALSE,
                                              zscore = TRUE,
                                              min_size = 5,
                                              show_progress = TRUE,
                                              test_top_n = NULL,
                                              cores = 1) {
  
  # Validate inputs
  if (!inherits(dromaset_object, c("DromaSet", "MultiDromaSet"))) {
    stop("Input must be a DromaSet or MultiDromaSet object from DROMA.Set package")
  }
  
  if (!is.character(gene_list) || length(gene_list) < min_size) {
    stop(paste0("gene_list must be a character vector with at least ", min_size, " genes"))
  }
  
  if (!is.character(pathway_name) || length(pathway_name) != 1) {
    stop("pathway_name must be a single character string")
  }
  
  # Validate feature_type
  valid_feature_types <- c("mRNA", "cnv", "meth", "proteinrppa", "proteinms",
                           "drug", "mutation_gene", "mutation_site", "fusion")
  if (!feature_type %in% valid_feature_types) {
    stop(paste0("Invalid feature_type. Please choose from: ",
                paste(valid_feature_types, collapse = ", ")))
  }
  
  # Validate gsva_method
  gsva_method <- match.arg(gsva_method)
  
  # Validate cores parameter
  max_cores <- ifelse(requireNamespace("parallel", quietly = TRUE),
                    parallel::detectCores() - 1,
                    1)
  if (!is.numeric(cores) || cores < 1 || cores > max_cores) {
    stop(paste0("cores must be a positive integer between 1 and ", max_cores))
  }
  
  # Determine if feature is continuous or discrete
  continuous_types <- c("drug", "cnv", "proteinrppa", "proteinms", "meth", "mRNA")
  is_continuous <- feature_type %in% continuous_types
  
  # Step 1: Calculate GSVA scores for the gene list
  if (show_progress) {
    message("Step 1: Calculating GSVA scores for pathway '", pathway_name, "'...")
  }
  
  gene_sets <- list(gene_list)
  names(gene_sets) <- pathway_name
  
  gsva_scores <- tryCatch({
    calculateGSVAScores(dromaset_object, gene_sets,
                       method = gsva_method,
                       data_type = data_type,
                       tumor_type = tumor_type,
                       overlap_only = overlap_only,
                       min_size = min_size)
  }, error = function(e) {
    stop("Failed to calculate GSVA scores: ", e$message)
  })
  
  if (is.null(gsva_scores) || length(gsva_scores) == 0) {
    stop("No GSVA scores could be calculated. Check if gene_list genes are present in expression data.")
  }
  
  # Extract GSVA scores for the pathway (convert to named vector format per dataset)
  gsva_data <- lapply(gsva_scores, function(score_mat) {
    if (pathway_name %in% rownames(score_mat)) {
      score_vec <- as.numeric(score_mat[pathway_name, ])
      names(score_vec) <- colnames(score_mat)
      return(score_vec[!is.na(score_vec)])
    }
    return(NULL)
  })
  gsva_data <- gsva_data[!sapply(gsva_data, is.null)]
  
  if (length(gsva_data) == 0) {
    stop("No valid GSVA scores calculated for pathway '", pathway_name, "'")
  }
  
  if (show_progress) {
    total_samples <- sum(sapply(gsva_data, length))
    message(sprintf("  Calculated GSVA scores for %d samples across %d dataset(s)", 
                   total_samples, length(gsva_data)))
  }
  
  # Step 2: Get feature list to test
  if (show_progress) {
    message("Step 2: Loading feature data for type '", feature_type, "'...")
  }
  
  # Get and validate feature list
  feature_list <- getAndValidateFeatureList(
    dromaset_object = dromaset_object,
    feature_type = feature_type,
    feature_names = feature_names,
    data_type = data_type,
    tumor_type = tumor_type
  )
  
  # Apply test_top_n filter
  if (!is.null(test_top_n) && is.numeric(test_top_n) && test_top_n > 0 && length(feature_list) > test_top_n) {
    feature_list <- feature_list[seq_len(min(test_top_n, length(feature_list)))]
  }
  
  if (length(feature_list) == 0) {
    stop("No features found for the selected feature type")
  }
  
  if (show_progress) {
    message(sprintf("  Found %d features to test", length(feature_list)))
  }
  
  # Step 3: Batch analysis
  if (show_progress) {
    message("Step 3: Running batch analysis...")
  }
  
  # Preload all feature data at once for efficiency
  preloaded_feature_data <- loadFeatureData(dromaset_object, feature_type, feature_list,
                                            data_type = data_type, 
                                            tumor_type = tumor_type,
                                            overlap_only = overlap_only, 
                                            is_continuous = is_continuous,
                                            return_all_samples = !is_continuous,
                                            zscore = zscore)
  
  # Filter preloaded data to ensure minimum sample size
  preloaded_feature_data <- filterFeatureData(preloaded_feature_data, min_samples = 3,
                                              is_discrete_with_all = !is_continuous)
  
  # Pre-compute sample intersection cache for continuous vs continuous
  sample_intersection_cache <- NULL
  if (is_continuous) {
    sample_intersection_cache <- list()
    for (i in seq_along(gsva_data)) {
      for (j in seq_along(preloaded_feature_data)) {
        cache_key <- paste(names(gsva_data)[i], names(preloaded_feature_data)[j], sep = "||")
        sample_intersection_cache[[cache_key]] <- intersect(
          names(gsva_data[[i]]),
          colnames(preloaded_feature_data[[j]])
        )
      }
    }
  }
  
  # Create worker function to process a single feature
  processPathwayFeature <- function(feature_idx) {
    tryCatch({
      feature_name <- feature_list[feature_idx]
      
      # Extract feature data from preloaded data
      if (is_continuous) {
        # Extract from preloaded continuous data
        feature_data <- lapply(preloaded_feature_data, function(data_matrix) {
          if (is.matrix(data_matrix) && feature_name %in% rownames(data_matrix)) {
            feature_vector <- as.numeric(data_matrix[feature_name, ])
            names(feature_vector) <- colnames(data_matrix)
            return(feature_vector[!is.na(feature_vector)])
          }
          return(NULL)
        })
        feature_data <- feature_data[!sapply(feature_data, is.null)]
      } else {
        # Extract from preloaded discrete data
        feature_data <- lapply(preloaded_feature_data, function(dataset) {
          if (is.list(dataset) && "present" %in% names(dataset)) {
            if (feature_name %in% names(dataset$present)) {
              return(list(present = dataset$present[[feature_name]], all = dataset$all))
            }
          }
          return(NULL)
        })
        feature_data <- feature_data[!sapply(feature_data, is.null)]
      }
      
      if (is.null(feature_data) || length(feature_data) == 0) {
        return(NULL)
      }
      
      # Pair and analyze based on feature type
      if (is_continuous) {
        # GSVA (continuous) vs feature (continuous): correlation
        myPairs <- pairContinuousFeatures(gsva_data, feature_data, merged = FALSE,
                                          intersection_cache = sample_intersection_cache)
        
        if (is.null(myPairs) || length(myPairs) == 0) {
          return(NULL)
        }
        
        if (length(myPairs) > 1) {
          meta_result <- metaCalcConCon(myPairs)
          
          if (!is.null(meta_result)) {
            p_val <- meta_result[["pval.random"]]
            eff_size <- meta_result[["TE.random"]]
            n_datasets <- length(meta_result[["studlab"]])
            
            if (is.na(p_val)) p_val <- 1
            if (is.na(eff_size)) eff_size <- 0
            
            return(data.frame(
              name = feature_name,
              p_value = p_val,
              effect_size = eff_size,
              n_datasets = n_datasets,
              stringsAsFactors = FALSE
            ))
          }
        } else {
          single_pair <- myPairs[[1]]
          if (length(single_pair$feature1) >= 3 && length(single_pair$feature2) >= 3) {
            cor_test <- cor.test(single_pair$feature1, single_pair$feature2, method = "spearman")
            
            return(data.frame(
              name = feature_name,
              p_value = cor_test$p.value,
              effect_size = cor_test$estimate,
              n_datasets = 1,
              stringsAsFactors = FALSE
            ))
          }
        }
      } else {
        # GSVA (continuous) vs feature (discrete): group comparison
        myPairs <- pairDiscreteFeatures(feature_data, gsva_data, merged = FALSE)
        
        if (is.null(myPairs) || length(myPairs) == 0) {
          return(NULL)
        }
        
        if (length(myPairs) > 1) {
          meta_result <- metaCalcConDis(myPairs)
          
          if (!is.null(meta_result)) {
            p_val <- meta_result[["pval.random"]]
            eff_size <- meta_result[["TE.random"]]
            n_datasets <- length(meta_result[["studlab"]])
            
            if (is.na(p_val)) p_val <- 1
            if (is.na(eff_size)) eff_size <- 0
            
            return(data.frame(
              name = feature_name,
              p_value = p_val,
              effect_size = eff_size,
              n_datasets = n_datasets,
              stringsAsFactors = FALSE
            ))
          }
        } else {
          single_pair <- myPairs[[1]]
          if (length(single_pair$no) >= 3 && length(single_pair$yes) >= 3) {
            wilcox_test <- wilcox.test(single_pair$no, single_pair$yes)
            pooled_sd <- sd(c(single_pair$no, single_pair$yes))
            effect_size <- (median(single_pair$yes) - median(single_pair$no)) / pooled_sd
            
            return(data.frame(
              name = feature_name,
              p_value = wilcox_test$p.value,
              effect_size = effect_size,
              n_datasets = 1,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      return(NULL)
    }, error = function(e) {
      # Silently skip on error
      return(NULL)
    })
  }
  
  # Initialize timing
  start_time <- Sys.time()
  total_features <- length(feature_list)
  
  message("Please be patient, it may take long time to run.")
  
  # Use parallel processing if cores > 1
  if (cores > 1) {
    # Setup future plan: use multicore (memory-efficient) on Unix/Mac, multisession on Windows
    # Increase memory limit for both backends
    old_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)  # 2GB
    
    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = cores)
      message(sprintf("Using fork-based parallelization with %d cores (memory-efficient)", cores))
    } else {
      future::plan(future::multisession, workers = cores)
      message(sprintf("Using multisession parallelization with %d cores (Windows)", cores))
    }
    
    # Unified cleanup to ensure proper execution order
    on.exit({
      future::plan(future::sequential)
      options(future.globals.maxSize = old_size)
    }, add = TRUE)
    
    # Adaptive task chunking based on feature count
    n_features <- length(feature_list)
    
    # Warn if over-parallelization
    if (cores > n_features / 100) {
      recommended_cores <- max(1, floor(n_features / 1000))
      message(sprintf("Note: Using %d cores for %d features may cause overhead. Consider cores=%d for optimal performance.",
                     cores, n_features, recommended_cores))
    }
    
    chunk_size <- if (n_features < 100) {
      max(10, ceiling(n_features / cores))  # Small: minimal chunking
    } else if (n_features < 1000) {
      ceiling(n_features / (cores * 4))  # Medium: 4 tasks/core
    } else {
      # Large: 8 tasks/core with minimum chunk_size to avoid excessive scheduling overhead
      max(50, min(200, ceiling(n_features / (cores * 8))))
    }
    chunks <- split(seq_along(feature_list),
                   ceiling(seq_along(feature_list) / chunk_size))
    
    # Process chunks in parallel with progressr support
    if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = n_features)
        chunk_results <- furrr::future_map(chunks, function(indices) {
          lapply(indices, function(idx) {
            result <- processPathwayFeature(idx)
            p()  # Update progress
            result
          })
        }, .options = furrr::furrr_options(seed = TRUE))
      })
    } else {
      if (show_progress) {
        message("Note: Install 'progressr' for progress tracking: install.packages('progressr')")
      }
      chunk_results <- furrr::future_map(chunks, function(indices) {
        lapply(indices, processPathwayFeature)
      }, .options = furrr::furrr_options(seed = TRUE))
    }
    
    # Flatten results
    results_list <- unlist(chunk_results, recursive = FALSE)
  } else {
    # Run sequential computation
    if (show_progress) {
      results_list <- lapply(seq_along(feature_list), function(i) {
        if (i %% 50 == 0 || i == 1 || i == total_features) {
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
          message(sprintf("  Processing feature %d/%d: %s (%.1fs elapsed)", 
                         i, total_features, feature_list[i], elapsed))
        }
        processPathwayFeature(i)
      })
    } else {
      results_list <- lapply(seq_along(feature_list), processPathwayFeature)
    }
  }
  
  # Filter out NULL results
  valid_results <- !sapply(results_list, is.null)
  results_list <- results_list[valid_results]
  
  # Combine results
  if (length(results_list) == 0) {
    return(data.frame(
      name = character(0),
      p_value = numeric(0),
      effect_size = numeric(0),
      n_datasets = integer(0),
      q_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  # Add q_value (FDR correction)
  results_df$q_value <- p.adjust(results_df$p_value, method = "BH")
  
  # Sort by p-value
  results_df <- results_df[order(results_df$p_value), ]
  
  if (show_progress) {
    total_time <- difftime(Sys.time(), start_time, units = "secs")
    message(sprintf("Analysis completed in %.1f seconds", as.numeric(total_time)))
    
    # Report filtered features if any
    n_filtered <- total_features - nrow(results_df)
    if (n_filtered > 0) {
      message(sprintf("Note: %d feature(s) were filtered out due to insufficient data or failed processing.", n_filtered))
    }
    
    message(sprintf("Found %d significant associations (q < 0.05) out of %d features",
                   sum(results_df$q_value < 0.05), nrow(results_df)))
  }
  
  return(results_df)
}
