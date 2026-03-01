# Functional Enrichment Analysis Functions ----

#' Perform GO enrichment analysis
#'
#' @description Performs Gene Ontology (GO) enrichment analysis using clusterProfiler.
#'              Automatically converts gene symbols to Entrez IDs and performs enrichment
#'              for the specified GO ontology (BP, CC, or MF).
#'
#' @param gene A character vector of gene symbols (e.g., c("TP53", "BRCA1", "BRCA2"))
#' @param use_padj Logical, whether to use adjusted p-value (p.adjust) for filtering.
#'                 If TRUE, filters by p.adjust; if FALSE, filters by pvalue. Default: TRUE
#' @param OrgDb Character string specifying the organism database (e.g., "org.Hs.eg.db").
#'              Default: "org.Hs.eg.db"
#' @param ont Character string specifying GO ontology: "BP" (Biological Process),
#'            "CC" (Cellular Component), or "MF" (Molecular Function). Default: "BP"
#' @param pvalueCutoff Numeric, p-value cutoff for enrichment analysis. Default: 0.05
#' @param pAdjustMethod Character string, method for p-value adjustment. Options include
#'                      "BH" (Benjamini-Hochberg), "bonferroni", "fdr", etc. Default: "BH"
#' @param qvalueCutoff Numeric, q-value cutoff for enrichment analysis. Default: 1
#' @param minGSSize Minimum gene set size. Default: 10
#' @param maxGSSize Maximum gene set size. Default: 500
#' @param readable Logical, whether to convert Entrez IDs to gene symbols in results.
#'                 Default: TRUE
#' @param keyType Character string, input gene ID type. Default: "SYMBOL"
#' @param ... Additional parameters passed to enrichGO()
#'
#' @return An enrichResult object from clusterProfiler containing GO enrichment results.
#'         Use as.data.frame() to convert to a data frame.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' go_result <- runGO(genes)
#'
#' # Use p-value instead of adjusted p-value
#' go_result <- runGO(genes, use_padj = FALSE)
#'
#' # Analyze Molecular Function ontology
#' go_result <- runGO(genes, ont = "MF")
#'
#' # Custom parameters
#' go_result <- runGO(genes, pvalueCutoff = 0.01, qvalueCutoff = 0.1)
#'
#' # Convert to data frame
#' go_df <- as.data.frame(go_result)
#' }
runGO <- function(gene,
                  use_padj = TRUE,
                  OrgDb = "org.Hs.eg.db",
                  ont = c("BP", "CC", "MF"),
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1,
                  minGSSize = 10,
                  maxGSSize = 500,
                  readable = TRUE,
                  keyType = "SYMBOL",
                  ...) {
  
  # Check if clusterProfiler is available
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Please install it using: ",
         "BiocManager::install('clusterProfiler')")
  }
  
  # Check if OrgDb package is available
  if (!requireNamespace(OrgDb, quietly = TRUE)) {
    stop("Organism database package '", OrgDb, "' is required. ",
         "Please install it using: BiocManager::install('", OrgDb, "')")
  }
  
  # Validate inputs
  if (missing(gene) || is.null(gene) || length(gene) == 0) {
    stop("gene parameter is required and must be a non-empty character vector")
  }
  
  if (!is.character(gene)) {
    stop("gene must be a character vector of gene symbols")
  }
  
  # Validate ont parameter
  ont <- match.arg(ont)
  
  # Validate use_padj
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  
  # Remove NA and empty strings
  gene <- unique(gene[!is.na(gene) & gene != ""])
  
  if (length(gene) == 0) {
    stop("No valid genes provided after removing NA and empty values")
  }
  
  # Load OrgDb object
  OrgDb_obj <- tryCatch({
    get(OrgDb, envir = asNamespace(OrgDb))
  }, error = function(e) {
    stop("Failed to load organism database '", OrgDb, "': ", e$message)
  })
  
  # Convert gene symbols to Entrez IDs if needed
  if (keyType == "SYMBOL") {
    gene_mapping <- tryCatch({
      clusterProfiler::bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb_obj)
    }, error = function(e) {
      stop("Failed to convert gene symbols to Entrez IDs: ", e$message)
    })
    
    if (is.null(gene_mapping) || nrow(gene_mapping) == 0) {
      stop("No genes could be converted to Entrez IDs. Check gene symbols and OrgDb.")
    }
    
    gene_ids <- unique(gene_mapping$ENTREZID)
    keyType_actual <- "ENTREZID"
  } else {
    gene_ids <- gene
    keyType_actual <- keyType
  }
  
  if (length(gene_ids) < minGSSize) {
    warning("Number of input genes (", length(gene_ids), ") is less than minGSSize (", 
            minGSSize, "). Results may be limited.")
  }
  
  # Perform GO enrichment analysis
  # Set pvalueCutoff and qvalueCutoff to 1 initially to get all results
  # Then filter based on use_padj parameter
  enrich_result <- tryCatch({
    clusterProfiler::enrichGO(
      gene = gene_ids,
      OrgDb = OrgDb_obj,
      keyType = keyType_actual,
      ont = ont,
      pvalueCutoff = 1,  # Get all results first
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = 1,  # Get all results first
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = readable,
      ...
    )
  }, error = function(e) {
    stop("GO enrichment analysis failed: ", e$message)
  })
  
  # Check if enrichment was successful
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    warning("No enriched GO terms found. Try adjusting pvalueCutoff, qvalueCutoff, ",
            "or gene set size parameters.")
    return(enrich_result)
  }
  
  # Filter results based on use_padj parameter
  if (use_padj) {
    # Filter by adjusted p-value
    if ("p.adjust" %in% colnames(enrich_result@result)) {
      enrich_result@result <- enrich_result@result[
        enrich_result@result$p.adjust <= pvalueCutoff, , drop = FALSE
      ]
    } else {
      warning("p.adjust column not found in results. Using pvalue instead.")
      enrich_result@result <- enrich_result@result[
        enrich_result@result$pvalue <= pvalueCutoff, , drop = FALSE
      ]
    }
  } else {
    # Filter by p-value
    enrich_result@result <- enrich_result@result[
      enrich_result@result$pvalue <= pvalueCutoff, , drop = FALSE
    ]
  }
  
  # Additional filtering by qvalue if specified
  if (qvalueCutoff < 1 && "qvalue" %in% colnames(enrich_result@result)) {
    enrich_result@result <- enrich_result@result[
      enrich_result@result$qvalue <= qvalueCutoff, , drop = FALSE
    ]
  }
  
  # Sort by p-value or adjusted p-value
  if (use_padj && "p.adjust" %in% colnames(enrich_result@result)) {
    enrich_result@result <- enrich_result@result[
      order(enrich_result@result$p.adjust), , drop = FALSE
    ]
  } else {
    enrich_result@result <- enrich_result@result[
      order(enrich_result@result$pvalue), , drop = FALSE
    ]
  }
  
  return(enrich_result)
}

#' Perform KEGG pathway enrichment analysis
#'
#' @description Performs KEGG pathway enrichment analysis using clusterProfiler.
#'              Automatically converts gene symbols to Entrez IDs and performs enrichment
#'              for KEGG pathways.
#'
#' @param gene A character vector of gene symbols (e.g., c("TP53", "BRCA1", "BRCA2"))
#' @param use_padj Logical, whether to use adjusted p-value (p.adjust) for filtering.
#'                 If TRUE, filters by p.adjust; if FALSE, filters by pvalue. Default: TRUE
#' @param organism Character string specifying KEGG organism code (e.g., "hsa" for human,
#'                 "mmu" for mouse, "rno" for rat). Default: "hsa"
#' @param keyType Character string, input gene ID type. Default: "SYMBOL"
#' @param OrgDb Character string specifying the organism database for ID conversion
#'              (e.g., "org.Hs.eg.db"). Only used when keyType is "SYMBOL".
#'              Default: "org.Hs.eg.db"
#' @param pvalueCutoff Numeric, p-value cutoff for enrichment analysis. Default: 0.05
#' @param pAdjustMethod Character string, method for p-value adjustment. Options include
#'                      "BH" (Benjamini-Hochberg), "bonferroni", "fdr", etc. Default: "BH"
#' @param qvalueCutoff Numeric, q-value cutoff for enrichment analysis. Default: 1
#' @param minGSSize Minimum gene set size. Default: 10
#' @param maxGSSize Maximum gene set size. Default: 500
#' @param use_internal_data Logical, whether to use KEGG.db (deprecated). Default: FALSE
#' @param ... Additional parameters passed to enrichKEGG()
#'
#' @return An enrichResult object from clusterProfiler containing KEGG enrichment results.
#'         Use as.data.frame() to convert to a data frame.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' kegg_result <- runKEGG(genes)
#'
#' # Use p-value instead of adjusted p-value
#' kegg_result <- runKEGG(genes, use_padj = FALSE)
#'
#' # Analyze mouse genes
#' kegg_result <- runKEGG(genes, organism = "mmu", OrgDb = "org.Mm.eg.db")
#'
#' # Custom parameters
#' kegg_result <- runKEGG(genes, pvalueCutoff = 0.01, qvalueCutoff = 0.1)
#'
#' # Convert to data frame
#' kegg_df <- as.data.frame(kegg_result)
#' }
runKEGG <- function(gene,
                     use_padj = TRUE,
                     organism = "hsa",
                     keyType = "SYMBOL",
                     OrgDb = "org.Hs.eg.db",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 1,
                     minGSSize = 10,
                     maxGSSize = 500,
                     use_internal_data = FALSE,
                     ...) {
  
  # Check if clusterProfiler is available
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler package is required. Please install it using: ",
         "BiocManager::install('clusterProfiler')")
  }
  
  # Validate inputs
  if (missing(gene) || is.null(gene) || length(gene) == 0) {
    stop("gene parameter is required and must be a non-empty character vector")
  }
  
  if (!is.character(gene)) {
    stop("gene must be a character vector of gene symbols")
  }
  
  # Validate use_padj
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  
  # Remove NA and empty strings
  gene <- unique(gene[!is.na(gene) & gene != ""])
  
  if (length(gene) == 0) {
    stop("No valid genes provided after removing NA and empty values")
  }
  
  # Convert gene symbols to Entrez IDs if needed
  if (keyType == "SYMBOL") {
    # Check if OrgDb package is available
    if (!requireNamespace(OrgDb, quietly = TRUE)) {
      stop("Organism database package '", OrgDb, "' is required for symbol conversion. ",
           "Please install it using: BiocManager::install('", OrgDb, "')")
    }
    
    # Load OrgDb object
    OrgDb_obj <- tryCatch({
      get(OrgDb, envir = asNamespace(OrgDb))
    }, error = function(e) {
      stop("Failed to load organism database '", OrgDb, "': ", e$message)
    })
    
    gene_mapping <- tryCatch({
      clusterProfiler::bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb_obj)
    }, error = function(e) {
      stop("Failed to convert gene symbols to Entrez IDs: ", e$message)
    })
    
    if (is.null(gene_mapping) || nrow(gene_mapping) == 0) {
      stop("No genes could be converted to Entrez IDs. Check gene symbols and OrgDb.")
    }
    
    gene_ids <- unique(gene_mapping$ENTREZID)
    keyType_actual <- "kegg"  # enrichKEGG expects "kegg" as keyType for Entrez IDs
  } else {
    gene_ids <- gene
    keyType_actual <- keyType
  }
  
  if (length(gene_ids) < minGSSize) {
    warning("Number of input genes (", length(gene_ids), ") is less than minGSSize (", 
            minGSSize, "). Results may be limited.")
  }
  
  # Perform KEGG enrichment analysis
  # Set pvalueCutoff and qvalueCutoff to 1 initially to get all results
  # Then filter based on use_padj parameter
  enrich_result <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene = gene_ids,
      organism = organism,
      keyType = keyType_actual,
      pvalueCutoff = 1,  # Get all results first
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = 1,  # Get all results first
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      use_internal_data = use_internal_data,
      ...
    )
  }, error = function(e) {
    stop("KEGG enrichment analysis failed: ", e$message)
  })
  
  # Check if enrichment was successful
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    warning("No enriched KEGG pathways found. Try adjusting pvalueCutoff, qvalueCutoff, ",
            "or gene set size parameters.")
    return(enrich_result)
  }
  
  # Filter results based on use_padj parameter
  if (use_padj) {
    # Filter by adjusted p-value
    if ("p.adjust" %in% colnames(enrich_result@result)) {
      enrich_result@result <- enrich_result@result[
        enrich_result@result$p.adjust <= pvalueCutoff, , drop = FALSE
      ]
    } else {
      warning("p.adjust column not found in results. Using pvalue instead.")
      enrich_result@result <- enrich_result@result[
        enrich_result@result$pvalue <= pvalueCutoff, , drop = FALSE
      ]
    }
  } else {
    # Filter by p-value
    enrich_result@result <- enrich_result@result[
      enrich_result@result$pvalue <= pvalueCutoff, , drop = FALSE
    ]
  }
  
  # Additional filtering by qvalue if specified
  if (qvalueCutoff < 1 && "qvalue" %in% colnames(enrich_result@result)) {
    enrich_result@result <- enrich_result@result[
      enrich_result@result$qvalue <= qvalueCutoff, , drop = FALSE
    ]
  }
  
  # Sort by p-value or adjusted p-value
  if (use_padj && "p.adjust" %in% colnames(enrich_result@result)) {
    enrich_result@result <- enrich_result@result[
      order(enrich_result@result$p.adjust), , drop = FALSE
    ]
  } else {
    enrich_result@result <- enrich_result@result[
      order(enrich_result@result$pvalue), , drop = FALSE
    ]
  }
  
  return(enrich_result)
}

#' Create a dotplot for enrichment results
#'
#' @description Creates a dotplot visualization for enrichResult objects, showing top pathways
#'              with p-value or adjusted p-value on x-axis and pathway names on y-axis.
#'              The dot size represents gene count and color represents significance.
#'
#' @param enrich_result An enrichResult object from clusterProfiler (from runGO or runKEGG)
#' @param use_padj Logical, whether to use adjusted p-value (p.adjust) for x-axis and coloring.
#'                 If TRUE, uses p.adjust; if FALSE, uses pvalue. Default: TRUE
#' @param showCategory Number of top pathways to display. Default: 10
#' @param title Plot title. If NULL, automatically generated based on enrichment type
#' @param x_label Label for x-axis. If NULL, automatically generated
#' @param size_by Character string, variable to control dot size: "Count" (default) or "GeneRatio"
#' @param color_low Color for low significance (high p-value). Default: "#BEBADAFF"
#' @param color_high Color for high significance (low p-value). Default: "#FB8072FF"
#' @param ... Additional parameters passed to ggplot2 functions
#'
#' @return A ggplot2 object with dotplot visualization
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform GO enrichment
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' go_result <- runGO(genes)
#'
#' # Create dotplot with adjusted p-value
#' p <- plotEnrichDotplot(go_result, use_padj = TRUE, showCategory = 15)
#' print(p)
#'
#' # Create dotplot with p-value
#' p <- plotEnrichDotplot(go_result, use_padj = FALSE, showCategory = 10)
#' print(p)
#' }
plotEnrichDotplot <- function(enrich_result,
                              use_padj = TRUE,
                              showCategory = 10,
                              title = NULL,
                              x_label = NULL,
                              size_by = c("Count", "GeneRatio"),
                              color_low = "#BEBADAFF",
                              color_high = "#FB8072FF",
                              ...) {
  
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it using: install.packages('ggplot2')")
  }
  
  # Validate input
  if (missing(enrich_result) || is.null(enrich_result)) {
    stop("enrich_result is required and must be an enrichResult object")
  }
  
  # Check if it's an enrichResult object
  if (!inherits(enrich_result, "enrichResult")) {
    stop("enrich_result must be an enrichResult object from clusterProfiler")
  }
  
  # Validate use_padj
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  
  # Validate size_by
  size_by <- match.arg(size_by)
  
  # Extract result data frame
  result_df <- enrich_result@result
  
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("No enrichment results found in enrich_result object")
  }
  
  # Determine which p-value column to use
  pval_col <- if (use_padj && "p.adjust" %in% colnames(result_df)) {
    "p.adjust"
  } else if ("pvalue" %in% colnames(result_df)) {
    "pvalue"
  } else {
    stop("Neither p.adjust nor pvalue column found in enrichment results")
  }
  
  # Select top pathways
  if (is.numeric(showCategory) && showCategory > 0) {
    # Sort by p-value (ascending = most significant first)
    result_df <- result_df[order(result_df[[pval_col]]), , drop = FALSE]
    result_df <- head(result_df, min(showCategory, nrow(result_df)))
  } else if (is.character(showCategory)) {
    # Show specific pathways
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
  
  # Prepare data for plotting
  plot_df <- result_df
  
  # Calculate -log10(p-value) for x-axis
  plot_df$x_value <- -log10(plot_df[[pval_col]])
  
  # Determine size variable
  if (size_by == "Count" && "Count" %in% colnames(plot_df)) {
    plot_df$size_value <- plot_df$Count
    size_label <- "Gene Count"
  } else if (size_by == "GeneRatio" && "GeneRatio" %in% colnames(plot_df)) {
    # Convert GeneRatio from "a/b" format to numeric
    plot_df$size_value <- sapply(plot_df$GeneRatio, function(x) {
      parts <- strsplit(as.character(x), "/")[[1]]
      if (length(parts) == 2) {
        as.numeric(parts[1]) / as.numeric(parts[2])
      } else {
        NA
      }
    })
    size_label <- "Gene Ratio"
  } else {
    # Fallback to Count if available, otherwise use a constant
    if ("Count" %in% colnames(plot_df)) {
      plot_df$size_value <- plot_df$Count
      size_label <- "Gene Count"
    } else {
      plot_df$size_value <- 1
      size_label <- "Size"
      warning("Neither Count nor GeneRatio found. Using constant size.")
    }
  }
  
  # Ensure Description column exists for y-axis labels
  if (!"Description" %in% colnames(plot_df)) {
    stop("Description column not found in enrichment results")
  }
  
  # Order pathways by significance (most significant at top)
  plot_df$Description <- factor(plot_df$Description, 
                                levels = rev(plot_df$Description[order(plot_df[[pval_col]])]))
  
  # Create plot
  p <- ggplot2::ggplot(plot_df, 
                       ggplot2::aes(x = .data$x_value, 
                                    y = .data$Description,
                                    size = .data$size_value,
                                    color = .data[[pval_col]])) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_gradient(low = color_high, 
                                 high = color_low,
                                 name = if (use_padj) "p.adjust" else "pvalue",
                                 trans = "log10") +
    ggplot2::scale_size_continuous(name = size_label, range = c(3, 8)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      title = ggplot2::element_text(size = 15, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right"
    ) +
    ggplot2::xlab(if (!is.null(x_label)) x_label else 
                 paste0("-log10(", if (use_padj) "p.adjust" else "pvalue", ")")) +
    ggplot2::ylab("Pathway")
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  } else {
    # Auto-generate title based on enrichment type
    enrich_type <- if (inherits(enrich_result, "enrichResult")) {
      if ("ONTOLOGY" %in% colnames(result_df)) {
        paste0("GO ", result_df$ONTOLOGY[1], " Enrichment")
      } else {
        "KEGG Pathway Enrichment"
      }
    } else {
      "Enrichment Analysis"
    }
    p <- p + ggplot2::ggtitle(enrich_type)
  }
  
  return(p)
}

#' Create a barplot for enrichment results
#'
#' @description Creates a barplot visualization for enrichResult objects, showing top pathways
#'              with p-value or adjusted p-value on x-axis and pathway names on y-axis.
#'              The bar color represents significance level.
#'
#' @param enrich_result An enrichResult object from clusterProfiler (from runGO or runKEGG)
#' @param use_padj Logical, whether to use adjusted p-value (p.adjust) for x-axis and coloring.
#'                 If TRUE, uses p.adjust; if FALSE, uses pvalue. Default: TRUE
#' @param showCategory Number of top pathways to display. Default: 10
#' @param title Plot title. If NULL, automatically generated based on enrichment type
#' @param x_label Label for x-axis. If NULL, automatically generated
#' @param color_low Color for low significance (high p-value). Default: "#BEBADAFF"
#' @param color_high Color for high significance (low p-value). Default: "#FB8072FF"
#' @param ... Additional parameters passed to ggplot2 functions
#'
#' @return A ggplot2 object with barplot visualization
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform KEGG enrichment
#' genes <- c("TP53", "BRCA1", "BRCA2", "ATM", "CHEK2")
#' kegg_result <- runKEGG(genes)
#'
#' # Create barplot with adjusted p-value
#' p <- plotEnrichBarplot(kegg_result, use_padj = TRUE, showCategory = 15)
#' print(p)
#'
#' # Create barplot with p-value
#' p <- plotEnrichBarplot(kegg_result, use_padj = FALSE, showCategory = 10)
#' print(p)
#' }
plotEnrichBarplot <- function(enrich_result,
                              use_padj = TRUE,
                              showCategory = 10,
                              title = NULL,
                              x_label = NULL,
                              color_low = "#BEBADAFF",
                              color_high = "#FB8072FF",
                              ...) {
  
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required. Please install it using: install.packages('ggplot2')")
  }
  
  # Validate input
  if (missing(enrich_result) || is.null(enrich_result)) {
    stop("enrich_result is required and must be an enrichResult object")
  }
  
  # Check if it's an enrichResult object
  if (!inherits(enrich_result, "enrichResult")) {
    stop("enrich_result must be an enrichResult object from clusterProfiler")
  }
  
  # Validate use_padj
  if (!is.logical(use_padj) || length(use_padj) != 1) {
    stop("use_padj must be a single logical value (TRUE or FALSE)")
  }
  
  # Extract result data frame
  result_df <- enrich_result@result
  
  if (is.null(result_df) || nrow(result_df) == 0) {
    stop("No enrichment results found in enrich_result object")
  }
  
  # Determine which p-value column to use
  pval_col <- if (use_padj && "p.adjust" %in% colnames(result_df)) {
    "p.adjust"
  } else if ("pvalue" %in% colnames(result_df)) {
    "pvalue"
  } else {
    stop("Neither p.adjust nor pvalue column found in enrichment results")
  }
  
  # Select top pathways
  if (is.numeric(showCategory) && showCategory > 0) {
    # Sort by p-value (ascending = most significant first)
    result_df <- result_df[order(result_df[[pval_col]]), , drop = FALSE]
    result_df <- head(result_df, min(showCategory, nrow(result_df)))
  } else if (is.character(showCategory)) {
    # Show specific pathways
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
  
  # Prepare data for plotting
  plot_df <- result_df
  
  # Calculate -log10(p-value) for x-axis
  plot_df$x_value <- -log10(plot_df[[pval_col]])
  
  # Ensure Description column exists for y-axis labels
  if (!"Description" %in% colnames(plot_df)) {
    stop("Description column not found in enrichment results")
  }
  
  # Order pathways by significance (most significant at top)
  plot_df$Description <- factor(plot_df$Description, 
                                levels = rev(plot_df$Description[order(plot_df[[pval_col]])]))
  
  # Create plot
  p <- ggplot2::ggplot(plot_df, 
                       ggplot2::aes(x = .data$x_value, 
                                    y = .data$Description,
                                    fill = .data[[pval_col]])) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.8) +
    ggplot2::scale_fill_gradient(low = color_high, 
                                high = color_low,
                                name = if (use_padj) "p.adjust" else "pvalue",
                                trans = "log10") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      title = ggplot2::element_text(size = 15, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.position = "right"
    ) +
    ggplot2::xlab(if (!is.null(x_label)) x_label else 
                 paste0("-log10(", if (use_padj) "p.adjust" else "pvalue", ")")) +
    ggplot2::ylab("Pathway")
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  } else {
    # Auto-generate title based on enrichment type
    enrich_type <- if (inherits(enrich_result, "enrichResult")) {
      if ("ONTOLOGY" %in% colnames(result_df)) {
        paste0("GO ", result_df$ONTOLOGY[1], " Enrichment")
      } else {
        "KEGG Pathway Enrichment"
      }
    } else {
      "Enrichment Analysis"
    }
    p <- p + ggplot2::ggtitle(enrich_type)
  }
  
  return(p)
}
