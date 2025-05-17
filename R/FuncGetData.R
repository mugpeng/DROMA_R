#' Select features from omics or drug data
#'
#' @description Function to select feature data including omic data and drug data
#' @param select_feas_type The type of feature to select (e.g., "mRNA", "cnv", "drug")
#' @param select_feas The specific feature to select within the feature type
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO" (patient-derived organoids), "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or any specific tumor type (e.g., "lung cancer", "breast cancer")
#' @return A list of selected features filtered by the specified data type and tumor type
#' @export
selectFeatures <- function(select_feas_type, select_feas,
                        data_type = "all", tumor_type = "all") {
  if(!select_feas_type %in% c("mRNA","cnv",
                               "meth", "proteinrppa", "proteinms", # continuous
                               "mutation_gene", "mutation_site", "fusion", # discrete
                               "drug", "drug_raw"
  )){
    stop("The select feature type doesn't exist. Please choose from drug, mRNA, meth, cnv, proteinms, proteinrppa, mutation_gene, mutation_site, or fusion.")
  }

  # Validate data_type parameter
  if(!data_type %in% c("all", "CellLine", "PDO", "PDC", "PDX")) {
    stop("Invalid data_type. Please choose from 'all', 'CellLine', 'PDO', 'PDC', or 'PDX'.")
  }

  # Get cells based on data_type and tumor_type filter
  filtered_cells <- NULL
  if(data_type != "all" || tumor_type != "all") {
    # Use "SampleID" as the sample ID column
    if(!"SampleID" %in% colnames(sample_anno)) {
      warning("'SampleID' column not found in sample_anno. Using all samples.")
    } else {
      # Create a logical vector for filtering
      filter_idx <- rep(TRUE, nrow(sample_anno))

      # Apply data_type filter if specified
      if(data_type != "all") {
        filter_idx <- filter_idx & (sample_anno$DataType == data_type)
      }

      # Apply tumor_type filter if specified
      if(tumor_type != "all") {
        if(!"TumorType" %in% colnames(sample_anno)) {
          warning("'TumorType' column not found in sample_anno. Ignoring tumor_type filter.")
        } else {
          filter_idx <- filter_idx & (sample_anno$TumorType == tumor_type)
        }
      }

      # Get filtered sample IDs
      filtered_cells <- sample_anno$SampleID[filter_idx]

      if(length(filtered_cells) == 0) {
        # Construct appropriate error message
        if(data_type != "all" && tumor_type != "all") {
          stop(paste0("No samples found with data type '", data_type, "' and tumor type '", tumor_type, "'"))
        } else if(data_type != "all") {
          stop(paste0("No samples found with data type: ", data_type))
        } else {
          stop(paste0("No samples found with tumor type: ", tumor_type))
        }
      }
    }
  }

  feas_sel_vec <- fea_list[[select_feas_type]]
  if(select_feas_type %in% c("mRNA","cnv",
                              "meth", "proteinrppa", "proteinms",
                              "drug", "drug_raw")){
    fea_sel_list <- lapply(feas_sel_vec, function(x){
      # x = feas_sel_vec[1]
      fea <- as.data.frame(base::get(paste0(x, "_", select_feas_type), envir = .GlobalEnv))
      fea_sel <- as.numeric(fea[select_feas,])
      if(all(is.na(fea_sel))){
        return(NULL)
      }
      names(fea_sel) <- colnames(fea)

      # Filter by data_type and tumor_type if specified
      if(!is.null(filtered_cells)) {
        common_cells <- intersect(names(fea_sel), filtered_cells)
        if(length(common_cells) == 0) {
          return(NULL)
        }
        fea_sel <- fea_sel[common_cells]
      }

      fea_sel
    })
  } else{
    fea_sel_list <- lapply(feas_sel_vec, function(x){
      # x = feas_sel_vec[1]
      fea <- as.data.frame(base::get(paste0(x, "_", select_feas_type), envir = .GlobalEnv))
      fea_sel <- fea$cells[fea[[1]] %in% select_feas]
      if(all(is.na(fea_sel))){
        return(NULL)
      }

      # Filter by data_type and tumor_type if specified
      if(!is.null(filtered_cells)) {
        fea_sel <- intersect(fea_sel, filtered_cells)
        if(length(fea_sel) == 0) {
          return(NULL)
        }
      }

      fea_sel
    })
  }
  names(fea_sel_list) <- feas_sel_vec
  fea_sel_list <- fea_sel_list[!sapply(fea_sel_list, is.null)]
  if(length(fea_sel_list) < 1){
    if(data_type == "all" && tumor_type == "all") {
      stop("The select feature doesn't exist, check if it is in your selected feature types. \nsuch as, 'YM-155' for 'drug'.")
    } else {
      # Construct appropriate error message
      filter_msg <- ""
      if(data_type != "all" && tumor_type != "all") {
        filter_msg <- paste0(" with data type '", data_type, "' and tumor type '", tumor_type, "'")
      } else if(data_type != "all") {
        filter_msg <- paste0(" with data type '", data_type, "'")
      } else {
        filter_msg <- paste0(" with tumor type '", tumor_type, "'")
      }
      stop(paste0("No data found for feature '", select_feas, "'", filter_msg, "."))
    }
  }
  fea_sel_list
}
