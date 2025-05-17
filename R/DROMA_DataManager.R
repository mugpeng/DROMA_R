#!/usr/bin/env Rscript

#' Setup DROMA Environment
#'
#' @description Initializes the DROMA environment by loading search vectors and setting up
#' search dataframes for drugs and omics features.
#' @return Invisibly returns TRUE if successful
#' @export
#' @examples
#' \dontrun{
#' setupDROMA()
#' # Now search vectors are loaded and ready for use
#' # drugs_search and omics_search are available in the global environment
#' }
setupDROMA <- function() {
  # Load search vectors
  if(file.exists(system.file("data", "search_vec.Rda", package = "DROMA"))) {
    load(system.file("data", "search_vec.Rda", package = "DROMA"))
    assign("fea_list", fea_list, envir = .GlobalEnv)
    assign("feas_search", feas_search, envir = .GlobalEnv)
    assign("samples_search", samples_search, envir = .GlobalEnv)

    # Create drug and omics search vectors
    drugs_search <- feas_search[feas_search$type %in% "drug",]
    omics_search <- feas_search[!feas_search$type %in% "drug",]

    # Assign to global environment
    assign("drugs_search", drugs_search, envir = .GlobalEnv)
    assign("omics_search", omics_search, envir = .GlobalEnv)
  } else {
    warning("search_vec.Rda not found. Some functionality may be limited.")
  }

  # Load sample annotation by default
  if(file.exists(system.file("data", "anno.Rda", package = "DROMA"))) {
    load(system.file("data", "anno.Rda", package = "DROMA"))
    assign("sample_anno", sample_anno, envir = .GlobalEnv)
    assign("drug_anno", drug_anno, envir = .GlobalEnv)
  } else {
    warning("anno.Rda not found. Sample and drug annotations will not be available.")
  }

  message("Environment variables are successfully loaded.")
  invisible(TRUE)
}

#' Load DROMA Data
#'
#' @description Loads specified data type from the DROMA package and optionally applies Z-score normalization.
#' @param data_type Character string specifying the data type to load. Options: "drug", "mRNA", "mut" (mutation),
#'        "meth" (methylation), "protein", "cnv" (copy number variation), "fusion", or "all".
#' @param apply_normalization Logical, whether to apply Z-score normalization to the loaded data (default: TRUE).
#'        Only applicable for continuous data types.
#' @return Invisibly returns TRUE if successful
#' @export
#' @examples
#' \dontrun{
#' setupDROMA()  # Load search vectors first
#' loadDROMA("drug")  # Load drug data
#' loadDROMA("mRNA")  # Load mRNA data
#' loadDROMA("mut")   # Load mutation data
#' }
loadDROMA <- function(data_type, apply_normalization = TRUE) {
  if (!is.character(data_type)) {
    stop("data_type must be a character string")
  }

  data_type <- tolower(data_type)
  valid_data_types <- c("drug", "mrna", "mut", "meth", "protein", "cnv", "fusion", "all")

  if (!data_type %in% valid_data_types) {
    stop("Invalid data_type. Must be one of: ", paste(valid_data_types, collapse = ", "))
  }

  # Helper function to load a specific data file
  load_data_file <- function(file_name) {
    data_path <- system.file("data", paste0(file_name, ".Rda"), package = "DROMA")
    if (file.exists(data_path)) {
      load(data_path, envir = .GlobalEnv)
      message(paste0(file_name, ".Rda loaded successfully"))
      return(TRUE)
    } else {
      warning(paste0(file_name, ".Rda not found. Skipping."))
      return(FALSE)
    }
  }

  # Load the specified data type
  if (data_type == "all") {
    # Load all data types
    load_success <- c(
      load_data_file("drug"),
      load_data_file("mRNA"),
      load_data_file("mut"),
      load_data_file("meth"),
      load_data_file("protein"),
      load_data_file("cnv"),
      load_data_file("fusion")
    )
  } else {
    # Load just the specified data type
    load_success <- load_data_file(data_type)
  }

  # Apply Z-score normalization if requested and applicable
  if (apply_normalization && any(load_success)) {
    continuous_data_types <- c("drug", "mrna", "meth", "protein", "cnv", "all")
    if (data_type %in% continuous_data_types) {
      # Source the Z-score normalization function if needed
      if (!exists("zscoreNormalize", mode = "function")) {
        if (requireNamespace("DROMA", quietly = TRUE)) {
          zscoreNormalize <- DROMA::zscoreNormalize
          applyZscoreNormalization <- DROMA::applyZscoreNormalization
        } else {
          # If package not installed, try to source from current directory
          if (file.exists("R/FuncZscoreWhole.R")) {
            source("R/FuncZscoreWhole.R")
          } else {
            warning("Could not find normalization functions. Skipping normalization.")
            return(invisible(TRUE))
          }
        }
      }

      # Apply normalization
      if (exists("applyZscoreNormalization", mode = "function")) {
        applyZscoreNormalization()
        message("Z-score normalization applied")
      }
    }
  } else if (apply_normalization == FALSE && data_type %in% c("drug", "all")) {
    tavor_drug <- tavor_drug_raw
    PDTXBreast_drug <- PDTXBreast_drug_raw
    ccle_drug <- ccle_drug_raw
    ctrp1_drug <- ctrp1_drug_raw
    ctrp2_drug <- ctrp2_drug_raw
    gdsc1_drug <- gdsc1_drug_raw
    gdsc2_drug <- gdsc2_drug_raw
    gCSI_drug <- gCSI_drug_raw
    prism_drug <- prism_drug_raw
    FIMM_drug <- FIMM_drug_raw
    UHNBreast_drug <- UHNBreast_drug_raw
    GRAY_drug <- GRAY_drug_raw
    NCI60_drug <- NCI60_drug_raw
    UMPDO1_drug <- UMPDO1_drug_raw
    UMPDO2_drug <- UMPDO2_drug_raw
    UMPDO3_drug <- UMPDO3_drug_raw
    Xeva_drug <- Xeva_drug_raw
    message("Z-score normalization not applied")
  }

  invisible(TRUE)
}
