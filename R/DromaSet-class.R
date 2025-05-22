#!/usr/bin/env Rscript

#' DromaSet Class
#'
#' @description A class to represent a DROMA project dataset with drug response and omics data.
#' Each DromaSet contains treatment response data (drug sensitivity) and multiple molecular profiles
#' (omics data) for a specific project.
#'
#' @slot name Character string, the name of the dataset (e.g., "gCSI", "CCLE")
#' @slot treatmentResponse List containing drug response data matrices
#' @slot molecularProfiles List containing different types of molecular data (omics)
#' @slot sampleMetadata Data frame with sample annotations
#' @slot treatmentMetadata Data frame with drug annotations
#' @slot datasetType Character string, the type of dataset (e.g., "CellLine", "PDX", "PDO")
#' @slot db_info List containing database connection information
#' @export
setClass("DromaSet",
  slots = c(
    name = "character",
    treatmentResponse = "list",
    molecularProfiles = "list",
    sampleMetadata = "data.frame",
    treatmentMetadata = "data.frame",
    datasetType = "character",
    db_info = "list"
  ),
  prototype = list(
    name = NA_character_,
    treatmentResponse = list(),
    molecularProfiles = list(),
    sampleMetadata = data.frame(),
    treatmentMetadata = data.frame(),
    datasetType = NA_character_,
    db_info = list()
  )
)

#' Create a DromaSet Object
#'
#' @description Creates a DromaSet object to store drug response and omics data for a specific project
#' @param name Character string, the name of the dataset (e.g., "gCSI", "CCLE")
#' @param treatmentResponse List containing drug response data matrices
#' @param molecularProfiles List containing different types of molecular data (omics)
#' @param sampleMetadata Data frame with sample annotations
#' @param treatmentMetadata Data frame with drug annotations
#' @param datasetType Character string, the type of dataset (e.g., "CellLine", "PDX", "PDO")
#' @param db_info List containing database connection information
#' @return A DromaSet object
#' @export
#' @examples
#' \dontrun{
#' # Create a DromaSet object for gCSI data
#' gCSI <- DromaSet(name = "gCSI",
#'                  datasetType = "CellLine",
#'                  db_info = list(db_path = "~/droma.sqlite", db_group = "gCSI"))
#' }
DromaSet <- function(name,
                   treatmentResponse = list(),
                   molecularProfiles = list(),
                   sampleMetadata = data.frame(),
                   treatmentMetadata = data.frame(),
                   datasetType = NA_character_,
                   db_info = list()) {

  # Create new DromaSet object
  object <- new("DromaSet",
              name = name,
              treatmentResponse = treatmentResponse,
              molecularProfiles = molecularProfiles,
              sampleMetadata = sampleMetadata,
              treatmentMetadata = treatmentMetadata,
              datasetType = datasetType,
              db_info = db_info)

  # Return the object
  return(object)
}

#' Show Method for DromaSet objects
#'
#' @description Displays information about a DromaSet object
#' @param object A DromaSet object
#' @return NULL, prints information to console
#' @export
setMethod("show", "DromaSet", function(object) {
  cat("DromaSet Object:", object@name, "\n")
  cat("Dataset Type:", object@datasetType, "\n\n")

  # Treatment Response
  cat("Treatment Response Data:\n")
  if (length(object@treatmentResponse) > 0) {
    tr_info <- lapply(object@treatmentResponse, function(x) {
      if(is.matrix(x) || is.data.frame(x)) {
        paste0("  (", nrow(x), " drugs x ", ncol(x), " samples)")
      } else {
        "  (data in database)"
      }
    })
    cat(paste(names(object@treatmentResponse), tr_info, sep = ": ", collapse = "\n"), "\n")
  } else {
    cat("  None loaded (may be available in database)\n")
  }

  # Molecular Profiles
  cat("\nMolecular Profiles:\n")
  if (length(object@molecularProfiles) > 0) {
    mp_info <- lapply(object@molecularProfiles, function(x) {
      if(is.matrix(x) || is.data.frame(x)) {
        paste0("  (", nrow(x), " features x ", ncol(x), " samples)")
      } else {
        "  (data in database)"
      }
    })
    cat(paste(names(object@molecularProfiles), mp_info, sep = ": ", collapse = "\n"), "\n")
  } else {
    cat("  None loaded (may be available in database)\n")
  }

  # Database information
  cat("\nDatabase Connection Information:\n")
  if (length(object@db_info) > 0) {
    cat("  Path:", ifelse(is.null(object@db_info$db_path), "Not specified", object@db_info$db_path), "\n")
    cat("  Group:", ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group), "\n")
  } else {
    cat("  No database information available\n")
  }

  invisible(NULL)
})

#' Get Available Molecular Profile Types
#'
#' @description Returns the types of molecular profiles available for a DromaSet
#' @param object A DromaSet object
#' @param include_db Logical, whether to include profiles available in the database but not loaded (default: TRUE)
#' @return Character vector of available molecular profile types
#' @export
setGeneric("availableMolecularProfiles", function(object, include_db = TRUE) standardGeneric("availableMolecularProfiles"))

#' @rdname availableMolecularProfiles
#' @export
setMethod("availableMolecularProfiles", "DromaSet", function(object, include_db = TRUE) {
  # Get profiles already loaded in the object
  loaded_profiles <- names(object@molecularProfiles)

  # If requested, check database for additional profiles
  if (include_db && length(object@db_info) > 0 && !is.null(object@db_info$db_path)) {
    if (file.exists(object@db_info$db_path)) {
      # Connect to database
      con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      # Get group prefix (if specified) or use dataset name
      group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

      # List all tables with this prefix
      all_tables <- DBI::dbListTables(con)
      project_tables <- grep(paste0("^", group_prefix, "_"), all_tables, value = TRUE)

      # Extract molecular profile types from table names
      if (length(project_tables) > 0) {
        # Remove prefix to get profile types
        db_profiles <- unique(sub(paste0("^", group_prefix, "_"), "", project_tables))
        # Remove 'drug' and 'drug_raw' which are treatment responses
        db_profiles <- db_profiles[!db_profiles %in% c("drug", "drug_raw")]

        # Combine with loaded profiles
        return(unique(c(loaded_profiles, db_profiles)))
      }
    }
  }

  return(loaded_profiles)
})

#' Get Available Treatment Response Types
#'
#' @description Returns the types of treatment response data available for a DromaSet
#' @param object A DromaSet object
#' @param include_db Logical, whether to include data available in the database but not loaded (default: TRUE)
#' @return Character vector of available treatment response types
#' @export
setGeneric("availableTreatmentResponses", function(object, include_db = TRUE) standardGeneric("availableTreatmentResponses"))

#' @rdname availableTreatmentResponses
#' @export
setMethod("availableTreatmentResponses", "DromaSet", function(object, include_db = TRUE) {
  # Get responses already loaded in the object
  loaded_responses <- names(object@treatmentResponse)

  # If requested, check database for additional responses
  if (include_db && length(object@db_info) > 0 && !is.null(object@db_info$db_path)) {
    if (file.exists(object@db_info$db_path)) {
      # Connect to database
      con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
      on.exit(DBI::dbDisconnect(con), add = TRUE)

      # Get group prefix (if specified) or use dataset name
      group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

      # Look for drug and drug_raw tables
      all_tables <- DBI::dbListTables(con)
      drug_tables <- grep(paste0("^", group_prefix, "_(drug|drug_raw)$"), all_tables, value = TRUE)

      # Extract response types from table names
      if (length(drug_tables) > 0) {
        # Remove prefix to get response types
        db_responses <- unique(sub(paste0("^", group_prefix, "_"), "", drug_tables))

        # Combine with loaded responses
        return(unique(c(loaded_responses, db_responses)))
      }
    }
  }

  return(loaded_responses)
})

