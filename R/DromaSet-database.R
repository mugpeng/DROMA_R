#!/usr/bin/env Rscript

#' Load Molecular Profiles from Database
#'
#' @description Loads specific molecular profile data from the database into a DromaSet object
#' @param object A DromaSet object
#' @param molecular_type The type of molecular data to load (e.g., "mRNA", "cnv", "mutation_gene")
#' @param features Optional vector of feature names to load. If NULL, loads all features.
#' @param samples Optional vector of sample IDs to load. If NULL, loads all samples.
#' @return Updated DromaSet object with loaded molecular data
#' @export
setGeneric("loadMolecularProfiles", function(object, molecular_type, features = NULL, samples = NULL)
           standardGeneric("loadMolecularProfiles"))

#' @rdname loadMolecularProfiles
#' @export
setMethod("loadMolecularProfiles", "DromaSet", function(object, molecular_type, features = NULL, samples = NULL) {
  # Verify we have database connection info
  if (length(object@db_info) == 0 || is.null(object@db_info$db_path)) {
    stop("No database connection information available")
  }

  if (!file.exists(object@db_info$db_path)) {
    stop("Database file not found: ", object@db_info$db_path)
  }

  # Get group prefix (if specified) or use dataset name
  group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if the table exists
  table_name <- paste0(group_prefix, "_", molecular_type)
  all_tables <- DBI::dbListTables(con)

  if (!table_name %in% all_tables) {
    stop("Molecular profile type '", molecular_type, "' not found for dataset '", object@name, "'")
  }

  # Construct query based on features and samples
  if (molecular_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms")) {
    # For continuous data matrices
    query <- paste0("SELECT * FROM ", table_name)

    # Add feature filter if specified
    if (!is.null(features)) {
      features_str <- paste0("'", features, "'", collapse = ",")
      query <- paste0(query, " WHERE feature_id IN (", features_str, ")")
    }

    # Execute query
    data <- DBI::dbGetQuery(con, query)

    # Reshape to matrix format
    if (nrow(data) > 0) {
      # Extract feature_id column
      feature_ids <- data$feature_id
      data$feature_id <- NULL

      # Convert to matrix
      mat <- as.matrix(data)
      rownames(mat) <- feature_ids

      # Filter by samples if needed
      if (!is.null(samples)) {
        common_samples <- intersect(colnames(mat), samples)
        if (length(common_samples) == 0) {
          warning("No samples match the specified filter")
          mat <- matrix(nrow = 0, ncol = 0)
        } else {
          mat <- mat[, common_samples, drop = FALSE]
        }
      }

      # Store in object
      object@molecularProfiles[[molecular_type]] <- mat
    } else {
      warning("No data found for molecular profile type: ", molecular_type)
      object@molecularProfiles[[molecular_type]] <- matrix(nrow = 0, ncol = 0)
    }
  } else if (molecular_type %in% c("mutation_gene", "mutation_site", "fusion")) {
    # For discrete data (mutations, fusions)
    query <- paste0("SELECT * FROM ", table_name)

    # Add feature filter if specified
    if (!is.null(features)) {
      features_str <- paste0("'", features, "'", collapse = ",")
      query <- paste0(query, " WHERE gene IN (", features_str, ")")
    }

    # Execute query
    data <- DBI::dbGetQuery(con, query)

    # Filter by samples if needed
    if (!is.null(samples) && nrow(data) > 0) {
      data <- data[data$cells %in% samples, ]
    }

    # Store in object
    if (nrow(data) > 0) {
      object@molecularProfiles[[molecular_type]] <- data
    } else {
      warning("No data found for molecular profile type: ", molecular_type)
      object@molecularProfiles[[molecular_type]] <- data.frame()
    }
  } else {
    warning("Unrecognized molecular profile type: ", molecular_type)
  }

  return(object)
})

#' Load Treatment Response from Database
#'
#' @description Loads drug response data from the database into a DromaSet object
#' @param object A DromaSet object
#' @param response_type The type of response data to load (e.g., "drug" or "drug_raw")
#' @param drugs Optional vector of drug names to load. If NULL, loads all drugs.
#' @param samples Optional vector of sample IDs to load. If NULL, loads all samples.
#' @return Updated DromaSet object with loaded drug response data
#' @export
setGeneric("loadTreatmentResponse", function(object, response_type = "drug", drugs = NULL, samples = NULL)
           standardGeneric("loadTreatmentResponse"))

#' @rdname loadTreatmentResponse
#' @export
setMethod("loadTreatmentResponse", "DromaSet", function(object, response_type = "drug", drugs = NULL, samples = NULL) {
  # Verify we have database connection info
  if (length(object@db_info) == 0 || is.null(object@db_info$db_path)) {
    stop("No database connection information available")
  }

  if (!file.exists(object@db_info$db_path)) {
    stop("Database file not found: ", object@db_info$db_path)
  }

  # Get group prefix (if specified) or use dataset name
  group_prefix <- ifelse(is.null(object@db_info$db_group), object@name, object@db_info$db_group)

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), object@db_info$db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if the table exists
  table_name <- paste0(group_prefix, "_", response_type)
  all_tables <- DBI::dbListTables(con)

  if (!table_name %in% all_tables) {
    stop("Treatment response type '", response_type, "' not found for dataset '", object@name, "'")
  }

  # Construct query
  query <- paste0("SELECT * FROM ", table_name)

  # Add drug filter if specified
  if (!is.null(drugs)) {
    drugs_str <- paste0("'", drugs, "'", collapse = ",")
    query <- paste0(query, " WHERE feature_id IN (", drugs_str, ")")
  }

  # Execute query
  data <- DBI::dbGetQuery(con, query)

  # Reshape to matrix format
  if (nrow(data) > 0) {
    # Extract feature_id column (drug names)
    drug_ids <- data$feature_id
    data$feature_id <- NULL

    # Convert to matrix
    mat <- as.matrix(data)
    rownames(mat) <- drug_ids

    # Filter by samples if needed
    if (!is.null(samples)) {
      common_samples <- intersect(colnames(mat), samples)
      if (length(common_samples) == 0) {
        warning("No samples match the specified filter")
        mat <- matrix(nrow = 0, ncol = 0)
      } else {
        mat <- mat[, common_samples, drop = FALSE]
      }
    }

    # Store in object
    object@treatmentResponse[[response_type]] <- mat
  } else {
    warning("No data found for treatment response type: ", response_type)
    object@treatmentResponse[[response_type]] <- matrix(nrow = 0, ncol = 0)
  }

  return(object)
})

#' Create a DromaSet from Database
#'
#' @description Creates a DromaSet object linked to data in a SQLite database
#' @param project_name The name of the project/dataset (e.g., "gCSI", "CCLE")
#' @param db_path Path to the SQLite database
#' @param db_group Optional group name in the database, if different from project_name
#' @param load_metadata Logical, whether to load sample and treatment metadata (default: TRUE)
#' @param dataset_type Optional dataset type (e.g., "CellLine", "PDX", "PDO")
#' @return A DromaSet object linked to the database
#' @export
#' @examples
#' \dontrun{
#' # Create a DromaSet for gCSI data from database
#' gCSI <- createDromaSetFromDatabase("gCSI", "~/droma.sqlite")
#' }
createDromaSetFromDatabase <- function(project_name, db_path = file.path(path.expand("~"), "droma.sqlite"),
                                    db_group = NULL, load_metadata = TRUE, dataset_type = NULL) {

  if (!file.exists(db_path)) {
    stop("Database file not found: ", db_path)
  }

  # Set db_group to project_name if not specified
  if (is.null(db_group)) {
    db_group <- project_name
  }

  # Connect to database
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Check if tables exist for this project
  all_tables <- DBI::dbListTables(con)
  project_tables <- grep(paste0("^", db_group, "_"), all_tables, value = TRUE)

  if (length(project_tables) == 0) {
    stop("No tables found for project '", project_name, "' with group prefix '", db_group, "'")
  }

  # Load metadata if requested
  sample_metadata <- data.frame()
  treatment_metadata <- data.frame()

  if (load_metadata) {
    # Try to load sample metadata
    if ("sample_anno" %in% all_tables) {
      # Construct query to get samples for this project
      sample_query <- paste0(
        "SELECT * FROM sample_anno WHERE SampleID IN (",
        "SELECT DISTINCT name FROM (",
        paste(lapply(project_tables, function(t) {
          paste0("SELECT name FROM pragma_table_info('", t, "') ",
                "WHERE name != 'feature_id'")
        }), collapse = " UNION "),
        "))"
      )

      tryCatch({
        sample_metadata <- DBI::dbGetQuery(con, sample_query)
        sample_metadata <- sample_metadata[sample_metadata$ProjectID %in% project_name,]
        sample_metadata <- unique(sample_metadata)
      }, error = function(e) {
        warning("Error loading sample metadata: ", e$message)
      })
    }

    # Try to load treatment metadata
    if ("drug_anno" %in% all_tables) {
      # Construct query to get drugs for this project
      drug_query <- paste0(
        "SELECT * FROM drug_anno WHERE DrugName IN (",
        "SELECT feature_id FROM ", db_group, "_drug)"
      )

      tryCatch({
        treatment_metadata <- DBI::dbGetQuery(con, drug_query)
        treatment_metadata <- treatment_metadata[treatment_metadata$ProjectID %in% project_name,]
        treatment_metadata <- unique(treatment_metadata)
      }, error = function(e) {
        warning("Error loading treatment metadata: ", e$message)
      })
    }
  }

  # If dataset_type is not provided, try to infer from sample metadata
  if (is.null(dataset_type) && nrow(sample_metadata) > 0 && "DataType" %in% colnames(sample_metadata)) {
    dataset_type <- unique(sample_metadata$DataType)[1]
  }

  # Create the DromaSet object
  object <- DromaSet(
    name = project_name,
    sampleMetadata = sample_metadata,
    treatmentMetadata = treatment_metadata,
    datasetType = ifelse(is.null(dataset_type), NA_character_, dataset_type),
    db_info = list(
      db_path = db_path,
      db_group = db_group
    )
  )

  return(object)
}

#' Convert Rda Data to DromaSet Objects in Database
#'
#' @description Converts original DROMA Rda files to a structure of DromaSet objects in a SQLite database
#' @param db_path Path where the SQLite database file should be created
#' @param rda_dir Directory containing the Rda files to convert
#' @param projects Optional vector of project names to include. If NULL, includes all.
#' @return Invisibly returns the path to the created database
#' @export
#' @examples
#' \dontrun{
#' # Convert all Rda files to DromaSet database
#' convertRdaToDromaSetDatabase()
#' }
convertRdaToDromaSetDatabase <- function(db_path = file.path(path.expand("~"), "droma.sqlite"),
                                      rda_dir = system.file("data", package = "DROMA"),
                                      projects = NULL) {

  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Please install them with install.packages(c('RSQLite', 'DBI'))")
  }

  # Create database connection
  message("Creating database at ", db_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Get all Rda files
  rda_files <- list.files(rda_dir, pattern = "\\.Rda$", full.names = TRUE)

  # First load annotation files
  anno_file <- file.path(rda_dir, "anno.Rda")
  if (file.exists(anno_file)) {
    message("Processing annotations...")
    e <- new.env()
    load(anno_file, envir = e)

    # Write annotation tables
    if (exists("sample_anno", envir = e)) {
      DBI::dbWriteTable(con, "sample_anno", e$sample_anno, overwrite = TRUE)
      message("  - Wrote sample annotations")
    }

    if (exists("drug_anno", envir = e)) {
      DBI::dbWriteTable(con, "drug_anno", e$drug_anno, overwrite = TRUE)
      message("  - Wrote drug annotations")
    }
  }

  # Load search vectors if available
  search_file <- file.path(rda_dir, "search_vec.Rda")
  if (file.exists(search_file)) {
    message("Processing search vectors...")
    e <- new.env()
    load(search_file, envir = e)

    # Store search vectors as serialized objects
    for (obj_name in ls(e)) {
      obj <- get(obj_name, envir = e)
      if (is.list(obj)) {
        df <- data.frame(
          name = obj_name,
          value = I(list(serialize(obj, NULL))),
          stringsAsFactors = FALSE
        )
        DBI::dbWriteTable(con, "search_vectors", df, append = TRUE)
        message("  - Stored search vector: ", obj_name)
      }
    }
  }

  # Process each data file
  for (rda_file in rda_files) {
    file_name <- basename(rda_file)

    # Skip already processed files
    if (file_name %in% c("anno.Rda", "search_vec.Rda")) {
      next
    }

    message("Processing ", file_name)
    data_type <- tools::file_path_sans_ext(file_name)

    # Load the Rda file in a new environment
    e <- new.env()
    load(rda_file, envir = e)

    # Process each object in the environment
    for (obj_name in ls(e)) {
      obj <- get(obj_name, envir = e)

      # Extract project name from object name
      parts <- strsplit(obj_name, "_")[[1]]
      if (length(parts) < 2) {
        message("  Skipping ", obj_name, " (not in project_datatype format)")
        next
      }

      project_name <- parts[1]
      obj_type <- paste(parts[-1], collapse = "_")  # In case there are underscores in the type

      # Skip if projects is specified and this project is not in the list
      if (!is.null(projects) && !project_name %in% projects) {
        next
      }

      # Table name is project_datatype
      table_name <- paste0(project_name, "_", obj_type)
      message("  - Processing ", table_name)

      if (is.matrix(obj)) {
        # Convert matrix to data frame with feature_id column
        df <- as.data.frame(obj)
        df$feature_id <- rownames(df)

        # Write to database
        DBI::dbWriteTable(con, table_name, df, overwrite = TRUE)

        # Create index on feature_id for faster lookups
        DBI::dbExecute(con, paste0("CREATE INDEX idx_", table_name, "_feature_id ON ",
                                  table_name, " (feature_id)"))
      } else if (is.data.frame(obj)) {
        # Check if the data frame has rownames and preserve them
        if (!is.null(rownames(obj)) && !all(rownames(obj) == seq_len(nrow(obj)))) {
          # If data frame has meaningful rownames, add them as feature_id column
          df <- obj
          df$feature_id <- rownames(df)

          # Write to database
          DBI::dbWriteTable(con, table_name, df, overwrite = TRUE)

          # Create index on feature_id for faster lookups
          DBI::dbExecute(con, paste0("CREATE INDEX idx_", table_name, "_feature_id ON ",
                                    table_name, " (feature_id)"))
        } else {
          # If no meaningful rownames, write directly
          DBI::dbWriteTable(con, table_name, obj, overwrite = TRUE)
        }
      } else {
        message("    Skipping (unsupported type: ", class(obj)[1], ")")
      }
    }
  }

  # Create a projects table with metadata
  project_names <- unique(sapply(DBI::dbListTables(con), function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 2 && !t %in% c("sample_anno", "drug_anno", "search_vectors")) {
      return(parts[1])
    } else {
      return(NA)
    }
  }))
  project_names <- project_names[!is.na(project_names)]

  # Get data types for each project
  project_metadata <- lapply(project_names, function(proj) {
    tables <- grep(paste0("^", proj, "_"), DBI::dbListTables(con), value = TRUE)
    data_types <- unique(sapply(tables, function(t) {
      sub(paste0("^", proj, "_"), "", t)
    }))

    # Try to determine dataset type
    dataset_type <- NA_character_
    if ("sample_anno" %in% DBI::dbListTables(con)) {
      # Get sample IDs for this project
      sample_ids <- DBI::dbGetQuery(con, paste0(
        "SELECT name FROM (",
        paste(lapply(tables, function(t) {
          paste0("SELECT column_name AS name FROM pragma_table_info('", t, "') ",
                "WHERE column_name != 'feature_id'")
        }), collapse = " UNION "),
        ")"
      ))$name

      if (length(sample_ids) > 0) {
        # Look up dataset type in sample_anno
        types_query <- paste0(
          "SELECT DISTINCT DataType FROM sample_anno WHERE SampleID IN ('",
          paste(sample_ids, collapse = "','"),
          "')"
        )
        types <- DBI::dbGetQuery(con, types_query)$DataType

        if (length(types) > 0) {
          dataset_type <- types[1]
        }
      }
    }

    data.frame(
      project_name = proj,
      dataset_type = dataset_type,
      data_types = paste(data_types, collapse = ","),
      sample_count = length(unique(sample_ids)),
      stringsAsFactors = FALSE
    )
  })

  project_table <- do.call(rbind, project_metadata)
  DBI::dbWriteTable(con, "projects", project_table, overwrite = TRUE)

  message("Database creation complete. Database contains ", length(DBI::dbListTables(con)), " tables",
         " for ", length(project_names), " projects.")
  invisible(db_path)
}
