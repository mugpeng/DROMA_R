#!/usr/bin/env Rscript

#' Create DROMA SQLite Database
#'
#' @description Converts all DROMA data files to a SQLite database with project-oriented structure
#' @param db_path Path where the SQLite database file should be created. Default is "droma.sqlite" in the user's home directory.
#' @param rda_dir Directory containing the Rda files to convert (default is the data directory in the package)
#' @param projects Optional vector of project names to include. If NULL, includes all.
#' @return Invisibly returns the path to the created database
#' @export
#' @examples
#' \dontrun{
#' createDROMADatabase()
#' # Creates a SQLite database with all DROMA data organized by project
#' }
createDROMADatabase <- function(db_path = file.path(path.expand("~"), "droma.sqlite"),
                               rda_dir = system.file("data", package = "DROMA"),
                               projects = NULL) {
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Please install them with install.packages(c('RSQLite', 'DBI'))")
  }
  # db_path = "sql_db/droma.sqlite"
  # Create database connection
  message("Creating project-oriented database at ", db_path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Get all Rda files
  # rda_dir = "data"
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
  for (rda_file in rda_files[4:9]) {
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
        obj <- as.data.frame(obj)
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
          paste0("SELECT name FROM pragma_table_info('", t, "') ",
                "WHERE name != 'feature_id'")
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

#' Connect to DROMA Database
#'
#' @description Establishes a connection to the DROMA SQLite database
#' @param db_path Path to the SQLite database file
#' @return A database connection object
#' @export
connectDROMADatabase <- function(db_path = file.path(path.expand("~"), "droma.sqlite")) {
  if (!requireNamespace("RSQLite", quietly = TRUE) ||
      !requireNamespace("DBI", quietly = TRUE)) {
    stop("Packages 'RSQLite' and 'DBI' are required. Please install them with install.packages(c('RSQLite', 'DBI'))")
  }

  if (!file.exists(db_path)) {
    stop("Database file not found. Create the database first with createDROMADatabase()")
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  message("Connected to DROMA database at ", db_path)

  # Store the connection in the package environment
  assign("droma_db_connection", con, envir = .GlobalEnv)

  # Register an exit handler to close the connection when R exits
  reg.finalizer(.GlobalEnv, function(e) {
    if (exists("droma_db_connection", envir = e) &&
        inherits(get("droma_db_connection", envir = e), "DBIConnection")) {
      DBI::dbDisconnect(get("droma_db_connection", envir = e))
    }
  }, onexit = TRUE)

  return(con)
}

#' Retrieve Feature Data from DROMA Database
#'
#' @description Fetches specific feature data from the DROMA database based on selection criteria
#' @param select_feas_type The type of feature to select (e.g., "mRNA", "cnv", "drug")
#' @param select_feas The specific feature to select within the feature type
#' @param data_sources Vector of data sources to select from (e.g., c("ccle", "gdsc"))
#' @param data_type Filter by data type: "all" (default), "CellLine", "PDO", "PDC", or "PDX"
#' @param tumor_type Filter by tumor type: "all" (default) or specific tumor type
#' @param connection Optional database connection object. If NULL, uses global connection.
#' @return A list of selected features from specified data sources
#' @export
#' @note This function is provided for backward compatibility. For new code, consider
#' using the DromaSet object approach instead.
getFeatureFromDatabase <- function(select_feas_type, select_feas,
                                 data_sources = "all",
                                 data_type = "all", tumor_type = "all",
                                 connection = NULL) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required. Please install with install.packages('DBI')")
  }

  # Get connection from global environment if not provided
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Get data source tables that match the feature type
  all_tables <- DBI::dbListTables(connection)
  pattern <- paste0("_", select_feas_type, "$")
  feature_tables <- grep(pattern, all_tables, value = TRUE)

  if (length(feature_tables) == 0) {
    stop("No tables found for feature type: ", select_feas_type)
  }

  # Filter by data sources if specified
  if (!identical(data_sources, "all")) {
    feature_tables <- grep(paste0("^(", paste(data_sources, collapse = "|"), ")_"),
                           feature_tables, value = TRUE)
  }

  if (length(feature_tables) == 0) {
    stop("No matching tables found for the specified data sources")
  }

  # Get sample IDs from sample_anno based on data_type and tumor_type
  filtered_samples <- NULL
  if (data_type != "all" || tumor_type != "all") {
    # Construct SQL query for sample filtering
    sample_query <- "SELECT SampleID FROM sample_anno WHERE 1=1"

    if (data_type != "all") {
      sample_query <- paste0(sample_query, " AND DataType = '", data_type, "'")
    }

    if (tumor_type != "all") {
      sample_query <- paste0(sample_query, " AND TumorType = '", tumor_type, "'")
    }

    # Execute the query
    filtered_samples <- DBI::dbGetQuery(connection, sample_query)$SampleID

    if (length(filtered_samples) == 0) {
      stop("No samples match the specified data_type and tumor_type criteria")
    }
  }

  # Retrieve data for each table
  result_list <- list()

  for (table in feature_tables) {
    # Extract data source name from table name
    data_source <- sub(paste0("_", select_feas_type, "$"), "", table)

    # Query for the specified feature
    if (select_feas_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
      # For continuous data, get the row for the feature
      query <- paste0("SELECT * FROM ", table, " WHERE feature_id = '", select_feas, "'")
      feature_data <- DBI::dbGetQuery(connection, query)

      if (nrow(feature_data) == 0) {
        next  # Skip if feature not found
      }

      # Convert to vector format (excluding feature_id column)
      feature_vector <- as.numeric(as.vector(feature_data[1, -which(names(feature_data) == "feature_id")]))
      names(feature_vector) <- colnames(feature_data)[-which(names(feature_data) == "feature_id")]
    } else {
      # For discrete data like mutations, get sample IDs where feature is present
      query <- paste0("SELECT cells FROM ", table, " WHERE gene = '", select_feas, "'")
      feature_data <- DBI::dbGetQuery(connection, query)

      if (nrow(feature_data) == 0) {
        next  # Skip if feature not found
      }

      feature_vector <- feature_data$cells
    }

    # Filter by samples if needed
    if (!is.null(filtered_samples)) {
      if (select_feas_type %in% c("mRNA", "cnv", "meth", "proteinrppa", "proteinms", "drug", "drug_raw")) {
        common_samples <- intersect(names(feature_vector), filtered_samples)
        if (length(common_samples) == 0) {
          next  # Skip if no samples match the filter
        }
        feature_vector <- feature_vector[common_samples]
      } else {
        feature_vector <- intersect(feature_vector, filtered_samples)
        if (length(feature_vector) == 0) {
          next  # Skip if no samples match the filter
        }
      }
    }

    # Add to result list
    result_list[[data_source]] <- feature_vector
  }

  if (length(result_list) == 0) {
    stop("No data found for feature '", select_feas, "' with the specified criteria")
  }

  return(result_list)
}

#' Close DROMA Database Connection
#'
#' @description Closes the connection to the DROMA database
#' @param connection Optional database connection object. If NULL, uses global connection.
#' @return TRUE if successfully disconnected
#' @export
closeDROMADatabase <- function(connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      message("No open database connection found")
      return(invisible(FALSE))
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  DBI::dbDisconnect(connection)

  if (exists("droma_db_connection", envir = .GlobalEnv)) {
    rm("droma_db_connection", envir = .GlobalEnv)
  }

  message("Database connection closed")
  invisible(TRUE)
}

#' List Available Tables in DROMA Database
#'
#' @description Provides information about tables available in the DROMA database
#' @param pattern Optional regex pattern to filter table names
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return A data frame with table information
#' @export
listDROMADatabaseTables <- function(pattern = NULL, connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Get table list
  tables <- DBI::dbListTables(connection)

  # Filter by pattern if provided
  if (!is.null(pattern)) {
    tables <- grep(pattern, tables, value = TRUE)
  }

  # Get metadata for each table
  if ("droma_metadata" %in% DBI::dbListTables(connection)) {
    metadata <- DBI::dbReadTable(connection, "droma_metadata")
    result <- metadata[metadata$table_name %in% tables, ]
  } else {
    # If metadata table doesn't exist, create basic info
    result <- data.frame(
      table_name = tables,
      row_count = sapply(tables, function(t)
        DBI::dbGetQuery(connection, paste0("SELECT COUNT(*) FROM ", t))[1,1]),
      stringsAsFactors = FALSE
    )
  }

  # Add categorization by data type
  result$data_type <- sub("_.*$", "", result$table_name)
  result$feature_type <- sub("^.*_", "", result$table_name)

  return(result)
}

#' List Available Projects in DROMA Database
#'
#' @description Lists all projects available in the DROMA database
#' @param connection Optional database connection object. If NULL, uses global connection
#' @return A data frame with project information
#' @export
listDROMARojects <- function(connection = NULL) {
  if (is.null(connection)) {
    if (!exists("droma_db_connection", envir = .GlobalEnv)) {
      stop("No database connection found. Connect first with connectDROMADatabase()")
    }
    connection <- get("droma_db_connection", envir = .GlobalEnv)
  }

  # Check if projects table exists
  if ("projects" %in% DBI::dbListTables(connection)) {
    return(DBI::dbReadTable(connection, "projects"))
  }

  # Otherwise, try to infer projects from table names
  tables <- DBI::dbListTables(connection)

  # Extract project names from table prefixes
  project_names <- unique(sapply(tables, function(t) {
    parts <- strsplit(t, "_")[[1]]
    if (length(parts) >= 2 && !t %in% c("sample_anno", "drug_anno", "droma_metadata", "search_vectors")) {
      return(parts[1])
    } else {
      return(NA)
    }
  }))
  project_names <- project_names[!is.na(project_names)]

  if (length(project_names) == 0) {
    message("No projects found in database")
    return(data.frame())
  }

  # Create a data frame with project information
  result <- data.frame(
    project_name = project_names,
    stringsAsFactors = FALSE
  )

  return(result)
}

