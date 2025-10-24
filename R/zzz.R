#' Package Startup Functions
#'
#' @description Functions that run when the package is attached or loaded
#' @keywords internal

.onAttach <- function(libname, pkgname) {
  # List of packages to attach
  packages_to_load <- c(
    "meta",
    "metafor",
    "effsize",
    "parallel",
    "snowfall",
    "ggplot2",
    "ggpubr",
    "dplyr",
    "stats",
    "DT",
    "htmltools",
    "utils",
    "patchwork",
    "ggrepel",
    "gg.gap",
    "gridExtra",
    "grid",
    "RSQLite",
    "DBI",
    "methods"
  )
  
  packageStartupMessage("Loading DROMA.R (version ", 
                        utils::packageVersion("DROMA.R"), ")")
  
  # Check for missing packages
  missing_packages <- c()
  for (pkg in packages_to_load) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  # If there are missing packages, stop with informative message
  if (length(missing_packages) > 0) {
    stop_msg <- paste(
      "The following required packages are not installed:\n  ",
      paste(missing_packages, collapse = ", "),
      "\n\nPlease install the missing packages and try loading DROMA.R again.",
      sep = ""
    )
    stop(stop_msg, call. = FALSE)
  }
  
  # Attach all required packages
  for (pkg in packages_to_load) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
  
  packageStartupMessage("All dependency packages loaded successfully!")
}

