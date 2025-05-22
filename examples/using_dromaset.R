#!/usr/bin/env Rscript

# Using the DromaSet Class for DROMA
# This example demonstrates how to use the new project-oriented DromaSet class
# for more efficient data access and organization

library(DROMA)

# ---------- Setup ---------- #

# First-time setup: Create the database from .Rda files (only needed once)
if (!file.exists("~/droma.sqlite")) {
  cat("Creating DROMA database (this may take a few minutes)...\n")
  createDROMADatabase()  # Directly creates project-oriented database
}

# ---------- List Available Projects ---------- #

# Connect to the database
db_conn <- connectDROMADatabase("sql_db/droma.sqlite")

# View all projects in the database
projects <- listDROMARojects()
cat("Available projects in database:", nrow(projects), "\n")
print(projects)

# ---------- Create DromaSet Objects ---------- #

# Create a DromaSet object for gCSI data
gCSI <- createDromaSetFromDatabase("gCSI", db_path = "sql_db/droma.sqlite")

# Display information about the DromaSet
print(gCSI)

# Check what molecular profiles are available
profiles <- availableMolecularProfiles(gCSI)
cat("\nAvailable molecular profiles for gCSI:\n")
print(profiles)

# Check what treatment response data is available
responses <- availableTreatmentResponses(gCSI)
cat("\nAvailable treatment response data for gCSI:\n")
print(responses)

# ---------- Load Specific Data ---------- #

# Load drug response data
gCSI <- loadTreatmentResponse(gCSI, response_type = "drug")

# Load Copy Number data for specific genes
gCSI <- loadMolecularProfiles(
  gCSI,
  molecular_type = "cnv",
  features = c("BRCA1", "BRCA2", "TP53")
)

# ---------- Example Analysis ---------- #

# Analyze relationship between BRCA1 copy number and drug response
if ("cnv" %in% names(gCSI@molecularProfiles) && "drug" %in% names(gCSI@treatmentResponse)) {
  # Get BRCA1 CNV data
  brca1_cnv <- gCSI@molecularProfiles$cnv["BRCA1", ]

  # Get response to a specific drug (e.g., first drug in the matrix)
  drug_name <- rownames(gCSI@treatmentResponse$drug)[1]
  drug_response <- gCSI@treatmentResponse$drug[drug_name, ]

  # Find common samples
  common_samples <- intersect(names(brca1_cnv), names(drug_response))

  # Calculate correlation if there are enough samples
  if (length(common_samples) >= 10) {
    correlation <- cor.test(
      brca1_cnv[common_samples],
      drug_response[common_samples],
      method = "pearson"
    )

    cat("\nCorrelation between BRCA1 CNV and", drug_name, "response:\n")
    cat("r =", round(correlation$estimate, 3),
        ", p-value =", round(correlation$p.value, 5), "\n")
  } else {
    cat("\nNot enough common samples for correlation analysis\n")
  }
}

# ---------- Multiple Projects ---------- #

# Create DromaSet objects for multiple projects
ccle <- createDromaSetFromDatabase("ccle")
gdsc <- createDromaSetFromDatabase("gdsc")

# Load the same molecular data from each project
projects <- list(gCSI = gCSI, CCLE = ccle, GDSC = gdsc)

# For each project, load TP53 mutation data
for (proj_name in names(projects)) {
  tryCatch({
    projects[[proj_name]] <- loadMolecularProfiles(
      projects[[proj_name]],
      molecular_type = "mutation_gene",
      features = "TP53"
    )
    cat("\nLoaded TP53 mutation data for", proj_name, "\n")
  }, error = function(e) {
    cat("\nCould not load TP53 mutation data for", proj_name, ":", e$message, "\n")
  })
}

# ---------- Clean Up ---------- #

# Close the database connection when done
closeDROMADatabase()

cat("\nExample completed successfully\n")
