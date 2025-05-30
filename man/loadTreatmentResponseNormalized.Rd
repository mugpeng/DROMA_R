% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuncDataLoaders.R
\name{loadTreatmentResponseNormalized}
\alias{loadTreatmentResponseNormalized}
\title{Load treatment response data with optional z-score normalization}
\usage{
loadTreatmentResponseNormalized(
  dromaset_object,
  drugs = NULL,
  data_type = "all",
  tumor_type = "all",
  zscore = TRUE,
  return_data = TRUE
)
}
\arguments{
\item{dromaset_object}{A DromaSet object}

\item{drugs}{Optional vector of specific drugs to load. If NULL, loads all drugs}

\item{data_type}{Filter by data type ("all", "CellLine", "PDO", "PDC", "PDX")}

\item{tumor_type}{Filter by tumor type ("all" or specific tumor types)}

\item{zscore}{Logical, whether to apply z-score normalization (default: TRUE)}

\item{return_data}{Logical, whether to return the data directly (default: TRUE)}
}
\value{
Matrix with treatment response data, optionally z-score normalized
}
\description{
Wrapper function for loadTreatmentResponse that applies z-score normalization by default
}
\examples{
\dontrun{
# Load drug data with z-score normalization (default)
gCSI <- createDromaSetFromDatabase("gCSI", "path/to/droma.sqlite")
drug_data <- loadTreatmentResponseNormalized(gCSI, drugs = "Paclitaxel")

# Load without z-score normalization
drug_raw <- loadTreatmentResponseNormalized(gCSI, drugs = "Paclitaxel", zscore = FALSE)

# Load all drug data with normalization
all_drugs <- loadTreatmentResponseNormalized(gCSI)
}
}
