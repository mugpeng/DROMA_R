% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuncStatMultiDromaSet.R
\name{generateOverlapPlots}
\alias{generateOverlapPlots}
\title{Generate overlap plots for samples and drugs from database}
\usage{
generateOverlapPlots(selected_projects, connection)
}
\arguments{
\item{selected_projects}{Character vector of project names}

\item{connection}{Database connection object}
}
\value{
List of UpSet plot objects
}
\description{
Creates UpSet plots showing overlaps between projects using database data
}
