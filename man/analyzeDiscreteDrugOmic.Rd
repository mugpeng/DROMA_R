% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuncDrugOmicPair.R
\name{analyzeDiscreteDrugOmic}
\alias{analyzeDiscreteDrugOmic}
\title{Analyze discrete drug-omics pairs}
\usage{
analyzeDiscreteDrugOmic(myPairs)
}
\arguments{
\item{myPairs}{List of paired drug-omics datasets from pairDiscreteDrugOmic}
}
\value{
Meta-analysis results object or NULL if analysis couldn't be performed
}
\description{
Performs meta-analysis on discrete drug-omic pairs using Wilcoxon test and Cliff's Delta
}
