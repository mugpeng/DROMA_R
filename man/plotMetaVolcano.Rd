% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FuncBatchFeature.R
\name{plotMetaVolcano}
\alias{plotMetaVolcano}
\title{Create a volcano plot from meta-analysis results}
\usage{
plotMetaVolcano(
  meta_df,
  es_t = 0.4,
  P_t = 0.001,
  label = TRUE,
  top_label_each = 5,
  label_size = 5,
  point_size = 2.5,
  point_alpha = 0.6,
  title = NULL,
  p_adj_method = "none",
  custom_colors = NULL
)
}
\arguments{
\item{meta_df}{Data frame containing meta-analysis results with columns:
effect_size, p_value, and name}

\item{es_t}{Effect size threshold to consider significant}

\item{P_t}{P-value threshold to consider significant}

\item{label}{Whether to add labels to top points (TRUE/FALSE)}

\item{top_label_each}{Number of top points in each direction to label}

\item{label_size}{Size of text labels}

\item{point_size}{Size of points}

\item{point_alpha}{Alpha transparency of points}

\item{title}{Plot title (NULL for no title)}

\item{p_adj_method}{Method for p-value adjustment ("none", "BH", "bonferroni")}

\item{custom_colors}{Custom color vector for Up, NS, Down (NULL for defaults)}
}
\value{
ggplot object with volcano plot
}
\description{
Create a volcano plot from meta-analysis results
}
