# Clinical Trial Database (CTRDB) Visualization Functions ----

#' Plot clinical drug response for all datasets
#'
#' @description Creates and combines plots for all clinical datasets showing drug response differences
#' @param patient_data_list List of patient data
#' @param select_omics Name of the omics feature
#' @param select_drugs Name of the drug
#' @return A combined plot with all patient comparisons or NULL if no valid data
#' @export
plotAllClinicaldatasets <- function(patient_data_list, select_omics, select_drugs) {
  # Transform patient data to match expected format for plotMultipleGroupComparisons
  # Response = yes, Non_response = no
  transformed_pairs <- lapply(patient_data_list, function(patient_data) {
    list(
      yes = patient_data$response,      # Response samples
      no = patient_data$non_response    # Non-response samples
    )
  })
  names(transformed_pairs) <- names(patient_data_list)

  # Use existing plotMultipleGroupComparisons function from FuncPairVisualization
  combined_plot <- plotMultipleGroupComparisons(
    pairs_list = transformed_pairs,
    x_label = select_omics,
    y_label = "Expression",
    group_labels = c("NR", "R")
  )

  # Add overall title if plot is not NULL
  if (!is.null(combined_plot) && requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = paste(select_omics, "Expression vs", select_drugs, "Response (Clinical Data)"),
        theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
      )
  }

  return(combined_plot)
}

