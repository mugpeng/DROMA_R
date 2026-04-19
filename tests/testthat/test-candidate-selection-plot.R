source("../../R/FuncVisualization.R")

test_that("candidate selection plot add keeps specified genes", {
  candidate_df <- data.frame(
    name = c("A", "B", "C", "D"),
    effect_size_clinical = c(4, 3, 2, 1),
    effect_size_preclinical = c(4, 3, 2, 1),
    cell_supported = c(TRUE, TRUE, FALSE, FALSE),
    pdopdx_supported = c(TRUE, FALSE, TRUE, FALSE),
    clinical_supported = c(TRUE, TRUE, FALSE, FALSE),
    direction_concordant = c(TRUE, TRUE, TRUE, FALSE),
    retained = c(TRUE, FALSE, FALSE, FALSE)
  )

  plot_obj <- plotClinicallySupportedCandidateSelection(
    candidate_df,
    n_retained = 1,
    n_filtered = 1,
    add = "D"
  )

  expect_true("D" %in% as.character(plot_obj[[2]]$data$name_wrapped))
})

test_that("candidate selection effect plot uses effect size label", {
  candidate_df <- data.frame(
    name = c("A", "B"),
    effect_size_clinical = c(2, 1),
    effect_size_preclinical = c(2, 1),
    cell_supported = c(TRUE, FALSE),
    pdopdx_supported = c(TRUE, TRUE),
    clinical_supported = c(TRUE, FALSE),
    direction_concordant = c(TRUE, TRUE),
    retained = c(TRUE, FALSE)
  )

  plot_obj <- plotClinicallySupportedCandidateSelection(candidate_df)

  expect_equal(plot_obj[[2]]$labels$x, "effect size")
})
