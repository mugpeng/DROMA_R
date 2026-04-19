source("../../R/FuncClinical.R")

test_that("clinical merged values use raw data when ComBat fails", {
  patient_data <- list(
    GSE1 = list(
      response = c(S1 = 1, S2 = 2),
      non_response = c(S3 = 3, S4 = 4)
    ),
    GSE2 = list(
      response = c(S5 = 5, S6 = 6),
      non_response = c(S7 = 7, S8 = 8)
    )
  )

  expect_warning(
    merged <- .prepareClinicalMergedValues(
      patient_data,
      normalization = "combat",
      combat_fun = function(...) stop("mock ComBat failure")
    ),
    "ComBat normalization failed"
  )

  expect_equal(unname(merged$response), c(1, 2, 5, 6))
  expect_equal(unname(merged$non_response), c(3, 4, 7, 8))
  expect_equal(merged$normalization_used, "none")
})

test_that("clinical merged values can z-score pooled data", {
  patient_data <- list(
    GSE1 = list(
      response = c(S1 = 1, S2 = 2),
      non_response = c(S3 = 3, S4 = 4)
    ),
    GSE2 = list(
      response = c(S5 = 5, S6 = 6),
      non_response = c(S7 = 7, S8 = 8)
    )
  )

  merged <- .prepareClinicalMergedValues(patient_data, normalization = "zscore")
  expected <- as.numeric(scale(1:8))

  expect_equal(unname(c(merged$response, merged$non_response)), expected[c(1, 2, 5, 6, 3, 4, 7, 8)])
  expect_equal(merged$normalization_used, "zscore")
})
