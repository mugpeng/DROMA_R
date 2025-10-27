# Meta-analysis Functions Module ----

#' Meta-analysis for continuous vs continuous features
#'
#' @description Performs meta-analysis on correlation between two continuous features
#' @param selected_pair List of paired data with two continuous variables
#' @return Meta-analysis result object or NULL if insufficient data
metaCalcConCon <- function(selected_pair) {
  if (length(selected_pair) < 1) return(NULL)

  # Test pairs one by one
  cal_list <- lapply(1:length(selected_pair), function(y) {
    fea1_sel <- selected_pair[[y]][[1]]  # Feature 1 values (e.g., drug sensitivity)
    fea2_sel <- selected_pair[[y]][[2]]  # Feature 2 values (e.g., gene expression)

    # Check for minimum length
    if (length(fea1_sel) < 3 || length(fea2_sel) < 3) return(NULL)

    # Spearman correlation: positive means both increase together
    cor_re <- tryCatch(
      suppressWarnings(cor.test(fea1_sel, fea2_sel, method = "spearman")),
      error = function(x) { return(NULL) }
    )

    if (is.null(cor_re)) return(NULL)

    data.frame(
      study = names(selected_pair)[y],  # Study name
      p = cor_re$p.value,
      effect = cor_re$estimate,  # Correlation coefficient: positive means positive association
      N = length(fea2_sel)
    )
  })

  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)
  rownames(cal_re) <- cal_re$study
  
  # Handle extreme correlation values to prevent infinite z-scores
  # Clamp correlation coefficients to avoid numerical issues in Fisher's z transform
  cal_re$effect <- pmax(pmin(cal_re$effect, 0.999), -0.999)
  
  cal_re$se <- sqrt((1 - cal_re$effect^2) / (cal_re$N - 2))
  cal_re$z <- 0.5 * log((1 + cal_re$effect) / (1 - cal_re$effect))  # Fisher's z
  cal_re$se_z <- 1 / sqrt(cal_re$N - 3)

  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = z, seTE = se_z, studlab = study, data = cal_re, sm = "Z",
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )}),
    error = function(x) { return(NULL) }
  )
  cal_meta_re
}

#' Meta-analysis for continuous vs discrete features
#'
#' @description Performs meta-analysis comparing continuous values between discrete groups
#' @param selected_pair List of paired data with continuous and discrete variables
#' @return Meta-analysis result object or NULL if insufficient data
metaCalcConDis <- function(selected_pair) {
  if (length(selected_pair) < 1) return(NULL)

  cal_list <- lapply(1:length(selected_pair), function(y) {
    yes_drugs <- selected_pair[[y]][[1]]  # Samples with feature present (e.g., mutated)
    no_drugs <- selected_pair[[y]][[2]]   # Samples without feature (e.g., wild-type)

    # Check for minimum length
    if (length(yes_drugs) < 3 || length(no_drugs) < 3) return(NULL)

    # Wilcox test: yes vs no (positive effect means yes > no)
    wilcox_re <- tryCatch(
      suppressWarnings(wilcox.test(yes_drugs, no_drugs)),
      error = function(x) { return(NULL) }
    )
    if (is.null(wilcox_re)) return(NULL)

    # Cliff's Delta: yes vs no (positive effect means yes > no)
    # This ensures that when feature presence (e.g., mutation) increases values,
    # the effect size is positive
    cliff_delta <- tryCatch(
      suppressWarnings(cliff.delta(yes_drugs, no_drugs)),
      error = function(x) { return(NULL) }
    )
    if (is.null(cliff_delta)) return(NULL)

    data.frame(
      study = names(selected_pair)[y],  # Study name
      p = wilcox_re$p.value,
      effect = cliff_delta$estimate,
      N = length(yes_drugs) + length(no_drugs),
      n1 = length(yes_drugs),
      n2 = length(no_drugs)
    )
  })

  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)
  rownames(cal_re) <- cal_re$study
  # Calculate standard error for Cliff's Delta
  cal_re$se <- sqrt((1 - cal_re$effect^2) * (cal_re$n1 + cal_re$n2 + 1) /
                      (12 * cal_re$n1 * cal_re$n2))

  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = effect,
              seTE = se,
              studlab = study,
              data = cal_re,
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001),
              sm = "CMD"  # Custom Mean Difference (using Cliff's Delta)
      )
    }),
    error = function(x) { return(NULL) }
  )
  cal_meta_re
}

#' Meta-analysis for discrete vs discrete features
#'
#' @description Performs meta-analysis on two discrete features using odds ratio
#' @param selected_pair List of paired data with contingency tables
#' @return Meta-analysis result object or NULL if insufficient data
metaCalcDisDis <- function(selected_pair) {
  # Check if we have enough pairs for meta-analysis
  if (length(selected_pair) < 1) return(NULL)

  # Calculate statistics for each pair
  cal_list <- lapply(1:length(selected_pair), function(y) {
    cont_table <- selected_pair[[y]]$cont_table

    # Skip if any cell has too few observations (e.g., < 3)
    if (any(cont_table < 3)) return(NULL)

    # Calculate odds ratio and its standard error
    tryCatch({
      # Extract values from contingency table
      # Contingency table structure:
      #              Feature2
      #              Yes    No
      # Feature1 Yes  a     b
      #          No   c     d
      a <- cont_table[1,1] # both features present (e.g., both mutations)
      b <- cont_table[1,2] # feature1 yes, feature2 no
      c <- cont_table[2,1] # feature1 no, feature2 yes
      d <- cont_table[2,2] # both features absent

      # Calculate log odds ratio and its standard error
      # OR = (a*d)/(b*c): positive association → OR > 1 → log(OR) > 0
      log_or <- log((a * d)/(b * c))
      se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)

      # Calculate Fisher's exact test p-value
      fisher_test <- fisher.test(cont_table)

      data.frame(
        study = names(selected_pair)[y],  # Study name
        log_or = log_or,  # Positive value means positive association between features
        se = se_log_or,
        p = fisher_test$p.value,
        N = sum(cont_table)
      )
    }, error = function(x) NULL)
  })

  # Remove NULL results and combine
  cal_list <- cal_list[!sapply(cal_list, is.null)]
  if (length(cal_list) < 1) return(NULL)

  cal_re <- do.call(rbind, cal_list)
  rownames(cal_re) <- cal_re$study

  # Perform meta-analysis using random effects model
  cal_meta_re <- tryCatch(
    suppressWarnings({
      metagen(TE = log_or,
              seTE = se,
              studlab = study,
              data = cal_re,
              sm = "OR", # Specify odds ratio as summary measure
              control = list(maxiter = 2000,
                             stepadj = 0.1,
                             threshold = 0.000001)
      )
    }),
    error = function(x) NULL
  )

  cal_meta_re
}