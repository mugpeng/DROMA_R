# Priority Scoring Functions ----

.addPriorityDirection <- function(df, direction_col = "direction") {
  df <- as.data.frame(df)

  if (!direction_col %in% colnames(df)) {
    df[[direction_col]] <- if ("effect_size" %in% colnames(df)) {
      ifelse(df$effect_size >= 0, "Up", "Down")
    } else {
      rep(NA_character_, nrow(df))
    }
  }

  df
}

.renamePrioritySourceStats <- function(df, prefix) {
  if (is.null(df) || nrow(df) == 0) {
    out <- data.frame(name = character(0), stringsAsFactors = FALSE, row.names = integer(0))
    out[[paste0("p_", prefix)]] <- numeric(0)
    out[[paste0("es_", prefix)]] <- numeric(0)
    out[[paste0("n_ds_", prefix)]] <- numeric(0)
    out[[paste0("q_", prefix)]] <- numeric(0)
    out[[paste0("dir_", prefix)]] <- character(0)
    return(out)
  }

  df <- .addPriorityDirection(df)
  out <- data.frame(name = df$name, stringsAsFactors = FALSE)
  out[[paste0("p_", prefix)]] <- if ("p_value" %in% colnames(df)) df$p_value else NA_real_
  out[[paste0("es_", prefix)]] <- if ("effect_size" %in% colnames(df)) df$effect_size else NA_real_
  out[[paste0("n_ds_", prefix)]] <- if ("n_datasets" %in% colnames(df)) df$n_datasets else NA_real_
  out[[paste0("q_", prefix)]] <- if ("q_value" %in% colnames(df)) df$q_value else NA_real_
  out[[paste0("dir_", prefix)]] <- if ("direction" %in% colnames(df)) df$direction else NA_character_
  out
}

.scorePriorityPComponent <- function(p_value, p_floor = 1e-300) {
  out <- pmin(-log10(pmax(p_value, p_floor)), 10) / 10
  out[is.na(out)] <- 0
  out
}

.scorePriorityEffectComponent <- function(effect_size) {
  out <- pmin(abs(effect_size), 0.5) / 0.5
  out[is.na(out)] <- 0
  out
}

.scorePriorityDatasetComponent <- function(n_datasets) {
  out <- pmin(log1p(n_datasets) / log1p(10), 1)
  out[is.na(out)] <- 0
  out
}

.computePriorityDirectionScore <- function(dir_cell, dir_pdopdx, dir_clin) {
  dirs <- c(cell = dir_cell, pdopdx = dir_pdopdx, clin = dir_clin)
  dirs <- dirs[!is.na(dirs) & dirs != ""]

  n_valid <- length(dirs)
  if (n_valid <= 1) {
    return(0)
  }
  if (n_valid == 2) {
    return(ifelse(length(unique(dirs)) == 1, 0.5, -0.5))
  }
  if (length(unique(dirs)) == 1) {
    return(1)
  }
  if (!is.na(dir_cell) && !is.na(dir_pdopdx) &&
      dir_cell == dir_pdopdx &&
      !is.na(dir_clin) &&
      dir_clin != dir_cell) {
    return(-2 / 3)
  }

  -1 / 3
}

#' Build a unified priority table from cell, PDO/PDX, and clinical evidence
#'
#' @description Combines evidence from cell-set, PDO/PDX, and clinical result
#'   tables into a single ranking table with normalized component scores and a
#'   final priority score.
#' @param candidate_names Character vector of candidate feature or pathway names.
#' @param cell_raw Raw result table for cell-set evidence.
#' @param pdopdx_raw Raw result table for PDO/PDX evidence.
#' @param clinical_raw Raw result table for clinical evidence.
#' @param cell_sig Significant candidate table derived from \code{cell_raw}.
#' @param pdopdx_sig Significant candidate table derived from \code{pdopdx_raw}.
#' @param clinical_sig Significant candidate table derived from \code{clinical_raw}.
#' @return A ranked data frame containing per-source statistics, component
#'   scores, and \code{PriorityScore}.
#' @export
buildPriorityTable <- function(candidate_names,
                               cell_raw,
                               pdopdx_raw,
                               clinical_raw,
                               cell_sig,
                               pdopdx_sig,
                               clinical_sig) {
  candidate_names <- unique(candidate_names)
  if (length(candidate_names) == 0) {
    return(data.frame(name = character(0), stringsAsFactors = FALSE))
  }

  out <- data.frame(name = candidate_names, stringsAsFactors = FALSE)
  out <- merge(
    out,
    .renamePrioritySourceStats(cell_raw[cell_raw$name %in% candidate_names, , drop = FALSE], "cell"),
    by = "name",
    all.x = TRUE
  )
  out <- merge(
    out,
    .renamePrioritySourceStats(pdopdx_raw[pdopdx_raw$name %in% candidate_names, , drop = FALSE], "pdopdx"),
    by = "name",
    all.x = TRUE
  )
  out <- merge(
    out,
    .renamePrioritySourceStats(clinical_raw[clinical_raw$name %in% candidate_names, , drop = FALSE], "clin"),
    by = "name",
    all.x = TRUE
  )

  out$selected_cellsets <- out$name %in% cell_sig$name
  out$selected_pdopdx <- out$name %in% pdopdx_sig$name
  out$selected_clinical <- out$name %in% clinical_sig$name
  out$n_selected_sources <- rowSums(cbind(out$selected_cellsets, out$selected_pdopdx, out$selected_clinical))

  out$P_cell <- .scorePriorityPComponent(out$p_cell)
  out$E_cell <- .scorePriorityEffectComponent(out$es_cell)
  out$D_cell <- .scorePriorityDatasetComponent(out$n_ds_cell)

  out$P_pdopdx <- .scorePriorityPComponent(out$p_pdopdx)
  out$E_pdopdx <- .scorePriorityEffectComponent(out$es_pdopdx)
  out$D_pdopdx <- .scorePriorityDatasetComponent(out$n_ds_pdopdx)

  out$P_clin <- .scorePriorityPComponent(out$p_clin)
  out$E_clin <- .scorePriorityEffectComponent(out$es_clin)
  out$D_clin <- .scorePriorityDatasetComponent(out$n_ds_clin)

  out$S_cellsets <- 0.50 * out$P_cell + 0.35 * out$E_cell + 0.15 * out$D_cell
  out$S_pdopdx <- 0.50 * out$P_pdopdx + 0.35 * out$E_pdopdx + 0.15 * out$D_pdopdx
  out$S_clinical <- 0.55 * out$P_clin + 0.35 * out$E_clin + 0.10 * out$D_clin
  out$S_direction_adj <- mapply(.computePriorityDirectionScore, out$dir_cell, out$dir_pdopdx, out$dir_clin)

  out$PriorityScore <- 0.25 * out$S_cellsets +
    0.30 * out$S_pdopdx +
    0.35 * out$S_clinical +
    0.10 * out$S_direction_adj

  out <- out[order(
    -out$PriorityScore,
    -out$S_clinical,
    -out$S_pdopdx,
    -out$S_cellsets,
    -abs(ifelse(is.na(out$es_clin), 0, out$es_clin)),
    out$name
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}
