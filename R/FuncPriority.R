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

.hasPrioritySourceRecord <- function(p_value, effect_size, n_datasets, q_value) {
  !is.na(p_value) | !is.na(effect_size) | !is.na(n_datasets) | !is.na(q_value)
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
#' @param require_sig_in_all_available Logical, whether to retain only
#'   candidates that are significant in every source where they have a record.
#' @return A ranked data frame containing per-source statistics, component
#'   scores, and \code{PriorityScore}.
#' @export
buildPriorityTable <- function(candidate_names,
                               cell_raw,
                               pdopdx_raw,
                               clinical_raw,
                               cell_sig,
                               pdopdx_sig,
                               clinical_sig,
                               require_sig_in_all_available = TRUE) {
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

  out$available_cellsets <- .hasPrioritySourceRecord(out$p_cell, out$es_cell, out$n_ds_cell, out$q_cell)
  out$available_pdopdx <- .hasPrioritySourceRecord(out$p_pdopdx, out$es_pdopdx, out$n_ds_pdopdx, out$q_pdopdx)
  out$available_clinical <- .hasPrioritySourceRecord(out$p_clin, out$es_clin, out$n_ds_clin, out$q_clin)

  out$selected_cellsets <- ifelse(out$available_cellsets, out$name %in% cell_sig$name, NA)
  out$selected_pdopdx <- ifelse(out$available_pdopdx, out$name %in% pdopdx_sig$name, NA)
  out$selected_clinical <- ifelse(out$available_clinical, out$name %in% clinical_sig$name, NA)

  out$n_available_sources <- rowSums(cbind(out$available_cellsets, out$available_pdopdx, out$available_clinical), na.rm = TRUE)
  out$n_selected_sources <- rowSums(cbind(
    ifelse(is.na(out$selected_cellsets), FALSE, out$selected_cellsets),
    ifelse(is.na(out$selected_pdopdx), FALSE, out$selected_pdopdx),
    ifelse(is.na(out$selected_clinical), FALSE, out$selected_clinical)
  ))

  if (isTRUE(require_sig_in_all_available)) {
    keep_idx <- out$n_available_sources > 0 & out$n_selected_sources == out$n_available_sources
    out <- out[keep_idx, , drop = FALSE]
    if (nrow(out) == 0) {
      return(out)
    }
  }

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

#' Prepare merged priority table for visualization
#'
#' @description Merges gene and pathway priority tables, annotates candidate type,
#'   strips common pathway prefixes, and derives summary columns used by
#'   visualization helpers.
#' @param gene_priority Data frame produced by \code{buildPriorityTable()} for genes.
#' @param pathway_priority Data frame produced by \code{buildPriorityTable()} for pathways.
#' @param top_n Optional integer. If provided, keep only the top N candidates by
#'   \code{PriorityScore} after merging.
#' @param pathway_prefix Regular expression prefix to strip from pathway names for display.
#' @return A merged and annotated data frame ready for plotting.
#' @export
preparePriorityVisualizationData <- function(gene_priority,
                                             pathway_priority,
                                             top_n = NULL,
                                             pathway_prefix = "^(HALLMARK_|GOBP_|KEGG_)") {
  gene_priority <- as.data.frame(gene_priority)
  pathway_priority <- as.data.frame(pathway_priority)

  if (!"name" %in% colnames(gene_priority) || !"PriorityScore" %in% colnames(gene_priority)) {
    stop("gene_priority must contain at least 'name' and 'PriorityScore' columns")
  }
  if (!"name" %in% colnames(pathway_priority) || !"PriorityScore" %in% colnames(pathway_priority)) {
    stop("pathway_priority must contain at least 'name' and 'PriorityScore' columns")
  }

  gene_priority$candidate_class <- "Gene"
  pathway_priority$candidate_class <- "Pathway"

  merged_priority <- rbind(pathway_priority, gene_priority)
  merged_priority <- merged_priority[order(-merged_priority$PriorityScore, merged_priority$name), , drop = FALSE]
  rownames(merged_priority) <- NULL

  if (!is.null(top_n)) {
    top_n <- min(as.integer(top_n), nrow(merged_priority))
    merged_priority <- utils::head(merged_priority, top_n)
  }

  merged_priority$display_name <- ifelse(
    merged_priority$candidate_class == "Pathway",
    gsub(pathway_prefix, "", merged_priority$name),
    merged_priority$name
  )
  merged_priority$rank <- seq_len(nrow(merged_priority))
  merged_priority$n_selected_sources <- if ("n_selected_sources" %in% colnames(merged_priority)) {
    merged_priority$n_selected_sources
  } else {
    rowSums(cbind(
      if ("selected_cellsets" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_cellsets), FALSE, merged_priority$selected_cellsets) else FALSE,
      if ("selected_pdopdx" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_pdopdx), FALSE, merged_priority$selected_pdopdx) else FALSE,
      if ("selected_clinical" %in% colnames(merged_priority)) ifelse(is.na(merged_priority$selected_clinical), FALSE, merged_priority$selected_clinical) else FALSE
    ))
  }
  merged_priority$S_preclinical <- rowMeans(cbind(
    if ("S_cellsets" %in% colnames(merged_priority)) merged_priority$S_cellsets else NA_real_,
    if ("S_pdopdx" %in% colnames(merged_priority)) merged_priority$S_pdopdx else NA_real_
  ), na.rm = TRUE)
  merged_priority$S_preclinical[!is.finite(merged_priority$S_preclinical)] <- 0

  merged_priority$dir_cell_selected <- if ("selected_cellsets" %in% colnames(merged_priority) &&
                                           "dir_cell" %in% colnames(merged_priority)) {
    ifelse(!is.na(merged_priority$selected_cellsets) & merged_priority$selected_cellsets, merged_priority$dir_cell, NA_character_)
  } else {
    NA_character_
  }
  merged_priority$dir_pdopdx_selected <- if ("selected_pdopdx" %in% colnames(merged_priority) &&
                                             "dir_pdopdx" %in% colnames(merged_priority)) {
    ifelse(!is.na(merged_priority$selected_pdopdx) & merged_priority$selected_pdopdx, merged_priority$dir_pdopdx, NA_character_)
  } else {
    NA_character_
  }
  merged_priority$dir_clin_selected <- if ("selected_clinical" %in% colnames(merged_priority) &&
                                           "dir_clin" %in% colnames(merged_priority)) {
    ifelse(!is.na(merged_priority$selected_clinical) & merged_priority$selected_clinical, merged_priority$dir_clin, NA_character_)
  } else {
    NA_character_
  }

  merged_priority$direction_pattern <- apply(
    cbind(
      merged_priority$dir_cell_selected,
      merged_priority$dir_pdopdx_selected,
      merged_priority$dir_clin_selected
    ),
    1,
    function(x) paste(ifelse(is.na(x) | x == "", "NA", x), collapse = " / ")
  )

  merged_priority$direction_support_n <- rowSums(cbind(
    !is.na(merged_priority$dir_cell_selected),
    !is.na(merged_priority$dir_pdopdx_selected),
    !is.na(merged_priority$dir_clin_selected)
  ))

  merged_priority$direction_pattern_short <- apply(
    cbind(
      merged_priority$dir_cell_selected,
      merged_priority$dir_pdopdx_selected,
      merged_priority$dir_clin_selected
    ),
    1,
    function(x) {
      source_codes <- c("C", "P", "L")
      valid <- !is.na(x) & x != ""
      if (!any(valid)) {
        return("NA")
      }
      dir_codes <- ifelse(x[valid] == "Up", "+", ifelse(x[valid] == "Down", "-", "?"))
      paste0(source_codes[valid], dir_codes, collapse = " ")
    }
  )

  merged_priority$direction_consensus <- apply(
    cbind(
      merged_priority$dir_cell_selected,
      merged_priority$dir_pdopdx_selected,
      merged_priority$dir_clin_selected
    ),
    1,
    function(x) {
      dirs <- x[!is.na(x) & x != ""]
      if (length(dirs) == 0) {
        return("NA")
      }
      if (length(unique(dirs)) == 1) {
        return(dirs[1])
      }
      tab <- sort(table(dirs), decreasing = TRUE)
      if (length(tab) > 1 && tab[1] == tab[2]) {
        return("Mixed")
      }
      names(tab)[1]
    }
  )

  merged_priority$display_name <- factor(
    merged_priority$display_name,
    levels = rev(merged_priority$display_name)
  )

  merged_priority
}
