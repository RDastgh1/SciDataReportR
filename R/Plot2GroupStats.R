#' Plot & Summarize Group Stats via MakeComparisonTable
#'
#' @param Data data.frame
#' @param Variables character vector of variables to analyze
#' @param VariableCategories optional df with columns: Variable, Category
#' @param impClust,normalClust group labels (impClust is the "rightward" direction)
#' @param GroupVar grouping column name in Data
#' @param missing_threshold drop vars with > this fraction missing (default 0.80)
#' @param max_levels drop factors with > max_levels (default 10)
#' @param adjust_method p-adjust method for FDR (default "fdr")
#' @param label_q q-value threshold for labeling (default 0.05)
#' @param mct_args list of extra args forwarded to MakeComparisonTable()
#' @return list(plot, table, pvaltable, data_used)
#' @export
Plot2GroupStats <- function(
    Data, Variables, VariableCategories = NULL,
    impClust, normalClust, GroupVar,
    missing_threshold = 0.80, max_levels = 10,
    adjust_method = "fdr", label_q = 0.05,
    mct_args = list()
){
  stopifnot(is.data.frame(Data), is.character(Variables), length(Variables) >= 1)
  if (!GroupVar %in% names(Data)) stop("`GroupVar` not found in `Data`.")
  miss <- setdiff(Variables, names(Data))
  if (length(miss)) stop("Variables not in `Data`: ", paste(miss, collapse=", "))

  # prep & filter to the two groups of interest
  tData <- dplyr::mutate(Data, GroupVar = .data[[GroupVar]]) |>
    dplyr::select(GroupVar, dplyr::all_of(Variables)) |>
    dplyr::filter(.data$GroupVar %in% c(normalClust, impClust))
  tData$GroupVar <- factor(tData$GroupVar, levels = c(normalClust, impClust))

  # harmonize types
  tData <- dplyr::mutate(tData, dplyr::across(where(is.character) & !dplyr::all_of("GroupVar"), ~factor(.x)))

  # drop by missingness / factor levels / numeric nzv
  miss_prop <- colMeans(is.na(tData))
  drop_missing <- setdiff(names(miss_prop[miss_prop > missing_threshold]), "GroupVar")
  if (length(drop_missing)) tData <- dplyr::select(tData, -dplyr::all_of(drop_missing))

  fac_levels <- vapply(tData, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, 1L)
  drop_low  <- names(fac_levels[!is.na(fac_levels) & fac_levels < 2])
  drop_high <- names(fac_levels[!is.na(fac_levels) & fac_levels > max_levels])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(c(drop_low, drop_high), "GroupVar")))

  num_sd <- vapply(tData, function(x) if (is.numeric(x)) stats::sd(x, na.rm=TRUE) else NA_real_, 1.0)
  drop_nzv <- names(num_sd[is.na(num_sd) | num_sd == 0])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(drop_nzv, "GroupVar")))

  # ---- Call your wrapper ----------------------------------------------------
  # Default args for MakeComparisonTable; allow user to override via mct_args
  mct_defaults <- list(
    DataFrame    = tData,
    Variables    = setdiff(names(tData), "GroupVar"),
    CompVariable = "GroupVar",
    ValueDigits  = 2,
    pDigits      = 3
  )
  mct_call <- utils::modifyList(mct_defaults, mct_args, keep.null = TRUE)
  MCT <- do.call(SciDataReportR::MakeComparisonTable, mct_call)

  # Ensure p and q exist. If wrapper already computed p, we leave it;
  # if not, add Wilcoxon/chi-sq with exact=FALSE to avoid ties warnings.
  if (inherits(MCT, "gtsummary")) {
    tb <- MCT$table_body
    if (!"p.value" %in% names(tb)) {
      MCT <- MCT |>
        gtsummary::add_p(
          test = list(
            gtsummary::all_continuous()  ~ "wilcox.test",
            gtsummary::all_categorical() ~ "chisq.test"
          ),
          test.args = list(
            gtsummary::all_continuous() ~ list(exact = FALSE)
          )
        )
    }
    tb <- MCT$table_body
    if (!"q.value" %in% names(tb)) {
      MCT <- gtsummary::add_q(MCT, method = adjust_method)
    }
    tb <- MCT$table_body
  } else {
    stop("MakeComparisonTable did not return a gtsummary table.")
  }

  # ---- Extract effect size (robustly) --------------------------------------
  # Try common column names; otherwise compute direction from the data.
  tb <- MCT$table_body
  effect_cols <- intersect(
    c("std.diff","SMD","effect_size","estimate","diff","Difference","cohens_d","hedges_g","logOR","OR","RR"),
    names(tb)
  )
  effect <- suppressWarnings(sapply(tb[[effect_cols[1]]], as.numeric))
  if (length(effect) == 0 || all(is.na(effect))) effect <- NULL

  # build p-value table (one row per variable label row)
  pvaltable <- tb |>
    dplyr::filter(.data$row_type == "label") |>
    dplyr::select(variable, label, dplyr::any_of(c("p.value","q.value")))

  # Effect sign: use effect if numeric; otherwise compute from raw data
  eff_sign <- function(vname) {
    if (!is.null(effect)) {
      return(sign(effect[match(vname, tb$variable)]))
    } else {
      v <- tData[[vname]]
      if (is.factor(v)) v <- suppressWarnings(as.numeric(as.character(v)))
      if (!is.numeric(v)) return(NA_real_)
      m_imp  <- mean(v[tData$GroupVar == impClust],  na.rm=TRUE)
      m_norm <- mean(v[tData$GroupVar == normalClust], na.rm=TRUE)
      return(sign(m_imp - m_norm))
    }
  }
  pvaltable$effect_sign <- vapply(pvaltable$variable, eff_sign, 1.0)

  # magnitude: |effect| if available (for point size)
  if (!is.null(effect)) {
    pvaltable$effect_mag <- abs(effect[match(pvaltable$variable, tb$variable)])
  } else {
    pvaltable$effect_mag <- NA_real_
  }

  # signed -log10(p)
  tiny <- .Machine$double.xmin
  pvaltable$signed_log10p <- -log10(pmax(pvaltable$p.value, tiny)) * pvaltable$effect_sign

  # optional categories
  if (!is.null(VariableCategories) && all(c("Variable","Category") %in% names(VariableCategories))) {
    pvaltable$Category <- VariableCategories$Category[
      match(pvaltable$variable, VariableCategories$Variable)
    ]
  } else pvaltable$Category <- NA_character_

  # ---- Plot -----------------------------------------------------------------
  label_df <- dplyr::filter(pvaltable, !is.na(.data$q.value) & .data$q.value < label_q)
  p <- ggplot2::ggplot(
    pvaltable,
    ggplot2::aes(x = .data$signed_log10p,
                 y = .data$variable,
                 size = .data$effect_mag,
                 color = if (any(!is.na(pvaltable$Category))) .data$Category else NULL)
  ) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = rev(pvaltable$variable)) +
    ggplot2::scale_size_continuous(name = "Effect size", range = c(2,6), guide = "legend") +
    ggplot2::xlab(paste0("signed -log10(p)   (→ higher in ", impClust, ")")) +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(paste("---> Higher in", impClust)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

  if (rlang::is_installed("paletteer") && any(!is.na(pvaltable$Category)))
    p <- p + paletteer::scale_color_paletteer_d("fishualize::Scarus_quoyi")

  if (rlang::is_installed("ggrepel") && nrow(label_df)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df, ggplot2::aes(label = .data$variable),
      max.overlaps = 30, min.segment.length = 0
    )
  }

  list(
    plot       = p,
    table      = MCT,            # your MakeComparisonTable output (gtsummary)
    pvaltable  = pvaltable,
    data_used  = dplyr::as_tibble(tData)
  )
}
