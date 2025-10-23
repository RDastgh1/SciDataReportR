#' Plot & Summarize Group Stats via MakeComparisonTable
#'
#' Creates a volcano-style comparison plot and returns the full comparison table.
#' Uses SciDataReportR::MakeComparisonTable() for p-values (and effect sizes when available).
#'
#' @param Data data.frame with your data
#' @param Variables character vector of variable names to analyze
#' @param VariableCategories optional data.frame with columns: Variable, Category
#' @param impClust character/number; the "important" group label (plotted to the right)
#' @param normalClust character/number; the comparison group label
#' @param GroupVar character; column name in `Data` with the grouping variable
#' @param missing_threshold numeric in [0,1]; drop vars with > this fraction missing (default 0.80)
#' @param max_levels integer; drop factors with > this many levels (default 10)
#' @param adjust_method character; p.adjust method for FDR (default "fdr")
#' @param label_q numeric; q-value threshold for labeling points on the plot (default 0.05)
#' @param mct_args list; extra args forwarded to SciDataReportR::MakeComparisonTable()
#'
#' @return list(plot=ggplot, table=gtsummary object, pvaltable=data.frame, data_used=tibble)
#' @export
#' @import dplyr ggplot2 gtsummary
Plot2GroupStats <- function(
    Data, Variables, VariableCategories = NULL,
    impClust, normalClust, GroupVar,
    missing_threshold = 0.80, max_levels = 10,
    adjust_method = "fdr", label_q = 0.05,
    mct_args = list()
){
  # ---------- checks ----------
  stopifnot(is.data.frame(Data), is.character(Variables), length(Variables) >= 1)
  if (!GroupVar %in% names(Data)) stop("`GroupVar` not found in `Data`.")
  miss <- setdiff(Variables, names(Data))
  if (length(miss)) stop("Variables not in `Data`: ", paste(miss, collapse = ", "))

  # ---------- subset & harmonize ----------
  tData <- dplyr::mutate(Data, GroupVar = .data[[GroupVar]]) |>
    dplyr::select(GroupVar, dplyr::all_of(Variables)) |>
    dplyr::filter(.data$GroupVar %in% c(normalClust, impClust))
  # ensure both groups exist
  if (!all(c(normalClust, impClust) %in% unique(tData$GroupVar))) {
    stop("After filtering, one of the groups is missing. Check `impClust` / `normalClust` and data.")
  }
  tData$GroupVar <- factor(tData$GroupVar, levels = c(normalClust, impClust))
  # convert characters (not the GroupVar) to factors
  tData <- tData |>
    dplyr::mutate(dplyr::across(where(is.character), factor))

  # drop by missingness
  miss_prop <- colMeans(is.na(tData))
  drop_missing <- setdiff(names(miss_prop[miss_prop > missing_threshold]), "GroupVar")
  if (length(drop_missing)) {
    tData <- dplyr::select(tData, -dplyr::all_of(drop_missing))
  }
  # drop factor vars with <2 or >max_levels levels
  fac_levels <- vapply(tData, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, 1L)
  drop_low  <- names(fac_levels[!is.na(fac_levels) & fac_levels < 2])
  drop_high <- names(fac_levels[!is.na(fac_levels) & fac_levels > max_levels])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(c(drop_low, drop_high), "GroupVar")))
  # drop numeric zero-variance
  num_sd <- vapply(tData, function(x) if (is.numeric(x)) stats::sd(x, na.rm = TRUE) else NA_real_, 1.0)
  drop_nzv <- names(num_sd[is.na(num_sd) | num_sd == 0])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(drop_nzv, "GroupVar")))

  # if nothing left to analyze, stop early
  vars_in <- setdiff(names(tData), "GroupVar")
  if (length(vars_in) == 0) stop("No analyzable variables remain after filtering/missingness checks.")

  # ---------- call your wrapper (robustly) ----------
  # Build defaults; user can override via mct_args (e.g., AddEffectSize=TRUE)
  mct_defaults <- list(
    DataFrame    = tData,
    Variables    = vars_in,
    CompVariable = "GroupVar",
    ValueDigits  = 2,
    pDigits      = 3
  )
  mct_call <- utils::modifyList(mct_defaults, mct_args, keep.null = TRUE)

  # Try MakeComparisonTable with provided args; if an unknown arg sneaks in, retry without it.
  MCT <- tryCatch(
    do.call(SciDataReportR::MakeComparisonTable, mct_call),
    error = function(e) {
      # strip args not in formals if possible
      fn <- tryCatch(get("MakeComparisonTable", asNamespace("SciDataReportR")), error = function(e2) NULL)
      if (!is.null(fn)) {
        valid <- names(formals(fn))
        cleaned <- mct_call[intersect(names(mct_call), valid)]
        do.call(SciDataReportR::MakeComparisonTable, cleaned)
      } else stop(e)
    }
  )

  if (!inherits(MCT, "gtsummary")) {
    stop("MakeComparisonTable did not return a gtsummary table.")
  }

  # ensure p and q exist; add quietly if missing (Wilcoxon exact=FALSE to avoid tie warnings)
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

  # ---------- extract a clean p/q/effect table ----------
  # label rows are one-per-variable
  pvaltable <- tb |>
    dplyr::filter(.data$row_type == "label") |>
    dplyr::select(variable, label, dplyr::any_of(c("p.value","q.value")))

  # Try common effect-size column names from your wrapper; may be absent.
  effect_cols <- intersect(
    c("std.diff","SMD","effect_size","estimate","diff","Difference",
      "cohens_d","hedges_g","logOR","OR","RR"),
    names(tb)
  )
  effect_vec <- NULL
  if (length(effect_cols) > 0) {
    tmp <- suppressWarnings(as.numeric(tb[[effect_cols[1]]]))
    if (!all(is.na(tmp))) {
      effect_vec <- tmp
      names(effect_vec) <- tb$variable
    }
  }

  # Direction function: use effect sign if available; else compute from data
  dir_from_data <- function(vname) {
    v <- tData[[vname]]
    g <- tData$GroupVar
    if (is.numeric(v)) {
      return(sign(mean(v[g == impClust], na.rm = TRUE) -
                    mean(v[g == normalClust], na.rm = TRUE)))
    }
    if (is.factor(v)) {
      lv <- levels(v)
      # pick the level with the largest absolute proportion difference
      pdiff <- sapply(lv, function(l) {
        mean(v[g == impClust] == l, na.rm = TRUE) -
          mean(v[g == normalClust] == l, na.rm = TRUE)
      })
      return(sign(pdiff[which.max(abs(pdiff))]))
    }
    return(NA_real_)
  }

  if (!is.null(effect_vec)) {
    es <- effect_vec[match(pvaltable$variable, names(effect_vec))]
    pvaltable$effect_sign <- sign(es)
    pvaltable$effect_mag  <- abs(es)
  } else {
    pvaltable$effect_sign <- vapply(pvaltable$variable, dir_from_data, 1.0)
    pvaltable$effect_mag  <- NA_real_
  }

  # signed -log10(p)
  tiny <- .Machine$double.xmin
  pvaltable$signed_log10p <- -log10(pmax(pvaltable$p.value, tiny)) * pvaltable$effect_sign

  # significance bins for coloring when no categories provided
  if ("q.value" %in% names(pvaltable)) {
    brks <- c(-Inf, 0.001, 0.01, 0.05, Inf)
    labs <- c("q<0.001","q<0.01","q<0.05","ns")
    pvaltable$Sig <- cut(pvaltable$q.value, breaks = brks, labels = labs, right = TRUE)
  } else {
    brks <- c(-Inf, 0.001, 0.01, 0.05, Inf)
    labs <- c("p<0.001","p<0.01","p<0.05","ns")
    pvaltable$Sig <- cut(pvaltable$p.value, breaks = brks, labels = labs, right = TRUE)
  }

  # optional category mapping
  if (!is.null(VariableCategories) &&
      all(c("Variable","Category") %in% names(VariableCategories))) {
    pvaltable$Category <- VariableCategories$Category[
      match(pvaltable$variable, VariableCategories$Variable)
    ]
  } else {
    pvaltable$Category <- NA_character_
  }

  # ---------- plot ----------
  have_cat <- any(!is.na(pvaltable$Category))
  label_df <- dplyr::filter(pvaltable, !is.na(q.value) & q.value < label_q)

  aes_color <- if (have_cat) ggplot2::aes(color = Category) else ggplot2::aes(color = Sig)
  p <- ggplot2::ggplot(
    pvaltable,
    ggplot2::aes(x = signed_log10p, y = variable, size = effect_mag)
  ) +
    aes_color +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = rev(pvaltable$variable)) +
    ggplot2::scale_size_continuous(name = "Effect size", range = c(2, 6), guide = "legend") +
    ggplot2::xlab(paste0("signed -log10(p)   (→ higher in ", impClust, ")")) +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(paste("---> Higher in", impClust)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

  if (requireNamespace("paletteer", quietly = TRUE) && have_cat) {
    p <- p + paletteer::scale_color_paletteer_d("fishualize::Scarus_quoyi")
  }

  if (nrow(label_df) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df, ggplot2::aes(label = variable),
        max.overlaps = 30, min.segment.length = 0
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_df, ggplot2::aes(label = variable),
        hjust = -0.05, size = 3
      )
    }
  }

  # ---------- return ----------
  list(
    plot      = p,
    table     = MCT,                   # your MakeComparisonTable output (gtsummary)
    pvaltable = pvaltable,
    data_used = dplyr::as_tibble(tData)
  )
}
