#' Plot & Summarize Group Stats via MakeComparisonTable
#' (BH q from p; SHAPE by p; COLOR by Category (vector or data frame); stable point size; palette via paletteer)
#'
#' @param Data data.frame
#' @param Variables character vector of variables to analyze
#' @param VariableCategories optional:
#'        - data frame with columns Variable, Category; OR
#'        - vector of categories (named by variable OR unnamed aligned to `Variables`)
#' @param impClust,normalClust labels for the two groups (impClust plotted to the RIGHT for signed axes)
#' @param GroupVar column name in `Data` holding the group labels
#' @param missing_threshold drop vars with > this fraction missing (default 0.80)
#' @param max_levels drop factors with > this many levels (default 10)
#' @param label_q label threshold using q (default 0.05)
#' @param x_axis one of c("signed_logp","signed_effect","effect","logp")
#' @param sort_by one of c("q","p","effect","signed_logp","signed_effect","none")
#' @param mct_args list of extra args to SciDataReportR::MakeComparisonTable(); e.g., AddEffectSize=TRUE
#' @param palette paletteer palette string for category colors (default "pals::alphabet")
#' @param point_size numeric constant for point size (default 3.5)
#' @param verbose logical; TRUE prints stage messages for debugging (default FALSE)
#'
#' @return list(plot=ggplot, table=gtsummary, pvaltable=data.frame, data_used=tibble)
#' @export
#' @import dplyr ggplot2 gtsummary
Plot2GroupStats <- function(
    Data, Variables, VariableCategories = NULL,
    impClust, normalClust, GroupVar,
    missing_threshold = 0.80, max_levels = 10,
    label_q = 0.05,
    x_axis = c("signed_logp","signed_effect","effect","logp"),
    sort_by = c("q","p","effect","signed_logp","signed_effect","none"),
    mct_args = list(),
    palette = "pals::alphabet",
    point_size = 3.5,
    verbose = FALSE
){
  x_axis  <- match.arg(x_axis)
  sort_by <- match.arg(sort_by)
  say <- function(...) if (isTRUE(verbose)) message(...)

  # ---- helpers --------------------------------------------------------------
  safe_num <- function(x) { if (is.character(x)) suppressWarnings(as.numeric(x)) else x }
  parse_pvals <- function(x) {
    if (is.numeric(x)) return(x)
    y <- gsub(",", "", x); y <- gsub("^\\s*[<≤]\\s*", "", y)
    suppressWarnings(as.numeric(y))
  }
  # Return a named character vector Category[variable]
  resolve_categories <- function(VariableCategories, variables_arg, vars_present, verbose = FALSE){
    if (is.null(VariableCategories))
      return(stats::setNames(rep(NA_character_, length(vars_present)), vars_present))

    if (is.data.frame(VariableCategories) &&
        all(c("Variable","Category") %in% names(VariableCategories))) {
      cat_map <- stats::setNames(as.character(VariableCategories$Category),
                                 as.character(VariableCategories$Variable))
      out <- unname(cat_map[vars_present])
      return(stats::setNames(ifelse(is.na(out), NA_character_, out), vars_present))
    }

    if (is.atomic(VariableCategories)) {
      v <- VariableCategories
      # Named vector -> match by names
      if (!is.null(names(v)) && any(nzchar(names(v)))) {
        cat_map <- stats::setNames(as.character(v), names(v))
        out <- unname(cat_map[vars_present])
        return(stats::setNames(ifelse(is.na(out), NA_character_, out), vars_present))
      }
      # Unnamed vector: align to Variables first, else to present vars
      if (length(v) == length(variables_arg)) {
        cat_map <- stats::setNames(as.character(v), variables_arg)
        out <- unname(cat_map[vars_present])
        return(stats::setNames(ifelse(is.na(out), NA_character_, out), vars_present))
      } else if (length(v) == length(vars_present)) {
        return(stats::setNames(as.character(v), vars_present))
      } else {
        if (verbose) message("VariableCategories length mismatch; categories set to NA.")
        return(stats::setNames(rep(NA_character_, length(vars_present)), vars_present))
      }
    }

    if (verbose) message("Unrecognized VariableCategories format; categories set to NA.")
    stats::setNames(rep(NA_character_, length(vars_present)), vars_present)
  }

  # ---- checks ---------------------------------------------------------------
  stopifnot(is.data.frame(Data), is.character(Variables), length(Variables) >= 1)
  if (!GroupVar %in% names(Data)) stop("`GroupVar` not found in `Data`.")
  miss <- setdiff(Variables, names(Data))
  if (length(miss)) stop("Variables not in `Data`: ", paste(miss, collapse = ", "))

  # ---- working data ---------------------------------------------------------
  say("stage: prepare data")
  tData <- Data
  tData$GroupVar <- Data[[GroupVar]]
  tData <- tData |>
    dplyr::select(GroupVar, dplyr::all_of(Variables)) |>
    dplyr::filter(GroupVar %in% c(normalClust, impClust))
  if (!all(c(normalClust, impClust) %in% unique(tData$GroupVar))) {
    stop("After filtering, one of the groups is missing. Check `impClust` / `normalClust` and data.")
  }
  tData$GroupVar <- factor(tData$GroupVar, levels = c(normalClust, impClust))
  tData <- tData |> dplyr::mutate(dplyr::across(where(is.character), factor))

  # ---- pruning --------------------------------------------------------------
  say("stage: prune variables")
  miss_prop <- colMeans(is.na(tData))
  drop_missing <- setdiff(names(miss_prop[miss_prop > missing_threshold]), "GroupVar")
  if (length(drop_missing)) tData <- dplyr::select(tData, -dplyr::all_of(drop_missing))

  fac_levels <- vapply(tData, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, 1L)
  drop_low  <- names(fac_levels[!is.na(fac_levels) & fac_levels < 2])
  drop_high <- names(fac_levels[!is.na(fac_levels) & fac_levels > max_levels])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(c(drop_low, drop_high), "GroupVar")))

  num_sd <- vapply(tData, function(x) if (is.numeric(x)) stats::sd(x, na.rm = TRUE) else NA_real_, 1.0)
  drop_nzv <- names(num_sd[is.na(num_sd) | num_sd == 0])
  tData <- dplyr::select(tData, -dplyr::all_of(setdiff(drop_nzv, "GroupVar")))

  vars_in <- setdiff(names(tData), "GroupVar")
  if (length(vars_in) == 0) stop("No analyzable variables remain after filtering/missingness checks.")

  # ---- MakeComparisonTable --------------------------------------------------
  say("stage: MakeComparisonTable")
  mct_defaults <- list(
    DataFrame    = tData,
    Variables    = vars_in,
    CompVariable = "GroupVar",
    ValueDigits  = 2,
    pDigits      = 3
  )
  # Ask for absolute effect sizes unless user overrode; direction comes from data
  if ("AddEffectSize" %in% names(formals(SciDataReportR::MakeComparisonTable)) &&
      !"AddEffectSize" %in% names(mct_args)) {
    mct_defaults$AddEffectSize <- TRUE
  }
  mct_call <- utils::modifyList(mct_defaults, mct_args, keep.null = TRUE)
  MCT <- tryCatch(
    do.call(SciDataReportR::MakeComparisonTable, mct_call),
    error = function(e) {
      fn <- try(get("MakeComparisonTable", asNamespace("SciDataReportR")), silent = TRUE)
      if (!inherits(fn, "try-error")) {
        valid <- names(formals(fn))
        cleaned <- mct_call[intersect(names(mct_call), valid)]
        do.call(SciDataReportR::MakeComparisonTable, cleaned)
      } else stop(e)
    }
  )
  if (!inherits(MCT, "gtsummary")) stop("MakeComparisonTable did not return a gtsummary table.")

  # ---- extract p, compute BH q ---------------------------------------------
  say("stage: build p/q table")
  tb <- MCT$table_body
  if (!all(c("row_type","variable") %in% names(tb))) {
    stop("Unexpected structure from MakeComparisonTable$table_body (needs `row_type` and `variable`).")
  }
  pvaltable <- tb |>
    dplyr::filter(.data$row_type == "label") |>
    dplyr::select(dplyr::any_of(c("variable","label","p.value","effect_size","es_method")))
  if (!"p.value" %in% names(pvaltable)) {
    stop("`p.value` is missing in MakeComparisonTable output.")
  }
  pnum <- pvaltable$p.value
  if (!is.numeric(pnum)) pnum <- parse_pvals(pnum)
  tiny <- .Machine$double.xmin
  pnum <- pmax(pnum, tiny)
  pvaltable$p.value <- pnum
  pvaltable$q.value <- p.adjust(pnum, method = "BH")

  # ---- absolute effect magnitude & direction -------------------------------
  pvaltable$effect_abs <- safe_num(pvaltable$effect_size)

  say("stage: direction from data")
  dir_from_data <- function(vname) {
    if (!vname %in% names(tData)) return(NA_real_)
    v <- tData[[vname]]; g <- tData$GroupVar
    if (is.numeric(v)) {
      return(sign(mean(v[g == impClust], na.rm = TRUE) -
                    mean(v[g == normalClust], na.rm = TRUE)))
    }
    v <- as.factor(v)
    lv <- levels(v); if (length(lv) < 1) return(NA_real_)
    pdiff <- sapply(lv, function(l) {
      mean(v[g == impClust] == l, na.rm = TRUE) -
        mean(v[g == normalClust] == l, na.rm = TRUE)
    })
    sign(pdiff[which.max(abs(pdiff))])
  }
  pvaltable$effect_sign <- vapply(pvaltable$variable, dir_from_data, 1.0)

  # ---- metrics & y-order ----------------------------------------------------
  say("stage: metrics & ordering")
  pvaltable$logp        <- -log10(pvaltable$p.value)
  pvaltable$signed_logp <- pvaltable$logp * pvaltable$effect_sign
  pvaltable$signed_es   <- pvaltable$effect_abs * pvaltable$effect_sign

  # Shapes by p-value bins; colors by category (or q as fallback)
  brks_p <- c(-Inf, 0.001, 0.01, 0.05, Inf); labs_p <- c("p<0.001","p<0.01","p<0.05","ns")
  pvaltable$SigP <- cut(pvaltable$p.value, breaks = brks_p, labels = labs_p, right = TRUE)

  brks_q <- c(-Inf, 0.001, 0.01, 0.05, Inf); labs_q <- c("q<0.001","q<0.01","q<0.05","ns")
  pvaltable$SigQ <- cut(pvaltable$q.value, breaks = brks_q, labels = labs_q, right = TRUE)

  # Categories (accept vector or data frame)
  vars_present <- as.character(pvaltable$variable)
  cat_vec <- resolve_categories(VariableCategories, variables_arg = Variables,
                                vars_present = vars_present, verbose = verbose)
  pvaltable$Category <- unname(cat_vec[vars_present])

  # Order y-axis
  order_key <- switch(
    sort_by,
    q             = pvaltable$q.value,
    p             = pvaltable$p.value,
    effect        = -pvaltable$effect_abs,
    signed_logp   = -pvaltable$signed_logp,
    signed_effect = -pvaltable$signed_es,
    none          = seq_len(nrow(pvaltable))
  )
  ord <- order(order_key, na.last = TRUE)
  pvaltable <- pvaltable[ord, , drop = FALSE]
  pvaltable$variable <- factor(pvaltable$variable, levels = rev(pvaltable$variable))

  # Choose x & precompute
  x_map <- switch(
    x_axis,
    signed_logp   = list(var = "signed_logp",  lab = paste0("signed -log10(p)  (→ higher in ", impClust, ")"), signed = TRUE),
    signed_effect = list(var = "signed_es",    lab = paste0("signed effect size  (→ higher in ", impClust, ")"), signed = TRUE),
    effect        = list(var = "effect_abs",   lab = "absolute effect size",                                      signed = FALSE),
    logp          = list(var = "logp",         lab = "-log10(p)",                                                 signed = FALSE)
  )
  pvaltable$xvar <- pvaltable[[x_map$var]]

  # Labels by q
  label_df <- pvaltable[!is.na(pvaltable$q.value) & pvaltable$q.value < label_q, , drop = FALSE]

  # ---- plot (stable point size) --------------------------------------------
  say("stage: plot")
  use_categories <- any(!is.na(pvaltable$Category))
  aes_color <- if (use_categories) ggplot2::aes(color = Category) else ggplot2::aes(color = SigQ)

  p <- ggplot2::ggplot(
    pvaltable,
    ggplot2::aes(x = xvar, y = variable, shape = SigP)
  ) +
    aes_color +
    { if (isTRUE(x_map$signed)) ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.08) else NULL } +
    { if (isTRUE(x_map$signed)) ggplot2::geom_vline(xintercept = 0) else NULL } +
    ggplot2::geom_point(size = point_size) +  # <-- constant size
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_shape_manual(
      name   = "Significance (p)",
      values = c("p<0.001" = 17, "p<0.01" = 15, "p<0.05" = 16, "ns" = 1)
    ) +
    ggplot2::xlab(x_map$lab) +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(paste("---> Higher in", impClust)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

  # Dynamic category palette via paletteer (default pals::alphabet); robust fallback
  if (use_categories) {
    ucat <- unique(stats::na.omit(pvaltable$Category))
    n_cat <- length(ucat)
    cols <- NULL

    if (requireNamespace("paletteer", quietly = TRUE)) {
      cols <- tryCatch(
        as.character(paletteer::paletteer_d(palette, n = n_cat)),
        error = function(e) NULL
      )
    }
    if (is.null(cols) || length(cols) < n_cat) {
      cols <- grDevices::hcl.colors(n_cat, palette = "Okabe-Ito")
    }
    names(cols) <- ucat
    p <- p + ggplot2::scale_color_manual(name = "Category", values = cols, drop = FALSE)
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

  list(
    plot      = p,
    table     = MCT,
    pvaltable = pvaltable,
    data_used = dplyr::as_tibble(tData)
  )
}
