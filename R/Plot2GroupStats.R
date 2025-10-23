#' Plot and Summarize Group Statistics (safer, faster)
#'
#' @param Data data.frame
#' @param Variables character vector of column names to analyze
#' @param VariableCategories optional data.frame with columns: Variable, Category
#' @param impClust character, name of the "important" group (will be plotted as positive direction)
#' @param normalClust character, name of the comparison group
#' @param GroupVar character, column name of grouping variable in `Data`
#' @param missing_threshold proportion in [0,1]; drop vars with > this proportion missing (default 0.80)
#' @param max_levels drop factors with > max_levels (default 10)
#' @param adjust_method p-adjust method for FDR control (default "fdr")
#' @param label_q q-value threshold for labeling points (default 0.05)
#' @return list(plot=gg, tbl_raw=gtsummary, tbl_dummy=gtsummary, pvaltable=data.frame, data_used=tibble)
#' @export
#' @import dplyr ggplot2 gtsummary
Plot2GroupStats <- function(
    Data, Variables, VariableCategories = NULL,
    impClust, normalClust, GroupVar,
    missing_threshold = 0.80, max_levels = 10,
    adjust_method = "fdr", label_q = 0.05
){
  # --- dependency guards for optional goodies
  has_paletteer <- rlang::is_installed("paletteer")
  has_ggrepel   <- rlang::is_installed("ggrepel")
  has_fastdum   <- rlang::is_installed("fastDummies")
  has_stringr   <- rlang::is_installed("stringr")
  has_gtools    <- rlang::is_installed("gtools")

  # --- input checks
  stopifnot(is.data.frame(Data), is.character(Variables), length(Variables) >= 1)
  if (!GroupVar %in% names(Data)) stop("`GroupVar` not found in `Data`.")
  missing_vars <- setdiff(Variables, names(Data))
  if (length(missing_vars)) stop("Variables not found in `Data`: ", paste(missing_vars, collapse = ", "))

  # --- prep
  tData <- dplyr::mutate(Data, GroupVar = .data[[GroupVar]]) |>
    dplyr::select(GroupVar, dplyr::all_of(Variables)) |>
    dplyr::filter(.data$GroupVar %in% c(normalClust, impClust))
  tData$GroupVar <- factor(tData$GroupVar, levels = c(normalClust, impClust))

  # convert characters to factors (excluding GroupVar)
  tData <- tData |>
    dplyr::mutate(dplyr::across(where(is.character) & !dplyr::all_of("GroupVar"), ~factor(.x)))

  # --- drop vars with > missing_threshold missing
  miss_prop <- colMeans(is.na(tData))
  drop_missing <- names(miss_prop[miss_prop > missing_threshold])
  drop_missing <- setdiff(drop_missing, "GroupVar")
  if (length(drop_missing)) tData <- dplyr::select(tData, -dplyr::all_of(drop_missing))

  # --- drop factors with <2 levels and >max_levels
  fac_levels <- vapply(tData, function(x) if (is.factor(x)) nlevels(x) else NA_integer_, 1L)
  drop_low   <- names(fac_levels[!is.na(fac_levels) & fac_levels < 2])
  drop_high  <- names(fac_levels[!is.na(fac_levels) & fac_levels > max_levels])
  tData <- tData |>
    dplyr::select(-dplyr::all_of(setdiff(c(drop_low, drop_high), "GroupVar")))

  # --- drop numeric zero-variance
  num_sd <- vapply(tData, function(x) if (is.numeric(x)) stats::sd(x, na.rm=TRUE) else NA_real_, 1.0)
  drop_nzv <- names(num_sd[is.na(num_sd) | num_sd == 0])
  tData <- tData |>
    dplyr::select(-dplyr::all_of(setdiff(drop_nzv, "GroupVar")))

  # keep a copy for the raw summary
  tRaw <- tData

  # --- gtsummary table on raw data
  gtabp <- tRaw |>
    gtsummary::tbl_summary(
      by = GroupVar,
      type = list(where(is.numeric) ~ "continuous"),
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)")
    ) |>
    gtsummary::add_n() |>
    gtsummary::add_p() |>
    gtsummary::add_q(method = adjust_method) |>
    gtsummary::bold_p() |>
    gtsummary::bold_labels()

  # --- dummy encode (if needed)
  d <- tData
  if (has_fastdum) {
    # dummy only factor columns except GroupVar
    fac_cols <- names(d)[vapply(d, is.factor, TRUE)]
    fac_cols <- setdiff(fac_cols, "GroupVar")
    if (length(fac_cols)) {
      d <- fastDummies::dummy_cols(
        d, select_columns = fac_cols,
        remove_selected_columns = TRUE, remove_first_dummy = TRUE
      )
      # make new dummies factors "0/1" for nice categorical summaries
      newVars <- setdiff(names(d), names(tData))
      d[newVars] <- lapply(d[newVars], function(x) factor(x, levels = c(0,1)))
    }
  }

  # --- gtsummary on (possibly) dummified data
  gtabd <- d |>
    gtsummary::tbl_summary(
      by = GroupVar,
      type = list(where(is.numeric) ~ "continuous",
                  where(is.factor)  ~ "categorical"),
      statistic = list(all_continuous() ~ "{mean} ({sd})",
                       all_categorical() ~ "{n} ({p}%)")
    ) |>
    gtsummary::add_n() |>
    gtsummary::add_p() |>
    gtsummary::add_q(method = adjust_method) |>
    gtsummary::bold_p()

  # --- extract a clean p-value table
  tb <- gtabd$table_body
  pvaltable <- tb |> dplyr::filter(.data$row_type == "label") |>
    dplyr::select(variable, label, dplyr::any_of(c("p.value","q.value")))
  pvaltable <- pvaltable[!is.na(pvaltable$p.value), , drop = FALSE]

  # --- compute effect direction from data (numeric means or dummy 0/1 proportions)
  eff_sign <- function(v, grp) {
    x <- v
    # try to coerce factors "0/1" to numeric
    if (is.factor(x)) {
      xnum <- suppressWarnings(as.numeric(as.character(x)))
      if (!any(is.na(xnum))) x <- xnum
    }
    if (!is.numeric(x)) return(NA_real_)
    m_imp  <- mean(x[grp == impClust], na.rm = TRUE)
    m_norm <- mean(x[grp == normalClust], na.rm = TRUE)
    sign(m_imp - m_norm)
  }

  # align variable names to column names in `d`
  vars_present <- intersect(pvaltable$variable, names(d))
  sgn <- vapply(vars_present, function(v) eff_sign(d[[v]], d$GroupVar), 1.0)
  pvaltable$effect_sign <- sgn[match(pvaltable$variable, vars_present)]

  # signed -log10(p)
  tiny <- .Machine$double.xmin
  pvaltable$signed_log10p <- -log10(pmax(pvaltable$p.value, tiny)) * pvaltable$effect_sign

  # significance labels
  if (has_gtools) {
    pvaltable$Sig <- gtools::stars.pval(pvaltable$p.value)
  } else {
    brks <- c(-Inf, 0.001, 0.01, 0.05, Inf)
    labs <- c("***","**","*","ns")
    pvaltable$Sig <- cut(pvaltable$p.value, breaks = brks, labels = labs, right = TRUE)
  }

  # category mapping (optional)
  if (!is.null(VariableCategories) && all(c("Variable","Category") %in% names(VariableCategories))) {
    pvaltable$Category <- VariableCategories$Category[
      match(pvaltable$variable, VariableCategories$Variable)
    ]
  } else {
    pvaltable$Category <- NA_character_
  }

  # --- plot
  base <- ggplot2::ggplot(
    pvaltable,
    ggplot2::aes(x = .data$signed_log10p,
                 y = .data$variable,
                 color = if (any(!is.na(pvaltable$Category))) .data$Category else .data$Sig,
                 shape = .data$Sig)
  ) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(limits = rev(pvaltable$variable)) +
    ggplot2::xlab("signed -log10(p)   (→ higher in {impClust})") +
    ggplot2::ylab(NULL) +
    ggplot2::ggtitle(paste("---> Higher in", impClust)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1))

  if (has_paletteer) {
    base <- base + paletteer::scale_color_paletteer_d("fishualize::Scarus_quoyi")
  }

  if (has_ggrepel) {
    base <- base + ggrepel::geom_text_repel(
      data = subset(pvaltable, !is.na(.data$q.value) & .data$q.value < label_q),
      ggplot2::aes(label = .data$variable),
      max.overlaps = 30, min.segment.length = 0
    )
  } else {
    base <- base + ggplot2::geom_text(
      data = subset(pvaltable, !is.na(.data$q.value) & .data$q.value < label_q),
      ggplot2::aes(label = .data$variable),
      hjust = -0.05, size = 3
    )
  }

  # --- return
  list(
    p         = base,
    tbl_raw   = gtabp,
    tbl_dummy = gtabd,
    pvaltable = pvaltable,
    data_used = dplyr::as_tibble(d)
  )
}
