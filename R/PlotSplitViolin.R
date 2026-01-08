#' Split violin with aligned half-boxplots and a centered significance label
#'
#' Draws a split (left/right) violin for up to two groups at a single x-position,
#' overlays per-group boxplots aligned with each half, and optionally annotates
#' the plot with a p-value significance label (e.g., "*", "**", "***", "ns").
#' The y-axis title is taken from the variable's label attribute when available
#' (via \pkg{labelled} or \pkg{Hmisc}); otherwise the column name is used.
#'
#' @param data A data frame containing `Var`, `Group`, and any covariates.
#' @param Var Column of the numeric outcome to plot. Tidy-eval friendly (unquoted or string).
#' @param Group Column of the grouping variable. Must have 1 or 2 unique values.
#'   No factor relabeling is performed inside; set levels externally if needed.
#' @param covars Character vector of covariate column names (default `"Age"`).
#'   Use `character(0)` for no covariates.
#' @param nonparametric Logical. If `FALSE` (default), compute the p-value using
#'   a linear model and an \pkg{emmeans} contrast on `Group`. If `TRUE`:
#'   with no covariates, use Wilcoxon rank-sum on `Var ~ Group`; with covariates,
#'   first residualize `Var` on the covariates (no `Group`), then apply Wilcoxon
#'   rank-sum to the residuals by `Group`. With only one group, no test is performed.
#' @param annotation_text Optional label to draw (e.g., `"*"`, `"**"`, `"***"`,
#'   `"ns"`, `"★"`). If `NULL`, the label is derived from the computed p-value.
#' @param show_ns Logical; if `TRUE`, show `"ns"` for non-significant results
#'   when `annotation_text` is `NULL`.
#' @param left_group Optional. The `Group` value to place on the **left** half.
#'   If omitted, the first factor level (or first observed value) is used.
#' @param color_palette Optional color specification. Supply a named vector keyed
#'   by `Group` values; if unnamed, values are interpreted as `(left, right)`.
#' @param box_offset Horizontal offset for the left/right boxplots from the center
#'   (single x = 1). Small values (e.g., `0.10-0.14`) sit the boxes inside each half.
#' @param box_width Boxplot width (default `0.15`). For clean separation,
#'   keep `box_width <= 2 * box_offset`.
#' @param star_from Where to anchor the vertical position of the significance label:
#'   one of `"quantile"` (default), `"data_max"`, or `"whisker"`.
#' @param star_quantile Quantile used when `star_from = "quantile"` (default `0.995`).
#' @param star_pad Fraction of the y-range added above the anchor for the label (default `0.03`).
#' @param headroom Fraction of the y-range added to the top limit to avoid clipping (default `0.08`).
#' @param star_size Text size for the significance label (default `6`).
#'
#' @return A \pkg{ggplot2} object.
#' @export
PlotSplitViolin <- function(
    data,
    Var,
    Group,
    covars = NULL,
    nonparametric = FALSE,
    annotation_text = NULL,
    show_ns = FALSE,
    # layout / colors
    left_group = NULL,
    color_palette = NULL,
    box_offset = 0.11,
    box_width  = 0.15,
    # star placement
    star_from = c("quantile","data_max","whisker"),
    star_quantile = 0.995,
    star_pad = 0.03,
    headroom = 0.08,
    star_size = 6
) {
  stopifnot(requireNamespace("gghalves", quietly = TRUE),
            requireNamespace("ggpubr",   quietly = TRUE),
            requireNamespace("emmeans",  quietly = TRUE))

  star_from <- match.arg(star_from)

  # normalize covars: allow NULL from caller
  if (is.null(covars)) covars <- character(0)

  # tidy-eval: allow Var/Group passed as symbols or strings
  Var   <- rlang::ensym(Var)
  Group <- rlang::ensym(Group)
  var_nm <- rlang::as_string(Var)
  grp_nm <- rlang::as_string(Group)

  need <- unique(c(var_nm, grp_nm, covars))
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Missing columns in `data`: ", paste(miss, collapse = ", "))

  df <- data |>
    dplyr::select(!!Group, !!Var, dplyr::all_of(covars)) |>
    tidyr::drop_na(!!Group, !!Var, dplyr::all_of(covars)) |>
    # KEY FIX: force grouping to be discrete for plotting and modeling
    dplyr::mutate(.Group = as.factor(!!Group))

  # use .Group going forward
  Group_plot   <- rlang::sym(".Group")
  grp_nm_plot  <- ".Group"

  # determine left/right using factor levels
  g_order <- levels(df[[grp_nm_plot]])
  g_order <- g_order[g_order %in% unique(df[[grp_nm_plot]])]

  if (length(g_order) == 0) stop("`Group` has no non-missing values.")
  if (length(g_order) > 2)  stop("`Group` must have <= 2 unique values for a split violin.")

  if (!is.null(left_group)) left_group <- as.character(left_group)

  if (!is.null(left_group)) {
    if (!left_group %in% g_order) stop("`left_group` not found in `Group`.")
    g_left  <- left_group
    g_right <- setdiff(g_order, left_group)
    g_right <- if (length(g_right)) g_right else NULL
  } else {
    g_left  <- g_order[1]
    g_right <- if (length(g_order) == 2) g_order[2] else NULL
  }

  # palette in left→right order (names must match factor levels)
  present_groups <- as.character(stats::na.omit(c(g_left, g_right)))

  if (is.null(color_palette)) {
    base_cols <- c("#E8A007", "#8C1A45")
    color_palette <- stats::setNames(base_cols[seq_len(length(present_groups))], present_groups)
  } else if (is.null(names(color_palette))) {
    color_palette <- stats::setNames(color_palette[seq_len(length(present_groups))], present_groups)
  } else {
    color_palette <- color_palette[present_groups]
  }

  # y label
  get_var_label <- function(v) {
    lab <- attr(v, "label", exact = TRUE)
    if (is.null(lab) && requireNamespace("labelled", quietly = TRUE))
      lab <- tryCatch(labelled::var_label(v, unlist = TRUE), error = function(e) NULL)
    if (is.null(lab) && requireNamespace("Hmisc", quietly = TRUE))
      lab <- tryCatch({
        out <- Hmisc::label(v)
        if (is.character(out) && length(out) == 1 && nzchar(out)) out else NULL
      }, error = function(e) NULL)
    if (is.null(lab) || !nzchar(lab)) lab <- var_nm
    as.character(lab)
  }

  y_lab  <- get_var_label(data[[var_nm]])
  y_vals <- df[[var_nm]]
  y_min  <- suppressWarnings(min(y_vals, na.rm = TRUE)); if (!is.finite(y_min))  y_min  <- 0
  y_maxv <- suppressWarnings(max(y_vals, na.rm = TRUE)); if (!is.finite(y_maxv)) y_maxv <- 1
  dr     <- y_maxv - y_min; if (!is.finite(dr) || dr <= 0) dr <- 1

  # p -> stars (only if not provided; skip if single group)
  stars_from_p <- function(p) {
    if (is.na(p)) return(if (show_ns) "ns" else "")
    if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < .05) "*" else if (show_ns) "ns" else ""
  }

  if (is.null(annotation_text) && !is.null(g_right)) {
    rhs <- paste(c(grp_nm_plot, covars), collapse = " + ")
    fml <- stats::as.formula(paste0("`", var_nm, "` ~ ", rhs))

    if (nonparametric) {
      if (length(covars) == 0) {
        p_val <- tryCatch(
          stats::wilcox.test(
            stats::as.formula(paste0("`", var_nm, "` ~ `", grp_nm_plot, "`")),
            data = df, exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
      } else {
        # residualize on covariates (exclude Group), then Wilcoxon on residuals by Group
        cov_fml <- stats::as.formula(paste0("`", var_nm, "` ~ ", paste(covars, collapse = " + ")))
        res <- tryCatch(stats::residuals(stats::lm(cov_fml, data = df)), error = function(e) NA_real_)
        p_val <- tryCatch(stats::wilcox.test(res ~ df[[grp_nm_plot]], exact = FALSE)$p.value,
                          error = function(e) NA_real_)
      }
    } else {
      p_val <- tryCatch({
        fit <- stats::lm(fml, data = df)
        as.data.frame(
          emmeans::contrast(emmeans::emmeans(fit, grp_nm_plot),
                            method = "revpairwise", adjust = "none")
        )$p.value[1]
      }, error = function(e) NA_real_)

      if (is.na(p_val)) {
        p_val <- tryCatch(
          stats::t.test(
            stats::as.formula(paste0("`", var_nm, "` ~ `", grp_nm_plot, "`")),
            data = df
          )$p.value,
          error = function(e) NA_real_
        )
      }
    }

    annotation_text <- stars_from_p(p_val)
  } else if (is.null(annotation_text)) {
    annotation_text <- ""
  }

  # star position
  safe_whisker <- function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    out <- tryCatch(grDevices::boxplot.stats(z)$stats[5], error = function(e) NA_real_)
    if (!is.finite(out)) max(z, na.rm = TRUE) else out
  }

  anchor <- switch(
    star_from,
    quantile = as.numeric(stats::quantile(y_vals, probs = star_quantile, na.rm = TRUE, names = FALSE)),
    data_max = y_maxv,
    whisker  = if (!is.null(g_right)) {
      max(
        safe_whisker(y_vals[df[[grp_nm_plot]] == g_left]),
        safe_whisker(y_vals[df[[grp_nm_plot]] == g_right]),
        na.rm = TRUE
      )
    } else safe_whisker(y_vals)
  )

  if (!is.finite(anchor)) anchor <- y_maxv
  y_star   <- anchor + star_pad * dr
  ylim_top <- max(y_star + headroom * dr, y_maxv + headroom * dr)
  ylim_bot <- y_min - 0.05 * dr

  # build plot (fixed x for boxes; fill = .Group for both)
  df_left <- df[df[[grp_nm_plot]] == g_left, , drop = FALSE]

  p <- ggplot2::ggplot() +
    gghalves::geom_half_violin(
      data = df_left,
      ggplot2::aes(x = 1, y = !!Var, fill = !!Group_plot),
      side = "l", color = NA, trim = FALSE, alpha = 0.8
    )

  if (!is.null(g_right)) {
    df_right <- df[df[[grp_nm_plot]] == g_right, , drop = FALSE]
    p <- p +
      gghalves::geom_half_violin(
        data = df_right,
        ggplot2::aes(x = 1, y = !!Var, fill = !!Group_plot),
        side = "r", color = NA, trim = FALSE, alpha = 0.8
      )
  }

  p <- p +
    ggplot2::geom_boxplot(
      data = dplyr::mutate(df_left, .x = 1 - if (is.null(g_right)) 0 else box_offset),
      ggplot2::aes(x = .x, y = !!Var, fill = !!Group_plot),
      width = box_width, alpha = 0.60, outlier.shape = NA,
      color = "black", notch = FALSE, linewidth = 0.6, position = "identity"
    )

  if (!is.null(g_right)) {
    df_right <- df[df[[grp_nm_plot]] == g_right, , drop = FALSE]
    p <- p +
      ggplot2::geom_boxplot(
        data = dplyr::mutate(df_right, .x = 1 + box_offset),
        ggplot2::aes(x = .x, y = !!Var, fill = !!Group_plot),
        width = box_width, alpha = 0.60, outlier.shape = NA,
        color = "black", notch = FALSE, linewidth = 0.6, position = "identity"
      )
  }

  p <- p +
    ggplot2::scale_fill_manual(
      values = color_palette,
      breaks = present_groups,
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = y_lab, title = NULL) +
    ggplot2::coord_cartesian(ylim = c(ylim_bot, ylim_top), clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 10)
    )

  if (!is.null(annotation_text) && nzchar(annotation_text)) {
    ann_df <- data.frame(
      label = annotation_text,
      y.position = y_star,
      xmin = 1, xmax = 1,
      group1 = 1, group2 = 1
    )

    p <- p + ggpubr::stat_pvalue_manual(
      data = ann_df, label = "label", xmin = "xmin", xmax = "xmax",
      y.position = "y.position", remove.bracket = TRUE, tip.length = 0,
      size = star_size, color = "black"
    )
  }

  p
}
