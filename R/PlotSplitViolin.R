#' Split violin with aligned half-boxplots, significance label, sample sizes, and label-aware title
#'
#' Draws a split (left/right) violin for up to two groups at a single x-position,
#' overlays per-group boxplots aligned with each half, and optionally annotates
#' the plot with a p-value significance label. Supports displaying sample sizes (n)
#' and automatically using variable labels for axis and title when available.
#'
#' @param data A data frame containing `Var`, `Group`, and optional covariates.
#' @param Var Numeric outcome variable (tidy-eval).
#' @param Group Grouping variable (<= 2 unique values).
#' @param covars Character vector of covariates (default `NULL`).
#' @param nonparametric Logical. If `FALSE`, uses linear model + emmeans contrast.
#' If `TRUE`, uses Wilcoxon (with residualization if covariates are present).
#' @param annotation_text Optional manual annotation (e.g., "*", "ns").
#' @param show_ns Logical; if `TRUE`, display "ns" for non-significant results.
#' @param plot_title Optional custom plot title.
#' @param use_var_label_as_title Logical; if `TRUE`, uses variable label as title.
#' @param show_n Logical; if `TRUE`, display sample size per group.
#' @param n_position Where to display n: "legend" or "top".
#' @param n_size Text size for n labels when shown on top.
#' @param left_group Optional group to force on left side.
#' @param color_palette Optional named vector of colors.
#' @param box_offset Horizontal offset for boxplots.
#' @param box_width Width of boxplots.
#' @param star_from Position method ("quantile","data_max","whisker").
#' @param star_quantile Quantile used for placement.
#' @param star_pad Padding above anchor.
#' @param headroom Extra space above data.
#' @param star_size Text size for annotation.
#'
#' @return A ggplot2 object.
#'
#' @export
PlotSplitViolin <- function(
    data,
    Var,
    Group,
    covars = NULL,
    nonparametric = FALSE,
    annotation_text = NULL,
    show_ns = FALSE,
    plot_title = NULL,
    use_var_label_as_title = FALSE,
    show_n = TRUE,
    n_position = c("legend","top"),
    n_size = 3.5,
    left_group = NULL,
    color_palette = NULL,
    box_offset = 0.11,
    box_width  = 0.15,
    star_from = c("quantile","data_max","whisker"),
    star_quantile = 0.995,
    star_pad = 0.03,
    headroom = 0.08,
    star_size = 6
) {

  # ---- validate inputs ----
  stopifnot(requireNamespace("gghalves", quietly = TRUE),
            requireNamespace("emmeans", quietly = TRUE))

  star_from  <- match.arg(star_from)
  n_position <- match.arg(n_position)

  if (is.null(covars)) covars <- character(0)

  Var   <- rlang::ensym(Var)
  Group <- rlang::ensym(Group)
  var_nm <- rlang::as_string(Var)
  grp_nm <- rlang::as_string(Group)

  need <- unique(c(var_nm, grp_nm, covars))
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Missing columns in `data`: ", paste(miss, collapse = ", "))
  }

  # ---- prepare data ----
  df <- data %>%
    dplyr::select(!!Group, !!Var, dplyr::all_of(covars)) %>%
    tidyr::drop_na(!!Group, !!Var, dplyr::all_of(covars)) %>%
    dplyr::mutate(.Group = as.factor(!!Group))

  if (nlevels(df$.Group) > 2) {
    stop("`Group` must have <= 2 unique values.")
  }

  # ---- sample sizes ----
  n_df <- df %>% dplyr::count(.Group, name = "n")

  g_order <- levels(df$.Group)

  if (!is.null(left_group)) {
    g_left  <- as.character(left_group)
    g_right <- setdiff(g_order, g_left)
  } else {
    g_left  <- g_order[1]
    g_right <- if (length(g_order) == 2) g_order[2] else NULL
  }

  present_groups <- na.omit(c(g_left, g_right))

  # ---- colors ----
  if (is.null(color_palette)) {
    base_cols <- c("#E8A007", "#8C1A45")
    color_palette <- stats::setNames(base_cols[seq_along(present_groups)], present_groups)
  }

  # ---- legend labels ----
  if (show_n && n_position == "legend") {
    label_map <- n_df %>%
      dplyr::mutate(label = paste0(.Group, " (n=", n, ")"))
    legend_labels <- stats::setNames(label_map$label, label_map$.Group)
  } else {
    legend_labels <- present_groups
  }

  # ---- variable label ----
  get_var_label <- function(v) {
    lab <- attr(v, "label", exact = TRUE)
    if (is.null(lab) && requireNamespace("labelled", quietly = TRUE)) {
      lab <- tryCatch(labelled::var_label(v, unlist = TRUE), error = function(e) NULL)
    }
    if (is.null(lab) && requireNamespace("Hmisc", quietly = TRUE)) {
      lab <- tryCatch(Hmisc::label(v), error = function(e) NULL)
    }
    if (is.null(lab) || !nzchar(lab)) lab <- var_nm
    as.character(lab)
  }

  y_lab <- get_var_label(data[[var_nm]])

  # ---- title logic ----
  if (!is.null(plot_title)) {
    final_title <- plot_title
  } else if (use_var_label_as_title) {
    final_title <- y_lab
  } else {
    final_title <- NULL
  }

  y_vals <- df[[var_nm]]
  y_min  <- min(y_vals, na.rm = TRUE)
  y_maxv <- max(y_vals, na.rm = TRUE)
  dr <- y_maxv - y_min
  if (!is.finite(dr) || dr <= 0) dr <- 1

  # ---- p-value ----
  stars_from_p <- function(p) {
    if (is.na(p)) return(if (show_ns) "ns" else "")
    if (p < 1e-3) "***"
    else if (p < 1e-2) "**"
    else if (p < .05) "*"
    else if (show_ns) "ns"
    else ""
  }

  if (is.null(annotation_text) && !is.null(g_right)) {

    if (nonparametric) {
      if (length(covars) == 0) {
        p_val <- tryCatch(
          wilcox.test(df[[var_nm]] ~ df$.Group)$p.value,
          error = function(e) NA_real_
        )
      } else {
        cov_fml <- as.formula(paste0(var_nm, " ~ ", paste(covars, collapse = "+")))
        res <- residuals(lm(cov_fml, data = df))
        p_val <- tryCatch(
          wilcox.test(res ~ df$.Group)$p.value,
          error = function(e) NA_real_
        )
      }
    } else {
      rhs <- paste(c(".Group", covars), collapse = " + ")
      fml <- as.formula(paste0(var_nm, " ~ ", rhs))

      p_val <- tryCatch({
        fit <- lm(fml, data = df)
        as.data.frame(
          emmeans::contrast(emmeans::emmeans(fit, ".Group"),
                            method = "revpairwise")
        )$p.value[1]
      }, error = function(e) NA_real_)

      if (is.na(p_val)) {
        p_val <- tryCatch(
          t.test(df[[var_nm]] ~ df$.Group)$p.value,
          error = function(e) NA_real_
        )
      }
    }

    annotation_text <- stars_from_p(p_val)
  }

  # ---- star placement ----
  safe_whisker <- function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    tryCatch(boxplot.stats(z)$stats[5], error = function(e) max(z))
  }

  anchor <- switch(
    star_from,
    quantile = quantile(y_vals, star_quantile, na.rm = TRUE),
    data_max = y_maxv,
    whisker  = if (!is.null(g_right)) {
      max(
        safe_whisker(y_vals[df$.Group == g_left]),
        safe_whisker(y_vals[df$.Group == g_right]),
        na.rm = TRUE
      )
    } else safe_whisker(y_vals)
  )

  y_star   <- anchor + star_pad * dr
  ylim_top <- max(y_star + headroom * dr, y_maxv + headroom * dr)
  ylim_bot <- y_min - 0.05 * dr

  # ---- build plot ----
  df_left <- df[df$.Group == g_left, ]

  p <- ggplot2::ggplot() +
    gghalves::geom_half_violin(
      data = df_left,
      ggplot2::aes(x = 1, y = !!Var, fill = .Group),
      side = "l", alpha = 0.8
    )

  if (!is.null(g_right)) {
    df_right <- df[df$.Group == g_right, ]
    p <- p +
      gghalves::geom_half_violin(
        data = df_right,
        ggplot2::aes(x = 1, y = !!Var, fill = .Group),
        side = "r", alpha = 0.8
      )
  }

  p <- p +
    ggplot2::geom_boxplot(
      data = dplyr::mutate(df_left, .x = 1 - box_offset),
      ggplot2::aes(x = .x, y = !!Var, fill = .Group),
      width = box_width, outlier.shape = NA
    )

  if (!is.null(g_right)) {
    df_right <- df[df$.Group == g_right, ]
    p <- p +
      ggplot2::geom_boxplot(
        data = dplyr::mutate(df_right, .x = 1 + box_offset),
        ggplot2::aes(x = .x, y = !!Var, fill = .Group),
        width = box_width, outlier.shape = NA
      )
  }

  # ---- n labels ----
  if (show_n && n_position == "top") {
    n_pos_df <- n_df %>%
      dplyr::mutate(
        x = ifelse(.Group == g_left, 1 - box_offset, 1 + box_offset),
        y = y_maxv + 0.08 * dr,
        label = paste0("n=", n)
      )

    p <- p +
      ggplot2::geom_text(
        data = n_pos_df,
        ggplot2::aes(x = x, y = y, label = label),
        size = n_size,
        inherit.aes = FALSE
      )
  }

  # ---- annotation (fixed, no ggpubr dependency) ----
  if (!is.null(annotation_text) && nzchar(annotation_text)) {
    p <- p +
      ggplot2::annotate(
        "text",
        x = 1,
        y = y_star,
        label = annotation_text,
        size = star_size
      )
  }

  # ---- finalize ----
  p +
    ggplot2::scale_fill_manual(values = color_palette, labels = legend_labels) +
    ggplot2::labs(y = y_lab, title = final_title) +
    ggplot2::coord_cartesian(ylim = c(ylim_bot, ylim_top), clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}
