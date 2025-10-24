#' Add r-values and exactly one set of significance stars to a correlations heatmap
#' built by PlotCorrelationsHeatmap().
#'
#' @param res  List returned by PlotCorrelationsHeatmap()
#' @param star_from One of: "raw","fdr","P","P_adj","column"
#'                  ("raw" uses `stars` or computes from `P`;
#'                   "fdr" uses `stars_FDR` or computes from `P_adj`;
#'                   "P"/"P_adj" compute from those columns;
#'                   "column" uses `star_col`.)
#' @param star_col Column to use when star_from = "column"
#' @param r_var    Correlation column (default "R")
#' @param r_digits Decimal places for r labels
#' @param r_size   Text size for r labels
#' @param r_color  Color for r labels
#' @param r_nudge_y Vertical nudge for r (negative = below center)
#' @param star_size Text size for star labels (default 3 to match originals)
#' @param star_color Color for star labels
#' @param star_nudge_y Vertical nudge for stars (positive = above center)
#' @param p_breaks Cutpoints for p→stars
#' @param p_labels Labels for p→stars
#' @return ggplot
#' @export
add_r_and_stars <- function(res,
                            star_from = c("raw","fdr","P","P_adj","column"),
                            star_col = NULL,
                            r_var = "R",
                            r_digits = 2,
                            r_size = 3,
                            r_color = "black",
                            r_nudge_y = -0.28,
                            star_size = 3,             # small, like original
                            star_color = "black",
                            star_nudge_y = 0.28,
                            p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            p_labels = c("***","**","*","")) {
  star_from <- match.arg(star_from)

  # ----- helpers -------------------------------------------------------------
  .is_text <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")
  .label_name <- function(expr) {
    if (is.null(expr)) return(NA_character_)
    s <- paste(deparse(expr), collapse = "")
    s <- tolower(s)
    s <- sub("^after_stat\\((.*)\\)$", "\\1", s)  # after_stat(stars) -> stars
    s <- sub("^\\.?data\\$", "", s)               # .data$stars -> stars
    s <- gsub("`", "", s)                         # `stars` -> stars
    gsub("[^a-z0-9_]", "", s)                     # strip punctuation/spaces
  }
  .stars_from_p <- function(pvec, n) {
    if (!length(pvec)) return(rep("", n))
    as.character(cut(pvec, breaks = p_breaks, labels = p_labels))
  }
  .has_plot <- function(x) is.list(x) && inherits(x$plot, "ggplot")

  # ----- choose base plot STRICTLY by requested source -----------------------
  if (star_from %in% c("raw","P")) {
    if (!.has_plot(res$Unadjusted)) stop("Unadjusted plot not found in `res`.")
    p <- res$Unadjusted$plot
  } else if (star_from %in% c("fdr","P_adj")) {
    if (!.has_plot(res$FDRCorrected)) stop("FDRCorrected plot not found in `res`.")
    p <- res$FDRCorrected$plot
  } else { # column
    if (is.null(star_col)) stop("star_from='column' requires `star_col`.")
    if (.has_plot(res$Unadjusted) && star_col %in% names(res$Unadjusted$plot$data)) {
      p <- res$Unadjusted$plot
    } else if (.has_plot(res$FDRCorrected) && star_col %in% names(res$FDRCorrected$plot$data)) {
      p <- res$FDRCorrected$plot
    } else stop("`star_col` not found in either plot's data.")
  }

  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # ----- compute EXACTLY the stars you asked for (no cross-fallbacks) -------
  label_star <- switch(
    star_from,
    raw = {
      if ("stars" %in% names(d)) as.character(d$stars)
      else if ("P" %in% names(d)) .stars_from_p(d$P, nrow(d))
      else stop("For RAW stars, need 'stars' or 'P' in plot data.")
    },
    fdr = {
      if ("stars_FDR" %in% names(d)) as.character(d$stars_FDR)
      else if ("P_adj" %in% names(d)) .stars_from_p(d$P_adj, nrow(d))
      else stop("For FDR stars, need 'stars_FDR' or 'P_adj' in plot data.")
    },
    P = {
      if (!"P" %in% names(d)) stop("Column 'P' not in plot data.")
      .stars_from_p(d$P, nrow(d))
    },
    P_adj = {
      if (!"P_adj" %in% names(d)) stop("Column 'P_adj' not in plot data.")
      .stars_from_p(d$P_adj, nrow(d))
    },
    column = {
      if (!star_col %in% names(d)) stop("`star_col` not found in plot data.")
      as.character(d[[star_col]])
    }
  )

  # ----- label data for drawing ---------------------------------------------
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- label_star

  # ----- robust cleanup: drop prior star/r text layers on the tiles ---------
  if (length(p$layers)) {
    keep <- vapply(p$layers, function(ly) {
      if (!.is_text(ly)) return(TRUE)
      nm <- .label_name(ly$mapping$label)
      # remove any text layer that was a stars layer or a prior helper layer
      is.na(nm) || !(nm %in% c("stars","starsfdr","labelstar","labelr","labelcomb"))
    }, logical(1))
    p$layers <- p$layers[keep]
  }

  # ----- compose: r below, stars above --------------------------------------
  p +
    ggplot2::geom_text(
      data = d, ggplot2::aes(label = label_r),
      inherit.aes = TRUE, size = r_size, color = r_color,
      position = ggplot2::position_nudge(y = r_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = d, ggplot2::aes(label = label_star),
      inherit.aes = TRUE, size = star_size, color = star_color,
      position = ggplot2::position_nudge(y = star_nudge_y),
      check_overlap = TRUE, na.rm = TRUE
    )
}
