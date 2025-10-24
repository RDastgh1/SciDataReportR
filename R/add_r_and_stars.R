#' Add r-values and exactly one set of significance stars to a correlations heatmap
#' (pass the full `res` from PlotCorrelationsHeatmap()).
#'
#' @param res  List returned by PlotCorrelationsHeatmap()
#' @param star_from One of: "raw","fdr","P","P_adj","column"
#'                  ("raw" uses `stars` or computes from `P`;
#'                   "fdr" uses `stars_FDR` or computes from `P_adj`;
#'                   "P"/"P_adj" compute from those columns;
#'                   "column" uses `star_col`.)
#' @param star_col Column name to use when star_from = "column"
#' @param r_var    Correlation column (default "R")
#' @param r_digits Decimal places for r labels
#' @param r_size   Text size for r labels
#' @param r_color  Color for r labels
#' @param r_nudge_y Vertical nudge for r (negative = below center)
#' @param star_size Text size for star labels
#' @param star_color Color for star labels
#' @param star_nudge_y Vertical nudge for stars (positive = above center)
#' @param remove_existing_stars If TRUE, drop any pre-existing star/r text layers
#' @param p_breaks Cutpoints for p→stars
#' @param p_labels Labels for p→stars
#' @return ggplot
#' @export
add_r_and_stars <- function(res,
                            star_from = c("existing","raw","fdr","P","P_adj","column"),
                            star_col = NULL,
                            r_var = "R",
                            r_digits = 2,
                            r_size = 3,
                            r_color = "black",
                            r_nudge_y = -0.28,
                            star_size = 6,
                            star_color = "black",
                            star_nudge_y = 0.28,
                            remove_existing_stars = TRUE,
                            p_breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            p_labels = c("***","**","*","")) {
  star_from <- match.arg(star_from)

  has_fdr <- is.list(res$FDRCorrected) && inherits(res$FDRCorrected$plot, "ggplot")
  has_raw <- is.list(res$Unadjusted)   && inherits(res$Unadjusted$plot,   "ggplot")

  safe_label <- function(x) if (is.null(x)) NA_character_ else paste(deparse(x), collapse = "")
  is_text    <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")

  plot_has_star_map <- function(p) {
    if (!inherits(p, "ggplot")) return(FALSE)
    any(vapply(p$layers, function(ly) {
      if (!is_text(ly)) return(FALSE)
      lab <- safe_label(ly$mapping$label)
      !is.na(lab) && grepl("^stars", lab, ignore.case = TRUE)
    }, logical(1)))
  }

  pick_plot <- function() {
    if (star_from %in% c("fdr","P_adj")) {
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from %in% c("raw","P")) {
      if (has_raw) return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
    } else if (star_from == "existing") {
      if (has_fdr && plot_has_star_map(res$FDRCorrected$plot)) return(res$FDRCorrected$plot)
      if (has_raw && plot_has_star_map(res$Unadjusted$plot))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from == "column") {
      if (is.null(star_col)) stop("star_from='column' needs `star_col`.")
      if (has_fdr && star_col %in% names(res$FDRCorrected$plot$data)) return(res$FDRCorrected$plot)
      if (has_raw && star_col %in% names(res$Unadjusted$plot$data))   return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    }
    stop("Couldn't find a valid ggplot in `res`.")
  }

  p <- pick_plot()
  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # which stars column did the chosen plot originally use?
  existing_star_col <- {
    labs <- vapply(p$layers, function(ly) if (is_text(ly)) safe_label(ly$mapping$label) else NA_character_, character(1))
    labs <- labs[!is.na(labs)]
    labs <- labs[grepl("^stars", labs, ignore.case = TRUE)]
    labs[labs %in% names(d)][1]
  }

  stars_from_p <- function(pvec) as.character(cut(pvec, breaks = p_breaks, labels = p_labels))

  pick_stars <- function() {
    cols <- names(d)
    if (star_from == "existing") {
      if (!is.na(existing_star_col)) return(as.character(d[[existing_star_col]]))
      if ("stars_FDR" %in% cols) return(as.character(d$stars_FDR))
      if ("P_adj"    %in% cols)  return(stars_from_p(d$P_adj))
      if ("stars"    %in% cols)  return(as.character(d$stars))
      if ("P"        %in% cols)  return(stars_from_p(d$P))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "raw") {
      if ("stars" %in% cols) return(as.character(d$stars))
      if ("P"     %in% cols) return(stars_from_p(d$P))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "fdr") {
      if ("stars_FDR" %in% cols) return(as.character(d$stars_FDR))
      if ("P_adj"    %in% cols)  return(stars_from_p(d$P_adj))
      return(rep(NA_character_, nrow(d)))
    }
    if (star_from == "P")     return(stars_from_p(d$P))
    if (star_from == "P_adj") return(stars_from_p(d$P_adj))
    if (star_from == "column") {
      if (!star_col %in% names(d)) stop("`star_col` not found in plot data.")
      return(as.character(d[[star_col]]))
    }
    rep(NA_character_, nrow(d))
  }

  # build labels
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- pick_stars()

  # 🔧 STRONG, IDEMPOTENT CLEANUP: always drop any prior star/r text layers
  if (remove_existing_stars && length(p$layers)) {
    drop_re <- "(^stars$|^stars_FDR$|^label_star$|^label_r$|^label_comb$)"
    keep <- vapply(p$layers, function(ly) {
      if (!is_text(ly)) return(TRUE)
      lab <- safe_label(ly$mapping$label)
      if (is.na(lab)) return(TRUE)
      !grepl(drop_re, lab)
    }, logical(1))
    p$layers <- p$layers[keep]
  }

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
