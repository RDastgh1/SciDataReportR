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
                            star_from = c("raw","fdr","P","P_adj","column"),
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

  # --- pick a base plot (decoupled from star choice): prefer Unadjusted -----
  base_plot <-
    if (is.list(res$Unadjusted) && inherits(res$Unadjusted$plot, "ggplot")) {
      res$Unadjusted$plot
    } else if (is.list(res$FDRCorrected) && inherits(res$FDRCorrected$plot, "ggplot")) {
      res$FDRCorrected$plot
    } else {
      stop("`res` must contain a ggplot in $Unadjusted$plot or $FDRCorrected$plot.")
    }

  d <- base_plot$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the plot data.", r_var))
  }

  # --- compute exactly the stars you asked for ------------------------------
  stars_from_p <- function(pvec) {
    if (!length(pvec)) return(rep("", nrow(d)))
    as.character(cut(pvec, breaks = p_breaks, labels = p_labels))
  }

  label_star <- switch(
    star_from,
    raw   = if ("stars" %in% names(d)) as.character(d$stars)
    else if ("P" %in% names(d)) stars_from_p(d$P)
    else rep("", nrow(d)),
    fdr   = if ("stars_FDR" %in% names(d)) as.character(d$stars_FDR)
    else if ("P_adj" %in% names(d)) stars_from_p(d$P_adj)
    else rep("", nrow(d)),
    P     = { if (!"P" %in% names(d)) stop("Column 'P' not in plot data."); stars_from_p(d$P) },
    P_adj = { if (!"P_adj" %in% names(d)) stop("Column 'P_adj' not in plot data."); stars_from_p(d$P_adj) },
    column= { if (is.null(star_col) || !star_col %in% names(d)) stop("Provide a valid `star_col`."); as.character(d[[star_col]]) }
  )

  # --- label data ------------------------------------------------------------
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- label_star

  # --- strong, idempotent cleanup of any prior text star/r layers -----------
  if (remove_existing_stars && length(base_plot$layers)) {
    # normalize label expressions to names: stars, stars_FDR, label_star, label_r, label_comb
    norm <- function(expr) {
      if (is.null(expr)) return(NA_character_)
      s <- paste(deparse(expr), collapse = "")
      s <- tolower(s)
      s <- gsub("^after_stat\\(|\\)$", "", s)
      s <- gsub("^\\.data\\$", "", s)
      s <- gsub("[^a-z0-9_]", "", s)  # remove backticks/punct
      s
    }
    is_text <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")
    drop_names <- c("stars","starsfdr","labelstar","labelr","labelcomb")

    keep <- vapply(base_plot$layers, function(ly) {
      if (!is_text(ly)) return(TRUE)
      nm <- norm(ly$mapping$label)
      if (is.na(nm)) return(TRUE)
      !(nm %in% drop_names)
    }, logical(1))
    base_plot$layers <- base_plot$layers[keep]
  }

  # --- compose: r below, chosen stars above ---------------------------------
  base_plot +
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
