#' Add r-values and exactly one set of significance stars to a correlations heatmap
#' (pass the full `res` from PlotCorrelationsHeatmap()).
#'
#' @param res  List returned by PlotCorrelationsHeatmap()
#' @param star_from One of: "existing","raw","fdr","P","P_adj","column"
#'                  ("existing" repositions what the chosen plot already used;
#'                   "raw" uses `stars` or computes from `P`;
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

  # --- helpers ---------------------------------------------------------------
  .is_text <- function(ly) inherits(ly$geom, "GeomText") || inherits(ly$geom, "GeomLabel")
  .label_name <- function(expr) {
    if (is.null(expr)) return(NA_character_)
    s <- paste(deparse(expr), collapse = "")
    s <- tolower(s)
    s <- sub("^after_stat\\((.*)\\)$", "\\1", s)   # after_stat(stars) -> stars
    s <- sub("^\\.?data\\$", "", s)                # .data$stars -> stars
    s <- gsub("`", "", s)                          # `stars` -> stars
    gsub("[^a-z0-9_]", "", s)                      # drop other punctuation/spaces
  }
  .same_data <- function(ly, p) is.null(ly$data) || isTRUE(identical(ly$data, p$data))
  .stars_from_p <- function(pvec, n) {
    if (!length(pvec)) return(rep("", n))
    as.character(cut(pvec, breaks = p_breaks, labels = p_labels))
  }

  # pick base plot (keeps semantics intuitive)
  pick_plot <- function() {
    if (star_from %in% c("fdr","P_adj")) {
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from %in% c("raw","P")) {
      if (has_raw) return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
    } else if (star_from == "existing") {
      if (has_fdr) {
        p <- res$FDRCorrected$plot
        used <- any(vapply(p$layers, function(ly) .is_text(ly) && .label_name(ly$mapping$label) %in% c("stars","starsfdr"), logical(1)))
        if (used) return(p)
      }
      if (has_raw) {
        p <- res$Unadjusted$plot
        used <- any(vapply(p$layers, function(ly) .is_text(ly) && .label_name(ly$mapping$label) %in% c("stars","starsfdr"), logical(1)))
        if (used) return(p)
      }
      if (has_fdr) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
    } else if (star_from == "column") {
      if (is.null(star_col)) stop("star_from='column' needs `star_col`.")
      if (has_raw && star_col %in% names(res$Unadjusted$plot$data))   return(res$Unadjusted$plot)
      if (has_fdr && star_col %in% names(res$FDRCorrected$plot$data)) return(res$FDRCorrected$plot)
      if (has_raw) return(res$Unadjusted$plot)
      if (has_fdr) return(res$FDRCorrected$plot)
    }
    stop("Couldn't find a valid ggplot in `res`.")
  }

  p <- pick_plot()
  d <- p$data
  if (!is.data.frame(d) || !r_var %in% names(d)) {
    stop(sprintf("Column '%s' not found in the selected plot's data.", r_var))
  }

  # which stars column did the chosen plot originally map?
  existing_star_col <- {
    nm <- vapply(p$layers, function(ly)
      if (.is_text(ly)) .label_name(ly$mapping$label) else NA_character_, character(1))
    nm <- nm[!is.na(nm)]
    if ("starsfdr" %in% nm) "stars_FDR" else if ("stars" %in% nm) "stars" else NA_character_
  }

  # pick stars exactly as requested
  label_star <- switch(
    star_from,
    existing = {
      if (!is.na(existing_star_col)) as.character(d[[existing_star_col]])
      else if ("stars_FDR" %in% names(d)) as.character(d$stars_FDR)
      else if ("P_adj" %in% names(d))     .stars_from_p(d$P_adj, nrow(d))
      else if ("stars" %in% names(d))     as.character(d$stars)
      else if ("P" %in% names(d))         .stars_from_p(d$P, nrow(d))
      else rep("", nrow(d))
    },
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
      if (is.null(star_col) || !star_col %in% names(d)) stop("Provide a valid `star_col` in plot data.")
      as.character(d[[star_col]])
    }
  )

  # label data
  d$label_r    <- sprintf(paste0("%.", r_digits, "f"), d[[r_var]])
  d$label_star <- label_star

  # idempotent cleanup: drop prior star/r text layers that annotate the tiles
  if (remove_existing_stars && length(p$layers)) {
    drop_names <- c("stars","starsfdr","labelstar","labelr","labelcomb")
    keep <- vapply(p$layers, function(ly) {
      if (!.is_text(ly)) return(TRUE)
      if (!.same_data(ly, p)) return(TRUE)  # keep unrelated annotations
      nm <- .label_name(ly$mapping$label)
      if (is.na(nm)) return(TRUE)
      !(nm %in% drop_names)
    }, logical(1))
    p$layers <- p$layers[keep]
  }

  # compose: r below center, chosen stars above center
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
