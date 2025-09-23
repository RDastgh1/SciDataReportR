#' Plot Phi Correlations Between Binary Variables
#'
#' Computes pairwise phi coefficients (Φ) between binary categorical variables
#' with explicit 0/1 coding (1 == PositiveLevel from `createBinaryMapping()`),
#' then renders heatmap-style plots with raw and FDR-adjusted significance.
#'
#' @param Data A dataframe.
#' @param CatVars Character vector of binary categorical variable names.
#' @param Relabel Logical; if TRUE, uses sjlabelled variable labels for axes.
#' @param binary_map Optional mapping as returned by `createBinaryMapping()`.
#'   If NULL, a mapping is created internally for `CatVars`.
#' @return A list with:
#'   - `Unadjusted`: list(PvalTable, plot)
#'   - `FDRCorrected`: list(PvalTable, plot)
#'   - `method` = "Phi"
#'   - `Relabel`
#'   - `BinaryMapping` (used)
#' @export
PlotPhiHeatmap <- function(Data, CatVars, Relabel = TRUE, binary_map = NULL) {

  # ---- helpers -------------------------------------------------------------
  # encode a binary vector to {0,1} using mapping (1 == PositiveLevel)
  .to01 <- function(x, var, map) {
    # strip haven labels; then coerce to character for matching
    if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)
    pos <- map$PositiveLevel[match(var, map$Variable)]
    if (is.na(pos)) {
      stop("No PositiveLevel found in mapping for variable: ", var)
    }
    as.integer(as.character(x) == pos)
  }

  # is a vector truly binary (2 unique non-NA values)?
  .is_binary <- function(x) length(unique(stats::na.omit(x))) == 2

  # safe chisq p-value for 2x2; fallback to cor.test if needed
  .pval_binbin <- function(x01, y01) {
    keep <- stats::complete.cases(x01, y01)
    x <- x01[keep]; y <- y01[keep]
    if (length(x) < 3L) return(NA_real_)
    tbl <- table(x, y)
    if (all(dim(tbl) == c(2, 2))) {
      # chisq test without continuity correction to match correlation-style tests
      suppressWarnings(stats::chisq.test(tbl, correct = FALSE)$p.value)
    } else {
      # degenerate (one level missing) -> NA
      NA_real_
    }
  }

  # ---- data prep -----------------------------------------------------------
  if (length(CatVars) < 2) {
    stop("Need at least two binary variables for a Phi heatmap.")
  }

  # keep only requested vars, preserve labels
  DataSubset <- Data[CatVars]
  DataSubset <- ReplaceMissingLabels(DataSubset)

  # check all are binary
  non_binary <- vapply(DataSubset, .is_binary, logical(1))
  if (!all(non_binary)) {
    bad <- names(non_binary)[!non_binary]
    stop("The following variables are not binary (exactly 2 unique non-NA values required): ",
         paste(bad, collapse = ", "))
  }

  # mapping: use supplied or create
  if (is.null(binary_map)) {
    binary_map <- createBinaryMapping(DataSubset, CatVars)
  }

  # pre-encode each variable once to 0/1 with 1 == PositiveLevel
  encoded <- lapply(CatVars, function(v) .to01(DataSubset[[v]], v, binary_map))
  names(encoded) <- CatVars

  # ---- pairwise phi and p-values ------------------------------------------
  pairs <- expand.grid(X = CatVars, Y = CatVars, stringsAsFactors = FALSE)

  stat.list <- lapply(seq_len(nrow(pairs)), function(i) {
    vx <- pairs$X[i]; vy <- pairs$Y[i]
    x01 <- encoded[[vx]]
    y01 <- encoded[[vy]]

    keep <- stats::complete.cases(x01, y01)
    x <- x01[keep]; y <- y01[keep]
    nPairs <- length(x)

    if (nPairs < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
      Phi <- NA_real_
      p   <- NA_real_
    } else {
      # Pearson correlation of 0/1 vectors equals Φ
      Phi <- suppressWarnings(stats::cor(x, y))
      # Use χ²-based p-value; fallback if necessary
      p <- .pval_binbin(x, y)
      if (is.na(p)) {
        p <- suppressWarnings(stats::cor.test(x, y)$p.value)
      }
    }

    data.frame(
      XVar   = vx,
      YVar   = vy,
      Phi    = Phi,
      p_value = p,
      nPairs = nPairs,
      stringsAsFactors = FALSE
    )
  })

  stat.test <- do.call(rbind, stat.list)

  # multiple testing (FDR); keep NA as NA
  stat.test$p.adj <- stats::p.adjust(stat.test$p_value, method = "fdr")

  # significance stars
  stat.test <- stat.test %>%
    rstatix::add_significance("p_value", output.col = "p<.05") %>%
    rstatix::add_significance("p.adj",   output.col = "p.adj.signif")

  stat.test$`p<.05`[stat.test$`p<.05` == "ns"] <- ""
  stat.test$p.adj.signif[stat.test$p.adj.signif == "ns"] <- ""
  stat.test$`p<.05` <- factor(stat.test$`p<.05`,
                              levels = c("ns","*","**","***","****"),
                              labels = c("ns","*","**","***","***"))
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif,
                                   levels = c("ns","*","**","***","****"),
                                   labels = c("ns","*","**","***","***"))

  stat.test$test <- "Phi"

  # ---- labels --------------------------------------------------------------
  if (Relabel) {
    XLabel <- vapply(stat.test$XVar, function(v) {
      lab <- sjlabelled::get_label(DataSubset[[v]])
      if (is.null(lab) || is.na(lab) || lab == "") v else lab
    }, character(1))
    YLabel <- vapply(stat.test$YVar, function(v) {
      lab <- sjlabelled::get_label(DataSubset[[v]])
      if (is.null(lab) || is.na(lab) || lab == "") v else lab
    }, character(1))
  } else {
    XLabel <- stat.test$XVar
    YLabel <- stat.test$YVar
  }
  stat.test$XLabel <- XLabel
  stat.test$YLabel <- YLabel

  # tooltip
  PlotText <- paste0(
    "</br>X Label: ", stat.test$XLabel,
    "</br>X Var: ",   stat.test$XVar,
    "</br>Y Label: ", stat.test$YLabel,
    "</br>Y Var: ",   stat.test$YVar,
    "</br> Φ: ",      signif(stat.test$Phi, 3),
    "</br> P-Value: ", signif(stat.test$p_value, 3), " ", stat.test$`p<.05`,
    "</br> FDR P: ",   signif(stat.test$p.adj, 3),   " ", stat.test$p.adj.signif,
    "</br> nPairs: ",  stat.test$nPairs
  )

  # ---- plots ---------------------------------------------------------------
  p <- ggplot2::ggplot(
    stat.test,
    ggplot2::aes(x = XLabel, y = YLabel, fill = Phi, text = PlotText)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = `p<.05`), color = "black") +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      name = "\u03A6",                       # Phi
      low  = scales::muted("purple"),
      high = scales::muted("green")
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  p_FDR <- ggplot2::ggplot(
    stat.test,
    ggplot2::aes(x = XLabel, y = YLabel, fill = Phi, text = PlotText)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = p.adj.signif), color = "black") +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      name = "\u03A6",
      low  = scales::muted("purple"),
      high = scales::muted("green")
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # ---- return --------------------------------------------------------------
  M     <- list(PvalTable = stat.test, plot = p)
  M_FDR <- list(PvalTable = stat.test, plot = p_FDR)

  list(
    Unadjusted    = M,
    FDRCorrected  = M_FDR,
    method        = "Phi",
    Relabel       = Relabel,
    BinaryMapping = binary_map
  )
}
