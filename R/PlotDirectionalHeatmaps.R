#' Create directional heatmaps across continuous & binary variables
#'
#' Combines:
#'  - Continuous~Continuous (Pearson/Spearman)
#'  - Binary~Binary (Phi; 1 == PositiveLevel)
#'  - Binary~Continuous (r_pb; 1 == PositiveLevel)
#' into a single square heatmap with raw and FDR-star overlays.
#'
#' Constant variables (no variation in the current data) are
#' automatically excluded before computing any tiles.
#'
#' @param Data A dataframe.
#' @param xVars Character vector of variables for x-axis (subset of Data cols).
#'              If NULL, uses all detected continuous + binary vars.
#' @param yVars Character vector for y-axis (defaults to xVars if NULL).
#' @param Relabel Logical; use sjlabelled variable labels if present.
#' @param Ordinal Logical; reserved for future use.
#' @return list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)
#' @export
PlotDirectionalHeatmaps <- function (Data, xVars = NULL, yVars = NULL, Relabel = TRUE, Ordinal = TRUE) {

  # ---- Validate inputs ----------------------------------------------------
  if (!is.data.frame(Data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  # ---- Prepare data: protect against make.names() behavior downstream -----
  orig_names <- colnames(Data)
  safe_names <- make.names(orig_names, unique = TRUE)

  # Map original -> safe and safe -> original
  map_orig_to_safe <- stats::setNames(safe_names, orig_names)
  map_safe_to_orig <- stats::setNames(orig_names, safe_names)

  Data_safe <- Data
  colnames(Data_safe) <- safe_names

  # Helper to map vectors while preserving unknowns
  to_safe <- function(v) {
    v <- as.character(v)
    out <- unname(map_orig_to_safe[v])
    out[is.na(out)] <- v[is.na(out)]
    out
  }
  to_orig <- function(v) {
    v <- as.character(v)
    out <- unname(map_safe_to_orig[v])
    out[is.na(out)] <- v[is.na(out)]
    out
  }

  # ---- helpers ------------------------------------------------------------
  uniq2 <- function(x) length(unique(stats::na.omit(x)))

  has_two_levels <- function(v_safe) {
    uniq2(Data_safe[[v_safe]]) == 2
  }

  has_variance <- function(v_safe) {
    x <- Data_safe[[v_safe]]
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    length(x) >= 2 && stats::var(x) > 0
  }

  # ---- 1) Determine candidate vars ---------------------------------------
  if (is.null(xVars)) {
    Cont_all_safe <- getNumVars(Data_safe)
    Cat_all_safe  <- getBinaryVars(Data_safe)
    xVars_safe <- colnames(Data_safe)[colnames(Data_safe) %in% c(Cat_all_safe, Cont_all_safe)]
    xVars_orig <- to_orig(xVars_safe)
  } else {
    xVars_orig <- as.character(xVars)
    xVars_safe <- to_safe(xVars_orig)
  }

  if (is.null(yVars)) {
    yVars_safe <- xVars_safe
    yVars_orig <- xVars_orig
  } else {
    yVars_orig <- as.character(yVars)
    yVars_safe <- to_safe(yVars_orig)
  }

  cand_safe <- unique(c(xVars_safe, yVars_safe))

  # Guard against truly missing vars (in either naming universe)
  missing_safe <- setdiff(cand_safe, colnames(Data_safe))
  if (length(missing_safe) > 0) {
    missing_orig <- to_orig(missing_safe)
    stop(
      paste0(
        "These variables were not found in Data: ",
        paste(missing_orig, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # initial typing from your helpers, but we'll enforce non-constant below
  Cat0_safe  <- intersect(getBinaryVars(dplyr::select(Data_safe, dplyr::all_of(cand_safe))), cand_safe)
  Cont0_safe <- intersect(getNumVars(   dplyr::select(Data_safe, dplyr::all_of(cand_safe))), cand_safe)

  # ---- 2) Drop constants BEFORE any subcalls ------------------------------
  CatVars_safe  <- Cat0_safe[vapply(Cat0_safe, has_two_levels, logical(1))]
  ContVars_safe <- Cont0_safe[vapply(Cont0_safe, has_variance, logical(1))]

  dropped_safe <- setdiff(cand_safe, c(CatVars_safe, ContVars_safe))
  if (length(dropped_safe)) {
    reasons <- vapply(dropped_safe, function(vs) {
      if (vs %in% Cont0_safe && !has_variance(vs))    return("no variance")
      if (vs %in% Cat0_safe  && !has_two_levels(vs))  return("not two levels")
      if (!(vs %in% c(Cat0_safe, Cont0_safe)))        return("not binary/continuous")
      "excluded"
    }, character(1))

    Excluded <- tibble::tibble(
      Variable = to_orig(dropped_safe),
      Reason   = reasons
    )
  } else {
    Excluded <- tibble::tibble(Variable = character(), Reason = character())
  }

  # refresh axes after exclusions
  xVars_safe <- xVars_safe[xVars_safe %in% c(CatVars_safe, ContVars_safe)]
  yVars_safe <- yVars_safe[yVars_safe %in% c(CatVars_safe, ContVars_safe)]

  xVars_orig <- to_orig(xVars_safe)
  yVars_orig <- to_orig(yVars_safe)

  # Build a single binary mapping (stable 1 == PositiveLevel)
  BinaryMapping <- NULL
  if (length(CatVars_safe) > 0) {
    BinaryMapping <- createBinaryMapping(Data_safe, CatVars_safe)

    # Best-effort: if mapping is a data frame with a Variable column, map back to original names
    if (is.data.frame(BinaryMapping) && "Variable" %in% names(BinaryMapping)) {
      BinaryMapping$Variable <- to_orig(BinaryMapping$Variable)
    }

    # Best-effort: if mapping is a named list keyed by variables, rename to original
    if (is.list(BinaryMapping) && !is.data.frame(BinaryMapping) && !is.null(names(BinaryMapping))) {
      nm <- names(BinaryMapping)
      names(BinaryMapping) <- to_orig(nm)
    }
  }

  # ---- 3) Collect tiles from each sub-plotter -----------------------------
  df_Combined <- tibble::tibble(
    XVar        = character(),
    YVar        = character(),
    correlation = numeric(),
    test        = character()
  )

  # Continuous ~ Continuous
  if (length(ContVars_safe) > 0) {
    O_ContCont <- PlotCorrelationsHeatmap(Data_safe, ContVars_safe, Relabel = Relabel)
    df_ContCont <- O_ContCont$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = R, test = O_ContCont$method)

    df_Combined <- dplyr::bind_rows(df_Combined, df_ContCont)
  }

  # Binary ~ Binary (Phi)
  if (length(CatVars_safe) > 1) {
    O_CatCat <- PlotPhiHeatmap(Data_safe, CatVars_safe, Relabel = Relabel, binary_map = BinaryMapping)
    df_CatCat <- O_CatCat$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = Phi, test = "Phi")

    df_Combined <- dplyr::bind_rows(df_Combined, df_CatCat)
  }

  # Binary ~ Continuous (point-biserial)
  if (length(CatVars_safe) > 0 && length(ContVars_safe) > 0) {
    O_CatCont <- PlotPointCorrelationsHeatmap(
      Data_safe, CatVars_safe, ContVars_safe,
      Relabel = Relabel, Ordinal = Ordinal, binary_map = BinaryMapping
    )

    df_CatCont <- O_CatCont$Unadjusted$plot$data %>%
      dplyr::mutate(stars = `p<.05`, stars_FDR = p.adj.signif) %>%
      dplyr::select(-`p<.05`, -p.adj.signif) %>%
      dplyr::mutate(xVar = CategoricalVariable, yVar = ContinuousVariable)

    # make square by duplicating with swapped axes/labels
    df_CatCont1 <- df_CatCont %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable) %>%
      dplyr::mutate(XVar = xVar, YVar = yVar) %>%
      dplyr::select(-xVar, -yVar)

    df_CatCont2 <- df_CatCont %>%
      dplyr::mutate(XVar = yVar, YVar = xVar) %>%
      dplyr::select(-xVar, -yVar) %>%
      dplyr::mutate(XL = XLabel, YL = YLabel) %>%
      dplyr::mutate(YLabel = XL, XLabel = YL) %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable, -XL, -YL)

    df_CatContSquare <- dplyr::bind_rows(df_CatCont1, df_CatCont2)
    df_CatContSquare$test <- "Point Correlation"

    df_Combined <- dplyr::bind_rows(df_Combined, df_CatContSquare)
  }

  # ---- 4) Map back to original var names and filter to requested axes -----
  # Ensure we have XVar/YVar columns and convert them back to original names
  if ("XVar" %in% names(df_Combined)) df_Combined$XVar <- to_orig(df_Combined$XVar)
  if ("YVar" %in% names(df_Combined)) df_Combined$YVar <- to_orig(df_Combined$YVar)

  df_Combined_plot <- df_Combined %>%
    dplyr::filter(XVar %in% xVars_orig, YVar %in% yVars_orig)

  # handle fully-empty safely
  if (nrow(df_Combined_plot) == 0L) {
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::ggtitle("No valid variable pairs after excluding constants")

    Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)
    FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)

    return(list(
      Unadjusted   = Unadjusted,
      FDRCorrected = FDRCorrected,
      Relabel      = Relabel,
      BinaryMapping= BinaryMapping,
      Excluded     = Excluded
    ))
  }

  # Lock label order using the original variable order
  ordered_xlabels <- vapply(xVars_orig, function(var) {
    df_Combined_plot$XLabel[df_Combined_plot$XVar == var][1]
  }, character(1))
  ordered_xlabels <- unique(ordered_xlabels[!is.na(ordered_xlabels)])

  ordered_ylabels <- vapply(yVars_orig, function(var) {
    df_Combined_plot$YLabel[df_Combined_plot$YVar == var][1]
  }, character(1))
  ordered_ylabels <- unique(ordered_ylabels[!is.na(ordered_ylabels)])

  df_Combined_plot$XLabel <- factor(df_Combined_plot$XLabel, levels = ordered_xlabels)
  df_Combined_plot$YLabel <- factor(df_Combined_plot$YLabel, levels = ordered_ylabels)

  # ---- 5) Build layered heatmap ------------------------------------------
  p <- ggplot2::ggplot()

  if (nrow(dplyr::filter(df_Combined_plot, test == "pearson")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot, test == "pearson"),
        ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = "r") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "spearman")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot, test == "spearman"),
        ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = "\u03C1") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "Phi")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot, test == "Phi"),
        ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(
        limits = c(-1, 1),
        name   = "\u03A6",
        low    = scales::muted("purple"),
        high   = scales::muted("green")
      ) +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot, test == "Point Correlation")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot, test == "Point Correlation"),
        ggplot2::aes(x = XLabel, y = YLabel, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(
        limits = c(-1, 1),
        name   = expression(r[pb]),
        low    = scales::muted("#FFA500"),
        high   = scales::muted("#008080")
      ) +
      ggnewscale::new_scale_fill()
  }

  # reset limits to preserve ordering
  p <- p +
    ggplot2::scale_x_discrete(limits = levels(df_Combined_plot$XLabel)) +
    ggplot2::scale_y_discrete(limits = levels(df_Combined_plot$YLabel))

  # annotate with raw and FDR stars
  p_raw <- p +
    ggplot2::geom_text(
      data = df_Combined_plot,
      ggplot2::aes(x = XLabel, y = YLabel, label = stars),
      color = "black"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  p_FDR <- p +
    ggplot2::geom_text(
      data = df_Combined_plot,
      ggplot2::aes(x = XLabel, y = YLabel, label = stars_FDR),
      color = "black"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  # ---- 6) Return ----------------------------------------------------------
  Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_raw)
  FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_FDR)

  list(
    Unadjusted    = Unadjusted,
    FDRCorrected  = FDRCorrected,
    Relabel       = Relabel,
    BinaryMapping = BinaryMapping,
    Excluded      = Excluded
  )
}
