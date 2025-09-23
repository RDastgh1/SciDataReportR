#' Plot Point-Biserial Correlations Between Binary and Continuous Variables
#'
#' Calculates point-biserial correlations (binary vs continuous) with explicit 0/1 coding
#' where 1 == PositiveLevel from `createBinaryMapping()`, and renders heatmap-style tiles.
#'
#' @param Data A dataframe.
#' @param CatVars Character vector of binary categorical variables.
#' @param ContVars Character vector of continuous variables.
#' @param Covariates Optional covariates (reserved).
#' @param Relabel Logical; use sjlabelled variable labels for axes.
#' @param Ordinal Logical; reserved for future use.
#' @param binary_map Optional mapping as returned by `createBinaryMapping()`.
#'   If NULL, a mapping is created internally for `CatVars`.
#' @return A list with Unadjusted, FDRCorrected, method ("R_pb"), Relabel, Covariates, BinaryMapping.
#' @export
PlotPointCorrelationsHeatmap <- function (
    Data, CatVars, ContVars, Covariates = NULL, Relabel = TRUE, Ordinal = TRUE, binary_map = NULL
) {
  # ---- helpers ----
  .is_binary <- function(x) length(unique(stats::na.omit(x))) == 2

  .to01 <- function(x, var, map) {
    if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)
    pos <- map$PositiveLevel[match(var, map$Variable)]
    if (is.na(pos)) stop("No PositiveLevel found in mapping for variable: ", var)
    as.integer(as.character(x) == pos)
  }

  .label_or_name <- function(vec, nm) {
    lab <- sjlabelled::get_label(vec)
    if (is.null(lab) || is.na(lab) || lab == "") nm else lab
  }

  # ---- data prep ----
  DataSubset <- Data[c(CatVars, ContVars, if (!is.null(Covariates)) Covariates)]
  DataSubset <- ReplaceMissingLabels(DataSubset)

  # sanity: binaries truly binary
  nb_chk <- vapply(DataSubset[CatVars], .is_binary, logical(1))
  if (!all(nb_chk)) {
    bad <- names(nb_chk)[!nb_chk]
    stop("The following variables are not binary (exactly 2 unique non-NA values required): ",
         paste(bad, collapse = ", "))
  }

  # mapping (reuse if supplied)
  if (is.null(binary_map)) {
    binary_map <- createBinaryMapping(DataSubset, CatVars)
  }

  # long format
  mData <- DataSubset %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(CatVars),
                                ~ as.character(haven::zap_labels(.)))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(CatVars),
                        names_to = "CategoricalVariable",
                        values_to = "CategoricalValue") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(ContVars),
                        names_to = "ContinuousVariable",
                        values_to = "ContinuousValue") %>%
    dplyr::mutate(
      ContinuousVariable = as.factor(ContinuousVariable),
      CategoricalValue   = as.factor(CategoricalValue)
    ) %>%
    stats::na.omit()

  # drop degenerate pairs (single level / zero variance)
  keep_pairs <- mData %>%
    dplyr::group_by(ContinuousVariable, CategoricalVariable) %>%
    dplyr::summarise(
      has_two_levels = length(unique(CategoricalValue)) > 1,
      var_ok         = stats::var(ContinuousValue) > 0,
      .groups = "drop"
    )
  mData <- mData %>%
    dplyr::inner_join(keep_pairs %>% dplyr::filter(has_two_levels, var_ok),
                      by = c("ContinuousVariable","CategoricalVariable"))

  # ---- stats: r_pb & p ----
  stat.test <- mData %>%
    dplyr::group_by(ContinuousVariable, CategoricalVariable) %>%
    dplyr::group_modify(~{
      df <- .x
      cat_var <- as.character(.y$CategoricalVariable)
      x01 <- .to01(df$CategoricalValue, cat_var, binary_map)
      r <- suppressWarnings(stats::cor(df$ContinuousValue, x01, use = "complete.obs"))
      p <- suppressWarnings(stats::cor.test(df$ContinuousValue, x01)$p.value)
      tibble::tibble(
        nPairs = sum(stats::complete.cases(df$CategoricalValue, df$ContinuousValue)),
        correlation = r,
        p_value = p
      )
    }) %>%
    dplyr::ungroup()

  # FDR & stars
  stat.test$p.adj <- stats::p.adjust(stat.test$p_value, method = "fdr")
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

  stat.test$test <- "point biserial correlation"

  # labels for axes
  if (Relabel) {
    stat.test$XLabel <- vapply(
      stat.test$CategoricalVariable,
      function(v) .label_or_name(DataSubset[[v]], v), character(1)
    )
    stat.test$YLabel <- vapply(
      stat.test$ContinuousVariable,
      function(v) .label_or_name(DataSubset[[v]], v), character(1)
    )
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  # tooltip
  PlotText <- paste0(
    "</br>Cat Var Label: ", stat.test$XLabel,
    "</br> Cat Var: ",      stat.test$CategoricalVariable,
    "</br> Cont Var Label: ", stat.test$YLabel,
    "</br> Cont Var: ",       stat.test$ContinuousVariable,
    "</br> P-Value: ", signif(stat.test$p_value, 3), " ", stat.test$`p<.05`,
    "</br> FDR-corrected P: ", signif(stat.test$p.adj, 3), " ", stat.test$p.adj.signif,
    "</br> nPairs: ", stat.test$nPairs
  )

  # ---- plots ----
  p <- ggplot2::ggplot(
    stat.test,
    ggplot2::aes(y = YLabel, x = XLabel, fill = correlation, text = PlotText)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = `p<.05`), color = "black") +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1),
                                  low = scales::muted("#FFA500"),
                                  high = scales::muted("#008080")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggplot2::labs(fill = expression(r[pb]))

  p_FDR <- ggplot2::ggplot(
    stat.test,
    ggplot2::aes(y = YLabel, x = XLabel, fill = correlation, text = PlotText)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = p.adj.signif), color = "black") +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1),
                                  low = scales::muted("#FFA500"),
                                  high = scales::muted("#008080")) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    ggplot2::labs(fill = expression(r[pb]))

  # ---- return ----
  M     <- list(PvalTable = stat.test, plot = p)
  M_FDR <- list(PvalTable = stat.test, plot = p_FDR)

  return(list(
    Unadjusted    = M,
    FDRCorrected  = M_FDR,
    method        = "R_pb",
    Relabel       = Relabel,
    Covariates    = Covariates,
    BinaryMapping = binary_map
  ))
}
