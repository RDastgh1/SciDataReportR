#' Plot Point-Biserial Correlations (Binary vs Continuous)
#'
#' Computes point-biserial correlations between binary categorical variables
#' and continuous variables using explicit 0/1 coding where
#' \code{1 == PositiveLevel} from \code{createBinaryMapping()}.
#' Produces tidy tables and a dot-plot visualization.
#'
#' @param Data A dataframe.
#' @param CatVars Character vector of binary categorical variables.
#' @param ContVars Character vector of continuous variables.
#' @param Covariates Optional covariates (reserved; not used in the current implementation).
#' @param Relabel Logical; if TRUE, uses sjlabelled variable labels for axes.
#' @param Ordinal Logical; reserved for future use.
#' @param binary_map Optional mapping as returned by \code{createBinaryMapping()}.
#'   If NULL, a mapping is created internally for \code{CatVars}.
#' @param sort_by One of \code{"abs"}, \code{"value"}, or \code{"p"} to control y-axis ordering.
#' @return A list with:
#' \itemize{
#'   \item \code{Unadjusted}: list(PvalTable, plot)
#'   \item \code{FDRCorrected}: list(PvalTable, plot)
#'   \item \code{method} = "R_pb"
#'   \item \code{Relabel}
#'   \item \code{Covariates}
#'   \item \code{BinaryMapping} (used)
#' }
#' @export
PlotPointCorrelations <- function(
    Data, CatVars, ContVars, Covariates = NULL, Relabel = TRUE, Ordinal = TRUE,
    binary_map = NULL, sort_by = c("abs","value","p")
) {

  sort_by <- match.arg(sort_by)

  # ---- helpers -------------------------------------------------------------
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

  # ---- data prep -----------------------------------------------------------
  DataSubset <- Data[c(CatVars, ContVars, if (!is.null(Covariates)) Covariates)]
  DataSubset <- ReplaceMissingLabels(DataSubset)

  # sanity: binaries are truly binary
  nb_chk <- vapply(DataSubset[CatVars], .is_binary, logical(1))
  if (!all(nb_chk)) {
    bad <- names(nb_chk)[!nb_chk]
    stop("The following variables are not binary (exactly 2 unique non-NA values required): ",
         paste(bad, collapse = ", "))
  }

  # mapping
  if (is.null(binary_map)) {
    binary_map <- createBinaryMapping(DataSubset, CatVars)
  }

  # long format (pair every Cat with every Cont)
  long_df <- DataSubset %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(CatVars),
                                ~ as.character(haven::zap_labels(.)))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(CatVars),
                        names_to = "CategoricalVariable",
                        values_to = "CategoricalValue") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(ContVars),
                        names_to = "ContinuousVariable",
                        values_to = "ContinuousValue") %>%
    dplyr::mutate(
      CategoricalValue   = as.factor(CategoricalValue)
    ) %>%
    stats::na.omit()

  # drop degenerate rows (single level or zero variance)
  keep_pairs <- long_df %>%
    dplyr::group_by(ContinuousVariable, CategoricalVariable) %>%
    dplyr::summarise(
      has_two_levels = length(unique(CategoricalValue)) > 1,
      var_ok         = stats::var(ContinuousValue) > 0,
      .groups = "drop"
    )
  long_df <- long_df %>%
    dplyr::inner_join(keep_pairs %>% dplyr::filter(has_two_levels, var_ok),
                      by = c("ContinuousVariable","CategoricalVariable"))

  # ---- stats: r_pb & p for each pair --------------------------------------
  stat.test <- long_df %>%
    dplyr::group_by(ContinuousVariable, CategoricalVariable) %>%
    dplyr::group_modify(~{
      df <- .x
      cat_var <- as.character(.y$CategoricalVariable)
      x01 <- .to01(df$CategoricalValue, cat_var, binary_map)
      # r and p
      r <- suppressWarnings(stats::cor(df$ContinuousValue, x01, use = "complete.obs"))
      p <- suppressWarnings(stats::cor.test(df$ContinuousValue, x01)$p.value)
      tibble::tibble(
        nPairs = sum(stats::complete.cases(df$CategoricalValue, df$ContinuousValue)),
        correlation = r,
        p_value = p
      )
    }) %>%
    dplyr::ungroup()

  # FDR
  stat.test$p.adj <- stats::p.adjust(stat.test$p_value, method = "fdr")

  # stars
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

  # labels
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

  # order Y per facet using requested sorting
  stat.test <- dplyr::group_by(stat.test, CategoricalVariable, XLabel)
  stat.test <- dplyr::mutate(
    stat.test,
    order_key = dplyr::case_when(
      sort_by == "abs"   ~ -abs(correlation),
      sort_by == "value" ~ -correlation,
      sort_by == "p"     ~ stat.test$p_value,
      TRUE               ~ -abs(correlation)
    )
  )
  stat.test <- dplyr::ungroup(stat.test)

  # per-facet y ordering
  stat.test <- stat.test %>%
    dplyr::group_by(CategoricalVariable, XLabel) %>%
    dplyr::arrange(order_key, .by_group = TRUE) %>%
    dplyr::mutate(YLabel = factor(YLabel, levels = unique(YLabel))) %>%
    dplyr::ungroup()

  # tooltip
  PlotText <- paste0(
    "</br>Binary (facet): ", stat.test$XLabel,
    "</br>Binary var: ",     stat.test$CategoricalVariable,
    "</br>Continuous: ",     stat.test$YLabel,
    "</br>Cont var: ",       stat.test$ContinuousVariable,
    "</br> r_pb: ",          signif(stat.test$correlation, 3),
    "</br> p: ",             signif(stat.test$p_value, 3), " ", stat.test$`p<.05`,
    "</br> FDR p: ",         signif(stat.test$p.adj, 3),   " ", stat.test$p.adj.signif,
    "</br> nPairs: ",        stat.test$nPairs
  )

  # ---- plots: dot plot faceted by binary var -------------------------------
  base_dot <- ggplot2::ggplot(
    stat.test,
    ggplot2::aes(x = correlation, y = YLabel, text = PlotText)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.text  = ggplot2::element_text(size = 15)
    ) +
    ggplot2::labs(x = expression(r[pb]))

  p <- base_dot +
    ggplot2::geom_text(ggplot2::aes(label = `p<.05`), nudge_x = 0.02, color = "black") +
    ggplot2::facet_wrap(~ XLabel, scales = "free_y")

  p_FDR <- base_dot +
    ggplot2::geom_text(ggplot2::aes(label = p.adj.signif), nudge_x = 0.02, color = "black") +
    ggplot2::facet_wrap(~ XLabel, scales = "free_y")

  # ---- return --------------------------------------------------------------
  M     <- list(PvalTable = stat.test, plot = p)
  M_FDR <- list(PvalTable = stat.test, plot = p_FDR)

  list(
    Unadjusted    = M,
    FDRCorrected  = M_FDR,
    method        = "R_pb",
    Relabel       = Relabel,
    Covariates    = Covariates,
    BinaryMapping = binary_map
  )
}
