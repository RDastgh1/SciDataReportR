#' Plot ANOVA Relationships Matrix
#'
#' This function plots the relationship between continuous and categorical variables
#' using ANOVA or Kruskal-Wallis tests. It generates a "heatmap" with points
#' colored and shaped based on statistical significance and effect size.
#'
#' @param Data The data frame containing the variables of interest.
#' @param CatVars Character vector of categorical variable names.
#' @param ContVars Character vector of continuous variable names.
#' @param Covariates Optional character vector of covariate names for ANCOVA analysis.
#' @param Relabel Logical indicating whether to relabel variables with their labels (default is TRUE).
#' @param Ordinal Logical, indicating whether ordinal variables should be considered.
#' @param Parametric Logical indicating whether to use parametric (ANOVA) or non-parametric (Kruskal-Wallis) tests (default is TRUE).
#' @return A list containing three ggplot objects: p (scatter plot without multiple comparison correction), p_FDR (scatter plot with FDR correction), and pvaltable (data frame of p-values and significance).
#' @import dplyr
#' @importFrom ggplot2 aes geom_point labs scale_shape_manual scale_color_gradientn guides theme
#' @importFrom sjlabelled get_label
#' @importFrom tidyr pivot_longer
#' @importFrom rstatix anova_test kruskal_test get_summary_stats add_significance adjust_pvalue
#' @importFrom stats var
#' @importFrom utils na.omit
#' @importFrom rstatix add_significance adjust_pvalue anova_test kruskal_test get_summary_stats
#' @export
PlotAnovaRelationshipsMatrix <- function(
    Data,
    CatVars,
    ContVars,
    Covariates = NULL,
    Relabel = TRUE,
    Parametric = TRUE,
    Ordinal = FALSE,
    min_n = 4,
    eps = 1e-8
) {

  method <- if (isTRUE(Parametric)) "Anova" else "Kruskal"

  CatVars <- unique(intersect(as.character(CatVars), names(Data)))
  ContVars <- unique(intersect(as.character(ContVars), names(Data)))
  Covariates <- if (!is.null(Covariates) && length(Covariates) > 0) unique(intersect(as.character(Covariates), names(Data))) else character(0)

  empty_return <- function() {
    empty <- data.frame()
    p0 <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_blank() + ggplot2::theme_void()
    list(
      Unadjusted = list(PvalTable = empty, plot = p0),
      FDRCorrected = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    )
  }

  if (length(CatVars) == 0 || length(ContVars) == 0) return(empty_return())

  # drop categorical vars with < 2 observed levels
  CatVars <- CatVars[vapply(CatVars, function(v) {
    x <- Data[[v]]
    x <- x[!is.na(x)]
    length(unique(x)) >= 2
  }, logical(1))]
  if (length(CatVars) == 0) return(empty_return())

  # continuous vars: must be numeric or coercible with enough non-NA + variance
  ContVars <- ContVars[vapply(ContVars, function(v) {
    x <- Data[[v]]
    if (!(is.numeric(x) || is.integer(x))) x <- suppressWarnings(as.numeric(as.character(x)))
    x <- x[is.finite(x)]
    if (length(x) < min_n) return(FALSE)
    vv <- stats::var(x, na.rm = TRUE)
    is.finite(vv) && vv > eps
  }, logical(1))]
  if (length(ContVars) == 0) return(empty_return())

  # covars: keep only non-constant and usable
  if (length(Covariates) > 0) {
    Covariates <- Covariates[vapply(Covariates, function(v) {
      z <- Data[[v]]
      z <- z[!is.na(z)]
      if (length(z) < min_n) return(FALSE)
      if (is.factor(z) || is.character(z)) return(length(unique(z)) >= 2)
      z2 <- suppressWarnings(as.numeric(as.character(z)))
      z2 <- z2[is.finite(z2)]
      if (length(z2) < min_n) return(FALSE)
      vv <- stats::var(z2, na.rm = TRUE)
      is.finite(vv) && vv > eps
    }, logical(1))]
  }

  vars_keep <- unique(c(CatVars, ContVars, Covariates))

  df0 <- Data |>
    dplyr::select(dplyr::all_of(vars_keep)) |>
    dplyr::mutate(.rowID = dplyr::row_number()) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(CatVars), ~ as.factor(.x)))

  # Keep covars attached by rowID (avoid pivot mixing types)
  cov_df <- if (length(Covariates) > 0) {
    df0 |> dplyr::select(.rowID, dplyr::all_of(Covariates))
  } else {
    df0 |> dplyr::select(.rowID)
  }

  mData <- df0 |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(CatVars),
      names_to = "CategoricalVariable",
      values_to = "CategoricalValue"
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(ContVars),
      names_to = "ContinuousVariable",
      values_to = "ContinuousValue"
    ) |>
    dplyr::select(-dplyr::any_of(Covariates)) |>
    dplyr::left_join(cov_df, by = ".rowID") |>
    dplyr::mutate(
      CategoricalValue = as.factor(CategoricalValue),
      ContinuousVariable = as.factor(ContinuousVariable),
      ContinuousValue = suppressWarnings(as.numeric(as.character(ContinuousValue)))
    ) |>
    tidyr::drop_na(CategoricalValue, ContinuousValue)

  if (length(Covariates) > 0) {
    mData <- tidyr::drop_na(mData, dplyr::all_of(Covariates))
  }

  if (nrow(mData) == 0) return(empty_return())

  # pair-level guards
  mData <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
    dplyr::filter(
      dplyr::n() >= min_n,
      dplyr::n_distinct(CategoricalValue) > 1,
      dplyr::n_distinct(ContinuousValue) > 2,
      stats::sd(ContinuousValue, na.rm = TRUE) > eps
    ) |>
    dplyr::ungroup()

  if (nrow(mData) == 0) return(empty_return())

  rhs <- "CategoricalValue"
  if (length(Covariates) > 0) rhs <- paste(rhs, paste(Covariates, collapse = " + "), sep = " + ")
  fml <- stats::as.formula(paste("ContinuousValue ~", rhs))

  safe_test_df <- function(df) {
    fit <- tryCatch(stats::lm(fml, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    rss <- sum(stats::resid(fit)^2)
    if (!is.finite(rss) || rss <= eps) return(NULL)

    out <- tryCatch({
      if (isTRUE(Parametric)) rstatix::anova_test(df, fml) else rstatix::kruskal_test(df, fml)
    }, error = function(e) NULL)

    if (is.null(out)) return(NULL)
    as.data.frame(out)  # strip classes
  }

  stat.test <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
    dplyr::group_modify(~{
      res <- safe_test_df(.x)
      if (is.null(res)) return(tibble::tibble())
      tibble::as_tibble(res)
    }) |>
    dplyr::ungroup()

  if (nrow(stat.test) == 0) return(empty_return())

  if (!("p" %in% names(stat.test)) && ("p.value" %in% names(stat.test))) {
    stat.test <- dplyr::rename(stat.test, p = p.value)
  }

  stat.test <- stat.test |>
    rstatix::add_significance(p.col = "p", output.col = "p.signif") |>
    rstatix::adjust_pvalue(p.col = "p", output.col = "p.adj", method = "fdr") |>
    rstatix::add_significance(p.col = "p.adj", output.col = "p.adj.signif")

  stat.test$logp <- -log10(stat.test$p)
  stat.test$logp_FDR <- -log10(stat.test$p.adj)
  stat.test$`p<.05` <- factor(stat.test$p.signif, levels = c("ns","*","**","***","****"))
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif, levels = c("ns","*","**","***","****"))

  # Keep only main effect row if present
  if ("Effect" %in% names(stat.test)) {
    stat.test <- stat.test |>
      dplyr::filter(Effect == "CategoricalValue")
  }

  # labels
  if (isTRUE(Relabel)) {
    Data <- ReplaceMissingLabels(Data)

    xlabels <- sjlabelled::get_label(Data[as.character(stat.test$CategoricalVariable)],
                                     def.value = stat.test$CategoricalVariable) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable", "label")

    ylabels <- sjlabelled::get_label(Data[as.character(stat.test$ContinuousVariable)],
                                     def.value = stat.test$ContinuousVariable) |>
      as.data.frame() |>
      tibble::rownames_to_column()
    colnames(ylabels) <- c("Variable", "label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  # ordering
  if (isTRUE(Relabel)) {
    stat.test$XLabel <- factor(stat.test$XLabel,
                               levels = sjlabelled::get_label(Data[as.character(CatVars)], def.value = CatVars))
    stat.test$YLabel <- factor(stat.test$YLabel,
                               levels = sjlabelled::get_label(Data[as.character(ContVars)], def.value = ContVars))
  } else {
    stat.test$XLabel <- factor(stat.test$XLabel, levels = CatVars)
    stat.test$YLabel <- factor(stat.test$YLabel, levels = ContVars)
  }

  stat.test$PlotText <- paste(
    "</br>Cat Var Label:", stat.test$XLabel,
    "</br>Cat Var:", stat.test$CategoricalVariable,
    "</br>Cont Var Label:", stat.test$YLabel,
    "</br>Cont Var:", stat.test$ContinuousVariable,
    "</br>P-Value:", stat.test$p, stat.test$p.signif,
    "</br>FDR P:", stat.test$p.adj, stat.test$p.adj.signif,
    if ("ges" %in% names(stat.test)) paste0("</br>GES:", stat.test$ges) else ""
  )

  stat.test <- stat.test |>
    dplyr::mutate(Test = method) |>
    as.data.frame()

  # plots: use ges if available, else color by p
  colour_var <- if ("ges" %in% names(stat.test)) "ges" else "p"

  p <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = `p<.05`, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = `p<.05`, colour = .data[[colour_var]])) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  p_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = p.adj.signif, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = p.adj.signif, colour = .data[[colour_var]])) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  list(
    Unadjusted = list(PvalTable = stat.test, plot = p),
    FDRCorrected = list(PvalTable = stat.test, plot = p_FDR),
    method = method, Relabel = Relabel, Covariates = Covariates
  )
}
