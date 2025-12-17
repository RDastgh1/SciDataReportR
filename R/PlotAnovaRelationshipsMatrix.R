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

  # ---- 0) Input hygiene ----
  CatVars <- unique(as.character(CatVars))
  ContVars <- unique(as.character(ContVars))
  Covariates <- if (!is.null(Covariates) && length(Covariates) > 0) unique(as.character(Covariates)) else character(0)

  # Keep only variables that exist
  CatVars <- intersect(CatVars, names(Data))
  ContVars <- intersect(ContVars, names(Data))
  Covariates <- intersect(Covariates, names(Data))

  # Drop categorical vars with <2 observed levels
  if (length(CatVars) > 0) {
    CatVars <- CatVars[vapply(CatVars, function(v) {
      x <- Data[[v]]
      x <- x[!is.na(x)]
      length(unique(x)) >= 2
    }, logical(1))]
  }

  # Drop continuous vars that are not numeric/integer OR have near-zero variance
  if (length(ContVars) > 0) {
    ContVars <- ContVars[vapply(ContVars, function(v) {
      x <- Data[[v]]
      if (!(is.numeric(x) || is.integer(x))) return(FALSE)
      x <- x[is.finite(x)]
      if (length(x) < min_n) return(FALSE)
      vv <- stats::var(x, na.rm = TRUE)
      is.finite(vv) && vv > eps
    }, logical(1))]
  }

  # Drop covariates that are missing / all NA / constant
  if (length(Covariates) > 0) {
    Covariates <- Covariates[vapply(Covariates, function(v) {
      x <- Data[[v]]
      x <- x[!is.na(x)]
      if (length(x) < min_n) return(FALSE)
      if (is.factor(x) || is.character(x)) {
        length(unique(x)) >= 2
      } else if (is.numeric(x) || is.integer(x)) {
        vv <- stats::var(x, na.rm = TRUE)
        is.finite(vv) && vv > eps
      } else {
        TRUE
      }
    }, logical(1))]
  }

  method <- if (Parametric) "Anova" else "Kruskal"

  # If nothing left, return safely
  if (length(CatVars) == 0 || length(ContVars) == 0) {
    empty <- data.frame()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = ggplot2::ggplot()),
      FDRCorrected = list(PvalTable = empty, plot = ggplot2::ggplot()),
      method = method, Relabel = Relabel, Covariates = Covariates
    ))
  }

  # ---- 1) Build stable long data with covariates preserved ----
  vars_keep <- unique(c(CatVars, ContVars, Covariates))

  df0 <- Data |>
    dplyr::select(dplyr::all_of(vars_keep)) |>
    dplyr::mutate(.rowID = dplyr::row_number()) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(CatVars), ~ as.factor(.x)))

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

  # ---- 2) Pair-level guards (prevents df=0 and RSS≈0 issues) ----
  mData <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
    dplyr::filter(
      dplyr::n() >= min_n,
      dplyr::n_distinct(CategoricalValue) > 1,
      dplyr::n_distinct(ContinuousValue) > 2,
      stats::sd(ContinuousValue, na.rm = TRUE) > eps
    ) |>
    dplyr::ungroup()

  if (nrow(mData) == 0) {
    empty <- data.frame()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = ggplot2::ggplot()),
      FDRCorrected = list(PvalTable = empty, plot = ggplot2::ggplot()),
      method = method, Relabel = Relabel, Covariates = Covariates
    ))
  }

  # Model formula
  rhs <- "CategoricalValue"
  if (length(Covariates) > 0) rhs <- paste(rhs, paste(Covariates, collapse = " + "), sep = " + ")
  fml <- stats::as.formula(paste("ContinuousValue ~", rhs))

  # ---- 3) Safe per-group test: skip models that error or have ~0 RSS ----
  safe_test <- function(df) {
    # Fit lm first so we can detect perfect-fit / constant residuals before Anova()
    fit <- tryCatch(stats::lm(fml, data = df), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    rss <- sum(stats::resid(fit)^2)
    if (!is.finite(rss) || rss <= eps) return(NULL)

    out <- tryCatch({
      if (Parametric) rstatix::anova_test(df, fml) else rstatix::kruskal_test(df, fml)
    }, error = function(e) NULL)

    out
  }

  stat.test <- mData |>
    dplyr::group_by(ContinuousVariable, CategoricalVariable) |>
    rstatix::doo(~ safe_test(.x)) |>
    tidyr::unnest(cols = c(.results))

  if (nrow(stat.test) == 0) {
    empty <- data.frame()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = ggplot2::ggplot()),
      FDRCorrected = list(PvalTable = empty, plot = ggplot2::ggplot()),
      method = method, Relabel = Relabel, Covariates = Covariates
    ))
  }

  # ---- 4) Significance / FDR ----
  stat.test <- stat.test |>
    rstatix::add_significance() |>
    rstatix::adjust_pvalue(method = "fdr") |>
    rstatix::add_significance()

  # ---- 5) Keep your original downstream expectations ----
  stat.test$logp <- -log10(stat.test$p)
  stat.test$`p<.05` <- factor(stat.test$p.signif, levels = c("ns","*","**","***","****"))
  stat.test$logp_FDR <- -log10(stat.test$p.adj)
  stat.test$p.adj.signif <- factor(stat.test$p.adj.signif, levels = c("ns","*","**","***","****"))

  if (Relabel) {
    Data <- ReplaceMissingLabels(Data)

    xlabels <- sjlabelled::get_label(
      Data[as.character(stat.test$CategoricalVariable)],
      def.value = stat.test$CategoricalVariable
    ) |> as.data.frame() |> tibble::rownames_to_column()
    colnames(xlabels) <- c("Variable","label")

    ylabels <- sjlabelled::get_label(
      Data[as.character(stat.test$ContinuousVariable)],
      def.value = stat.test$ContinuousVariable
    ) |> as.data.frame() |> tibble::rownames_to_column()
    colnames(ylabels) <- c("Variable","label")

    stat.test$XLabel <- xlabels$label
    stat.test$YLabel <- ylabels$label
  } else {
    stat.test$XLabel <- stat.test$CategoricalVariable
    stat.test$YLabel <- stat.test$ContinuousVariable
  }

  stat.test$PlotText <- paste(
    "</br>Cat Var Label:", stat.test$XLabel,
    "</br> Cat Var:", stat.test$CategoricalVariable,
    "</br> Cont Var Label:", stat.test$YLabel,
    "</br> Cont Var:", stat.test$ContinuousVariable,
    "</br> P-Value: ", stat.test$p, stat.test$p.signif,
    "</br> FDR-corrected P: ", stat.test$p.adj, stat.test$p.adj.signif,
    "</br> GES Effect size: ", stat.test$ges,
    "</br> npairs: ", stat.test$DFd
  )

  stat.test <- stat.test |>
    dplyr::ungroup() |>
    as.data.frame() |>
    dplyr::filter(Effect == "CategoricalValue") |>
    dplyr::mutate(Test = method)

  # Plot (kept close to your original)
  p <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = `p<.05`, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = `p<.05`, colour = ges)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  p_FDR <- ggplot2::ggplot(stat.test, ggplot2::aes(y = YLabel, x = XLabel, shape = p.adj.signif, text = PlotText)) +
    ggplot2::geom_point(ggplot2::aes(size = p.adj.signif, colour = ges)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_shape_manual(values = c(7, 16, 17, 15, 18), drop = FALSE) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::xlab("") + ggplot2::ylab("")

  M <- list(PvalTable = stat.test, plot = p)
  M_FDR <- list(PvalTable = stat.test, plot = p_FDR)

  list(Unadjusted = M, FDRCorrected = M_FDR, method = method, Relabel = Relabel, Covariates = Covariates)
}
