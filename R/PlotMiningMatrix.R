#' PlotMiningMatrix
#'
#' Generate a matrix of statistical relationships between variables.
#'
#' @param Data A data frame.
#' @param OutcomeVars Outcome variables.
#' @param PredictorVars Predictor variables. If NULL, uses OutcomeVars.
#' @param Covariates Optional covariates (reserved for future use).
#' @param Relabel Use labels instead of names.
#' @param Parametric Use parametric tests.
#'
#' @return List with tables and plots.
#' @export
PlotMiningMatrix <- function(
    Data,
    OutcomeVars,
    PredictorVars = NULL,
    Covariates = NULL,
    Relabel = TRUE,
    Parametric = TRUE
) {

  method <- ifelse(Parametric, "pearson", "spearman")

  OutcomeVars <- unique(intersect(as.character(OutcomeVars), names(Data)))

  if (is.null(PredictorVars)) {
    PredictorVars <- OutcomeVars
  } else {
    PredictorVars <- unique(intersect(as.character(PredictorVars), names(Data)))
  }

  if (length(OutcomeVars) == 0 || length(PredictorVars) == 0) {
    empty <- data.frame()
    p0 <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = p0),
      FDRCorrected = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    ))
  }

  labels <- sjlabelled::get_label(Data, def.value = names(Data))
  names(labels) <- names(Data)

  num_vars <- SciDataReportR::getNumVars(Data)
  cat_vars <- SciDataReportR::getCatVars(Data)

  get_type <- function(v) {
    if (v %in% num_vars) return("num")
    if (v %in% cat_vars) return("cat")
    return("other")
  }

  pairs <- expand.grid(
    XVar = OutcomeVars,
    YVar = PredictorVars,
    stringsAsFactors = FALSE
  )

  pairs <- pairs[pairs$XVar != pairs$YVar, ]

  results <- purrr::map_dfr(seq_len(nrow(pairs)), function(i) {

    x <- pairs$XVar[i]
    y <- pairs$YVar[i]

    type_x <- get_type(x)
    type_y <- get_type(y)

    df <- Data %>%
      dplyr::select(dplyr::all_of(c(x, y))) %>%
      tidyr::drop_na()

    if (nrow(df) < 3) return(NULL)

    out <- NULL

    # numeric vs numeric
    if (type_x == "num" && type_y == "num") {

      if (stats::sd(df[[x]]) == 0 || stats::sd(df[[y]]) == 0) return(NULL)

      test <- suppressWarnings(stats::cor.test(df[[x]], df[[y]], method = method))

      out <- data.frame(
        XVar = x,
        YVar = y,
        p = test$p.value,
        EffectSize = unname(test$estimate),
        Test = method
      )
    }

    # numeric vs categorical
    if ((type_x == "num" && type_y == "cat") || (type_x == "cat" && type_y == "num")) {

      num_var <- ifelse(type_x == "num", x, y)
      cat_var <- ifelse(type_x == "cat", x, y)

      y_vec <- df[[num_var]]
      g_vec <- droplevels(as.factor(df[[cat_var]]))

      if (length(levels(g_vec)) < 2) return(NULL)
      if (stats::sd(y_vec) == 0) return(NULL)

      if (Parametric) {

        fit <- stats::lm(y_vec ~ g_vec)
        sm <- stats::anova(fit)

        pval <- sm[["Pr(>F)"]][1]
        ss <- sm[["Sum Sq"]]
        eta2 <- ss[1] / sum(ss)

        out <- data.frame(
          XVar = x,
          YVar = y,
          p = pval,
          EffectSize = eta2,
          Test = "ANOVA"
        )

      } else {

        test <- stats::kruskal.test(y_vec, g_vec)

        out <- data.frame(
          XVar = x,
          YVar = y,
          p = test$p.value,
          EffectSize = NA_real_,
          Test = "Kruskal"
        )
      }
    }

    # categorical vs categorical
    if (type_x == "cat" && type_y == "cat") {

      x_fac <- droplevels(as.factor(df[[x]]))
      y_fac <- droplevels(as.factor(df[[y]]))

      if (length(levels(x_fac)) < 2 || length(levels(y_fac)) < 2) return(NULL)

      tbl <- table(x_fac, y_fac)

      if (all(dim(tbl) > 1)) {

        test <- suppressWarnings(stats::chisq.test(tbl))

        n <- sum(tbl)
        cramer_v <- sqrt(test$statistic / (n * (min(dim(tbl)) - 1)))

        out <- data.frame(
          XVar = x,
          YVar = y,
          p = test$p.value,
          EffectSize = cramer_v,
          Test = "ChiSq"
        )
      }
    }

    out
  })

  if (nrow(results) == 0) {
    empty <- data.frame()
    p0 <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(
      Unadjusted = list(PvalTable = empty, plot = p0),
      FDRCorrected = list(PvalTable = empty, plot = p0),
      method = method, Relabel = Relabel, Covariates = Covariates
    ))
  }

  # labels
  if (Relabel) {
    results$XLabel <- labels[results$XVar]
    results$YLabel <- labels[results$YVar]
  } else {
    results$XLabel <- results$XVar
    results$YLabel <- results$YVar
  }

  # enforce ordering
  x_order <- if (Relabel) labels[OutcomeVars] else OutcomeVars
  y_order <- if (Relabel) labels[PredictorVars] else PredictorVars

  results$XLabel <- factor(results$XLabel, levels = x_order)
  results$YLabel <- factor(results$YLabel, levels = rev(y_order))

  # adjust p-values
  results$p_adj <- stats::p.adjust(results$p, method = "fdr")

  results <- results %>%
    rstatix::add_significance(p.col = "p", output.col = "stars") %>%
    rstatix::add_significance(p.col = "p_adj", output.col = "stars_fdr")

  # effect size for plotting
  results$EffectSizeAbs <- abs(results$EffectSize)

  # size mapping
  size_map <- c("ns" = 2, "*" = 3, "**" = 4, "***" = 5)
  results$size_val <- unname(size_map[as.character(results$stars)])
  results$size_val_fdr <- unname(size_map[as.character(results$stars_fdr)])

  # shape mapping
  shape_map <- c("ns" = 16, "*" = 17, "**" = 15, "***" = 18)

  col_scale <- ggplot2::scale_colour_gradientn(
    colours = as.character(paletteer::paletteer_c("grDevices::Purple-Yellow", n = 100, direction = -1)),
    limits = c(0, 1),
    name = "Effect Size"
  )

  p <- ggplot2::ggplot(results, ggplot2::aes(x = XLabel, y = YLabel)) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = EffectSizeAbs,
        size = size_val,
        shape = stars
      )
    ) +
    col_scale +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE, name = "Significance") +
    ggplot2::labs(subtitle = "No Multiple Comparison Correction") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank()
    )

  p_fdr <- ggplot2::ggplot(results, ggplot2::aes(x = XLabel, y = YLabel)) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = EffectSizeAbs,
        size = size_val_fdr,
        shape = stars_fdr
      )
    ) +
    col_scale +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_shape_manual(values = shape_map, drop = FALSE, name = "Significance") +
    ggplot2::labs(subtitle = "FDR Correction") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_blank()
    )

  list(
    Unadjusted = list(PvalTable = results, plot = p),
    FDRCorrected = list(PvalTable = results, plot = p_fdr),
    method = method,
    Relabel = Relabel,
    Covariates = Covariates
  )
}
