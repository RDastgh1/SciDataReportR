#' Plot correlations within two strata as split triangles
#'
#' Each cell represents the same predictor-outcome pair in two strata of long
#' data. The first selected stratum occupies the top-left triangle and the
#' second occupies the bottom-right triangle. Correlations are calculated by
#' [PlotCorrelationsHeatmap()] independently within each stratum, including
#' covariate adjustment and FDR correction.
#'
#' @param data A long data frame with one column per measured variable and a
#'   column identifying the stratum of each record.
#' @param predictor_vars,outcome_vars Character vectors passed to
#'   [PlotCorrelationsHeatmap()]. Predictors are rows and outcomes are columns.
#'   Defaults use numeric predictors (excluding `split_var`) and the same
#'   variables as outcomes.
#' @param split_var Name of the column identifying strata.
#' @param split_levels Two observed values of `split_var`, in triangle order.
#'   If omitted, exactly two observed values are required; factor order is
#'   respected, otherwise first-observed order is used.
#' @param covariates,method,Relabel,TreatOrdinalAs,min_n,eps,fdr_scope,IncludeTrendSymbol
#'   Passed to [PlotCorrelationsHeatmap()]. FDR correction is performed
#'   separately within each stratum, before any display alignment.
#' @param interactive Either `"none"` (default) or `"plotly"`. Static plots
#'   are always returned; Plotly widgets require the suggested plotly package.
#'
#' @details
#' Stars test the correlation within each stratum, not the difference between
#' strata. This function performs no between-stratum comparison test. Missing
#' split values are excluded with a warning. Unavailable correlations remain
#' visible as gray triangles, with counts and status in `Results`.
#' `Unadjusted` refers to raw p-values; requested covariate adjustment is
#' retained in both output branches.
#'
#' Labels are inherited from the source correlation results (normally attached
#' from a codebook using [RevalueData()]). Variable identifiers, rather than
#' labels, determine cell identity, so duplicate display labels are supported.
#'
#' @return A list with `Unadjusted` and `FDRCorrected`, each containing a
#'   static `plot` and named lists of matrices `r`, `p`, and `npairs`, indexed
#'   by split level. `p` and `p_fdr` alias these branches. `Correlations` holds
#'   the source results; `Results` contains one row per stratum and cell;
#'   `Metadata` records settings, axis order, counts, and unavailable results.
#'   Optional widgets are in `Interactive$Plotly$Unadjusted` and
#'   `Interactive$Plotly$FDRCorrected`.
#' @author Aparna Bhattacharyya (original split-triangle wrapper).
#' @examples
#' set.seed(42)
#' df_Example <- data.frame(
#'   Compartment = factor(rep(c("Plasma", "CSF"), each = 30)),
#'   IL6 = rnorm(60), TNF = rnorm(60), Cognition = rnorm(60),
#'   Age = runif(60, 25, 75)
#' )
#' df_Example$IL6 <- sjlabelled::set_label(df_Example$IL6, "Interleukin-6")
#' df_Example$TNF <- sjlabelled::set_label(df_Example$TNF, "Tumor necrosis factor")
#' df_Example$Cognition <- sjlabelled::set_label(
#'   df_Example$Cognition, "Cognitive score"
#' )
#' result <- PlotSplitCorrelationHeatmap(
#'   df_Example, predictor_vars = c("IL6", "TNF"),
#'   outcome_vars = "Cognition", split_var = "Compartment",
#'   split_levels = c("Plasma", "CSF"), covariates = "Age"
#' )
#' result$FDRCorrected$plot
#' # Optional interactive version:
#' # result <- PlotSplitCorrelationHeatmap(
#' #   df_Example, c("IL6", "TNF"), "Cognition", "Compartment",
#' #   covariates = "Age", interactive = "plotly"
#' # )
#' # result$Interactive$Plotly$FDRCorrected
#' @details
#' Categorical input safeguards apply to `split_var`.
#' Numeric columns are converted locally to factors, with a warning, only when
#' all observed values are finite whole-number codes and there are at most 10
#' distinct non-missing values. Other numeric inputs error: explicitly convert
#' the column to a factor first if categorical treatment is intended (for
#' example, `df_All$cohort <- factor(df_All$cohort)`). Few unique values alone
#' do not establish categorical meaning. Factors, character columns, and logical
#' columns are accepted; more than 10 observed categories warn and continue.
#' Entirely missing or unsupported columns error. Caller data are unchanged,
#' and existing restrictions on the required number of groups still apply.
#' @export
PlotSplitCorrelationHeatmap <- function(
    data,
    predictor_vars = NULL,
    outcome_vars = NULL,
    split_var,
    split_levels = NULL,
    covariates = NULL,
    method = "pearson",
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    min_n = 3,
    eps = 1e-12,
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    IncludeTrendSymbol = FALSE,
    interactive = c("none", "plotly")) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(split_var) || length(split_var) != 1L ||
      is.na(split_var) || !split_var %in% names(data)) {
    stop("`split_var` must name one column in `data`.", call. = FALSE)
  }
  if (!is.atomic(data[[split_var]]) || !is.null(dim(data[[split_var]]))) {
    stop("`split_var` must be an atomic vector.", call. = FALSE)
  }
  data <- ScidrPrepareCategorical(data, split_var, "split_var")
  interactive <- match.arg(interactive)
  fdr_scope <- match.arg(fdr_scope)
  if (interactive == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    stop("Install the suggested package 'plotly' for interactive output.", call. = FALSE)
  }

  split_values <- as.character(data[[split_var]])
  missing_split <- is.na(data[[split_var]])
  observed <- unique(split_values[!missing_split])
  if (is.factor(data[[split_var]])) {
    observed <- levels(data[[split_var]])[levels(data[[split_var]]) %in% observed]
  }
  if (is.null(split_levels)) {
    if (length(observed) != 2L) {
      stop("Specify `split_levels` when `split_var` does not have exactly two observed levels.",
        call. = FALSE)
    }
    split_levels <- observed
  }
  split_levels <- as.character(split_levels)
  if (length(split_levels) != 2L || anyNA(split_levels) ||
      anyDuplicated(split_levels) || !all(split_levels %in% observed)) {
    stop("`split_levels` must select two distinct observed levels with records.", call. = FALSE)
  }
  if (any(missing_split)) {
    warning(sum(missing_split), " rows with missing split values were excluded.", call. = FALSE)
  }

  # Select variables once so a numeric stratum column never becomes a default
  # predictor and both source calls receive the same variable specification.
  if (is.null(predictor_vars)) {
    predictor_vars <- setdiff(
      getNumVars(data, Ordinal = identical(TreatOrdinalAs, "Continuous")),
      split_var
    )
  }
  display_labels <- ScidrDisplayLabels(data, names(data), Relabel)
  correlations <- stats::setNames(lapply(split_levels, function(split_level) {
    df_Stratum <- data[!missing_split & split_values %in% split_level, , drop = FALSE]
    # Base-vector row subsetting can drop label attributes. Restore the
    # original codebook labels before the correlation function reads them.
    for (variable in names(df_Stratum)) {
      df_Stratum[[variable]] <- sjlabelled::set_label(
        df_Stratum[[variable]], display_labels[[variable]]
      )
    }
    PlotCorrelationsHeatmap(
      data = df_Stratum, predictor_vars = predictor_vars,
      outcome_vars = outcome_vars, covariates = covariates, method = method,
      Relabel = Relabel, TreatOrdinalAs = TreatOrdinalAs, min_n = min_n,
      eps = eps, fdr_scope = fdr_scope, IncludeTrendSymbol = IncludeTrendSymbol
    )
  }), split_levels)

  vars_Rows <- unique(unlist(lapply(correlations, function(result) {
    rownames(result$Unadjusted$r)
  }), use.names = FALSE))
  vars_Columns <- unique(unlist(lapply(correlations, function(result) {
    colnames(result$Unadjusted$r)
  }), use.names = FALSE))
  if (!length(vars_Rows) || !length(vars_Columns)) {
    stop("No valid correlation variables remain after variable handling.", call. = FALSE)
  }

  # Source calls use a full display, so plot data retains labels for every cell.
  df_Labels <- dplyr::bind_rows(lapply(correlations, function(result) {
    result$Unadjusted$plot$data
  }))
  row_labels <- stats::setNames(df_Labels$XLabel[match(vars_Rows, df_Labels$XVar)], vars_Rows)
  column_labels <- stats::setNames(df_Labels$YLabel[match(vars_Columns, df_Labels$YVar)], vars_Columns)

  df_Results <- dplyr::bind_rows(lapply(seq_along(split_levels), function(split_index) {
    result <- correlations[[split_index]]
    df_Cells <- expand.grid(
      XVar = vars_Rows, YVar = vars_Columns, stringsAsFactors = FALSE
    )
    matrix_index <- cbind(
      match(df_Cells$XVar, rownames(result$Unadjusted$r)),
      match(df_Cells$YVar, colnames(result$Unadjusted$r))
    )
    df_Cells$Split <- split_levels[split_index]
    df_Cells$Triangle <- c("top-left", "bottom-right")[split_index]
    df_Cells$XLabel <- unname(row_labels[df_Cells$XVar])
    df_Cells$YLabel <- unname(column_labels[df_Cells$YVar])
    df_Cells$R <- result$Unadjusted$r[matrix_index]
    df_Cells$P <- result$Unadjusted$p[matrix_index]
    df_Cells$P_adj <- result$FDRCorrected$p[matrix_index]
    df_Cells$nPairs <- result$Unadjusted$npairs[matrix_index]
    df_Cells$Status <- ifelse(is.finite(df_Cells$R), "available", "unavailable")
    df_Cells$stars <- .FormatPValueStars(df_Cells$P, IncludeTrendSymbol = IncludeTrendSymbol)
    df_Cells$stars_FDR <- .FormatPValueStars(df_Cells$P_adj, IncludeTrendSymbol = IncludeTrendSymbol)
    df_Cells$CellId <- match(df_Cells$XVar, vars_Rows) +
      (match(df_Cells$YVar, vars_Columns) - 1L) * length(vars_Rows)
    df_Cells$x <- match(df_Cells$YVar, vars_Columns)
    df_Cells$y <- length(vars_Rows) + 1L - match(df_Cells$XVar, vars_Rows)
    df_Cells
  }))

  EscapeHtml <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  FormatValue <- function(x) {
    ifelse(is.finite(x), formatC(x, digits = 3, format = "g"), "Unavailable")
  }
  split_label <- attr(data[[split_var]], "label", exact = TRUE)
  if (!isTRUE(Relabel) || is.null(split_label) || !nzchar(split_label)) split_label <- split_var
  df_Results$Tooltip <- paste0(
    "<b>", EscapeHtml(split_label), ":</b> ", EscapeHtml(df_Results$Split),
    "<br><b>Predictor:</b> ", EscapeHtml(df_Results$XLabel),
    "<br><b>Outcome:</b> ", EscapeHtml(df_Results$YLabel),
    "<br><b>Correlation:</b> ", FormatValue(df_Results$R),
    "<br><b>p-value:</b> ", FormatValue(df_Results$P),
    "<br><b>FDR p-value:</b> ", FormatValue(df_Results$P_adj),
    "<br><b>Complete pairs:</b> ", df_Results$nPairs,
    "<br><b>Status:</b> ", df_Results$Status
  )

  df_Polygons <- df_Results[rep(seq_len(nrow(df_Results)), each = 3L), , drop = FALSE]
  is_top <- df_Polygons$Triangle == "top-left"
  df_Polygons$x <- df_Polygons$x + rep(c(-0.48, -0.48, 0.48), nrow(df_Results))
  df_Polygons$x[!is_top] <- rep(df_Results$x[df_Results$Triangle == "bottom-right"], each = 3L) +
    rep(c(-0.48, 0.48, 0.48), sum(df_Results$Triangle == "bottom-right"))
  df_Polygons$y <- df_Polygons$y + ifelse(is_top,
    rep(c(-0.48, 0.48, 0.48), nrow(df_Results)),
    rep(c(-0.48, -0.48, 0.48), nrow(df_Results)))
  df_Polygons$PolygonId <- rep(seq_len(nrow(df_Results)), each = 3L)
  df_Stars <- df_Results
  df_Stars$x <- df_Stars$x + ifelse(df_Stars$Triangle == "top-left", -0.16, 0.16)
  df_Stars$y <- df_Stars$y + ifelse(df_Stars$Triangle == "top-left", 0.16, -0.16)

  # The invisible point layer supplies a separate orientation key, with actual
  # triangular glyphs rather than relying on a title or caption.
  DrawTriangleKey <- function(data, params, size) {
    top <- data$shape == 0
    grid::polygonGrob(
      x = if (top) c(0.1, 0.1, 0.9) else c(0.1, 0.9, 0.9),
      y = if (top) c(0.1, 0.9, 0.9) else c(0.1, 0.1, 0.9),
      gp = grid::gpar(fill = "grey65", col = "grey40")
    )
  }
  BuildPlot <- function(star_column) {
    ggplot2::ggplot(df_Polygons, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_polygon(
        ggplot2::aes(fill = .data$R, group = .data$PolygonId),
        colour = "grey85", linewidth = 0.2
      ) +
      ggplot2::geom_text(
        data = df_Stars,
        ggplot2::aes(label = .data[[star_column]]), size = 3.2, fontface = "bold"
      ) +
      ggplot2::geom_point(
        data = df_Results[c(1L, nrow(df_Results)), , drop = FALSE],
        ggplot2::aes(shape = .data$Split), alpha = 0, size = 0,
        key_glyph = DrawTriangleKey
      ) +
      ggplot2::scale_shape_manual(
        name = split_label, limits = split_levels,
        values = stats::setNames(c(0, 1), split_levels),
        labels = paste0(split_levels, c(" (top-left)", " (bottom-right)"))
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#B2182B", mid = "white", high = "#2166AC",
        midpoint = 0, limits = c(-1, 1), na.value = "grey95", name = "r"
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(vars_Columns), labels = unname(column_labels),
        limits = c(0.5, length(vars_Columns) + 0.5), expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq_along(vars_Rows), labels = rev(unname(row_labels)),
        limits = c(0.5, length(vars_Rows) + 0.5), expand = c(0, 0)
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
  BuildBranch <- function(branch, star_column) {
    list(
      r = lapply(correlations, function(result) result[[branch]]$r),
      p = lapply(correlations, function(result) result[[branch]]$p),
      npairs = lapply(correlations, function(result) result[[branch]]$npairs),
      plot = BuildPlot(star_column)
    )
  }
  out <- list(
    Unadjusted = BuildBranch("Unadjusted", "stars"),
    FDRCorrected = BuildBranch("FDRCorrected", "stars_FDR"),
    Correlations = correlations,
    Results = df_Results,
    Metadata = list(
      SplitVariable = split_var, SplitLevels = split_levels,
      TriangleOrder = stats::setNames(c("top-left", "bottom-right"), split_levels),
      AxisOrder = list(rows = vars_Rows, columns = vars_Columns),
      Method = correlations[[1]]$method, Covariates = covariates,
      CovariatesMissing = lapply(correlations, function(result) result$CovariatesMissing),
      Relabel = Relabel, TreatOrdinalAs = TreatOrdinalAs, MinN = min_n, Eps = eps,
      FDRScope = fdr_scope, FDRFamily = "within_stratum",
      IncludeTrendSymbol = IncludeTrendSymbol,
      MissingSplitRows = sum(missing_split),
      UnselectedSplitRows = sum(!missing_split & !split_values %in% split_levels),
      StratumRows = stats::setNames(vapply(split_levels, function(level) {
        sum(!missing_split & split_values %in% level)
      }, integer(1)), split_levels),
      Unavailable = df_Results[df_Results$Status == "unavailable",
        c("Split", "XVar", "YVar", "nPairs", "Status"), drop = FALSE]
    )
  )
  if (interactive == "plotly") {
    BuildWidget <- function(plot) {
      # Add the Plotly-specific aesthetic only during conversion, avoiding
      # unknown-aesthetic warnings and shared ggproto mutations in static plots.
      polygon_layer <- ggplot2::geom_polygon(
        ggplot2::aes(fill = .data$R, group = .data$PolygonId), colour = "grey85", linewidth = 0.2
      )
      polygon_layer$mapping$text <- ggplot2::aes(text = .data$Tooltip)$text
      plot$layers[[1]] <- polygon_layer
      widget <- plotly::ggplotly(plot, tooltip = "text")
      plotly::layout(widget, margin = list(b = 100, l = 100))
    }
    out$Interactive <- list(Plotly = list(
      Unadjusted = BuildWidget(out$Unadjusted$plot),
      FDRCorrected = BuildWidget(out$FDRCorrected$plot)
    ))
  }
  out$p <- out$Unadjusted
  out$p_fdr <- out$FDRCorrected
  out
}
