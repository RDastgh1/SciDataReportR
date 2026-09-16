#' Compare correlations between two independent groups
#'
#' Computes correlations or partial correlations separately within two
#' independent groups, compares corresponding correlations, and visualizes
#' the between-group difference as a heatmap. Tile color represents
#' `DeltaR = r_comparison - r_reference`, significance stars represent the
#' statistical test comparing the two correlations, and striped tiles indicate
#' correlations with opposite signs between groups.
#'
#' Correlations are calculated using [PlotCorrelationsHeatmap()] so variable
#' handling, covariate adjustment, ordinal handling, labels, and missing-data
#' behavior remain consistent with the SciDataReportR correlation workflow.
#'
#' @details
#' Pearson correlations without covariates use the usual independent-samples
#' Fisher-z comparison. Spearman and residualized partial-correlation
#' comparisons use an approximate Fisher-z calculation; inspect
#' `Results$InferenceStatus` or `Metadata$InferenceStatus` before interpreting
#' those p-values.
#'
#' When `comparison_group` and `reference_group` are both omitted for a
#' two-level factor, the first factor level is used as the reference group and
#' the second factor level as the comparison group. For non-factor grouping
#' variables, observed order is used. When more than two groups are present,
#' both groups must be specified explicitly.
#'
#' @param data A data frame.
#' @param predictor_vars Character vector of predictor variables. If `NULL`,
#'   variable selection is inherited from [PlotCorrelationsHeatmap()].
#' @param outcome_vars Optional character vector of outcome variables. If
#'   `NULL`, the same variables are used on both axes.
#' @param group_var Character string naming the grouping variable.
#' @param comparison_group Optional level of `group_var` used as the comparison
#'   group. Positive DeltaR values indicate a more positive correlation in this
#'   group relative to the reference group.
#' @param reference_group Optional level of `group_var` used as the reference
#'   group.
#' @param covariates Optional character vector of covariates used to calculate
#'   partial correlations within each group.
#' @param method Correlation method. Either `"pearson"` or `"spearman"`.
#' @param Relabel Logical indicating whether variable labels should be used
#'   when available.
#' @param TreatOrdinalAs Passed to [PlotCorrelationsHeatmap()].
#' @param min_n Minimum number of complete observations required for an
#'   individual correlation.
#' @param eps Variance tolerance passed to [PlotCorrelationsHeatmap()].
#' @param fdr_scope Scope for FDR correction of correlation-comparison tests.
#'   One of `"matrix"`, `"per_outcome"`, or `"per_predictor"`.
#' @param reversal_style How correlations with opposite signs should be shown.
#'   One of `"outline"` (default), `"stripe"`, or `"none"`. Stripes are a
#'   static-only display option and require `ggpattern`.
#' @param interactive Optional interactive output. One of `"none"` (default),
#'   `"plotly"`, `"girafe"`, or `"both"`. Static ggplots are always retained
#'   in `Unadjusted$plot` and `FDRCorrected$plot`; requested widgets are added
#'   under `Interactive`.
#' @param low_color Color representing negative DeltaR values.
#' @param mid_color Color representing DeltaR = 0.
#' @param high_color Color representing positive DeltaR values.
#' @param color_limits Limits for the DeltaR color scale. The theoretical
#'   range is -2 to 2.
#' @param triangle Display "full" (default), "upper", or "lower" half of a
#'   symmetric comparison matrix. This affects plots only; returned matrices
#'   and Results remain complete.
#'
#' @return A list containing:
#' \describe{
#'   \item{Correlations}{The original [PlotCorrelationsHeatmap()] objects for
#'     the comparison and reference groups.}
#'   \item{Unadjusted}{Matrices and heatmap using raw comparison p-values.}
#'   \item{FDRCorrected}{Matrices and heatmap using FDR-adjusted comparison
#'     p-values.}
#'   \item{Results}{A tibble with one row per correlation pair.}
#'   \item{DirectionReversal}{Logical matrix indicating opposite correlation
#'     signs between groups.}
#'   \item{Metadata}{Comparison settings, group information, and the
#'     inferential approximation used.}
#'   \item{Interactive}{Optional Plotly and/or ggiraph widgets.}
#' }
#'
#' @examples
#' data(SampleData)
#'
#' # If Sex is a factor with levels c("Male", "Female"),
#' # Male is automatically the reference and Female the comparison.
#' #
#' # res <- PlotCorrelationComparisons(
#' #   data = SampleData,
#' #   predictor_vars = c("age", "AXL", "Ferritin", "IL_6"),
#' #   outcome_vars = c("Cortisol", "Insulin"),
#' #   group_var = "Sex",
#' #   covariates = "education",
#' #   method = "spearman",
#' #   triangle = "upper"
#' # )
#' #
#' # res$FDRCorrected$plot
#'
#' @export
PlotCorrelationComparisons <- function(
    data,
    predictor_vars = NULL,
    outcome_vars = NULL,
    group_var,
    comparison_group = NULL,
    reference_group = NULL,
    covariates = NULL,
    method = "pearson",
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    min_n = 4,
    eps = 1e-12,
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    reversal_style = c("outline", "stripe", "none"),
    interactive = c("none", "plotly", "girafe", "both"),
    low_color = "#B2182B",
    mid_color = "white",
    high_color = "#2166AC",
    color_limits = c(-2, 2),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    triangle = c("full", "upper", "lower")) {

  # Validate inputs

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.logical(cluster_rows) || length(cluster_rows) != 1 || is.na(cluster_rows)) {
    stop("cluster_rows must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(cluster_columns) || length(cluster_columns) != 1 || is.na(cluster_columns)) {
    stop("cluster_columns must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    length(group_var) != 1 ||
    !is.character(group_var) ||
    !group_var %in% names(data)
  ) {
    stop(
      "`group_var` must be the name of one column in `data`.",
      call. = FALSE
    )
  }

  method <- match.arg(
    tolower(method),
    c("pearson", "spearman")
  )

  fdr_scope <- match.arg(fdr_scope)
  reversal_style <- match.arg(reversal_style)
  interactive <- match.arg(interactive)
  triangle <- match.arg(triangle)

  if (
    length(color_limits) != 2 ||
    !all(is.finite(color_limits)) ||
    color_limits[1] >= 0 ||
    color_limits[2] <= 0
  ) {
    stop(
      "`color_limits` must contain one negative and one positive value.",
      call. = FALSE
    )
  }

  if (color_limits[1] < -2 || color_limits[2] > 2) {
    stop(
      "`color_limits` cannot extend beyond the theoretical DeltaR range of -2 to 2.",
      call. = FALSE
    )
  }

  ScidrValidateVariables(
    data,
    group_var,
    "group_var"
  )

  ScidrValidateVariables(
    data,
    predictor_vars,
    "predictor_vars"
  )

  ScidrValidateVariables(
    data,
    outcome_vars,
    "outcome_vars"
  )

  ScidrValidateVariables(
    data,
    covariates,
    "covariates"
  )

  # Determine groups

  group_values <- data[[group_var]]

  if (is.factor(group_values)) {

    available_groups <- levels(
      droplevels(group_values)
    )

  } else {

    available_groups <- unique(
      as.character(
        group_values[!is.na(group_values)]
      )
    )
  }

  if (length(available_groups) < 2) {
    stop(
      sprintf(
        "`group_var` '%s' must contain at least two observed groups.",
        group_var
      ),
      call. = FALSE
    )
  }

  if (
    is.null(comparison_group) &&
    is.null(reference_group)
  ) {

    if (length(available_groups) != 2) {
      stop(
        paste0(
          "`group_var` '",
          group_var,
          "' contains ",
          length(available_groups),
          " observed groups: ",
          paste(available_groups, collapse = ", "),
          ". Specify `comparison_group` and `reference_group` explicitly ",
          "when more than two groups are present."
        ),
        call. = FALSE
      )
    }

    reference_group <- available_groups[1]
    comparison_group <- available_groups[2]
  }

  if (
    is.null(comparison_group) &&
    !is.null(reference_group)
  ) {

    if (!as.character(reference_group) %in% available_groups) {
      stop(
        sprintf(
          "`reference_group` '%s' was not found in `%s`. Available groups: %s.",
          reference_group,
          group_var,
          paste(available_groups, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    remaining_groups <- setdiff(
      available_groups,
      as.character(reference_group)
    )

    if (length(remaining_groups) != 1) {
      stop(
        paste0(
          "Could not uniquely determine `comparison_group` because `",
          group_var,
          "` contains more than two groups. Specify both groups explicitly."
        ),
        call. = FALSE
      )
    }

    comparison_group <- remaining_groups
  }

  if (
    !is.null(comparison_group) &&
    is.null(reference_group)
  ) {

    if (!as.character(comparison_group) %in% available_groups) {
      stop(
        sprintf(
          "`comparison_group` '%s' was not found in `%s`. Available groups: %s.",
          comparison_group,
          group_var,
          paste(available_groups, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    remaining_groups <- setdiff(
      available_groups,
      as.character(comparison_group)
    )

    if (length(remaining_groups) != 1) {
      stop(
        paste0(
          "Could not uniquely determine `reference_group` because `",
          group_var,
          "` contains more than two groups. Specify both groups explicitly."
        ),
        call. = FALSE
      )
    }

    reference_group <- remaining_groups
  }

  comparison_group <- as.character(comparison_group)
  reference_group <- as.character(reference_group)

  if (!comparison_group %in% available_groups) {
    stop(
      sprintf(
        "`comparison_group` '%s' was not found in `%s`. Available groups: %s.",
        comparison_group,
        group_var,
        paste(available_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!reference_group %in% available_groups) {
    stop(
      sprintf(
        "`reference_group` '%s' was not found in `%s`. Available groups: %s.",
        reference_group,
        group_var,
        paste(available_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (identical(comparison_group, reference_group)) {
    stop(
      "`comparison_group` and `reference_group` must be different.",
      call. = FALSE
    )
  }

  message(
    "Comparing ",
    comparison_group,
    " to ",
    reference_group,
    ": DeltaR = r(",
    comparison_group,
    ") - r(",
    reference_group,
    ")."
  )

  # Prepare data

  comparison_data <- data %>%
    dplyr::filter(
      as.character(.data[[group_var]]) == comparison_group
    )

  reference_data <- data %>%
    dplyr::filter(
      as.character(.data[[group_var]]) == reference_group
    )

  if (nrow(comparison_data) < min_n) {
    stop(
      sprintf(
        "Comparison group '%s' contains only %s observations; at least %s are required.",
        comparison_group,
        nrow(comparison_data),
        min_n
      ),
      call. = FALSE
    )
  }

  if (nrow(reference_data) < min_n) {
    stop(
      sprintf(
        "Reference group '%s' contains only %s observations; at least %s are required.",
        reference_group,
        nrow(reference_data),
        min_n
      ),
      call. = FALSE
    )
  }

  # Calculate group-specific correlations

  comparison_correlations <- PlotCorrelationsHeatmap(
    data = comparison_data,
    predictor_vars = predictor_vars,
    outcome_vars = outcome_vars,
    covariates = covariates,
    method = method,
    Relabel = Relabel,
    TreatOrdinalAs = TreatOrdinalAs,
    min_n = min_n,
    eps = eps,
    fdr_scope = fdr_scope
  )

  reference_correlations <- PlotCorrelationsHeatmap(
    data = reference_data,
    predictor_vars = predictor_vars,
    outcome_vars = outcome_vars,
    covariates = covariates,
    method = method,
    Relabel = Relabel,
    TreatOrdinalAs = TreatOrdinalAs,
    min_n = min_n,
    eps = eps,
    fdr_scope = fdr_scope
  )

  r_comparison <- comparison_correlations$Unadjusted$r
  r_reference <- reference_correlations$Unadjusted$r

  n_comparison <- comparison_correlations$Unadjusted$npairs
  n_reference <- reference_correlations$Unadjusted$npairs

  if (
    !identical(dim(r_comparison), dim(r_reference)) ||
    !identical(dimnames(r_comparison), dimnames(r_reference))
  ) {
    stop(
      paste(
        "The two group-specific correlation matrices do not have matching",
        "dimensions and variable names."
      ),
      call. = FALSE
    )
  }

  predictor_names <- rownames(r_comparison)
  outcome_names <- colnames(r_comparison)

  # Determine effective covariate dimensions

  k_comparison <- matrix(
    0,
    nrow = nrow(r_comparison),
    ncol = ncol(r_comparison),
    dimnames = dimnames(r_comparison)
  )

  k_reference <- k_comparison

  if (!is.null(covariates) && length(covariates) > 0) {

    for (i in seq_along(predictor_names)) {

      for (j in seq_along(outcome_names)) {

        x_name <- predictor_names[i]
        y_name <- outcome_names[j]

        vars_needed <- unique(
          c(x_name, y_name, covariates)
        )

        for (group_name in c("comparison", "reference")) {

          group_data <- if (group_name == "comparison") {
            comparison_data
          } else {
            reference_data
          }

          complete_rows <- stats::complete.cases(
            group_data[, vars_needed, drop = FALSE]
          )

          tmp <- group_data[
            complete_rows,
            vars_needed,
            drop = FALSE
          ]

          k_value <- 0

          if (nrow(tmp) >= min_n) {

            cov_df <- tmp[, covariates, drop = FALSE]

            for (nm in names(cov_df)) {

              if (
                is.character(cov_df[[nm]]) ||
                is.logical(cov_df[[nm]])
              ) {
                cov_df[[nm]] <- factor(cov_df[[nm]])
              }
            }

            mm <- tryCatch(
              stats::model.matrix(
                ~ .,
                data = cov_df
              ),
              error = function(e) NULL
            )

            if (!is.null(mm)) {

              mm <- mm[
                ,
                colnames(mm) != "(Intercept)",
                drop = FALSE
              ]

              if (ncol(mm) > 0) {

                keep <- vapply(
                  seq_len(ncol(mm)),
                  function(column_index) {

                    v <- suppressWarnings(
                      stats::var(mm[, column_index])
                    )

                    is.finite(v) &&
                      !is.na(v) &&
                      v > eps
                  },
                  logical(1)
                )

                mm <- mm[, keep, drop = FALSE]

                if (
                  ncol(mm) > 0 &&
                  nrow(tmp) > (ncol(mm) + 2)
                ) {
                  k_value <- qr(mm)$rank
                }
              }
            }
          }

          if (group_name == "comparison") {
            k_comparison[i, j] <- k_value
          } else {
            k_reference[i, j] <- k_value
          }
        }
      }
    }
  }

  # Compare correlations

  delta_r <- r_comparison - r_reference

  r_comparison_clipped <- pmin(
    pmax(r_comparison, -1 + 1e-12),
    1 - 1e-12
  )

  r_reference_clipped <- pmin(
    pmax(r_reference, -1 + 1e-12),
    1 - 1e-12
  )

  fisher_comparison <- atanh(
    r_comparison_clipped
  )

  fisher_reference <- atanh(
    r_reference_clipped
  )

  cohens_q <- fisher_comparison -
    fisher_reference

  denominator_comparison <-
    n_comparison - k_comparison - 3

  denominator_reference <-
    n_reference - k_reference - 3

  variance_multiplier <- if (method == "spearman") {
    1.06
  } else {
    1
  }

  valid_test <-
    is.finite(r_comparison) &
    is.finite(r_reference) &
    is.finite(denominator_comparison) &
    is.finite(denominator_reference) &
    denominator_comparison > 0 &
    denominator_reference > 0

  standard_error <- matrix(
    NA_real_,
    nrow = nrow(r_comparison),
    ncol = ncol(r_comparison),
    dimnames = dimnames(r_comparison)
  )

  standard_error[valid_test] <- sqrt(
    variance_multiplier /
      denominator_comparison[valid_test] +
      variance_multiplier /
      denominator_reference[valid_test]
  )

  fisher_z <- cohens_q /
    standard_error

  p_value <- 2 * stats::pnorm(
    -abs(fisher_z)
  )

  p_value[!valid_test] <- NA_real_

  direction_reversal <-
    sign(r_comparison) != sign(r_reference)

  direction_reversal[
    is.na(r_comparison) |
      is.na(r_reference) |
      r_comparison == 0 |
      r_reference == 0
  ] <- FALSE

  symmetric_matrix <-
    identical(predictor_names, outcome_names) &&
    nrow(p_value) == ncol(p_value)
  p_adjusted <- ApplyFDRCorrection(
    p_value,
    fdr_scope = fdr_scope,
    outcome_margin = 2,
    method = "fdr",
    symmetric = symmetric_matrix
  )

  # Apply labels

  # Reuse the display labels already prepared by PlotCorrelationsHeatmap().
  plot_data_labels <- comparison_correlations$Unadjusted$plot$data
  predictor_labels <- stats::setNames(
    as.character(plot_data_labels$XLabel),
    plot_data_labels$XVar
  )
  predictor_labels <- predictor_labels[!duplicated(names(predictor_labels))]
  predictor_labels <- predictor_labels[predictor_names]
  outcome_labels <- stats::setNames(
    as.character(plot_data_labels$YLabel),
    plot_data_labels$YVar
  )
  outcome_labels <- outcome_labels[!duplicated(names(outcome_labels))]
  outcome_labels <- outcome_labels[outcome_names]

  inference_status <- if (method == "pearson" && length(covariates) == 0) {
    "Fisher-z test for independent Pearson correlations"
  } else if (method == "spearman" && length(covariates) == 0) {
    "Approximate Fisher-z test for independent Spearman correlations"
  } else if (method == "pearson") {
    "Approximate Fisher-z test for residualized partial correlations"
  } else {
    "Approximate Fisher-z test for residualized Spearman correlations"
  }

  # Build tidy results

  results <- tidyr::expand_grid(
    Predictor = predictor_names,
    Outcome = outcome_names
  ) %>%
    dplyr::mutate(
      PredictorIndex = match(
        Predictor,
        predictor_names
      ),
      OutcomeIndex = match(
        Outcome,
        outcome_names
      ),
      RComparison = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ r_comparison[.x, .y]
      ),
      RReference = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ r_reference[.x, .y]
      ),
      DeltaR = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ delta_r[.x, .y]
      ),
      CohensQ = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ cohens_q[.x, .y]
      ),
      FisherZ = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ fisher_z[.x, .y]
      ),
      PValue = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ p_value[.x, .y]
      ),
      PAdjusted = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ p_adjusted[.x, .y]
      ),
      NComparison = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ n_comparison[.x, .y]
      ),
      NReference = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ n_reference[.x, .y]
      ),
      CovariateDfComparison = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ k_comparison[.x, .y]
      ),
      CovariateDfReference = purrr::map2_dbl(
        PredictorIndex,
        OutcomeIndex,
        ~ k_reference[.x, .y]
      ),
      DirectionReversal = purrr::map2_lgl(
        PredictorIndex,
        OutcomeIndex,
        ~ direction_reversal[.x, .y]
      ),
      TestAvailable = purrr::map2_lgl(
        PredictorIndex,
        OutcomeIndex,
        ~ valid_test[.x, .y]
      ),
      PredictorLabel = unname(
        predictor_labels[Predictor]
      ),
      OutcomeLabel = unname(
        outcome_labels[Outcome]
      ),
      Direction = dplyr::case_when(
        DirectionReversal ~ "Opposite signs",
        RComparison > 0 & RReference > 0 ~ "Both positive",
        RComparison < 0 & RReference < 0 ~ "Both negative",
        TRUE ~ "No reversal"
      ),
      InferenceStatus = inference_status,
      CellId = paste(Predictor, Outcome, sep = "__")
    ) %>%
    dplyr::select(
      Predictor,
      Outcome,
      PredictorLabel,
      OutcomeLabel,
      RComparison,
      RReference,
      DeltaR,
      CohensQ,
      FisherZ,
      PValue,
      PAdjusted,
      NComparison,
      NReference,
      CovariateDfComparison,
      CovariateDfReference,
      DirectionReversal,
      Direction,
      InferenceStatus,
      TestAvailable,
      CellId
    )

  if (symmetric_matrix) {

    results <- results %>%
      dplyr::filter(
        Predictor != Outcome
      )
  }

  results <- results %>%
    dplyr::mutate(
      HoverText = paste0(
        PredictorLabel,
        " × ",
        OutcomeLabel,
        "<br>",
        comparison_group,
        ": r = ",
        sprintf("%.2f", RComparison),
        ", n = ",
        NComparison,
        "<br>",
        reference_group,
        ": r = ",
        sprintf("%.2f", RReference),
        ", n = ",
        NReference,
        "<br>Δr (",
        comparison_group,
        " - ",
        reference_group,
        ") = ",
        sprintf("%.2f", DeltaR),
        "<br>Cohen's q = ",
        sprintf("%.2f", CohensQ),
        "<br>Z = ",
        sprintf("%.2f", FisherZ),
        "<br>p = ",
        format.pval(
          PValue,
          digits = 3,
          eps = 0.001
        ),
        "<br>FDR p = ",
        format.pval(
          PAdjusted,
          digits = 3,
          eps = 0.001
        ),
        "<br>",
        Direction,
        "<br>",
        InferenceStatus
      ),
      ComparisonStatus = dplyr::case_when(
        TestAvailable ~ "Testable",
        !is.finite(RComparison) | !is.finite(RReference) ~
          "Unavailable: correlation could not be estimated in one or both groups",
        TRUE ~
          "Unavailable: insufficient pairwise observations after covariate adjustment"
      ),
      InferenceApproximate = method == "spearman" || length(covariates) > 0,
      stars = .FormatPValueStars(PValue),
      stars_FDR = .FormatPValueStars(PAdjusted),
      ReversalPattern = dplyr::if_else(
        DirectionReversal,
        "Opposite signs",
        "Same sign"
      ),
      PredictorLabel = factor(
        PredictorLabel,
        levels = rev(unname(predictor_labels))
      ),
      OutcomeLabel = factor(
        OutcomeLabel,
        levels = unname(outcome_labels)
      )
    )

  triangle_display <- .ResolveHeatmapTriangle(
    data = results, row_id = "Predictor", column_id = "Outcome", value = "DeltaR",
    row_order = predictor_vars, column_order = outcome_vars,
    cluster_rows = cluster_rows, cluster_columns = cluster_columns,
    triangle = triangle
  )
  AxisOrder <- triangle_display$AxisOrder
  results_display <- triangle_display$PlotData
  predictor_axis_labels <- stats::setNames(
    as.character(predictor_labels[AxisOrder$rows]), AxisOrder$rows
  )
  outcome_axis_labels <- stats::setNames(
    as.character(outcome_labels[AxisOrder$columns]), AxisOrder$columns
  )
  results$PredictorOrder <- factor(results$Predictor, levels = rev(AxisOrder$rows))
  results$OutcomeOrder <- factor(results$Outcome, levels = AxisOrder$columns)
  results_display$PredictorOrder <- factor(results_display$Predictor, levels = rev(AxisOrder$rows))
  results_display$OutcomeOrder <- factor(results_display$Outcome, levels = AxisOrder$columns)

  # Build plots

  build_comparison_heatmap <- function(
      plot_data,
      star_column,
      backend = c("static", "plotly", "girafe"),
      style = reversal_style) {

    backend <- match.arg(backend)
    style <- match.arg(style, c("outline", "stripe", "none"))
    if (backend != "static" && style == "stripe") style <- "outline"
    if (style == "stripe" && !requireNamespace("ggpattern", quietly = TRUE)) {
      warning("Package 'ggpattern' is not installed; using reversal outlines.", call. = FALSE)
      style <- "outline"
    }

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = OutcomeOrder, y = PredictorOrder, fill = DeltaR)
    )

    if (backend == "girafe") {
      p <- p +
        ggiraph::geom_tile_interactive(
          ggplot2::aes(tooltip = HoverText, data_id = CellId),
          colour = "white",
          linewidth = 0.25
        )
    } else if (style == "stripe") {
      p <- p +
        ggpattern::geom_tile_pattern(
          ggplot2::aes(pattern = ReversalPattern),
          pattern_fill = NA,
          pattern_colour = "black",
          pattern_density = 0.08,
          pattern_spacing = 0.04,
          pattern_angle = 45,
          colour = "white",
          linewidth = 0.25
        ) +
        ggpattern::scale_pattern_manual(
          values = c("Same sign" = "none", "Opposite signs" = "stripe"),
          name = "Direction reversal"
        )
    } else if (backend == "plotly") {
      p <- p +
        suppressWarnings(
          ggplot2::geom_tile(
            ggplot2::aes(text = HoverText),
            colour = "white",
            linewidth = 0.25
          )
        )
    } else {
      p <- p + ggplot2::geom_tile(colour = "white", linewidth = 0.25)
    }

    if (style == "outline" && any(plot_data$DirectionReversal, na.rm = TRUE)) {
      p <- p +
        ggplot2::geom_tile(
          data = plot_data %>% dplyr::filter(DirectionReversal),
          colour = "black",
          fill = NA,
          linewidth = 0.7
        )
    }

    p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data[[star_column]]),
        size = 5,
        colour = "black",
        na.rm = TRUE
      ) +
      .GetHeatmapColorScale(
        low_color = low_color,
        mid_color = mid_color,
        high_color = high_color,
        fill_midpoint = 0,
        fill_limits = color_limits,
        fill_oob = scales::squish,
        name = paste0("Δr\n", comparison_group, " - ", reference_group)
      ) +
      ggplot2::scale_x_discrete(labels = outcome_axis_labels) +
      ggplot2::scale_y_discrete(labels = predictor_axis_labels) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        legend.position = "right"
      )
  }

  raw_plot <- build_comparison_heatmap(results_display, "stars")
  fdr_plot <- build_comparison_heatmap(results_display, "stars_FDR")

  # Return result

  out <- list(
    Correlations = list(
      Comparison = comparison_correlations,
      Reference = reference_correlations
    ),
    Unadjusted = list(
      delta_r = delta_r,
      cohens_q = cohens_q,
      z = fisher_z,
      p = p_value,
      r_comparison = r_comparison,
      r_reference = r_reference,
      n_comparison = n_comparison,
      n_reference = n_reference,
      plot = raw_plot
    ),
    FDRCorrected = list(
      delta_r = delta_r,
      cohens_q = cohens_q,
      z = fisher_z,
      p = p_adjusted,
      r_comparison = r_comparison,
      r_reference = r_reference,
      n_comparison = n_comparison,
      n_reference = n_reference,
      plot = fdr_plot
    ),
    Results = results,
    AxisOrder = AxisOrder,
    DirectionReversal = direction_reversal,
    Metadata = list(
      GroupVariable = group_var,
      ComparisonGroup = comparison_group,
      ReferenceGroup = reference_group,
      DeltaRDefinition = paste0(
        "r(",
        comparison_group,
        ") - r(",
        reference_group,
        ")"
      ),
      CorrelationMethod = method,
      Covariates = covariates,
      FDRScope = fdr_scope,
      InferenceStatus = inference_status,
      SpearmanVarianceMultiplier = if (method == "spearman") {
        1.06
      } else {
        NA_real_
      },
      ColorLimits = color_limits,
      ReversalStyle = reversal_style,
      Triangle = triangle,
      TriangleApplied = AxisOrder$triangle_applied
    )
  )

  if (interactive %in% c("plotly", "both")) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("Package 'plotly' is required for interactive = 'plotly' or 'both'.", call. = FALSE)
    }
    out$Interactive$Plotly <- list(
      Unadjusted = plotly::ggplotly(
        build_comparison_heatmap(results_display, "stars", backend = "plotly", style = "outline"),
        tooltip = "text"
      ),
      FDRCorrected = plotly::ggplotly(
        build_comparison_heatmap(results_display, "stars_FDR", backend = "plotly", style = "outline"),
        tooltip = "text"
      )
    )
  }

  if (interactive %in% c("girafe", "both")) {
    if (!requireNamespace("ggiraph", quietly = TRUE)) {
      stop("Package 'ggiraph' is required for interactive = 'girafe' or 'both'.", call. = FALSE)
    }
    out$Interactive$Girafe <- list(
      Unadjusted = ggiraph::girafe(
        ggobj = build_comparison_heatmap(results_display, "stars", backend = "girafe", style = "outline")
      ),
      FDRCorrected = ggiraph::girafe(
        ggobj = build_comparison_heatmap(results_display, "stars_FDR", backend = "girafe", style = "outline")
      )
    )
  }

  out$p <- out$Unadjusted
  out$p_fdr <- out$FDRCorrected

  out
}


#' Add DeltaR values and significance stars to a correlation comparison heatmap
#'
#' Adds correlation differences and significance stars to a heatmap returned
#' by [PlotCorrelationComparisons()]. This is the correlation-comparison
#' counterpart to [add_r_and_stars()].
#'
#' @param res An object returned by [PlotCorrelationComparisons()].
#' @param star_from Whether stars should use `"fdr"` or `"raw"` comparison
#'   p-values.
#' @param delta_digits Number of decimal places used for DeltaR.
#' @param delta_size Text size for DeltaR labels.
#' @param star_size Text size for significance stars.
#' @param delta_color Color for DeltaR labels.
#' @param star_color Color for significance stars.
#' @param delta_nudge_y Vertical position adjustment for DeltaR.
#' @param star_nudge_y Vertical position adjustment for stars.
#' @param remove_existing_stars Logical. Remove the original star-only layer
#'   before adding the combined annotations.
#'
#' @return A ggplot containing DeltaR values and significance stars.
#'
#' @examples
#' # res <- PlotCorrelationComparisons(...)
#' # add_delta_r_and_stars(res)
#'
#' @export
add_delta_r_and_stars <- function(
    res,
    star_from = c("fdr", "raw"),
    delta_digits = 2,
    delta_size = 3,
    star_size = 5,
    delta_color = "black",
    star_color = "black",
    delta_nudge_y = -0.18,
    star_nudge_y = 0.20,
    remove_existing_stars = TRUE) {

  star_from <- match.arg(star_from)

  if (
    is.null(res$Results) ||
    is.null(res$Unadjusted$plot) ||
    is.null(res$FDRCorrected$plot)
  ) {
    stop(
      paste(
        "`res` must be an object returned by",
        "`PlotCorrelationComparisons()`."
      ),
      call. = FALSE
    )
  }

  if (star_from == "fdr") {
    p <- res$FDRCorrected$plot
    p_column <- "PAdjusted"
  } else {
    p <- res$Unadjusted$plot
    p_column <- "PValue"
  }

  d <- p$data

  if (!all(c("DeltaR", p_column) %in% names(d))) {
    stop(
      "The selected comparison plot does not contain DeltaR and p-values.",
      call. = FALSE
    )
  }

  d <- d %>%
    dplyr::mutate(
      DeltaLabel = dplyr::if_else(
        is.na(DeltaR),
        "",
        sprintf(
          paste0("%+.", delta_digits, "f"),
          DeltaR
        )
      ),
      StarLabel = dplyr::case_when(
        is.na(.data[[p_column]]) ~ "",
        .data[[p_column]] < 0.001 ~ "***",
        .data[[p_column]] < 0.01 ~ "**",
        .data[[p_column]] < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  if (remove_existing_stars && length(p$layers) > 0) {

    keep_layer <- vapply(
      p$layers,
      function(layer) {

        if (!inherits(layer$geom, "GeomText")) {
          return(TRUE)
        }

        if (is.null(layer$mapping$label)) {
          return(TRUE)
        }

        label_expression <- paste(
          deparse(layer$mapping$label),
          collapse = ""
        )

        !grepl(
          "star",
          label_expression,
          ignore.case = TRUE
        )
      },
      logical(1)
    )

    p$layers <- p$layers[keep_layer]
  }

  p +
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(
        label = DeltaLabel
      ),
      inherit.aes = TRUE,
      colour = delta_color,
      size = delta_size,
      position = ggplot2::position_nudge(
        y = delta_nudge_y
      ),
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(
        label = StarLabel
      ),
      inherit.aes = TRUE,
      colour = star_color,
      size = star_size,
      position = ggplot2::position_nudge(
        y = star_nudge_y
      ),
      na.rm = TRUE
    )
}
