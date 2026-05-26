#' Plot volcano-style association effects
#'
#' Screen a set of predictor variables against one outcome and visualize the
#' association results as volcano plots. The function supports continuous
#' outcomes using standardized beta values, and two-group categorical outcomes
#' using either Cohen-style standardized effects or log2 fold change effects.
#' Optional covariates are included in the model for each predictor.
#'
#' Missingness is handled pairwise by model, so each predictor is analyzed using
#' all complete observations available for that predictor, the outcome, and any
#' covariates. Invalid predictor variables are removed with a warning.
#'
#' @param Data A data frame.
#' @param xVars Character vector of predictor variable names to screen.
#' @param yVar Character string naming the outcome variable.
#' @param Covariates Optional character vector of covariate variable names.
#' @param OutcomeType Outcome type. One of `"auto"`, `"continuous"`, or
#'   `"categorical"`. If `"auto"`, numeric outcomes are treated as continuous
#'   and nonnumeric two-level outcomes are treated as categorical.
#' @param EffectMetric Effect metric. One of `"auto"`, `"cohens_d"`, or
#'   `"log2fc"`. Used for two-group categorical outcomes. For continuous
#'   outcomes, the effect is always a standardized beta from a model using
#'   scaled predictor and scaled outcome.
#' @param AdjustMethod Multiple-comparison correction method passed to
#'   `stats::p.adjust()`. Default is `"fdr"`.
#' @param Alpha Significance threshold for raw and adjusted p-values. Default is
#'   `0.05`.
#' @param Format Color format. One of `"tiered"`, `"classic"`, `"fdr_only"`,
#'   `"directional"`, `"effect_gradient"`, `"minimal"`, or `"neon"`.
#' @param LabelMode Labeling mode. One of `"none"`, `"top_n"`,
#'   `"significant"`, `"fdr"`, or `"extreme"`.
#' @param TopN Number of variables to label when `LabelMode` is `"top_n"` or
#'   `"extreme"`. Default is `10`.
#' @param Relabel Logical. If `TRUE`, variable labels are used when available.
#'   Labels are pulled first from `Codebook` if supplied, then from variable
#'   label attributes. Default is `TRUE`.
#' @param Codebook Optional codebook data frame with columns `Variable` and
#'   `Label`.
#' @param InteractiveLabels Logical. If `TRUE`, a `text` aesthetic is added for
#'   compatibility with `plotly::ggplotly(tooltip = "text")`. Default is `TRUE`.
#' @return A named list with `RawPPlot`, `FDRPlot`, and `ResultsTable`.
#'   `RawPPlot` uses `-log10(PValue)` on the y-axis. `FDRPlot` uses
#'   `-log10(FDR)` on the y-axis. `ResultsTable` is a tibble with one row per
#'   analyzed predictor.
#' @examples
#' PlotVolcanoEffects(
#'   Data = mtcars,
#'   xVars = c("disp", "hp", "drat", "wt", "qsec"),
#'   yVar = "mpg",
#'   Covariates = "cyl",
#'   OutcomeType = "continuous",
#'   LabelMode = "top_n",
#'   TopN = 3
#' )
#'
#' PlotVolcanoEffects(
#'   Data = mtcars,
#'   xVars = c("mpg", "disp", "hp", "drat", "wt", "qsec"),
#'   yVar = "am",
#'   OutcomeType = "categorical",
#'   EffectMetric = "cohens_d",
#'   LabelMode = "top_n",
#'   TopN = 3
#' )
#' @export
PlotVolcanoEffects <- function(
    Data,
    xVars,
    yVar,
    Covariates = NULL,
    OutcomeType = c("auto", "continuous", "categorical"),
    EffectMetric = c("auto", "cohens_d", "log2fc"),
    AdjustMethod = "fdr",
    Alpha = 0.05,
    Format = c(
      "tiered",
      "classic",
      "fdr_only",
      "directional",
      "effect_gradient",
      "minimal",
      "neon"
    ),
    LabelMode = c(
      "none",
      "top_n",
      "significant",
      "fdr",
      "extreme"
    ),
    TopN = 10,
    Relabel = TRUE,
    Codebook = NULL,
    InteractiveLabels = TRUE
) {

  OutcomeType <- match.arg(OutcomeType)
  EffectMetric <- match.arg(EffectMetric)
  Format <- match.arg(Format)
  LabelMode <- match.arg(LabelMode)

  # Validate inputs

  if (!is.data.frame(Data)) {
    stop("Data must be a data frame.")
  }

  if (!is.character(xVars) || length(xVars) == 0) {
    stop("xVars must be a non-empty character vector of predictor variable names.")
  }

  if (!is.character(yVar) || length(yVar) != 1) {
    stop("yVar must be a single character string naming the outcome variable.")
  }

  if (!yVar %in% names(Data)) {
    stop(
      "yVar was not found in Data. Expected: ",
      yVar,
      ". Available columns include: ",
      paste(utils::head(names(Data), 20), collapse = ", ")
    )
  }

  if (!is.null(Covariates)) {
    missing_covars <- setdiff(Covariates, names(Data))

    if (length(missing_covars) > 0) {
      stop(
        "The following Covariates were not found in Data: ",
        paste(missing_covars, collapse = ", ")
      )
    }
  }

  missing_xvars <- setdiff(xVars, names(Data))

  if (length(missing_xvars) > 0) {
    warning(
      "The following xVars were not found in Data and were removed: ",
      paste(missing_xvars, collapse = ", "),
      call. = FALSE
    )
  }

  xVars <- intersect(xVars, names(Data))

  if (length(xVars) == 0) {
    stop("No valid xVars remain after checking against Data.")
  }

  if (!is.null(Codebook)) {
    if (!is.data.frame(Codebook)) {
      stop("Codebook must be a data frame when supplied.")
    }

    required_codebook_cols <- c("Variable", "Label")
    missing_codebook_cols <- setdiff(required_codebook_cols, names(Codebook))

    if (length(missing_codebook_cols) > 0) {
      stop(
        "Codebook must contain columns Variable and Label. Missing: ",
        paste(missing_codebook_cols, collapse = ", ")
      )
    }
  }

  # Prepare outcome

  outcome_vector <- Data[[yVar]]

  if (OutcomeType == "auto") {
    if (is.numeric(outcome_vector)) {
      OutcomeType <- "continuous"
    } else {
      OutcomeType <- "categorical"
    }
  }

  if (OutcomeType == "continuous" && !is.numeric(outcome_vector)) {
    stop("OutcomeType is continuous, but yVar is not numeric.")
  }

  if (OutcomeType == "categorical") {
    outcome_nonmissing <- outcome_vector[!is.na(outcome_vector)]
    outcome_levels <- unique(as.character(outcome_nonmissing))

    if (length(outcome_levels) != 2) {
      stop(
        "Categorical outcomes must have exactly two non-missing levels. ",
        "Received ",
        length(outcome_levels),
        " levels: ",
        paste(outcome_levels, collapse = ", ")
      )
    }

    outcome_levels <- sort(outcome_levels)
  } else {
    outcome_levels <- NULL
  }

  if (OutcomeType == "continuous") {
    EffectMetric <- "standardized_beta"
  } else if (EffectMetric == "auto") {
    EffectMetric <- "cohens_d"
  }

  # Prepare labels

  label_lookup <- stats::setNames(xVars, xVars)

  if (Relabel) {
    if (!is.null(Codebook)) {
      codebook_labels <- Codebook %>%
        dplyr::filter(.data$Variable %in% xVars) %>%
        dplyr::mutate(
          Label = dplyr::if_else(
            is.na(.data$Label) | .data$Label == "",
            .data$Variable,
            as.character(.data$Label)
          )
        ) %>%
        dplyr::select(.data$Variable, .data$Label)

      label_lookup[codebook_labels$Variable] <- codebook_labels$Label
    }

    for (this_var in xVars) {
      if (label_lookup[[this_var]] == this_var) {
        attr_label <- attr(Data[[this_var]], "label", exact = TRUE)

        if (!is.null(attr_label) && !is.na(attr_label) && attr_label != "") {
          label_lookup[[this_var]] <- as.character(attr_label)
        }
      }
    }
  }

  # Remove invalid predictors

  numeric_xvars <- xVars[vapply(Data[xVars], is.numeric, logical(1))]
  invalid_xvars <- setdiff(xVars, numeric_xvars)

  if (length(invalid_xvars) > 0) {
    warning(
      "The following xVars were removed because they are not numeric: ",
      paste(invalid_xvars, collapse = ", "),
      call. = FALSE
    )
  }

  xVars <- numeric_xvars

  if (length(xVars) == 0) {
    stop("No numeric xVars remain after removing invalid predictors.")
  }

  # Fit one model per predictor

  results <- purrr::map_dfr(xVars, function(this_var) {

    model_vars <- c(this_var, yVar, Covariates)

    model_data <- Data %>%
      dplyr::select(dplyr::all_of(model_vars)) %>%
      dplyr::filter(stats::complete.cases(.))

    final_n <- nrow(model_data)

    if (final_n < 3) {
      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = NA_real_,
          EffectType = ifelse(
            OutcomeType == "continuous",
            "Standardized beta",
            ifelse(EffectMetric == "log2fc", "Log2 fold change", "Cohen-style d")
          ),
          PValue = NA_real_,
          N = final_n,
          Note = "Too few complete observations"
        )
      )
    }

    if (OutcomeType == "continuous") {
      model_data <- model_data %>%
        dplyr::mutate(
          .volcano_y_scaled = as.numeric(scale(.data[[yVar]])),
          .volcano_x_scaled = as.numeric(scale(.data[[this_var]]))
        )

      if (
        all(is.na(model_data$.volcano_y_scaled)) ||
        all(is.na(model_data$.volcano_x_scaled)) ||
        stats::sd(model_data$.volcano_y_scaled, na.rm = TRUE) == 0 ||
        stats::sd(model_data$.volcano_x_scaled, na.rm = TRUE) == 0
      ) {
        return(
          tibble::tibble(
            Variable = this_var,
            Label = unname(label_lookup[[this_var]]),
            Outcome = yVar,
            OutcomeType = OutcomeType,
            Effect = NA_real_,
            EffectType = "Standardized beta",
            PValue = NA_real_,
            N = final_n,
            Note = "Zero variance in predictor or outcome"
          )
        )
      }

      model_formula <- stats::as.formula(
        paste(
          ".volcano_y_scaled ~ .volcano_x_scaled",
          if (!is.null(Covariates)) {
            paste("+", paste(Covariates, collapse = " + "))
          } else {
            ""
          }
        )
      )

      model_fit <- tryCatch(
        stats::lm(model_formula, data = model_data),
        error = function(e) NULL
      )

      if (is.null(model_fit)) {
        return(
          tibble::tibble(
            Variable = this_var,
            Label = unname(label_lookup[[this_var]]),
            Outcome = yVar,
            OutcomeType = OutcomeType,
            Effect = NA_real_,
            EffectType = "Standardized beta",
            PValue = NA_real_,
            N = final_n,
            Note = "Model failed"
          )
        )
      }

      coefficient_table <- summary(model_fit)$coefficients

      if (!".volcano_x_scaled" %in% rownames(coefficient_table)) {
        return(
          tibble::tibble(
            Variable = this_var,
            Label = unname(label_lookup[[this_var]]),
            Outcome = yVar,
            OutcomeType = OutcomeType,
            Effect = NA_real_,
            EffectType = "Standardized beta",
            PValue = NA_real_,
            N = final_n,
            Note = "Predictor coefficient not estimable"
          )
        )
      }

      effect <- unname(coefficient_table[".volcano_x_scaled", "Estimate"])
      p_value <- unname(coefficient_table[".volcano_x_scaled", "Pr(>|t|)"])

      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = effect,
          EffectType = "Standardized beta",
          PValue = p_value,
          N = final_n,
          Note = NA_character_
        )
      )
    }

    model_data <- model_data %>%
      dplyr::mutate(
        .volcano_group = factor(as.character(.data[[yVar]]), levels = outcome_levels)
      )

    group_counts <- table(model_data$.volcano_group)

    if (length(group_counts) != 2 || any(group_counts < 2)) {
      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = NA_real_,
          EffectType = ifelse(
            EffectMetric == "log2fc",
            ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
            ifelse(!is.null(Covariates), "Adjusted Cohen-style d", "Cohen-style d")
          ),
          PValue = NA_real_,
          N = final_n,
          Note = "Too few observations in one or both outcome groups"
        )
      )
    }

    if (EffectMetric == "log2fc") {
      if (any(model_data[[this_var]] <= 0, na.rm = TRUE)) {
        return(
          tibble::tibble(
            Variable = this_var,
            Label = unname(label_lookup[[this_var]]),
            Outcome = yVar,
            OutcomeType = OutcomeType,
            Effect = NA_real_,
            EffectType = ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
            PValue = NA_real_,
            N = final_n,
            Note = "Log2 fold change requires strictly positive predictor values"
          )
        )
      }

      model_data <- model_data %>%
        dplyr::mutate(.volcano_x_effect = log2(.data[[this_var]]))
    } else {
      model_data <- model_data %>%
        dplyr::mutate(.volcano_x_effect = as.numeric(scale(.data[[this_var]])))
    }

    if (
      all(is.na(model_data$.volcano_x_effect)) ||
      stats::sd(model_data$.volcano_x_effect, na.rm = TRUE) == 0
    ) {
      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = NA_real_,
          EffectType = ifelse(
            EffectMetric == "log2fc",
            ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
            ifelse(!is.null(Covariates), "Adjusted Cohen-style d", "Cohen-style d")
          ),
          PValue = NA_real_,
          N = final_n,
          Note = "Zero variance in predictor"
        )
      )
    }

    model_formula <- stats::as.formula(
      paste(
        ".volcano_x_effect ~ .volcano_group",
        if (!is.null(Covariates)) {
          paste("+", paste(Covariates, collapse = " + "))
        } else {
          ""
        }
      )
    )

    model_fit <- tryCatch(
      stats::lm(model_formula, data = model_data),
      error = function(e) NULL
    )

    if (is.null(model_fit)) {
      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = NA_real_,
          EffectType = ifelse(
            EffectMetric == "log2fc",
            ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
            ifelse(!is.null(Covariates), "Adjusted Cohen-style d", "Cohen-style d")
          ),
          PValue = NA_real_,
          N = final_n,
          Note = "Model failed"
        )
      )
    }

    coefficient_table <- summary(model_fit)$coefficients
    group_coef <- grep("^\\.volcano_group", rownames(coefficient_table), value = TRUE)

    if (length(group_coef) != 1) {
      return(
        tibble::tibble(
          Variable = this_var,
          Label = unname(label_lookup[[this_var]]),
          Outcome = yVar,
          OutcomeType = OutcomeType,
          Effect = NA_real_,
          EffectType = ifelse(
            EffectMetric == "log2fc",
            ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
            ifelse(!is.null(Covariates), "Adjusted Cohen-style d", "Cohen-style d")
          ),
          PValue = NA_real_,
          N = final_n,
          Note = "Group coefficient not estimable"
        )
      )
    }

    effect <- unname(coefficient_table[group_coef, "Estimate"])
    p_value <- unname(coefficient_table[group_coef, "Pr(>|t|)"])

    tibble::tibble(
      Variable = this_var,
      Label = unname(label_lookup[[this_var]]),
      Outcome = yVar,
      OutcomeType = OutcomeType,
      Effect = effect,
      EffectType = ifelse(
        EffectMetric == "log2fc",
        ifelse(!is.null(Covariates), "Adjusted log2 fold change", "Log2 fold change"),
        ifelse(!is.null(Covariates), "Adjusted Cohen-style d", "Cohen-style d")
      ),
      PValue = p_value,
      N = final_n,
      Note = NA_character_
    )
  })

  # Build outputs

  results <- results %>%
    dplyr::mutate(
      FDR = stats::p.adjust(.data$PValue, method = AdjustMethod),
      NegLog10P = -log10(.data$PValue),
      NegLog10FDR = -log10(.data$FDR),
      Significant = !is.na(.data$PValue) & .data$PValue < Alpha,
      FDRSignificant = !is.na(.data$FDR) & .data$FDR < Alpha,
      Direction = dplyr::case_when(
        is.na(.data$Effect) ~ "Missing",
        .data$Effect > 0 ~ "Positive",
        .data$Effect < 0 ~ "Negative",
        TRUE ~ "Zero"
      ),
      SignificanceTier = dplyr::case_when(
        .data$FDRSignificant ~ "FDR significant",
        .data$Significant ~ "Raw p significant",
        TRUE ~ "Not significant"
      ),
      Tooltip = paste0(
        "<b>", .data$Label, "</b>",
        "<br>Variable: ", .data$Variable,
        "<br>Outcome: ", .data$Outcome,
        "<br>Effect type: ", .data$EffectType,
        "<br>Effect: ", round(.data$Effect, 3),
        "<br>P: ", signif(.data$PValue, 3),
        "<br>FDR: ", signif(.data$FDR, 3),
        "<br>N: ", .data$N
      )
    )

  label_vars <- character(0)

  if (LabelMode == "top_n") {
    label_vars <- results %>%
      dplyr::filter(!is.na(.data$PValue)) %>%
      dplyr::arrange(.data$PValue) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::pull(.data$Variable)
  }

  if (LabelMode == "significant") {
    label_vars <- results %>%
      dplyr::filter(.data$Significant) %>%
      dplyr::pull(.data$Variable)
  }

  if (LabelMode == "fdr") {
    label_vars <- results %>%
      dplyr::filter(.data$FDRSignificant) %>%
      dplyr::pull(.data$Variable)
  }

  if (LabelMode == "extreme") {
    top_p <- results %>%
      dplyr::filter(!is.na(.data$PValue)) %>%
      dplyr::arrange(.data$PValue) %>%
      dplyr::slice_head(n = ceiling(TopN / 2)) %>%
      dplyr::pull(.data$Variable)

    top_effect <- results %>%
      dplyr::filter(!is.na(.data$Effect)) %>%
      dplyr::arrange(dplyr::desc(abs(.data$Effect))) %>%
      dplyr::slice_head(n = floor(TopN / 2)) %>%
      dplyr::pull(.data$Variable)

    label_vars <- unique(c(top_p, top_effect))
  }

  results <- results %>%
    dplyr::mutate(
      PlotLabel = dplyr::if_else(
        .data$Variable %in% label_vars,
        .data$Label,
        NA_character_
      )
    )

  x_axis_label <- unique(results$EffectType)
  x_axis_label <- x_axis_label[!is.na(x_axis_label)]

  if (length(x_axis_label) != 1) {
    x_axis_label <- "Effect"
  }

  make_plot <- function(plot_data, y_col, y_label, plot_title) {

    plot_data <- plot_data %>%
      dplyr::mutate(
        PlotY = .data[[y_col]],
        PlotY = dplyr::if_else(is.infinite(.data$PlotY), NA_real_, .data$PlotY)
      )

    if (Format == "effect_gradient") {
      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$Effect,
          y = .data$PlotY,
          color = .data$Effect
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(text = if (InteractiveLabels) .data$Tooltip else NULL),
          alpha = 0.85,
          size = 2.4
        ) +
        ggplot2::scale_color_gradient2(
          low = "#3B4CC0",
          mid = "grey80",
          high = "#B40426",
          midpoint = 0,
          na.value = "grey80"
        )
    } else {
      color_var <- dplyr::case_when(
        Format %in% c("tiered", "classic", "neon") ~ plot_data$SignificanceTier,
        Format == "fdr_only" & plot_data$FDRSignificant ~ "FDR significant",
        Format == "fdr_only" ~ "Not significant",
        Format == "directional" & plot_data$FDRSignificant & plot_data$Effect > 0 ~ "FDR positive",
        Format == "directional" & plot_data$FDRSignificant & plot_data$Effect < 0 ~ "FDR negative",
        Format == "directional" & plot_data$Significant & plot_data$Effect > 0 ~ "Raw positive",
        Format == "directional" & plot_data$Significant & plot_data$Effect < 0 ~ "Raw negative",
        Format == "directional" ~ "Not significant",
        Format == "minimal" & plot_data$FDRSignificant ~ "FDR significant",
        Format == "minimal" & plot_data$Significant ~ "Raw p significant",
        Format == "minimal" ~ "Not significant",
        TRUE ~ plot_data$SignificanceTier
      )

      plot_data <- plot_data %>%
        dplyr::mutate(PlotColor = color_var)

      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$Effect,
          y = .data$PlotY,
          color = .data$PlotColor
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(text = if (InteractiveLabels) .data$Tooltip else NULL),
          alpha = 0.85,
          size = 2.4
        )

      if (Format == "tiered") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey75",
              "Raw p significant" = "#F9C74F",
              "FDR significant" = "#D62828"
            ),
            na.value = "grey75"
          )
      }

      if (Format == "classic") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey70",
              "Raw p significant" = "#FDB515",
              "FDR significant" = "#B00020"
            ),
            na.value = "grey70"
          )
      }

      if (Format == "fdr_only") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey78",
              "FDR significant" = "#D62828"
            ),
            na.value = "grey78"
          )
      }

      if (Format == "directional") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey78",
              "Raw positive" = "#F9C74F",
              "Raw negative" = "#90BE6D",
              "FDR positive" = "#D62828",
              "FDR negative" = "#277DA1"
            ),
            na.value = "grey78"
          )
      }

      if (Format == "minimal") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey55",
              "Raw p significant" = "grey20",
              "FDR significant" = "black"
            ),
            na.value = "grey70"
          )
      }

      if (Format == "neon") {
        p <- p +
          ggplot2::scale_color_manual(
            values = c(
              "Not significant" = "grey25",
              "Raw p significant" = "#F9F871",
              "FDR significant" = "#FF00A8"
            ),
            na.value = "grey40"
          )
      }
    }

    p <- p +
      ggplot2::geom_hline(
        yintercept = -log10(Alpha),
        linetype = "dashed",
        linewidth = 0.4,
        color = "grey45"
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "solid",
        linewidth = 0.35,
        color = "grey70"
      ) +
      ggplot2::labs(
        title = plot_title,
        x = x_axis_label,
        y = y_label,
        color = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "right",
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )

    label_data <- plot_data %>%
      dplyr::filter(!is.na(.data$PlotLabel), !is.na(.data$Effect), !is.na(.data$PlotY))

    if (nrow(label_data) > 0) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p +
          ggrepel::geom_text_repel(
            data = label_data,
            ggplot2::aes(label = .data$PlotLabel),
            size = 3.3,
            max.overlaps = Inf,
            show.legend = FALSE
          )
      } else {
        p <- p +
          ggplot2::geom_text(
            data = label_data,
            ggplot2::aes(label = .data$PlotLabel),
            size = 3,
            check_overlap = TRUE,
            vjust = -0.5,
            show.legend = FALSE
          )
      }
    }

    p
  }

  raw_p_plot <- make_plot(
    plot_data = results,
    y_col = "NegLog10P",
    y_label = expression(-log[10](italic(p))),
    plot_title = paste0("Volcano plot using raw p-values: ", yVar)
  )

  fdr_plot <- make_plot(
    plot_data = results,
    y_col = "NegLog10FDR",
    y_label = expression(-log[10]("FDR")),
    plot_title = paste0("Volcano plot using FDR-adjusted p-values: ", yVar)
  )

  # Return result

  list(
    RawPPlot = raw_p_plot,
    FDRPlot = fdr_plot,
    ResultsTable = results
  )
}
