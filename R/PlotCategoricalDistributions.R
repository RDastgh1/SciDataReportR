#' Plot categorical distributions
#'
#' This function creates plots to visualize the distributions of categorical
#' variables in a dataframe.
#'
#' @param DataFrame The dataframe containing the variables to be plotted.
#' @param Variables Optional. A character vector specifying the names of the
#'   categorical variables to be plotted. If NULL, categorical variables are
#'   automatically detected.
#' @param Relabel Logical. If TRUE, missing labels in the dataframe are replaced
#'   with column names as labels for plotting.
#' @param Ordinal Logical, indicating whether ordinal variables should be included.
#' @param LabelType Character. Either "percent" or "count", indicating what
#'   should be shown on the x-axis and inside the bars.
#' @param MissingLabel Character label to use for missing values.
#'
#' @return A ggplot object visualizing the distributions of categorical variables.
#' @importFrom sjlabelled get_label
#' @importFrom magrittr %>%
#' @export
PlotCategoricalDistributions <- function(DataFrame,
                                         Variables = NULL,
                                         Relabel = TRUE,
                                         Ordinal = TRUE,
                                         LabelType = "percent",
                                         MissingLabel = "Missing") {

  # Validate inputs

  if (!is.data.frame(DataFrame)) {
    stop("`DataFrame` must be a data frame.")
  }

  if (is.null(Variables)) {
    Variables <- getCatVars(DataFrame)
  }

  if (!is.character(Variables)) {
    stop("`Variables` must be a character vector.")
  }

  missing_vars <- setdiff(Variables, colnames(DataFrame))
  if (length(missing_vars) > 0) {
    stop("Variables not found in DataFrame: ", paste(missing_vars, collapse = ", "))
  }

  if (!is.logical(Relabel) || length(Relabel) != 1 || is.na(Relabel)) {
    stop("`Relabel` must be TRUE or FALSE.")
  }

  if (!is.logical(Ordinal) || length(Ordinal) != 1 || is.na(Ordinal)) {
    stop("`Ordinal` must be TRUE or FALSE.")
  }

  if (!LabelType %in% c("percent", "count")) {
    stop("`LabelType` must be either 'percent' or 'count'.")
  }

  if (!is.character(MissingLabel) || length(MissingLabel) != 1 || is.na(MissingLabel)) {
    stop("`MissingLabel` must be a single character value.")
  }

  # Prepare data

  if (!Ordinal) {
    ordinal_vars <- Variables[vapply(DataFrame[Variables], is.ordered, logical(1))]
    Variables <- setdiff(Variables, ordinal_vars)
  }

  if (length(Variables) == 0) {
    stop("No variables available to plot.")
  }

  if (Relabel) {
    DataFrame <- ReplaceMissingLabels(DataFrame)

    facetlabels <- sjlabelled::get_label(
      DataFrame %>% dplyr::select(dplyr::all_of(Variables))
    )

    names(facetlabels) <- NULL
    facetlabels[is.na(facetlabels) | facetlabels == ""] <- Variables[is.na(facetlabels) | facetlabels == ""]
  } else {
    facetlabels <- Variables
  }

  plot_df <- purrr::map_dfr(Variables, function(var) {
    x <- DataFrame[[var]]

    if (is.factor(x) || is.integer(x) || is.logical(x) ||
        inherits(x, "Date") || inherits(x, "POSIXt")) {
      x <- as.character(x)
    }

    tibble::tibble(
      VariableRaw = var,
      Level = x
    ) %>%
      dplyr::mutate(
        Level = as.character(Level),
        Level = tidyr::replace_na(Level, MissingLabel)
      ) %>%
      dplyr::count(VariableRaw, Level, name = "n") %>%
      dplyr::group_by(VariableRaw) %>%
      dplyr::mutate(
        prop = n / sum(n)
      ) %>%
      dplyr::ungroup()
  })

  # Order levels within each variable so Missing is semantically last
  plot_df <- plot_df %>%
    dplyr::mutate(
      IsMissing = Level == MissingLabel
    ) %>%
    dplyr::group_by(VariableRaw) %>%
    dplyr::arrange(IsMissing, dplyr::desc(prop), .by_group = TRUE) %>%
    dplyr::mutate(
      stack_id = dplyr::row_number()
    ) %>%
    dplyr::ungroup()

  # Create unique fill keys so each variable can have its own within-bar order
  plot_df <- plot_df %>%
    dplyr::mutate(
      Variable = factor(VariableRaw, levels = Variables, labels = facetlabels),
      FillKey = paste0(VariableRaw, "___", Level)
    )

  fill_key_order <- plot_df %>%
    dplyr::arrange(Variable, stack_id) %>%
    dplyr::pull(FillKey)

  plot_df <- plot_df %>%
    dplyr::mutate(
      FillKey = factor(FillKey, levels = fill_key_order),
      Value = if (LabelType == "percent") prop else n,
      Label = dplyr::case_when(
        LabelType == "percent" & prop >= 0.04 ~ paste0(Level, "\n", scales::percent(prop, accuracy = 1)),
        LabelType == "count" & prop >= 0.04 ~ paste0(Level, "\n", n),
        TRUE ~ ""
      ),
      Tooltip = paste0(
        "Variable: ", as.character(Variable),
        "<br>Level: ", Level,
        "<br>Count: ", n,
        "<br>Percent: ", scales::percent(prop, accuracy = 0.1)
      )
    )

  # Build fill palette, with Missing always grey
  non_missing_keys <- plot_df %>%
    dplyr::filter(!IsMissing) %>%
    dplyr::pull(FillKey)

  non_missing_keys <- unique(as.character(non_missing_keys))

  fill_values <- stats::setNames(
    scales::hue_pal()(length(non_missing_keys)),
    non_missing_keys
  )

  missing_keys <- plot_df %>%
    dplyr::filter(IsMissing) %>%
    dplyr::pull(FillKey)

  missing_keys <- unique(as.character(missing_keys))

  if (length(missing_keys) > 0) {
    fill_values <- c(
      fill_values,
      stats::setNames(rep("grey70", length(missing_keys)), missing_keys)
    )
  }

  x_lab <- if (LabelType == "percent") "Percent" else "Count"

  # Build outputs

  d <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = Variable,
      y = Value,
      fill = FillKey,
      text = Tooltip
    )
  ) +
    ggplot2::geom_col(
      position = if (LabelType == "percent") {
        ggplot2::position_fill(reverse = TRUE)
      } else {
        ggplot2::position_stack(reverse = TRUE)
      },
      width = 0.8,
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Label),
      position = if (LabelType == "percent") {
        ggplot2::position_fill(vjust = 0.5, reverse = TRUE)
      } else {
        ggplot2::position_stack(vjust = 0.5, reverse = TRUE)
      },
      size = 3,
      lineheight = 0.9
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = x_lab,
      fill = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (LabelType == "percent") {
    d <- d + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    )
  }

  # Return result

  return(d)
}
