#' Heatmap of binary feature prevalence by cluster
#'
#' @description
#' Shows, for each cluster, the percentage of participants with each binary
#' feature (for example, a marker beyond a clinical cutoff, a diagnosis, or a
#' medication), alongside the whole cohort for reference. Use it to describe
#' phenotypes in clinical terms ("share with albuminuria") rather than by
#' means or z-scores, or to compare outcome prevalence across phenotypes.
#'
#' @param data A data frame.
#' @param ClusterVar Name of the cluster column in `data`.
#' @param variables Binary variables to show (factor, character, logical, or
#'   0/1 numeric).
#' @param positive Value counted as present for factor and character
#'   variables. Default `"Yes"`. Logical `TRUE` and numeric `1` are always
#'   counted as present.
#' @param Relabel Logical; if `TRUE` (default), rows use variable labels.
#' @param show_cohort Logical; if `TRUE` (default), add a column for all
#'   participants with a cluster.
#' @param digits Decimal places for the cell percentages. Default `0`.
#' @param title Plot title.
#'
#' @details
#' Percentages use participants with a non-missing value for that feature
#' within each cluster, so features measured in only some participants or
#' cycles are summarized over those who have them. Cells are shaded by
#' prevalence; the cell text gives the percentage and, in parentheses, the
#' number of participants with a value.
#'
#' @return A ggplot object. The underlying summary is attached as the
#'   `"prevalence"` attribute.
#'
#' @examples
#' df <- data.frame(
#'   Cluster = rep(c("A", "B"), each = 50),
#'   high_bp = rep(c("Yes", "No", "Yes", "No"), times = c(40, 10, 5, 45)),
#'   anemia = rep(c(TRUE, FALSE, TRUE, FALSE), times = c(5, 45, 20, 30))
#' )
#' attr(df$high_bp, "label") <- "Systolic BP >= 130 mmHg"
#' attr(df$anemia, "label") <- "Anemia"
#' PlotClusterPrevalence(df, "Cluster", c("high_bp", "anemia"))
#' @details
#' Categorical input safeguards apply to `ClusterVar`.
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
PlotClusterPrevalence <- function(data, ClusterVar, variables, positive = "Yes",
    Relabel = TRUE, show_cohort = TRUE, digits = 0,
    title = "Prevalence by cluster") {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  missing_vars <- setdiff(c(ClusterVar, variables), names(data))
  if (length(missing_vars)) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  data <- ScidrPrepareCategorical(data, ClusterVar, "ClusterVar")

  IsPresent <- function(x) {
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(ifelse(is.na(x), NA, x == 1))
    x <- as.character(x)
    ifelse(is.na(x), NA, x == positive)
  }
  labels <- vapply(variables, function(v) {
    label <- attr(data[[v]], "label", exact = TRUE)
    if (Relabel && !is.null(label) && nzchar(label)) as.character(label)[1] else v
  }, character(1))

  clustered <- data[!is.na(data[[ClusterVar]]), , drop = FALSE]
  cluster <- as.character(clustered[[ClusterVar]])
  cluster_levels <- if (is.factor(clustered[[ClusterVar]])) {
    levels(droplevels(clustered[[ClusterVar]]))
  } else {
    sort(unique(cluster))
  }

  summary <- dplyr::bind_rows(lapply(variables, function(v) {
    present <- IsPresent(clustered[[v]])
    groups <- c(if (show_cohort) "All", cluster_levels)
    dplyr::bind_rows(lapply(groups, function(g) {
      keep <- if (g == "All") rep(TRUE, length(present)) else cluster == g
      observed <- present[keep & !is.na(present)]
      dplyr::tibble(Variable = v, Label = labels[[v]], Cluster = g,
        n = length(observed),
        Prevalence = if (length(observed)) mean(observed) else NA_real_)
    }))
  }))
  summary$Cluster <- factor(summary$Cluster,
    levels = c(if (show_cohort) "All", cluster_levels))
  summary$Label <- factor(summary$Label, levels = rev(unique(labels)))
  summary$Text <- ifelse(is.na(summary$Prevalence), "",
    paste0(formatC(100 * summary$Prevalence, format = "f", digits = digits),
      "%\n(", summary$n, ")"))

  low_high <- unname(.SciDataColorValues(1))
  plot <- ggplot2::ggplot(summary, ggplot2::aes(x = .data$Cluster, y = .data$Label,
      fill = 100 * .data$Prevalence)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Text,
      colour = .data$Prevalence > 0.55), size = 3, lineheight = 0.9,
      show.legend = FALSE) +
    ggplot2::scale_fill_gradient(low = "white", high = low_high, limits = c(0, 100),
      na.value = "grey92", name = "% of cluster") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"),
      na.value = "black") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = NULL, y = NULL)
  if (show_cohort) {
    plot <- plot + ggplot2::geom_vline(xintercept = 1.5, colour = "grey40")
  }
  attr(plot, "prevalence") <- summary
  plot
}
