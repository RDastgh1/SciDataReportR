#' Plot retained histories for longitudinal anomaly review
#'
#' Create ordinary ggplot objects without Plotly or a browser runtime. Highlighted
#' records are review candidates, not confirmed incorrect values. Full retained
#' histories are shown, with separate lines for strata.
#' @param AnomalyObj A result from [DetectLongitudinalAnomalies()].
#' @param participants Optional vector of participant IDs. NULL selects participants
#'   represented by the top-ranked evaluated records (native flags first).
#' @param variables Optional character vector of retained model variable names.
#'   NULL uses all model variables.
#' @param TopN Positive integer limiting automatic selection to the first TopN
#'   unique candidate records, not the number of visits in each plotted history.
#' @return A named list of ggplot objects, one per participant and variable. No
#'   plots are saved. Names use participant IDs and codebook-derived labels.
#' @family data anomalies
#' @seealso [DetectLongitudinalAnomalies()]
#' @export
PlotDataAnomalyTrajectories <- function(AnomalyObj, participants = NULL,
    variables = NULL, TopN = 10) {
  if (!inherits(AnomalyObj, "LongitudinalAnomalyResult")) stop("Provide a longitudinal anomaly result.")
  if (!is.numeric(TopN) || length(TopN) != 1L || !is.finite(TopN) || TopN < 1 || TopN != floor(TopN)) stop("TopN must be a positive integer.")
  if (is.null(variables)) variables <- AnomalyObj$Metadata$Variables
  if (!is.character(variables) || !length(variables) || anyNA(variables) || anyDuplicated(variables) || any(!variables %in% AnomalyObj$Metadata$Variables)) stop("variables must select retained model variables.")
  df_Data <- AnomalyObj$ReviewData
  id_var <- AnomalyObj$Metadata$IDVariable
  time_var <- AnomalyObj$Metadata$TimeVariable
  df_Results <- AnomalyObj$Results
  df_Results <- df_Results[df_Results$Status == "evaluated", , drop = FALSE]
  df_Results <- df_Results[order(!df_Results$Flag, df_Results$Rank, na.last = TRUE), , drop = FALSE]
  candidate_rows <- head(unique(df_Results$.Row), TopN)
  if (is.null(participants)) participants <- unique(df_Data[[id_var]][match(candidate_rows, df_Data$.Row)])
  if (anyNA(participants) || any(!participants %in% df_Data[[id_var]])) stop("participants must identify retained participants.")
  plots <- list()
  for (participant in unique(participants)) {
    df_History <- df_Data[!is.na(df_Data[[id_var]]) & df_Data[[id_var]] == participant & !is.na(df_Data[[time_var]]), , drop = FALSE]
    df_History <- df_History[order(df_History[[time_var]]), , drop = FALSE]
    vars_Strata <- AnomalyObj$Metadata$StrataVariables
    df_History$.Series <- if (length(vars_Strata)) as.character(interaction(df_History[vars_Strata], drop = TRUE)) else "History"
    df_History$.Review <- ifelse(df_History$.Row %in% candidate_rows, "Review record", "Other visit")
    for (variable in variables) {
      p <- ggplot2::ggplot(df_History, ggplot2::aes(x = .data[[time_var]], y = .data[[variable]], group = .data$.Series)) +
        ggplot2::geom_line(na.rm = TRUE, color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(color = .data$.Review), size = 2.5, na.rm = TRUE) +
        ggplot2::scale_color_manual(values = c("Review record" = "#D55E00", "Other visit" = "#0072B2"), name = NULL) +
        ggplot2::labs(x = unname(AnomalyObj$Metadata$Labels[time_var]),
          y = unname(AnomalyObj$Metadata$Labels[variable])) + ggplot2::theme_bw()
      plots[[paste(participant, AnomalyObj$Metadata$Labels[[variable]], sep = ": ")]] <- p
    }
  }
  plots
}
