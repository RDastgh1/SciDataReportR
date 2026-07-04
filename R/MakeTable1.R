#' Create Summary Table using gtsummary
#'
#' This function is a wrapper around `gtsummary::tbl_summary` that ensures continuous variables are treated as continuous.
#'
#' @param data The dataframe to create the summary table from.
#' @param variables Optional. A character vector specifying the names of variables to include in the summary table. If NULL, all variables are included.
#' @param TreatOrdinalAs Character. Specifies how ordinal variables should be treated. Can be "Continuous", "Categorical", or "Both".
#' @param AutoDetectDistribution Logical. If TRUE, the function will attempt to automatically detect the distribution of variables. Default is FALSE.
#' @param IncludeMissing Character matching gtsummary criteria. Can be "no", "ifany", or "always". Default is "ifany"
#' @return A summary table created using gtsummary.
#' @importFrom gtsummary tbl_summary
#' @importFrom dplyr select all_of
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @export
MakeTable1 <- function(data,
    variables = NULL,
    TreatOrdinalAs = "Continuous",
    AutoDetectDistribution = FALSE,
    IncludeMissing = "ifany",
    DataFrame = lifecycle::deprecated(),
    Variables = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "MakeTable1(DataFrame)", "MakeTable1(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "MakeTable1(Variables)", "MakeTable1(variables)")
    variables <- Variables
  }
  Variables <- variables

  if (is.null(Variables)) {
    Variables <- colnames(DataFrame)
  }

  # Treat ordinal variables as specified
  if (TreatOrdinalAs %in% c("Continuous", "Both")) {
    Table1 <- DataFrame %>%
      select(all_of(Variables)) %>%
      gtsummary::tbl_summary(
        type = list(
          where(is.numeric) ~ "continuous"
        ),
        statistic = list(
          gtsummary::all_continuous() ~ "{mean} ({sd})"
        ),
        missing = IncludeMissing
      )
  } else {
    Table1 <- DataFrame %>%
      select(all_of(Variables)) %>%
      gtsummary::tbl_summary(
        type = list(
          where(is.numeric) ~ "continuous"
        ),
        statistic = list(
          gtsummary::all_continuous() ~ "{mean} ({sd})"
        ),
        missing = IncludeMissing
      )
  }

  return(Table1)
}
