#' Create Summary Table using gtsummary
#'
#' This function is a wrapper around `gtsummary::tbl_summary` that ensures continuous variables are treated as continuous.
#'
#' @param DataFrame The dataframe to create the summary table from.
#' @param Variables Optional. A character vector specifying the names of variables to include in the summary table. If NULL, all variables are included.
#' @param TreatOrdinalAs Character. Specifies how ordinal variables should be treated. Can be "Continuous", "Categorical", or "Both".
#' @param AutoDetectDistribution Logical. If TRUE, the function will attempt to automatically detect the distribution of variables. Default is FALSE.
#' @param IncludeMissing Character matching gtsummary criteria. Can be "no", "ifany", or "always". Default is "ifany"
#' @return A summary table created using gtsummary.
#' @importFrom gtsummary tbl_summary
#' @importFrom dplyr select all_of
#' @export
MakeTable1 <- function(DataFrame, Variables = NULL, TreatOrdinalAs = "Continuous", AutoDetectDistribution = FALSE, IncludeMissing = "ifany"){
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
