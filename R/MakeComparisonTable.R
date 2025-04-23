#' Make a Group Comparison Table (gtsummary wrapper)
#'
#' @description
#' This function creates a summary comparison table from a given data frame. It summarizes continuous variables with their mean (and standard deviation) and categorical variables with counts and percentages. Comparisons between groups are performed based on a specified grouping (comparison) variable.
#'
#' @param DataFrame A data frame containing the data for generating the comparison table.
#' @param Variables An optional character vector of variable names (besides the grouping variable) to include in the comparison table. If \code{NULL}, all variables in \code{DataFrame} are used.
#' @param CompVariable A character string specifying the name of the variable used to divide the data into comparison groups.
#' @param ValueDigits An integer specifying the number of digits to display for continuous variable statistics (mean and standard deviation). Default is 2.
#' @param pDigits An integer specifying the number of digits to display for p-values. Default is 3.
#'
#' @return A \code{tbl_summary} object containing the comparison table.
#'
#' @details
#' The function first determines the set of variables to analyze. If \code{Variables} is \code{NULL}, all variables from \code{DataFrame} are used; otherwise, it ensures that \code{CompVariable} is included once in the analysis. Next, it filters the data frame to include only the selected variables and excludes any factor variables with more than 15 unique levels. The summary is then created using \code{\link[gtsummary]{tbl_summary}} with continuous variables summarized as "mean (sd)" and categorical variables summarized as "n (p)". Additional information, such as sample sizes and p-values, is appended to the table.
#' @note This wrapper function is adapted from code written by Aparna Bhattacharyya.
#' @export
MakeComparisonTable <- function(DataFrame, Variables = NULL, CompVariable, ValueDigits = 2, pDigits = 3) {

  # If variables of interest are not provided then use all variables from provided data frame
  if (is.null(Variables)) {
    Variables <- colnames(DataFrame)
  } else {
    # Ensure the comparison variable is included once
    Variables <- unique(c(CompVariable, Variables))
  }

  # Filter the data frame to only include the selected variables.
  # Exclude factor variables that have more than 15 unique levels.
  Data_filtered <- DataFrame %>%
    dplyr::select(dplyr::all_of(Variables)) %>%
    dplyr::select_if(function(col) !(is.factor(col) && nlevels(col) > 15))

  # Create the comparison table using gtsummary functions and provided arguments.
  CompTable <- gtsummary::tbl_summary(
    data = Data_filtered,
    by = CompVariable,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(
      all_continuous() ~ ValueDigits
    ),
    type = list(
      where(is.double) ~ "continuous",
      where(~ is.factor(.) & nlevels(.) == 2) ~ "categorical"
    )
  ) %>%
    gtsummary::add_n() %>%
    gtsummary::add_p() %>%
    gtsummary::bold_p() %>%
    gtsummary::modify_fmt_fun(
      p.value ~ function(x) gtsummary::style_pvalue(x, digits = pDigits)
    )

  return(CompTable)
}
