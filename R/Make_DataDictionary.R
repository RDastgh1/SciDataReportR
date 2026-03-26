#' Create a data dictionary for a data frame
#'
#' This function generates a data dictionary for a given data frame, including
#' variable names, labels, types, missing value statistics, summary statistics,
#' and other relevant information. Use this when you want a stable, label-aware
#' summary table across mixed variable types.
#'
#' @param DataFrame A data frame.
#' @param numdecimals Number of decimals to display for numeric variables.
#' @return A tibble with one row per variable and stable summary columns.
#' @examples
#' df <- tibble::tibble(
#'   age = c(20, 25, 30),
#'   group = factor(c("A", "B", "A"))
#' )
#' Make_DataDictionary(df)
#' @export
Make_DataDictionary <- function(DataFrame, numdecimals = 2) {

  #Validate inputs

  if (!is.data.frame(DataFrame)) {
    stop("`DataFrame` must be a data frame.")
  }

  if (!is.numeric(numdecimals) || length(numdecimals) != 1 || is.na(numdecimals)) {
    stop("`numdecimals` must be a single numeric value.")
  }

 # Prepare data

  expected_skim_cols <- c(
    "skim_variable",
    "skim_type",
    "factor.ordered",
    "n_missing",
    "complete_rate",
    "factor.n_unique",
    "character.n_unique",
    "factor.top_counts",
    "numeric.mean",
    "numeric.median",
    "numeric.sd",
    "numeric.min",
    "numeric.max",
    "numeric.hist"
  )

  df_var <- sjlabelled::get_label(DataFrame)

  c_skim <- codebook::skim_codebook(DataFrame) %>%
    tibble::as_tibble()

  missing_cols <- setdiff(expected_skim_cols, colnames(c_skim))
  if (length(missing_cols) > 0) {
    c_skim[missing_cols] <- NA
  }

  c_skim <- c_skim %>%
    dplyr::rename(Variable = skim_variable) %>%
    dplyr::mutate(
      n_unique = dplyr::coalesce(factor.n_unique, character.n_unique)
    )

  cb <- tibble::tibble(
    Variable = colnames(DataFrame),
    Label = df_var
  ) %>%
    dplyr::left_join(
      c_skim %>%
        dplyr::select(
          dplyr::any_of(c(
            "Variable",
            "skim_type",
            "factor.ordered",
            "n_missing",
            "complete_rate",
            "n_unique",
            "factor.top_counts",
            "numeric.mean",
            "numeric.median",
            "numeric.sd",
            "numeric.min",
            "numeric.max",
            "numeric.hist"
          ))
        ),
      by = "Variable"
    )

  numeric_vars <- names(DataFrame)[vapply(DataFrame, is.numeric, logical(1))]

  if (length(numeric_vars) > 0) {
    cb <- cb %>%
      dplyr::mutate(
        n_unique = dplyr::if_else(
          Variable %in% numeric_vars,
          vapply(DataFrame[Variable], function(x) dplyr::n_distinct(x, na.rm = FALSE), numeric(1)),
          n_unique
        )
      )
  }

  #Apply labels

  is_num_col <- vapply(cb, is.numeric, logical(1))
  cb[is_num_col] <- lapply(cb[is_num_col], round, digits = numdecimals)

  cb <- cb %>%
    dplyr::rename(
      Type = skim_type,
      `Ordered Factor` = factor.ordered,
      `Top Counts` = factor.top_counts,
      Mean = numeric.mean,
      Median = numeric.median,
      SD = numeric.sd,
      Min = numeric.min,
      Max = numeric.max,
      Histogram = numeric.hist
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::everything(), ~ as.character(.x)),
      dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, " "))
    )

 # Build outputs

  cb <- cb %>%
    dplyr::select(
      Variable,
      Label,
      Type,
      `Ordered Factor`,
      n_missing,
      complete_rate,
      n_unique,
      `Top Counts`,
      Mean,
      Median,
      SD,
      Min,
      Max,
      Histogram
    )

 # Return result

  cb
}
