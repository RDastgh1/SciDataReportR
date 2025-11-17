#' MakeFacetCatComparisonTable
#'
#' @title Create a merged gtsummary table by faceting comparisons across multiple categorical variables
#'
#' @description
#' Generates a series of comparison tables using `MakeComparisonTable()` for each categorical
#' variable (facet) in the provided list and merges them side-by-side using `gtsummary::tbl_merge()`.
#' This function extends the functionality of `MakeComparisonTable()` by automatically detecting
#' which facet variables are categorical (factor or character) and producing a faceted summary
#' of how the main comparison variable (e.g., Cluster, TreatmentArm) differs across multiple
#' categorical dimensions such as Race, Sex, or HIV status.
#'
#' @param DataFrame A data frame containing all variables to be analyzed.
#' @param FacetVariables A character vector of variable names to facet by.
#'   The function automatically selects those that are categorical (`factor` or `character`).
#' @param Variables A character string naming the variable(s) being compared (e.g., "Cluster").
#' @param Covariates Optional character vector of covariate names to adjust for.
#' @param ValueDigits Number of decimal digits to display for numeric values (default = 2).
#' @param pDigits Number of decimal digits to display for p-values (default = 3).
#' @param AddEffectSize Logical; if TRUE, include effect sizes (default = FALSE).
#' @param EffectSizeDigits Decimal digits for effect size values (default = 2).
#' @param AddPairwise Logical; if TRUE, include pairwise comparisons (default = FALSE).
#' @param PairwiseMethod Method for pairwise comparison p-value adjustment (default = "bonferroni").
#' @param Parametric Logical; if TRUE, use parametric tests (default = TRUE).
#' @param ParametricDisplay Optional vector specifying which statistics to display for parametric tests.
#' @param IncludeOverallN Logical; if TRUE, adds overall N to the table (default = FALSE).
#' @param IncludeMissing Logical; if TRUE, includes missing categories (default = FALSE).
#' @param suppress_warnings Logical; suppress internal warnings (default = FALSE).
#' @param Referent Optional string specifying the referent category for binary or categorical comparisons.
#' @param IncludeOverallStats Logical; if TRUE, adds overall descriptive statistics (default = FALSE).
#' @param ShowPositiveBinaryOnLabel Logical; if TRUE, labels binary variables with positive outcome (default = TRUE).
#' @param CompFun Comparison function to apply; defaults to `MakeComparisonTable`.
#' @param ... Additional arguments passed to the comparison function.
#'
#' @return A `gtsummary` table created by merging each facet’s `MakeComparisonTable()` output
#'   side-by-side using `gtsummary::tbl_merge()`. Each facet variable is labeled with its own
#'   tab spanner header for clarity.
#' @importFrom gtsummary tbl_merge
#' @export
MakeFacetCatComparisonTable <- function(
    DataFrame,
    FacetVariables,
    Variables,
    Covariates = NULL,
    ValueDigits = 2,
    pDigits = 3,
    AddEffectSize = FALSE,
    EffectSizeDigits = 2,
    AddPairwise = FALSE,
    PairwiseMethod = "bonferroni",
    Parametric = TRUE,
    ParametricDisplay = NULL,
    IncludeOverallN = FALSE,
    IncludeMissing = FALSE,
    suppress_warnings = FALSE,
    Referent = NULL,
    IncludeOverallStats = FALSE,
    ShowPositiveBinaryOnLabel = TRUE,
    CompFun = MakeComparisonTable,
    ...
) {
  # Identify categorical facet variables
  FacetVariables <- FacetVariables[
    sapply(DataFrame[FacetVariables], function(x) is.character(x) || is.factor(x))
  ]

  if (length(FacetVariables) == 0) {
    stop("No categorical variables found in FacetVariables. Please provide factor or character columns.")
  }

  # Generate comparison subtables for each facet
  subtables <- lapply(FacetVariables, function(facet_var) {
    CompFun(
      DataFrame = DataFrame,
      CompVariable = facet_var,
      Variables = Variables,
      Covariates = Covariates,
      ValueDigits = ValueDigits,
      pDigits = pDigits,
      AddEffectSize = AddEffectSize,
      EffectSizeDigits = EffectSizeDigits,
      AddPairwise = AddPairwise,
      PairwiseMethod = PairwiseMethod,
      Parametric = Parametric,
      ParametricDisplay = ParametricDisplay,
      IncludeOverallN = IncludeOverallN,
      IncludeMissing = IncludeMissing,
      suppress_warnings = suppress_warnings,
      Referent = Referent,
      IncludeOverallStats = IncludeOverallStats,
      ShowPositiveBinaryOnLabel = ShowPositiveBinaryOnLabel,
      ...
    )
  })

  # Merge all subtables side-by-side with facet labels
  merged_tbl <- gtsummary::tbl_merge(
    tbls = subtables,
    tab_spanner = FacetVariables
  )

  return(merged_tbl)
}
