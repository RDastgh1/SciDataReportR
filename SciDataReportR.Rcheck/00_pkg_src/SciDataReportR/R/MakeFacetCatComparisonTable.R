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
#' @param data A data frame containing all variables to be analyzed.
#' @param FacetVariables A character vector of variable names to facet by.
#'   The function automatically selects those that are categorical (`factor` or `character`).
#' @param variables A character string naming the variable(s) being compared (e.g., "Cluster").
#' @param covariates Optional character vector of covariate names to adjust for.
#' @param value_digits Number of decimal digits to display for numeric values (default = 2).
#' @param p_digits Number of decimal digits to display for p-values (default = 3).
#' @param AddEffectSize Logical; if TRUE, include effect sizes (default = FALSE).
#' @param effect_size_digits Decimal digits for effect size values (default = 2).
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
#' @param Relabel Logical; if TRUE (default), use attached variable labels.
#' @param TreatOrdinalAs How ordinal variables are treated in each table.
#' @param ... Additional arguments passed to the comparison function.
#'
#' @return A `gtsummary` table created by merging each facet's `MakeComparisonTable()` output
#'   side-by-side using `gtsummary::tbl_merge()`. Each facet variable is labeled with its own
#'   tab spanner header for clarity.
#'
#' @details
#' Use this when the same set of measures needs to be compared across several
#' different groupings at once - by diagnosis, and by sex, and by APOE status -
#' and those comparisons belong in one table rather than three.
#'
#' Running [MakeComparisonTable()] once per grouping produces tables that are
#' each correct but cannot be read against each other: the variables are listed
#' three times, and spotting that a biomarker separates diagnosis groups but not
#' sexes means looking back and forth between pages. Merging them puts one row
#' per variable and one column block per grouping, so that comparison is a
#' left-to-right read.
#'
#' Facet variables that are not categorical are dropped automatically, so a
#' vector of candidate groupings can be passed without pre-filtering it.
#'
#' Every argument that controls the individual tables - `AddPairwise`,
#' `covariates`, `Parametric`, and the rest - is passed through to each facet,
#' so the blocks stay consistent with one another.
#'
#' @seealso [MakeComparisonTable()] for a single grouping, and [MakeTable1()]
#'   for a plain descriptive table with no grouping at all.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' vars_Compare <- c("age", "AXL", "Adiponectin", "tau", "p_tau")
#'
#' # The same five measures compared across two groupings at once
#' MakeFacetCatComparisonTable(
#'   data = Labelled,
#'   FacetVariables = c("Diagnosis", "sex"),
#'   variables = vars_Compare
#' )
#'
#' # Non-categorical facets are ignored, so a mixed vector is safe
#' MakeFacetCatComparisonTable(
#'   data = Labelled,
#'   FacetVariables = c("Diagnosis", "age", "sex"),
#'   variables = c("AXL", "tau")
#' )
#'
#' # Options are forwarded to every facet
#' MakeFacetCatComparisonTable(
#'   data = Labelled,
#'   FacetVariables = c("Diagnosis", "sex"),
#'   variables = vars_Compare,
#'   AddEffectSize = TRUE
#' )
#'
#' # Covariate adjustment, applied within each grouping
#' MakeFacetCatComparisonTable(
#'   data = Labelled,
#'   FacetVariables = c("Diagnosis", "sex"),
#'   variables = c("AXL", "Adiponectin", "tau", "p_tau"),
#'   covariates = "age"
#' )
#'
#' # A three-level grouping alongside a two-level one, with pairwise contrasts
#' Labelled$ApoeGroup <- factor(
#'   dplyr::case_when(
#'     Labelled$Genotype %in% c("E2E2", "E2E3") ~ "E2 carrier",
#'     Labelled$Genotype == "E3E3" ~ "E3/E3",
#'     TRUE ~ "E4 carrier"
#'   ),
#'   levels = c("E2 carrier", "E3/E3", "E4 carrier")
#' )
#'
#' MakeFacetCatComparisonTable(
#'   data = Labelled,
#'   FacetVariables = c("Diagnosis", "ApoeGroup"),
#'   variables = c("age", "tau", "p_tau"),
#'   AddPairwise = TRUE
#' )
#' }
#'
#' @importFrom gtsummary tbl_merge
#' @param DataFrame \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param Variables \strong{Deprecated} (since 19.15.0). Use \code{variables} instead.
#' @param Covariates \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @param ValueDigits \strong{Deprecated} (since 19.15.0). Use \code{value_digits} instead.
#' @param pDigits \strong{Deprecated} (since 19.15.0). Use \code{p_digits} instead.
#' @param EffectSizeDigits \strong{Deprecated} (since 19.15.0). Use \code{effect_size_digits} instead.
#' @export
MakeFacetCatComparisonTable <- function(data,
    FacetVariables,
    variables,
    covariates = NULL,
    value_digits = 2,
    p_digits = 3,
    AddEffectSize = FALSE,
    effect_size_digits = 2,
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
    Relabel = TRUE,
    TreatOrdinalAs = "Categorical",
    ...,
    DataFrame = lifecycle::deprecated(),
    Variables = lifecycle::deprecated(),
    Covariates = lifecycle::deprecated(),
    ValueDigits = lifecycle::deprecated(),
    pDigits = lifecycle::deprecated(),
    EffectSizeDigits = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DataFrame)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(DataFrame)", "MakeFacetCatComparisonTable(data)")
    data <- DataFrame
  }
  if (!missing(data)) DataFrame <- data
  if (lifecycle::is_present(Variables)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(Variables)", "MakeFacetCatComparisonTable(variables)")
    variables <- Variables
  }
  if (!missing(variables)) Variables <- variables
  if (lifecycle::is_present(Covariates)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(Covariates)", "MakeFacetCatComparisonTable(covariates)")
    covariates <- Covariates
  }
  Covariates <- covariates
  if (lifecycle::is_present(ValueDigits)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(ValueDigits)", "MakeFacetCatComparisonTable(value_digits)")
    value_digits <- ValueDigits
  }
  ValueDigits <- value_digits
  if (lifecycle::is_present(pDigits)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(pDigits)", "MakeFacetCatComparisonTable(p_digits)")
    p_digits <- pDigits
  }
  pDigits <- p_digits
  if (lifecycle::is_present(EffectSizeDigits)) {
    lifecycle::deprecate_warn("19.15.0", "MakeFacetCatComparisonTable(EffectSizeDigits)", "MakeFacetCatComparisonTable(effect_size_digits)")
    effect_size_digits <- EffectSizeDigits
  }
  EffectSizeDigits <- effect_size_digits

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
      data = DataFrame,
      group_var = facet_var,
      variables = Variables,
      covariates = Covariates,
      value_digits = ValueDigits,
      p_digits = pDigits,
      AddEffectSize = AddEffectSize,
      effect_size_digits = EffectSizeDigits,
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
      Relabel = Relabel,
      TreatOrdinalAs = TreatOrdinalAs,
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
