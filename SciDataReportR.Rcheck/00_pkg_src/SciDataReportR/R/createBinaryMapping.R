#' Create a Mapping Table for Binary Variables
#'
#' Identifies binary variables and returns a deterministic mapping for 0/1 coding.
#' - Factors with \emph{explicit order} (ordered = TRUE) do NOT use heuristics; the highest (last) level is Positive.
#' - Logicals map to Negative = "FALSE", Positive = "TRUE".
#' - Numeric 0/1 (or any 2-value numeric) maps Positive to the numeric maximum.
#' - Characters / unordered factors use minimal heuristics (no race/PWH/sex terms).
#'
#' @param data A dataframe.
#' @param CatVars Character vector of candidate binary variables.
#' @param prefer Optional named character vector of explicit positive levels,
#'   e.g., c(STATUS = "PWH", Smoker = "Yes"). This overrides other rules.
#' @return A data.frame with columns: Variable, Label, PositiveLevel, NegativeLevel.
#' @details
#' Modelling a two-level variable means scoring one level against the other,
#' and which level counts as "positive" decides the sign of every coefficient
#' and odds ratio that follows. This resolves that choice once, explicitly, and
#' returns it as a table so it can be checked rather than assumed.
#'
#' The rules are applied in order: an explicit `prefer` entry always wins; an
#' ordered factor uses its highest level; otherwise a short list of
#' conventional affirmative labels ("Yes", "Present", "Case", and similar) is
#' consulted, then the larger of two numbers, and finally the second level in
#' sorted order. The heuristics deliberately exclude race, sex, and serostatus
#' terms, because there is no defensible default "positive" level for those -
#' name them through `prefer`.
#'
#' @seealso [getBinaryVars()] to find the candidates.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#' vars_Binary <- getBinaryVars(Labelled)
#'
#' # One row per variable, naming the level scored as positive
#' mapping <- createBinaryMapping(Labelled, vars_Binary)
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(mapping, full_width = TRUE)
#' )))
#'
#' # `prefer` overrides the choice, which is how to set a direction the
#' # heuristics will not guess at - `sex` is exactly that case.
#' createBinaryMapping(
#'   Labelled, vars_Binary,
#'   prefer = c(Diagnosis = "Control", sex = "Female")
#' )
#'
#' # The rules on constructed variables: an ordered factor takes its highest
#' # level, a logical takes TRUE, a 0/1 numeric takes the larger number, and a
#' # conventional affirmative label is recognised.
#' df_Rules <- data.frame(
#'   Severity = factor(c("Mild", "Severe"), levels = c("Mild", "Severe"),
#'                     ordered = TRUE),
#'   Responded = c(TRUE, FALSE),
#'   Coded01 = c(0, 1),
#'   Smoker = c("Yes", "No")
#' )
#'
#' createBinaryMapping(
#'   df_Rules, c("Severity", "Responded", "Coded01", "Smoker")
#' )
#' }
#'
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @export
createBinaryMapping <- function(data,
    CatVars,
    prefer = NULL,
    Data = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "createBinaryMapping(Data)", "createBinaryMapping(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data


  # choose PositiveLevel for var given observed levels (character)
  choose_positive <- function(var, lvls_chr, is_ordered_factor = FALSE) {
    # 0) explicit user preference wins, even for ordered factors
    if (!is.null(prefer) && !is.null(prefer[[var]]) && prefer[[var]] %in% lvls_chr) {
      return(prefer[[var]])
    }

    # 1) ordered factor: take the highest level, no heuristics
    if (is_ordered_factor) {
      return(tail(lvls_chr, 1))
    }

    # 2) minimal heuristics (exclude race/PWH/sex)
    heuristics <- c("Yes", "TRUE", "True", "1", "Present", "Current",
                    "Ever", "High", "Positive", "Case", "Detected", "Success")
    hit <- heuristics[heuristics %in% lvls_chr]
    if (length(hit)) return(hit[1])

    # 3) numeric pair? choose the larger number
    num_try <- suppressWarnings(as.numeric(lvls_chr))
    if (!anyNA(num_try)) return(as.character(max(num_try)))

    # 4) stable fallback: second element in sorted order
    if (length(lvls_chr) == 2) return(sort(lvls_chr)[2])

    NA_character_
  }

  rows <- lapply(CatVars, function(var) {
    x <- Data[[var]]
    if (is.null(x)) return(NULL)

    # drop labels but keep values
    if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)

    # Collect two unique levels as characters, preserving appropriate order
    is_ord <- is.factor(x) && is.ordered(x)

    if (is.logical(x)) {
      lv_chr <- c("FALSE", "TRUE")
    } else if (is.factor(x)) {
      # factor keeps its defined level order
      lv_chr <- levels(x)
    } else {
      # character/numeric: unique non-NA values (convert to character)
      ux <- unique(stats::na.omit(x))
      # For numeric, sort to be stable; for character, sort too (stable fallback)
      if (is.numeric(ux)) ux <- sort(ux)
      lv_chr <- as.character(ux)
    }

    # ensure unique and exactly two levels
    lv_chr <- unique(lv_chr)
    if (length(lv_chr) != 2) return(NULL)

    # variable label or default to var name
    label <- sjlabelled::get_label(Data[[var]], def.value = var)

    pos <- choose_positive(var, lv_chr, is_ordered_factor = is_ord)
    if (is.na(pos)) return(NULL)
    neg <- setdiff(lv_chr, pos)[1]

    data.frame(
      Variable      = var,
      Label         = label,
      PositiveLevel = pos,
      NegativeLevel = neg,
      stringsAsFactors = FALSE
    )
  })

  mapping <- do.call(rbind, rows)
  if (is.null(mapping) || nrow(mapping) == 0) {
    stop("No binary variables found in the provided CatVars.")
  }

  rownames(mapping) <- NULL
  mapping
}
