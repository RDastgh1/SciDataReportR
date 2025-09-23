#' Create a Mapping Table for Binary Variables
#'
#' Identifies binary variables and returns a deterministic mapping for 0/1 coding.
#' - Factors with \emph{explicit order} (ordered = TRUE) do NOT use heuristics; the highest (last) level is Positive.
#' - Logicals map to Negative = "FALSE", Positive = "TRUE".
#' - Numeric 0/1 (or any 2-value numeric) maps Positive to the numeric maximum.
#' - Characters / unordered factors use minimal heuristics (no race/PWH/sex terms).
#'
#' @param Data A dataframe.
#' @param CatVars Character vector of candidate binary variables.
#' @param prefer Optional named character vector of explicit positive levels,
#'   e.g., c(STATUS = "PWH", Smoker = "Yes"). This overrides other rules.
#' @return A data.frame with columns: Variable, Label, PositiveLevel, NegativeLevel.
#' @export
createBinaryMapping <- function(Data, CatVars, prefer = NULL) {

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
