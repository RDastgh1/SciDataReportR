#' Revalue Data
#'
#' Revalues variables in a dataset using a VarTypes codebook.
#'
#' @param data A data.frame or tibble to be revalued.
#' @param codebook A data.frame with columns:
#'   Variable, Recode, Code, Type, Label, MissingCode. Only Variable is required.
#'   (Backward compatible: if MissingCode is absent/NA, will fall back to Missing.)
#' @param missingVal Default value to treat as missing when VarTypes$MissingCode is absent or NA.
#' @param splitchar Separator used in VarTypes$Code between pairs (default ";").
#' @param on_error Whether to stop at the first variable-level error (the default)
#'   or continue and record errors in the returned object.
#'
#' @return A list with:
#'   RevaluedData (data), warninglist (character), recodedvars (character),
#'   not_in_data (character), and errors (data frame with `Variable` and
#'   `Error` columns). In the default `on_error = "stop"` mode, an error names
#'   the offending variable and preserves the underlying message.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Before: no labels, and sex is stored as 0/1
#' sjlabelled::get_label(SampleData$age)   # NULL
#' class(SampleData$sex)                    # "integer"
#'
#' # Revalue using the codebook
#' revalued <- RevalueData(SampleData, SampleVariableTypes)
#' Labelled <- revalued$RevaluedData
#'
#' # After: labels attached, sex is a labelled factor
#' sjlabelled::get_label(Labelled$age)     # "Age"
#' levels(Labelled$sex)                     # "Female" "Male"
#'
#' # Recoded variables, and codebook entries not found in the data
#' revalued$recodedvars
#' revalued$not_in_data
#'
#' \donttest{
#' # A side-by-side view of what changed
#' vars_Show <- c("Diagnosis", "age", "sex", "Genotype", "AXL")
#'
#' df_Before <- utils::head(SampleData[, vars_Show], 6)
#' df_After <- utils::head(Labelled[, vars_Show], 6)
#'
#' ShowTable <- function(x, caption) {
#'   htmltools::browsable(htmltools::HTML(as.character(
#'     kableExtra::kable_styling(
#'       knitr::kable(x, format = "html", caption = caption, row.names = FALSE),
#'       bootstrap_options = c("striped", "hover", "condensed"),
#'       full_width = FALSE
#'     )
#'   )))
#' }
#'
#' ShowTable(df_Before, "Before: raw codes as imported")
#' ShowTable(df_After, "After: recoded and labelled")
#'
#' # The labels the codebook attached
#' df_Labels <- data.frame(
#'   Variable = vars_Show,
#'   Label = vapply(
#'     vars_Show,
#'     function(v) {
#'       lab <- sjlabelled::get_label(Labelled[[v]])
#'       if (is.null(lab) || is.na(lab)) v else lab
#'     },
#'     character(1)
#'   ),
#'   Class = vapply(vars_Show, function(v) class(Labelled[[v]])[1], character(1)),
#'   Levels = vapply(
#'     vars_Show,
#'     function(v) {
#'       lv <- levels(Labelled[[v]])
#'       if (is.null(lv)) "" else paste(lv, collapse = ", ")
#'     },
#'     character(1)
#'   ),
#'   row.names = NULL
#' )
#'
#' htmltools::browsable(htmltools::HTML(as.character(
#'   FreezeTableHeader(df_Labels, full_width = TRUE)
#' )))
#'
#' # Anything the codebook could not act on
#' revalued$warninglist
#' }
#'
#' @section What changes:
#' In the raw extract `sex` is a bare 0/1 column and nothing carries a label.
#' Afterwards it is a factor with real levels, and the labels are attached for
#' every downstream table and plot to pick up automatically - which is what
#' makes output readable without renaming anything by hand.
#'
#' Anything the codebook could not act on is reported in `warninglist` rather
#' than silently skipped, so a mistyped variable name or an unparseable `Code`
#' string is visible instead of quietly leaving a variable unrecoded.
#'
#' @param DatatoRevalue \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param VarTypes \strong{Deprecated} (since 19.15.0). Use \code{codebook} instead.
#' @export
RevalueData <- function(data,
    codebook,
    missingVal = -999,
    splitchar = ";",
    on_error = c("stop", "warn"),
    DatatoRevalue = lifecycle::deprecated(),
    VarTypes = lifecycle::deprecated()) {
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(DatatoRevalue)) {
    lifecycle::deprecate_warn("19.15.0", "RevalueData(DatatoRevalue)", "RevalueData(data)")
    data <- DatatoRevalue
  }
  if (!missing(data)) DatatoRevalue <- data
  if (lifecycle::is_present(VarTypes)) {
    lifecycle::deprecate_warn("19.15.0", "RevalueData(VarTypes)", "RevalueData(codebook)")
    codebook <- VarTypes
  }
  if (!missing(codebook)) VarTypes <- codebook
  on_error <- match.arg(on_error)

  if (!requireNamespace("sjlabelled", quietly = TRUE))
    stop("Package 'sjlabelled' is required.")

  ## --- Normalize VarTypes --------------------------------------------------
  if (is.null(VarTypes) || !NROW(VarTypes)) stop("VarTypes is empty.")

  norm_cols <- intersect(
    c("Variable","Recode","Code","Type","Label","MissingCode","Missing"),
    names(VarTypes)
  )
  for (cc in norm_cols) {
    VarTypes[[cc]] <- as.character(VarTypes[[cc]])
    VarTypes[[cc]] <- gsub("\u2018|\u2019", "'", VarTypes[[cc]], perl = TRUE)
    VarTypes[[cc]] <- gsub("\u201C|\u201D", "\"", VarTypes[[cc]], perl = TRUE)
    VarTypes[[cc]] <- trimws(VarTypes[[cc]])
  }

  if (!"Variable" %in% names(VarTypes))
    stop("VarTypes must have a 'Variable' column.")
  VarTypes$Variable <- trimws(as.character(VarTypes$Variable))
  VarTypes <- VarTypes[!is.na(VarTypes$Variable) & nzchar(VarTypes$Variable), , drop = FALSE]

  # Back-compat: fill MissingCode from Missing when absent/blank
  if (!"MissingCode" %in% names(VarTypes)) VarTypes$MissingCode <- NA_character_
  if ("Missing" %in% names(VarTypes)) {
    use_missing <- is.na(VarTypes$MissingCode) | VarTypes$MissingCode == ""
    VarTypes$MissingCode[use_missing] <- VarTypes$Missing[use_missing]
  }

  # Collapse duplicate Variable rows
  if (any(duplicated(VarTypes$Variable))) {
    split_list <- split(VarTypes, VarTypes$Variable)
    VarTypes <- do.call(rbind, lapply(split_list, function(df) {
      data.frame(
        Variable    = df$Variable[1],
        Recode      = { v <- df$Recode; v <- v[nzchar(v)]; if (length(v)) v[1] else NA_character_ },
        Code        = { v <- df$Code; v <- v[nzchar(v)]; if (length(v)) paste(unique(v), collapse = splitchar) else NA_character_ },
        Type        = { v <- df$Type; v <- v[nzchar(v)]; if (length(v)) v[1] else NA_character_ },
        Label       = { v <- df$Label; v <- v[nzchar(v)]; if (length(v)) v[1] else NA_character_ },
        MissingCode = { v <- df$MissingCode; v <- v[nzchar(v)]; if (length(v)) paste(unique(v), collapse = ",") else NA_character_ },
        stringsAsFactors = FALSE
      )
    }))
    rownames(VarTypes) <- NULL
  }

  ## --- Variables present in data ------------------------------------------
  RevaluedData        <- DatatoRevalue
  VariablesinData     <- colnames(DatatoRevalue)
  VariablesinVarTypes <- unique(VarTypes$Variable)
  not_in_data         <- setdiff(VariablesinVarTypes, VariablesinData)
  vt_in               <- VarTypes[VarTypes$Variable %in% VariablesinData, , drop = FALSE]
  vars                <- vt_in$Variable

  warninglist <- character(0)
  recodedvars <- character(0)
  errorlist <- data.frame(
    Variable = character(0),
    Error = character(0),
    stringsAsFactors = FALSE
  )

  # treat these as numeric/double storage
  num_types <- c("double","numeric","numerical","integer","continuous")

  ## --- Main loop -----------------------------------------------------------
  for (i in seq_along(vars)) {
    var <- vars[i]
    idx <- which(vt_in$Variable == var)[1]
    tryCatch({
      x <- RevaluedData[[var]]

    ## MissingCode -> NA
    mchr <- vt_in$MissingCode[idx]
    if (is.na(mchr) || mchr == "") mchr <- as.character(missingVal)
    mchr <- gsub("\\s+", "", mchr)
    mtok <- unlist(strsplit(mchr, "[,;|]", perl = TRUE))
    mtok <- mtok[nzchar(mtok)]

    if (is.numeric(x)) {
      suppressWarnings(mnum <- as.numeric(mtok))
      mnum <- mnum[!is.na(mnum)]
      if (length(mnum)) x[!is.na(x) & x %in% mnum] <- NA
    } else {
      x_chr <- trimws(as.character(x))
      if (length(mtok)) x[ x_chr %in% trimws(mtok) ] <- NA
    }
    RevaluedData[[var]] <- x

    ## Recoding & labeling
    rc_flag <- tolower(trimws(vt_in$Recode[idx]))
    valid_recode_flags <- c("yes", "y", "1", "true", "t", "no", "n", "0", "false", "f")
    if (!is.na(rc_flag) && nzchar(rc_flag) && !rc_flag %in% valid_recode_flags) {
      warninglist <- c(
        warninglist,
        paste0(var, ": Unrecognized Recode value '", vt_in$Recode[idx], "'. Expected Yes or No; recoding was skipped.")
      )
    }
    do_reco <- isTRUE(rc_flag %in% c("yes","y","1","true","t"))
    vartype <- tolower(trimws(vt_in$Type[idx]))

    if (do_reco) {
      codestr <- vt_in$Code[idx]
      if (is.na(codestr) || !nzchar(codestr)) {
        warninglist <- c(warninglist, paste0(var, ": Recode requested but Code is empty."))
      } else {
        parts <- strsplit(codestr, splitchar, fixed = TRUE)[[1]]
        parts <- parts[nzchar(trimws(parts))]

        from <- character(0); to <- character(0)
        for (p in parts) {
          p <- trimws(p)
          sep <- if (grepl("=>", p, fixed = TRUE)) "=>"
          else if (grepl("=", p,  fixed = TRUE)) "="
          else if (grepl(":", p,  fixed = TRUE)) ":"
          else NA_character_
          if (is.na(sep)) next
          ab <- strsplit(p, sep, fixed = TRUE)[[1]]
          ab <- trimws(ab)
          if (length(ab) >= 2 && nzchar(ab[1]) && nzchar(ab[2])) {
            from <- c(from, ab[1])
            to   <- c(to,   ab[2])
          }
        }

        if (!length(from)) {
          warninglist <- c(warninglist, paste0(var, ": Code parsing produced no valid key/value pairs."))
        } else {
          numeric_code <- all(grepl("^[-+]?[0-9]+(\\.[0-9]+)?$", from))
          if (numeric_code) {
            RevaluedData[[var]] <- suppressWarnings(as.numeric(RevaluedData[[var]]))
            labs_vec <- suppressWarnings(as.numeric(from))
            names(labs_vec) <- to  # names = labels, values = codes
          } else {
            labs_vec <- from
            names(labs_vec) <- to
          }

          tmp <- sjlabelled::set_labels(
            RevaluedData[[var]],
            labels = labs_vec,
            force.labels = TRUE,
            force.values = TRUE
          )

          if (vartype %in% c("categorical","factor","ordinal","ordered factor","ordered")) {
            RevaluedData[[var]] <- sjlabelled::as_label(tmp)
            if (vartype %in% c("ordinal","ordered factor","ordered"))
              RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
          } else if (vartype %in% num_types) {
            RevaluedData[[var]] <- tmp  # numeric with labels kept
          } else {
            RevaluedData[[var]] <- tmp
          }

          recodedvars <- c(recodedvars, var)
        }
      }
    } else {
      # No recode map; still coerce by Type if requested
      if (nzchar(vartype)) {
        if (vartype %in% c("categorical","factor")) {
          RevaluedData[[var]] <- sjlabelled::to_factor(RevaluedData[[var]])
        } else if (vartype %in% num_types) {
          RevaluedData[[var]] <- suppressWarnings(as.numeric(RevaluedData[[var]]))
        } else if (vartype %in% c("ordinal","ordered factor","ordered")) {
          RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        }
      }
    }

    # Variable label
    vlab <- vt_in$Label[idx]
    if (length(vlab) > 0 && !is.na(vlab) && nzchar(vlab)) {
      RevaluedData[[var]] <- sjlabelled::set_label(RevaluedData[[var]], label = vlab)
    }

    # Preserve the measurement-level decision from the codebook. Ordered
    # factors remain the categorical representation; the score metadata lets
    # downstream functions create a continuous representation without trying
    # to recover the original codes from display labels.
    if (vartype %in% c("ordinal", "ordered factor", "ordered")) {
      x_ordinal <- RevaluedData[[var]]
      if (!is.ordered(x_ordinal)) x_ordinal <- as.ordered(x_ordinal)
      attr(x_ordinal, "scidr_type") <- "ordinal"

      codestr <- vt_in$Code[idx]
      score_map <- numeric(0)
      if (!is.na(codestr) && nzchar(codestr)) {
        parts <- strsplit(codestr, splitchar, fixed = TRUE)[[1]]
        for (part in parts) {
          sep <- if (grepl("=>", part, fixed = TRUE)) "=>" else if (grepl("=", part, fixed = TRUE)) "=" else if (grepl(":", part, fixed = TRUE)) ":" else NA_character_
          if (is.na(sep)) next
          pair <- trimws(strsplit(part, sep, fixed = TRUE)[[1]])
          if (length(pair) < 2) next
          score <- suppressWarnings(as.numeric(pair[1]))
          if (!is.na(score) && nzchar(pair[2])) score_map[pair[2]] <- score
        }
      }

      ordinal_levels <- levels(x_ordinal)
      if (length(score_map) && all(ordinal_levels %in% names(score_map))) {
        attr(x_ordinal, "scidr_ordinal_scores") <- score_map[ordinal_levels]
        attr(x_ordinal, "scidr_ordinal_score_source") <- "codebook"
      } else {
        attr(x_ordinal, "scidr_ordinal_scores") <- stats::setNames(seq_along(ordinal_levels), ordinal_levels)
        attr(x_ordinal, "scidr_ordinal_score_source") <- "rank"
        warninglist <- c(
          warninglist,
          paste0(var, ": Ordinal continuous scores will use ordered-level ranks because Code does not provide a complete numeric mapping.")
        )
      }
      RevaluedData[[var]] <- x_ordinal
    }
    }, error = function(e) {
      error_message <- paste0("Error revaluing variable '", var, "': ", conditionMessage(e))
      errorlist <<- rbind(
        errorlist,
        data.frame(Variable = var, Error = conditionMessage(e), stringsAsFactors = FALSE)
      )
      if (on_error == "stop") stop(error_message, call. = FALSE)
      warninglist <<- c(warninglist, error_message)
    })
  }

  if (length(not_in_data) > 0) {
    warninglist <- c(
      warninglist,
      paste0("Variables listed in VarTypes but not found in data (ignored): ",
             paste(not_in_data, collapse = ", "))
    )
  }

  list(
    RevaluedData = RevaluedData,
    warninglist  = unique(warninglist),
    recodedvars  = unique(recodedvars),
    not_in_data  = not_in_data,
    errors       = errorlist
  )
}
