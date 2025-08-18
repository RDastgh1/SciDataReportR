#' Revalue Data
#'
#' Revalues variables in a dataset according to the specifications in a codebook.
#'
#' @param DatatoRevalue A data.frame/tibble to be revalued.
#' @param VarTypes A data.frame with columns:
#'   - Variable (variable names)
#'   - Recode ("yes"/"no", 1/0; optional)
#'   - Code (revalue codes separated by `splitchar`, e.g., "1=A;2=B"; optional)
#'   - Type (Categorical/Double/Ordinal; optional)
#'   - Label (optional variable label)
#'   - Missing (optional; value to treat as missing for that variable)
#' @param missingVal The value to use for missing values when VarTypes$Missing is absent/NA (default -999).
#' @param splitchar The character used to split codes in the Code column (default ";").
#' @return A list with:
#'   - RevaluedData: the revalued dataset
#'   - warninglist: character vector of warnings
#'   - recodedvars: character vector of variables that were recoded
#'   - not_in_data: character vector of VarTypes variables not found in DatatoRevalue
#' @export
RevalueData <- function(DatatoRevalue, VarTypes, missingVal = -999, splitchar = ";") {
  # deps
  if (!requireNamespace("sjlabelled", quietly = TRUE)) {
    stop("Package 'sjlabelled' is required.")
  }
  # dplyr only for convenience; fall back if absent
  has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
  `%>%` <- if (has_dplyr) get("%>%", asNamespace("dplyr")) else NULL

  # sanitize VarTypes
  if (has_dplyr) {
    VarTypes <- VarTypes %>% dplyr::filter(!is.na(Variable))
  } else {
    VarTypes <- VarTypes[!is.na(VarTypes$Variable), , drop = FALSE]
  }

  RevaluedData <- DatatoRevalue
  VariablesinData <- colnames(DatatoRevalue)
  VariablesinVarTypes <- as.character(VarTypes$Variable)

  # NEW: capture variables listed in VarTypes but absent from the data
  not_in_data <- setdiff(unique(VariablesinVarTypes), VariablesinData)

  # Only operate on variables that exist in the data
  if (has_dplyr) {
    VarTypes_in_data <- VarTypes %>% dplyr::filter(.data$Variable %in% VariablesinData)
  } else {
    VarTypes_in_data <- VarTypes[VarTypes$Variable %in% VariablesinData, , drop = FALSE]
  }
  vars <- unique(as.character(VarTypes_in_data$Variable))

  warninglist <- character(0)
  recodedvars <- character(0)

  # Ensure helper columns exist/are filled
  if (!("Recode" %in% colnames(VarTypes_in_data))) {
    VarTypes_in_data$Recode <- NA
  }
  if (!("Missing" %in% colnames(VarTypes_in_data))) {
    VarTypes_in_data$Missing <- missingVal
  }
  VarTypes_in_data$Missing[is.na(VarTypes_in_data$Missing)] <- missingVal

  # normalize Recode flag
  recode_flag <- function(x) {
    x <- tolower(trimws(as.character(x)))
    isTRUE(x %in% c("yes", "y", "1", "true"))
  }

  for (var in vars) {
    mV <- VarTypes_in_data$Missing[match(var, VarTypes_in_data$Variable)]
    x <- RevaluedData[[var]]
    x[x == mV] <- NA
    RevaluedData[[var]] <- x

    rc_raw <- VarTypes_in_data$Recode[match(var, VarTypes_in_data$Variable)]
    do_recode <- recode_flag(rc_raw)

    if (isTRUE(do_recode)) {
      recodedvars <- c(recodedvars, var)
      newcode <- VarTypes_in_data$Code[match(var, VarTypes_in_data$Variable)]

      lookup <- tryCatch({
        # Split like "1=A; 2=B"
        parts <- strsplit(as.character(newcode), splitchar)[[1]]
        kv <- lapply(parts, function(p) strsplit(p, "=", fixed = TRUE)[[1]])
        df <- as.data.frame(do.call(rbind, kv), stringsAsFactors = FALSE)
        if (ncol(df) != 2) stop("Code parsing produced malformed key/value pairs.")
        names(df) <- c("X1", "X2")
        df$X1 <- trimws(df$X1)
        df$X2 <- trimws(df$X2)
        df
      }, error = function(e) {
        warninglist <<- c(warninglist, paste0(var, ": recoding failed during lookup creation - ", e$message))
        NULL
      })

      if (!is.null(lookup)) {
        freqs1 <- summary(as.factor(RevaluedData[[var]]))
        nameList <- lookup$X2
        names(nameList) <- lookup$X1

        # apply labels and then convert to labelled -> factor with value labels
        d <- sjlabelled::set_labels(RevaluedData[[var]],
                                    labels = nameList,
                                    force.labels = TRUE,
                                    force.values = TRUE)
        RevaluedData[[var]] <- sjlabelled::as_label(d)
        freqs2 <- summary(as.factor(RevaluedData[[var]]))

        # simple consistency check on totals
        if (!identical(sum(freqs1), sum(freqs2))) {
          warninglist <- c(warninglist, paste0(var, ": final frequencies not consistent"))
        }
      }
    } else {
      type <- VarTypes_in_data$Type[match(var, VarTypes_in_data$Variable)]
      if (!is.na(type)) {
        type <- tolower(trimws(as.character(type)))
        if (type %in% c("categorical", "factor")) {
          RevaluedData[[var]] <- sjlabelled::to_factor(RevaluedData[[var]])
        }
        if (type %in% c("double", "numeric", "numerical")) {
          RevaluedData[[var]] <- suppressWarnings(as.numeric(RevaluedData[[var]]))
        }
        if (type %in% c("ordinal", "ordered factor", "ordered")) {
          RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        }
      }
    }

    # variable label (if provided)
    newVarLabel <- VarTypes_in_data$Label[match(var, VarTypes_in_data$Variable)]
    if (length(newVarLabel) == 1 && !is.na(newVarLabel) && nzchar(newVarLabel)) {
      RevaluedData[[var]] <- sjlabelled::set_label(RevaluedData[[var]], label = newVarLabel)
    }
  }

  # If there were any absent variables, note them in warnings, too
  if (length(not_in_data) > 0) {
    warninglist <- c(
      warninglist,
      paste0("Variables listed in VarTypes but not found in data (ignored): ",
             paste(not_in_data, collapse = ", "))
    )
  }

  return(list(
    RevaluedData = RevaluedData,
    warninglist  = warninglist,
    recodedvars  = recodedvars,
    not_in_data  = not_in_data
  ))
}
