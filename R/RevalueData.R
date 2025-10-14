#' Revalue Data
#'
#' Revalues variables in a dataset using a VarTypes codebook.
#'
#' @param DatatoRevalue A data.frame or tibble to be revalued.
#' @param VarTypes A data.frame with columns:
#'   Variable, Recode, Code, Type, Label, MissingCode. Only Variable is required.
#'   (Backward compatible: if MissingCode is absent/NA, will fall back to Missing.)
#' @param missingVal Default value to treat as missing when VarTypes$MissingCode is absent or NA.
#' @param splitchar Separator used in VarTypes$Code between pairs (default ";").
#'
#' @return A list with:
#'   RevaluedData (data), warninglist (character), recodedvars (character),
#'   not_in_data (character: variables listed in VarTypes but not found in the data).
#' @export
RevalueData <- function(DatatoRevalue, VarTypes, missingVal = -999, splitchar = ";") {
  if (!requireNamespace("sjlabelled", quietly = TRUE)) {
    stop("Package 'sjlabelled' is required.")
  }
  has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
  `%>%` <- if (has_dplyr) get("%>%", asNamespace("dplyr")) else NULL

  if (has_dplyr) {
    VarTypes <- VarTypes %>% dplyr::filter(!is.na(.data$Variable))
  } else {
    VarTypes <- VarTypes[!is.na(VarTypes$Variable), , drop = FALSE]
  }

  RevaluedData <- DatatoRevalue
  VariablesinData <- colnames(DatatoRevalue)
  VariablesinVarTypes <- as.character(VarTypes$Variable)
  not_in_data <- setdiff(unique(VariablesinVarTypes), VariablesinData)

  if (has_dplyr) {
    VarTypes_in_data <- VarTypes %>% dplyr::filter(.data$Variable %in% VariablesinData)
  } else {
    VarTypes_in_data <- VarTypes[VarTypes$Variable %in% VariablesinData, , drop = FALSE]
  }
  vars <- unique(as.character(VarTypes_in_data$Variable))

  warninglist <- character(0)
  recodedvars <- character(0)

  # --- normalize MissingCode (with backward compatibility) ---
  if (!("MissingCode" %in% colnames(VarTypes_in_data))) {
    VarTypes_in_data$MissingCode <- NA
  }
  if ("Missing" %in% colnames(VarTypes_in_data)) {
    use_missing <- is.na(VarTypes_in_data$MissingCode) | VarTypes_in_data$MissingCode == ""
    VarTypes_in_data$MissingCode[use_missing] <- VarTypes_in_data$Missing[use_missing]
  }
  VarTypes_in_data$MissingCode[is.na(VarTypes_in_data$MissingCode) | VarTypes_in_data$MissingCode == ""] <- missingVal

  parse_missing_codes <- function(x) {
    if (is.null(x) || all(is.na(x))) return(character(0))
    s <- trimws(as.character(x)[1])
    if (s == "") return(character(0))
    s <- gsub("\\s+", "", s)
    out <- unlist(strsplit(s, "[,;|]"))
    out[nzchar(out)]
  }

  recode_flag <- function(x) {
    x <- tolower(trimws(as.character(x)))
    isTRUE(x %in% c("yes", "y", "1", "true"))
  }

  for (var in vars) {
    idx <- match(var, VarTypes_in_data$Variable)
    x <- RevaluedData[[var]]

    # replace any listed MissingCode values with NA
    m_codes_raw <- VarTypes_in_data$MissingCode[idx]
    m_codes_chr <- parse_missing_codes(m_codes_raw)
    if (length(m_codes_chr) == 0) m_codes_chr <- as.character(missingVal)

    if (is.numeric(x)) {
      suppressWarnings(m_codes_num <- as.numeric(m_codes_chr))
      m_codes_num <- m_codes_num[!is.na(m_codes_num)]
      if (length(m_codes_num) > 0) x[!is.na(x) & x %in% m_codes_num] <- NA
    } else {
      x_chr <- trimws(as.character(x))
      if (length(m_codes_chr) > 0) x[ x_chr %in% trimws(m_codes_chr) ] <- NA
    }
    RevaluedData[[var]] <- x

    # recoding
    rc_raw <- VarTypes_in_data$Recode[idx]
    do_recode <- recode_flag(rc_raw)

    if (isTRUE(do_recode)) {
      recodedvars <- c(recodedvars, var)
      newcode <- VarTypes_in_data$Code[idx]

      lookup <- tryCatch({
        parts <- strsplit(as.character(newcode), splitchar)[[1]]
        kv <- lapply(parts, function(p) strsplit(p, "=", fixed = TRUE)[[1]])
        df <- as.data.frame(do.call(rbind, kv), stringsAsFactors = FALSE)
        if (ncol(df) != 2) stop("Code parsing produced malformed key/value pairs.")
        names(df) <- c("X1", "X2")
        df$X1 <- trimws(df$X1); df$X2 <- trimws(df$X2)
        df
      }, error = function(e) {
        warninglist <<- c(warninglist, paste0(var, ": recoding failed during lookup creation - ", e$message))
        NULL
      })

      if (!is.null(lookup)) {
        freqs1 <- summary(as.factor(RevaluedData[[var]]))
        nameList <- lookup$X2; names(nameList) <- lookup$X1
        d <- sjlabelled::set_labels(RevaluedData[[var]],
                                    labels = nameList,
                                    force.labels = TRUE,
                                    force.values = TRUE)
        RevaluedData[[var]] <- sjlabelled::as_label(d)
        freqs2 <- summary(as.factor(RevaluedData[[var]]))
        if (!identical(sum(freqs1), sum(freqs2))) {
          warninglist <- c(warninglist, paste0(var, ": final frequencies not consistent"))
        }
      }
    } else {
      type <- VarTypes_in_data$Type[idx]
      if (!is.na(type)) {
        type <- tolower(trimws(as.character(type)))
        if (type %in% c("categorical", "factor")) {
          RevaluedData[[var]] <- sjlabelled::to_factor(RevaluedData[[var]])
        } else if (type %in% c("double", "numeric", "numerical")) {
          RevaluedData[[var]] <- suppressWarnings(as.numeric(RevaluedData[[var]]))
        } else if (type %in% c("ordinal", "ordered factor", "ordered")) {
          RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        }
      }
    }

    newVarLabel <- VarTypes_in_data$Label[idx]
    if (length(newVarLabel) == 1 && !is.na(newVarLabel) && nzchar(newVarLabel)) {
      RevaluedData[[var]] <- sjlabelled::set_label(RevaluedData[[var]], label = newVarLabel)
    }
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
    warninglist  = warninglist,
    recodedvars  = recodedvars,
    not_in_data  = not_in_data
  )
}
