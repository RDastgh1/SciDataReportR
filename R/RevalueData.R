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
  if (!requireNamespace("sjlabelled", quietly = TRUE))
    stop("Package 'sjlabelled' is required.")

  ## ---- helpers ------------------------------------------------------------
  .trim <- function(x) {
    x <- as.character(x)
    x <- gsub("\u2018|\u2019", "'", x, perl = TRUE)   # smart single quotes -> '
    x <- gsub("\u201C|\u201D", "\"", x, perl = TRUE)  # smart double quotes -> "
    trimws(x)
  }

  .recode_flag <- function(x) {
    x <- tolower(.trim(x))
    isTRUE(x %in% c("yes","y","1","true","t"))
  }

  .parse_missing_codes <- function(x) {
    if (is.null(x) || all(is.na(x))) return(character(0))
    s <- .trim(x[1])
    if (identical(s, "") || is.na(s)) return(character(0))
    s <- gsub("\\s+", "", s)
    out <- unlist(strsplit(s, "[,;|]", perl = TRUE))
    out[nzchar(out)]
  }

  .parse_code_map <- function(code_string, splitchar = ";") {
    cs <- .trim(code_string)
    if (length(cs) == 0L || is.na(cs) || cs == "") return(NULL)
    parts <- unlist(strsplit(cs, splitchar, fixed = TRUE))
    parts <- parts[nzchar(.trim(parts))]
    if (!length(parts)) return(NULL)

    kv <- lapply(parts, function(p) {
      # allow "a=b" or "a => b" or "a:b"
      if (grepl("=>", p, fixed = TRUE)) sp <- "=>"
      else if (grepl("=", p, fixed = TRUE)) sp <- "="
      else if (grepl(":", p, fixed = TRUE)) sp <- ":"
      else return(c(NA, NA))
      ab <- strsplit(p, sp, fixed = TRUE)[[1]]
      ab <- .trim(ab)
      if (length(ab) != 2) ab <- c(NA, NA)
      ab
    })
    kv <- do.call(rbind, kv)
    kv <- as.data.frame(kv, stringsAsFactors = FALSE)
    names(kv) <- c("from","to")
    kv <- kv[!is.na(kv$from) & !is.na(kv$to) & nzchar(kv$from), , drop = FALSE]
    if (!nrow(kv)) return(NULL)
    kv
  }

  .first_nonempty <- function(x) {
    x <- .trim(x)
    x[which(nzchar(x))[1]]
  }

  .unique_collapse <- function(x, sep) {
    x <- .trim(x)
    x <- x[nzchar(x) & !is.na(x)]
    if (!length(x)) return(NA_character_)
    paste(unique(x), collapse = sep)
  }

  ## ---- normalize VarTypes -------------------------------------------------
  if (is.null(VarTypes) || !NROW(VarTypes))
    stop("VarTypes is empty.")

  # normalize expected columns to character, trim
  need_cols <- c("Variable","Recode","Code","Type","Label","MissingCode","Missing")
  present <- intersect(need_cols, names(VarTypes))
  for (cc in present) VarTypes[[cc]] <- .trim(VarTypes[[cc]])

  # ensure Variable exists
  if (!"Variable" %in% names(VarTypes))
    stop("VarTypes must have a 'Variable' column.")
  VarTypes$Variable <- .trim(VarTypes$Variable)
  VarTypes <- VarTypes[!is.na(VarTypes$Variable) & nzchar(VarTypes$Variable), , drop = FALSE]

  # back-compat: MissingCode <- Missing if MissingCode absent or blank
  if (!"MissingCode" %in% names(VarTypes)) VarTypes$MissingCode <- NA_character_
  if ("Missing" %in% names(VarTypes)) {
    use_missing <- is.na(VarTypes$MissingCode) | VarTypes$MissingCode == ""
    VarTypes$MissingCode[use_missing] <- VarTypes$Missing[use_missing]
  }

  # de-dup vars: if multiple rows per variable, collapse fields sensibly
  if (any(duplicated(VarTypes$Variable))) {
    collapsed <- lapply(split(VarTypes, VarTypes$Variable), function(df) {
      data.frame(
        Variable    = df$Variable[1],
        Recode      = .first_nonempty(df$Recode),
        Code        = .unique_collapse(df$Code, splitchar),
        Type        = .first_nonempty(df$Type),
        Label       = .first_nonempty(df$Label),
        MissingCode = .unique_collapse(df$MissingCode, ","),
        stringsAsFactors = FALSE
      )
    })
    VarTypes <- do.call(rbind, collapsed)
  }

  ## ---- select variables present in data ----------------------------------
  RevaluedData <- DatatoRevalue
  VariablesinData <- colnames(DatatoRevalue)
  VariablesinVarTypes <- unique(VarTypes$Variable)
  not_in_data <- setdiff(VariablesinVarTypes, VariablesinData)
  vt_in <- VarTypes[VarTypes$Variable %in% VariablesinData, , drop = FALSE]
  vars <- vt_in$Variable

  warninglist <- character(0)
  recodedvars <- character(0)

  ## ---- main loop ----------------------------------------------------------
  for (i in seq_along(vars)) {
    var <- vars[i]
    idx <- which(vt_in$Variable == var)[1]
    x   <- RevaluedData[[var]]

    # ---- missing codes to NA ---------------------------------------------
    m_codes_chr <- .parse_missing_codes(vt_in$MissingCode[idx])
    if (!length(m_codes_chr) || all(is.na(m_codes_chr))) {
      # fallback to function-level default if nothing is specified
      m_codes_chr <- as.character(missingVal)
    }

    if (is.numeric(x)) {
      suppressWarnings(m_codes_num <- as.numeric(m_codes_chr))
      m_codes_num <- m_codes_num[!is.na(m_codes_num)]
      if (length(m_codes_num)) x[!is.na(x) & x %in% m_codes_num] <- NA
    } else {
      x_chr <- .trim(as.character(x))
      if (length(m_codes_chr)) x[ x_chr %in% .trim(m_codes_chr) ] <- NA
    }
    RevaluedData[[var]] <- x

    # ---- optional recoding / labeling ------------------------------------
    rc_raw    <- vt_in$Recode[idx]
    do_recode <- .recode_flag(rc_raw)
    vartype   <- tolower(.trim(vt_in$Type[idx]))
    code_map  <- NULL

    if (do_recode) {
      code_map <- tryCatch(.parse_code_map(vt_in$Code[idx], splitchar),
                           error = function(e) {
                             warninglist <<- c(warninglist,
                                               paste0(var, ": recoding failed to parse Code - ", e$message))
                             NULL
                           })
      if (is.null(code_map)) {
        warninglist <- c(warninglist, paste0(var, ": no valid Code pairs to recode."))
      } else {
        # Build named vector for labels: names = from, values = to
        nameList <- code_map$to
        names(nameList) <- code_map$from

        # Attach value labels; storage handling based on Type
        d <- sjlabelled::set_labels(RevaluedData[[var]],
                                    labels = nameList,
                                    force.labels = TRUE,
                                    force.values = TRUE)

        if (vartype %in% c("categorical","factor","ordinal","ordered factor","ordered")) {
          RevaluedData[[var]] <- sjlabelled::as_label(d)
          if (vartype %in% c("ordinal","ordered factor","ordered"))
            RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        } else if (vartype %in% c("double","numeric","numerical","integer")) {
          suppressWarnings(RevaluedData[[var]] <- as.numeric(RevaluedData[[var]]))
          RevaluedData[[var]] <- sjlabelled::set_labels(RevaluedData[[var]],
                                                        labels = nameList,
                                                        force.labels = TRUE,
                                                        force.values = TRUE)
        } else {
          # Unknown / unspecified type: keep original storage, just attach labels
          RevaluedData[[var]] <- d
        }

        recodedvars <- c(recodedvars, var)
      }
    } else {
      # Type-only coercion if requested without an explicit Code map
      if (nzchar(vartype)) {
        if (vartype %in% c("categorical","factor")) {
          RevaluedData[[var]] <- sjlabelled::to_factor(RevaluedData[[var]])
        } else if (vartype %in% c("double","numeric","numerical","integer")) {
          RevaluedData[[var]] <- suppressWarnings(as.numeric(RevaluedData[[var]]))
        } else if (vartype %in% c("ordinal","ordered factor","ordered")) {
          RevaluedData[[var]] <- as.ordered(RevaluedData[[var]])
        }
      }
    }

    # ---- variable label ---------------------------------------------------
    vlabel <- vt_in$Label[idx]
    if (!is.na(vlabel) && nzchar(vlabel)) {
      RevaluedData[[var]] <- sjlabelled::set_label(RevaluedData[[var]], label = vlabel)
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
    warninglist  = unique(warninglist),
    recodedvars  = unique(recodedvars),
    not_in_data  = not_in_data
  )
}
