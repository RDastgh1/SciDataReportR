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
#'   not_in_data (character).
#' @export
RevalueData <- function(DatatoRevalue, VarTypes, missingVal = -999, splitchar = ";") {
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
    # normalize smart quotes
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

  # Collapse duplicate Variable rows (first non-empty for singletons; union for lists)
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

  ## --- Main loop -----------------------------------------------------------
  for (i in seq_along(vars)) {
    var <- vars[i]
    idx <- which(vt_in$Variable == var)[1]
    x   <- RevaluedData[[var]]

    ## MissingCode -> NA (supports multiple tokens: ",", ";", "|")
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
            # Expect "code = label"
            from <- c(from, ab[1])   # codes
            to   <- c(to,   ab[2])   # labels
          }
        }

        if (!length(from)) {
          warninglist <- c(warninglist, paste0(var, ": Code parsing produced no valid key/value pairs."))
        } else {
          # If codes are numeric, coerce storage to numeric; "." etc already NA above if char
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
          } else if (vartype %in% c("double","numeric","numerical","integer")) {
            RevaluedData[[var]] <- tmp  # keep numeric, labels attached
          } else {
            RevaluedData[[var]] <- tmp  # leave storage as-is with labels
          }

          recodedvars <- c(recodedvars, var)
        }
      }
    } else {
      # No recode map; still coerce by Type if requested
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

    # Variable label
    vlab <- vt_in$Label[idx]
    if (!is.na(vlab) && nzchar(vlab)) {
      RevaluedData[[var]] <- sjlabelled::set_label(RevaluedData[[var]], label = vlab)
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
