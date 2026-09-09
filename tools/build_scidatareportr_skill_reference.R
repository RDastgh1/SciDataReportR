# Build the checked-in SciDataReportR skill API reference from NAMESPACE and Rd.

ExtractNamespaceExports <- function(fn_Namespace) {
  lines_Namespace <- readLines(fn_Namespace, warn = FALSE)
  matches_Exports <- regexec('^export\\("?([^"\\)]+)"?\\)$', lines_Namespace)
  values_Exports <- regmatches(lines_Namespace, matches_Exports)
  exports <- vapply(
    values_Exports[vapply(values_Exports, length, integer(1)) > 0L],
    function(match) match[[2]],
    character(1)
  )

  sort(unique(exports))
}

ReadMacroContents <- function(text_Rd, macro, start_at = 1L) {
  pattern_Macro <- if (identical(macro, "")) "{" else paste0("\\", macro, "{")
  start_Macro <- regexpr(pattern_Macro, substr(text_Rd, start_at, nchar(text_Rd)), fixed = TRUE)

  if (start_Macro[[1]] < 0L) {
    return(NULL)
  }

  open_Brace <- start_at + start_Macro[[1]] + attr(start_Macro, "match.length") - 2L
  position <- open_Brace + 1L
  depth <- 1L

  while (position <= nchar(text_Rd) && depth > 0L) {
    character_Current <- substr(text_Rd, position, position)
    character_Previous <- if (position > 1L) substr(text_Rd, position - 1L, position - 1L) else ""

    if (character_Current == "{" && character_Previous != "\\") {
      depth <- depth + 1L
    }
    if (character_Current == "}" && character_Previous != "\\") {
      depth <- depth - 1L
    }
    position <- position + 1L
  }

  if (depth != 0L) {
    stop("Unbalanced Rd macro: ", macro, call. = FALSE)
  }

  list(
    text = substr(text_Rd, open_Brace + 1L, position - 2L),
    next_at = position
  )
}

CleanRdText <- function(text_Rd) {
  text_Clean <- gsub("\\", "", text_Rd, fixed = TRUE)
  text_Clean <- gsub("[{}]", "", text_Clean)
  text_Clean <- gsub("[[:space:]]+", " ", text_Clean)
  trimws(text_Clean)
}

ExtractRdItems <- function(text_Rd) {
  items <- list()
  position <- 1L

  repeat {
    item_Name <- ReadMacroContents(text_Rd, "item", start_at = position)
    if (is.null(item_Name)) {
      break
    }

    item_Description <- ReadMacroContents(text_Rd, "", start_at = item_Name$next_at)
    if (is.null(item_Description)) {
      break
    }

    items[[CleanRdText(item_Name$text)]] <- CleanRdText(item_Description$text)
    position <- item_Description$next_at
  }

  items
}

ExtractAliases <- function(text_Rd) {
  aliases <- character()
  position <- 1L

  repeat {
    alias <- ReadMacroContents(text_Rd, "alias", start_at = position)
    if (is.null(alias)) {
      break
    }
    aliases <- c(aliases, CleanRdText(alias$text))
    position <- alias$next_at
  }

  unique(aliases)
}

ReadRdEntry <- function(fn_Rd) {
  text_Rd <- paste(readLines(fn_Rd, warn = FALSE), collapse = "\n")
  usage <- ReadMacroContents(text_Rd, "usage")
  description <- ReadMacroContents(text_Rd, "description")
  value <- ReadMacroContents(text_Rd, "value")
  seealso <- ReadMacroContents(text_Rd, "seealso")
  arguments <- ReadMacroContents(text_Rd, "arguments")

  usage_Current <- if (is.null(usage)) "Not documented." else {
    usage_Lines <- strsplit(usage$text, "\n", fixed = TRUE)[[1]]
    usage_Lines <- usage_Lines[!grepl("lifecycle::deprecated", usage_Lines, fixed = TRUE)]
    closing_Lines <- which(trimws(usage_Lines) == ")")
    for (closing_Line in closing_Lines) {
      previous_Line <- closing_Line - 1L
      if (previous_Line > 0L) {
        usage_Lines[[previous_Line]] <- sub(",[[:space:]]*$", "", usage_Lines[[previous_Line]])
      }
    }
    trimws(paste(usage_Lines, collapse = "\n"))
  }

  list(
    aliases = ExtractAliases(text_Rd),
    title = CleanRdText(ReadMacroContents(text_Rd, "title")$text),
    usage = usage_Current,
    description = if (is.null(description)) "Not documented." else CleanRdText(description$text),
    arguments = if (is.null(arguments)) list() else ExtractRdItems(arguments$text),
    value = if (is.null(value)) "Not documented." else CleanRdText(value$text),
    seealso = if (is.null(seealso)) "None documented." else CleanRdText(seealso$text),
    deprecated = grepl("deprecated", text_Rd, ignore.case = TRUE)
  )
}

BuildSciDataReportRSkillReference <- function(
  package_root = ".",
  output_file = file.path(package_root, "skills", "scidatareportr", "references", "api-reference.md")
) {
  package_root <- normalizePath(package_root, mustWork = TRUE)
  fn_Namespace <- file.path(package_root, "NAMESPACE")
  dir_Rd <- file.path(package_root, "man")
  fn_Description <- file.path(package_root, "DESCRIPTION")

  if (!file.exists(fn_Namespace) || !dir.exists(dir_Rd) || !file.exists(fn_Description)) {
    stop("package_root must contain DESCRIPTION, NAMESPACE, and man/.", call. = FALSE)
  }

  exports <- ExtractNamespaceExports(fn_Namespace)
  fn_Rd <- list.files(dir_Rd, pattern = "\\.Rd$", full.names = TRUE)
  entries_Rd <- lapply(fn_Rd, ReadRdEntry)
  names(entries_Rd) <- fn_Rd
  aliases_Rd <- unlist(lapply(entries_Rd, `[[`, "aliases"), use.names = FALSE)
  missing_Documentation <- setdiff(exports, aliases_Rd)

  if (length(missing_Documentation) > 0L) {
    stop(
      "Exported functions without Rd documentation: ",
      paste(missing_Documentation, collapse = ", "),
      call. = FALSE
    )
  }

  description <- read.dcf(fn_Description)
  package_Version <- description[1, "Version"]
  lines_Output <- c(
    "# SciDataReportR API reference",
    "",
    paste0("Generated from `NAMESPACE` and `man/*.Rd` for SciDataReportR ", package_Version, "."),
    "Run `Rscript tools/build_scidatareportr_skill_reference.R` after changing exports or Rd documentation.",
    "",
    "This catalog has one entry per public export. Usage blocks omit arguments explicitly marked `lifecycle::deprecated()`; consult the compatibility guide for migration help.",
    ""
  )

  for (export in exports) {
    index_Rd <- which(vapply(entries_Rd, function(entry) export %in% entry$aliases, logical(1)))[[1]]
    entry <- entries_Rd[[index_Rd]]
    aliases_Exported <- setdiff(intersect(entry$aliases, exports), export)
    status_Deprecated <- if (entry$deprecated) "Contains deprecated compatibility interface(s); use the current usage below." else "Current documented interface."

    lines_Output <- c(
      lines_Output,
      paste0("## `", export, "`"),
      "",
      paste0("**Purpose:** ", entry$title),
      "",
      "**Canonical usage**",
      "```r",
      entry$usage,
      "```",
      "",
      paste0("**Description:** ", entry$description),
      "",
      paste0("**Deprecation status:** ", status_Deprecated),
      "",
      paste0("**Related exported aliases:** ", if (length(aliases_Exported) == 0L) "None." else paste(paste0("`", aliases_Exported, "`"), collapse = ", ")),
      "",
      "**Arguments**"
    )

    if (length(entry$arguments) == 0L) {
      lines_Output <- c(lines_Output, "- None documented.")
    } else {
      lines_Output <- c(
        lines_Output,
        vapply(
          names(entry$arguments),
          function(argument) paste0("- `", argument, "`: ", entry$arguments[[argument]]),
          character(1)
        )
      )
    }

    lines_Output <- c(
      lines_Output,
      "",
      paste0("**Returns:** ", entry$value),
      "",
      paste0("**See also:** ", entry$seealso),
      ""
    )
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines_Output, output_file, useBytes = TRUE)

  invisible(list(
    exports = exports,
    output_file = normalizePath(output_file),
    package_version = package_Version
  ))
}

if (sys.nframe() == 0L) {
  BuildSciDataReportRSkillReference()
}
