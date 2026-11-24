# Shared validation, dependency reporting, and console display for anomaly objects.

ScidrAnomalyPercentile <- function(flag_percentile) {
  if (is.null(flag_percentile)) return(NULL)
  if (!is.numeric(flag_percentile) || length(flag_percentile) != 1L || !is.finite(flag_percentile) ||
      flag_percentile <= 0 || flag_percentile >= 1)
    stop("flag_percentile must be NULL or a single number strictly between 0 and 1.", call. = FALSE)
  as.numeric(flag_percentile)
}

# Report, rather than discover halfway through a long run, which optional packages are absent.
ScidrAnomalyDependencies <- function(methods = character(), explorer = FALSE, workbook = FALSE,
    plots = FALSE, explain = FALSE) {
  needs <- list()
  Need <- function(package, purpose, required) {
    if (requireNamespace(package, quietly = TRUE)) return(invisible(NULL))
    needs[[length(needs) + 1L]] <<- tibble::tibble(Package = package, Purpose = purpose, Required = required)
  }
  backends <- c(conditional = "outliertree", isolation = "isotree", local = "dbscan",
    distance = "robustbase", cellwise = "cellWise", trajectory = "nlme")
  for (method in intersect(methods, names(backends)))
    Need(backends[[method]], paste0("the '", method, "' method"), TRUE)
  if (explain) Need("fastshap", "SHAP explanations (sensitivity needs no package)", FALSE)
  if (explorer) for (package in c("DT", "htmltools")) Need(package, "the interactive explorer", TRUE)
  if (plots) Need("plotly", "interactive figures inside the explorer", FALSE)
  if (explorer) Need("rmarkdown", "inlining dependencies into one standalone HTML file", FALSE)
  if (workbook) for (package in c("openxlsx", "jsonlite")) Need(package, "the Excel review workbook", TRUE)
  dplyr::bind_rows(c(list(tibble::tibble(Package = character(), Purpose = character(),
    Required = logical())), needs))
}

ScidrAnomalyRequireDependencies <- function(df_Needs) {
  Describe <- function(df_Rows) paste0("  ", df_Rows$Package, " - ", df_Rows$Purpose, collapse = "\n")
  optional <- df_Needs[!df_Needs$Required, , drop = FALSE]
  required <- df_Needs[df_Needs$Required, , drop = FALSE]
  if (nrow(required)) stop("Install these packages first:\n", Describe(required), call. = FALSE)
  if (nrow(optional)) message("These optional packages are missing; those features are skipped:\n", Describe(optional))
  invisible(df_Needs)
}

ScidrAnomalyPrintFrame <- function(df_Data) {
  df_Data <- as.data.frame(df_Data)
  if (!nrow(df_Data)) return(invisible(NULL))
  print(df_Data, row.names = FALSE, right = FALSE)
  invisible(df_Data)
}

#' @export
print.DataAnomalyResult <- function(x, ...) {
  kind <- if (inherits(x, "LongitudinalAnomalyResult")) "Longitudinal" else "Cross-sectional"
  cat(kind, "data anomaly screen\n")
  cat("  Records screened: ", nrow(x$ReviewData), "\n", sep = "")
  cat("  Variables screened: ", length(x$Metadata$Variables), "\n", sep = "")
  cat("  Methods: ", paste(x$Metadata$Methods, collapse = ", "), "\n", sep = "")
  if (length(x$Metadata$IDVariable)) cat("  ID column: ", x$Metadata$IDVariable, "\n", sep = "")
  if (length(x$Metadata$StrataVariables)) cat("  Strata: ", paste(x$Metadata$StrataVariables, collapse = ", "), "\n", sep = "")
  cat("\nFlags per method and stratum:\n")
  if (is.null(x$Summary) || !nrow(x$Summary)) cat("  No records evaluated.\n") else ScidrAnomalyPrintFrame(x$Summary)
  cat("\nCutoffs: ", x$Metadata$ThresholdProvenance, "\n", sep = "")
  cat("Flags are review candidates, not confirmed errors.\n")
  diagnostics <- if (is.null(x$Diagnostics)) 0L else nrow(x$Diagnostics)
  if (diagnostics) cat("Diagnostics: ", diagnostics, " note(s) in $Diagnostics.\n", sep = "")
  cat("\nNext: CreateDataAnomalyReport() for files to send out, or ExploreDataAnomalies() to view here.\n")
  invisible(x)
}

#' @export
summary.DataAnomalyResult <- function(object, ...) {
  df_Summary <- object$Summary
  if (is.null(df_Summary) || !nrow(df_Summary)) return(df_Summary)
  df_Results <- object$Results
  df_Summary$Threshold <- vapply(seq_len(nrow(df_Summary)), function(i) {
    values <- df_Results$Threshold[df_Results$Method == df_Summary$Method[i] &
      (df_Results$Stratum == df_Summary$Stratum[i] | (is.na(df_Results$Stratum) & is.na(df_Summary$Stratum[i])))]
    values <- unique(values[is.finite(values)])
    if (length(values) == 1L) values else NA_real_
  }, numeric(1))
  # A native rule flags without a numeric threshold, so classify by whether the
  # method produced any flag at all rather than by the threshold column alone.
  classified <- vapply(seq_len(nrow(df_Summary)), function(i)
    any(!is.na(df_Results$Flag[df_Results$Method == df_Summary$Method[i] &
      (df_Results$Stratum == df_Summary$Stratum[i] | (is.na(df_Results$Stratum) & is.na(df_Summary$Stratum[i])))])),
    logical(1))
  df_Summary$CutoffSource <- ifelse(df_Summary$Method %in% names(object$Metadata$Thresholds), "User cutoff",
    ifelse(df_Summary$Method %in% object$Metadata$FlagPercentileMethods, "Default percentile cutoff",
      ifelse(classified, "Native cutoff", "No cutoff; scores only")))
  df_Summary$ScoreDefinition <- unname(object$Metadata$ScoreDefinitions[df_Summary$Method])
  df_Summary
}

#' @export
print.AnomalyReviewQueue <- function(x, ...) {
  cat("Anomaly review queue\n")
  cat("  Analyses: ", paste(names(x$Sources), collapse = ", "), "\n", sep = "")
  cat("  Issues: ", nrow(x$Issues), " (evidence rows: ", nrow(x$Evidence), ")\n", sep = "")
  if (nrow(x$Issues) && "ReviewStatus" %in% names(x$Issues)) {
    status <- table(factor(ifelse(is.na(x$Issues$ReviewStatus), "Unreviewed", x$Issues$ReviewStatus),
      levels = unique(c(ScidrReviewStatuses(), stats::na.omit(x$Issues$ReviewStatus)))))
    status <- status[status > 0]
    cat("\nReview status:\n")
    ScidrAnomalyPrintFrame(data.frame(Status = names(status), Issues = as.integer(status)))
  }
  decided <- if (!nrow(x$Issues) || !"ReviewStatus" %in% names(x$Issues)) 0L else
    sum(!is.na(x$Issues$ReviewStatus) & x$Issues$ReviewStatus != "Unreviewed")
  cat("\nDecisions imported: ", if (decided) paste0("yes (", decided, " issue(s))") else "no", "\n", sep = "")
  diagnostics <- if (is.null(x$Diagnostics)) 0L else nrow(x$Diagnostics)
  if (diagnostics) cat("Diagnostics: ", diagnostics, " note(s) in $Diagnostics.\n", sep = "")
  cat("\nNext: ExportAnomalyReviewQueue() to send out, SummarizeAnomalyDecisions() once decisions return.\n")
  invisible(x)
}

# RStudio and Quarto bundle pandoc without putting it on PATH, so a plain R session
# would otherwise fall back to a sidecar dependency directory unnecessarily.
ScidrAnomalyFindPandoc <- function() {
  if (rmarkdown::pandoc_available("1.17")) return(invisible(TRUE))
  candidates <- c(Sys.getenv("RSTUDIO_PANDOC"), Sys.getenv("QUARTO_PANDOC"))
  quarto <- Sys.which("quarto")
  roots <- c(if (nzchar(quarto)) file.path(dirname(dirname(ScidrResolveLink(quarto))), "share"),
    "/Applications/quarto", "/usr/local/opt/quarto",
    "/Applications/RStudio.app/Contents/Resources/app/quarto",
    file.path(Sys.getenv("HOME"), "Applications/RStudio.app/Contents/Resources/app/quarto"))
  for (root in roots[nzchar(roots)]) candidates <- c(candidates,
    Sys.glob(file.path(root, "bin", "tools", "*", "pandoc")), Sys.glob(file.path(root, "bin", "tools", "pandoc")))
  for (candidate in unique(candidates[nzchar(candidates)])) {
    directory <- if (dir.exists(candidate)) candidate else dirname(candidate)
    if (!file.exists(file.path(directory, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"))) next
    rmarkdown::find_pandoc(cache = FALSE, dir = directory)
    if (rmarkdown::pandoc_available("1.17")) return(invisible(TRUE))
  }
  invisible(FALSE)
}

# Sys.readlink() returns "" for a real file; resolve only genuine symlinks.
ScidrResolveLink <- function(path) {
  target <- Sys.readlink(path)
  if (!nzchar(target)) return(path)
  if (grepl("^/", target)) target else file.path(dirname(path), target)
}
