#' Use the Get Started Script Template
#'
#' Copies the get-started R script template to the working directory.
#' @param filename The name to save the R script as (default:
#'   "Scripts/script_GetStarted.R").
#' @return Invisibly returns `filename` after copying the template. Called
#'   for its side effect of creating the script file.
#' @export
use_GetStartedScript <- function(filename = "Scripts/script_GetStarted.R") {
  template_path <- system.file("rmarkdown/templates/script_GetStarted.R", package = "SciDataReportR")
  if (template_path == "") {
    stop("Template not found. Please reinstall the package.", call. = FALSE)
  }
  file.copy(template_path, filename, overwrite = FALSE)
  message("Template copied to ", filename)

  # Open the file in RStudio if available
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    rstudioapi::navigateToFile(filename)
  } else {
    message("Open the file manually: ", filename)
  }
  invisible(filename)
}
