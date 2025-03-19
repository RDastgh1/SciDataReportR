#' Use the EDATemplate Quarto Template
#'
#' Copies the EDATemplate template to the working directory.
#' @param filename The name to save the Quarto file as (default: "EDA_Report.qmd").
#' @export
use_EDATemplate <- function(filename = "Scripts/script_GetStarted.R") {
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
}
