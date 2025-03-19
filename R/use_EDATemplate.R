#' Use the EDATemplate Quarto Template
#'
#' Copies the EDATemplate template to the working directory.
#' @param filename The name to save the Quarto file as (default: "EDA_Report.qmd").
#' @export
use_EDATemplate <- function(filename = "EDA_Report.qmd") {
  template_path <- system.file("rmarkdown/templates/EDATemplate/skeleton/EDATemplate.qmd", package = "yourpackage")
  if (template_path == "") {
    stop("Template not found. Please reinstall the package.", call. = FALSE)
  }
  file.copy(template_path, filename, overwrite = FALSE)
  message("Template copied to ", filename)
}
