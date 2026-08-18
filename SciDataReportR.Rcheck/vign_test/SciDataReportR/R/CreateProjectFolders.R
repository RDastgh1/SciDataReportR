#' Create Project Folder Structure
#'
#' This function creates a project folder structure with the following directories:
#' - Data/
#'   - Raw/
#'   - Clean/
#' - Scripts/
#' - Reports/
#'
#' @param base_path A character string specifying the base directory where the folders should be created. Defaults to the current working directory.
#' @return A message indicating that the project structure has been created successfully.
#' @export
CreateProjectFolders <- function(base_path = ".") {
  # Define the folder paths
  data_folder <- file.path(base_path, "Data")
  raw_folder <- file.path(data_folder, "Raw")
  clean_folder <- file.path(data_folder, "Clean")
  scripts_folder <- file.path(base_path, "Scripts")
  reports_folder <- file.path(base_path, "Reports")

  # Create the folders
  dir.create(data_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(raw_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(clean_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(scripts_folder, showWarnings = FALSE, recursive = TRUE)
  dir.create(reports_folder, showWarnings = FALSE, recursive = TRUE)

  # Print a message indicating the folders have been created
  message("Project folders created successfully!")
}


