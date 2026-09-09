# Install the repository-owned SciDataReportR skill into a Codex skills directory.

GetRepositoryRoot <- function() {
  argument_File <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument_File) != 1L) {
    stop("Run this script with Rscript.", call. = FALSE)
  }
  fn_Script <- gsub("~\\+~", " ", sub("^--file=", "", argument_File))
  dirname(dirname(normalizePath(fn_Script)))
}

ValidateSciDataReportRSkill <- function(dir_Skill) {
  fn_Skill <- file.path(dir_Skill, "SKILL.md")
  fn_WorkflowMap <- file.path(dir_Skill, "references", "workflow-map.md")
  fn_Compatibility <- file.path(dir_Skill, "references", "compatibility.md")
  fn_Api <- file.path(dir_Skill, "references", "api-reference.md")
  fn_Agent <- file.path(dir_Skill, "agents", "openai.yaml")
  required_Files <- c(fn_Skill, fn_WorkflowMap, fn_Compatibility, fn_Api, fn_Agent)

  if (!all(file.exists(required_Files))) {
    stop("Skill is missing required files: ", paste(required_Files[!file.exists(required_Files)], collapse = ", "), call. = FALSE)
  }

  contents_Skill <- paste(readLines(fn_Skill, warn = FALSE), collapse = "\n")
  if (!grepl("(?m)^name: scidatareportr$", contents_Skill, perl = TRUE) ||
      !grepl("(?m)^description:", contents_Skill, perl = TRUE)) {
    stop("SKILL.md has invalid or incomplete frontmatter.", call. = FALSE)
  }

  contents_Api <- paste(readLines(fn_Api, warn = FALSE), collapse = "\n")
  if (!grepl("# SciDataReportR API reference", contents_Api, fixed = TRUE) ||
      !grepl("## `RevalueData`", contents_Api, fixed = TRUE)) {
    stop("Generated API reference is missing expected content.", call. = FALSE)
  }

  invisible(TRUE)
}

ParseInstallerArguments <- function(arguments) {
  destination <- NULL
  force <- FALSE
  check_Only <- FALSE
  position <- 1L

  while (position <= length(arguments)) {
    argument <- arguments[[position]]
    if (argument == "--destination") {
      position <- position + 1L
      if (position > length(arguments)) {
        stop("--destination requires a path.", call. = FALSE)
      }
      destination <- arguments[[position]]
    } else if (argument == "--force") {
      force <- TRUE
    } else if (argument == "--check-only") {
      check_Only <- TRUE
    } else if (argument == "--help") {
      cat("Usage: Rscript tools/install_scidatareportr_skill.R [--destination PATH] [--force] [--check-only]\n")
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", argument, call. = FALSE)
    }
    position <- position + 1L
  }

  list(destination = destination, force = force, check_only = check_Only)
}

InstallSciDataReportRSkill <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  options_Installer <- ParseInstallerArguments(arguments)
  dir_Repository <- GetRepositoryRoot()
  dir_Source <- file.path(dir_Repository, "skills", "scidatareportr")
  dir_Default <- file.path(path.expand("~"), ".codex", "skills", "scidatareportr")
  dir_Destination <- normalizePath(options_Installer$destination %||% dir_Default, mustWork = FALSE)

  if (basename(dir_Destination) != "scidatareportr") {
    stop("Destination must be the scidatareportr skill directory.", call. = FALSE)
  }

  ValidateSciDataReportRSkill(dir_Source)
  if (options_Installer$check_only) {
    cat("SciDataReportR skill source is valid.\n")
    return(invisible(dir_Source))
  }

  if (dir.exists(dir_Destination) && !options_Installer$force) {
    stop("Destination exists. Re-run with --force to replace it: ", dir_Destination, call. = FALSE)
  }

  dir.create(dirname(dir_Destination), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(dir_Destination)) {
    unlink(dir_Destination, recursive = TRUE, force = TRUE)
  }
  copied <- file.copy(dir_Source, dirname(dir_Destination), recursive = TRUE)
  if (!copied) {
    stop("Could not copy skill to destination.", call. = FALSE)
  }

  ValidateSciDataReportRSkill(dir_Destination)
  package_Version <- read.dcf(file.path(dir_Repository, "DESCRIPTION"))[1, "Version"]
  cat("Installed SciDataReportR skill from package version ", package_Version, " to ", dir_Destination, ".\n", sep = "")
  invisible(dir_Destination)
}

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}

if (sys.nframe() == 0L) {
  InstallSciDataReportRSkill()
}
