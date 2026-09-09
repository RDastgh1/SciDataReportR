GetSkillTestRoot <- function() {
  dirs_Candidates <- unique(c(
    testthat::test_path("..", ".."),
    getwd()
  ))
  dirs_Candidates <- dirs_Candidates[dir.exists(dirs_Candidates)]
  dirs_WithSkillSource <- dirs_Candidates[
    file.exists(file.path(dirs_Candidates, "tools", "build_scidatareportr_skill_reference.R")) &
      dir.exists(file.path(dirs_Candidates, "skills", "scidatareportr"))
  ]

  if (length(dirs_WithSkillSource) == 0L) {
    skip("Repository skill source is not included in this package build.")
  }

  normalizePath(dirs_WithSkillSource[[1]], mustWork = TRUE)
}

GetApiEntry <- function(lines_Api, export) {
  heading <- paste0("## `", export, "`")
  starts <- which(lines_Api == heading)
  if (length(starts) != 1L) {
    stop("Expected one API entry for ", export, call. = FALSE)
  }
  next_Heading <- which(grepl("^## `", lines_Api) & seq_along(lines_Api) > starts)[[1]]
  lines_Api[starts:(next_Heading - 1L)]
}

GetCanonicalUsage <- function(entry_Api) {
  start_Usage <- which(entry_Api == "```r")[[1]]
  end_Usage <- which(entry_Api == "```" & seq_along(entry_Api) > start_Usage)[[1]]
  entry_Api[(start_Usage + 1L):(end_Usage - 1L)]
}

test_that("skill generator catalogs every public export exactly once", {
  dir_Root <- GetSkillTestRoot()
  source(file.path(dir_Root, "tools", "build_scidatareportr_skill_reference.R"))
  fn_Api <- tempfile(fileext = ".md")

  result <- BuildSciDataReportRSkillReference(dir_Root, fn_Api)
  lines_Api <- readLines(fn_Api, warn = FALSE)
  headings <- sub("^## `([^`]+)`$", "\\1", grep("^## `", lines_Api, value = TRUE))

  expect_equal(sort(headings), result$exports)
  expect_equal(length(headings), length(unique(headings)))
  expect_true(all(nzchar(result$exports)))
  expect_match(
    paste(lines_Api, collapse = "\n"),
    paste0("SciDataReportR ", result$package_version),
    fixed = TRUE
  )
})

test_that("skill generator preserves current calls and exposes deprecation guidance", {
  dir_Root <- GetSkillTestRoot()
  source(file.path(dir_Root, "tools", "build_scidatareportr_skill_reference.R"))
  fn_Api <- tempfile(fileext = ".md")

  BuildSciDataReportRSkillReference(dir_Root, fn_Api)
  lines_Api <- readLines(fn_Api, warn = FALSE)
  entry_Comparison <- GetApiEntry(lines_Api, "MakeComparisonTable")
  entry_Correlations <- GetApiEntry(lines_Api, "PlotCorrelationsHeatmap")
  fn_Compatibility <- file.path(dir_Root, "skills", "scidatareportr", "references", "compatibility.md")
  compatibility <- paste(readLines(fn_Compatibility, warn = FALSE), collapse = "\n")

  usage_Comparison <- GetCanonicalUsage(entry_Comparison)
  expect_match(paste(usage_Comparison, collapse = "\n"), "group_var = NULL", fixed = TRUE)
  expect_match(paste(usage_Comparison, collapse = "\n"), "TreatOrdinalAs = \"Categorical\"\n)", fixed = TRUE)
  expect_false(any(grepl("DataFrame = lifecycle::deprecated", usage_Comparison, fixed = TRUE)))
  expect_match(paste(entry_Correlations, collapse = "\n"), "predictor_vars = NULL", fixed = TRUE)
  expect_match(compatibility, "`DataFrame`, `CompVariable`", fixed = TRUE)
  expect_match(compatibility, "`xVars`, `yVars`", fixed = TRUE)
})

test_that("skill installer validates and installs an isolated copy", {
  dir_Root <- GetSkillTestRoot()
  dir_Target <- file.path(tempdir(), "scidatareportr-skill-test", "scidatareportr")
  fn_Installer <- file.path(dir_Root, "tools", "install_scidatareportr_skill.R")
  fn_Rscript <- file.path(R.home("bin"), "Rscript")

  output <- system2(
    fn_Rscript,
    c(shQuote(fn_Installer), "--destination", shQuote(dir_Target)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")

  expect_true(is.null(status) || status == 0L)
  expect_true(file.exists(file.path(dir_Target, "SKILL.md")))
  expect_true(file.exists(file.path(dir_Target, "references", "api-reference.md")))
  expect_match(paste(output, collapse = "\n"), "Installed SciDataReportR skill from package version", fixed = TRUE)
})
