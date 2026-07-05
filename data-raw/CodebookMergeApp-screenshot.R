# =============================================================================
# Regenerate the CodebookMergeApp reference screenshot
# =============================================================================
#
# CodebookMergeApp() launches an interactive Shiny dashboard, so its reference
# page cannot render a live example. Instead we capture a static screenshot of
# the running app and commit it to man/figures/, where the roxygen `\figure{}`
# directive in R/MergeCodebooks.R embeds it on the reference page.
#
# Run this script whenever the app UI changes:
#
#   source("data-raw/CodebookMergeApp-screenshot.R")
#
# Requirements (all Suggests / dev-only): shinytest2, and a headless Chrome
# available to chromote/webshot2. Install with:
#
#   install.packages(c("shinytest2", "chromote"))
#
# =============================================================================

if (!requireNamespace("shinytest2", quietly = TRUE)) {
  stop("shinytest2 is required to regenerate the screenshot. ",
       "Install it with install.packages('shinytest2').")
}

library(SciDataReportR)

# Build two small, overlapping codebooks to give the app something to harmonize.
data("SampleVariableTypes", package = "SciDataReportR")

base_cb <- SampleVariableTypes[
  seq_len(min(12, nrow(SampleVariableTypes))),
  c("Variable", "Label", "Type")
]

# A second codebook that overlaps on some variables and uses type synonyms,
# so the "Variable Presence" and "Type Conflicts" panels are populated.
alt_cb <- base_cb
alt_cb$Type[alt_cb$Type == "Double"] <- "numeric"
alt_cb$Type[alt_cb$Type == "Categorical"] <- "factor"
alt_cb <- alt_cb[-c(1, 2), ]              # drop a couple to create presence gaps
alt_cb$Label <- paste0(alt_cb$Label, " (v2)")

codebooks <- list(CohortA = base_cb, CohortB = alt_cb)

app <- CodebookMergeApp(codebooks, VariableCol = "Variable")

# Capture the initial view of the dashboard.
out_dir <- "man/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_png <- file.path(out_dir, "CodebookMergeApp.png")

driver <- shinytest2::AppDriver$new(
  app,
  name           = "CodebookMergeApp",
  width          = 1200,
  height         = 800,
  load_timeout   = 30000,
  screenshot_args = list(selector = "viewport")
)
on.exit(driver$stop(), add = TRUE)

driver$get_screenshot(out_png)

message("Wrote ", out_png)
