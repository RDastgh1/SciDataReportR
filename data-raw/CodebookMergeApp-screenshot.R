# =============================================================================
# Regenerate the CodebookMergeApp reference screenshots
# =============================================================================
#
# CodebookMergeApp() launches an interactive Shiny dashboard, so its reference
# page cannot render a live example. Instead we capture each dashboard tab and
# commit the static screenshots to man/figures/, where the roxygen `\figure{}`
# directives in R/MergeCodebooks.R embed them on the reference page.
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

if (!exists("CodebookMergeApp", mode = "function")) {
  library(SciDataReportR)
}

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

# Capture the Overview, Harmonization, and Export views of the dashboard.
out_dir <- "man/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
overview_png <- file.path(out_dir, "CodebookMergeApp.png")
harmonization_png <- file.path(out_dir, "CodebookMergeApp-harmonization.png")
export_png <- file.path(out_dir, "CodebookMergeApp-export.png")

driver <- shinytest2::AppDriver$new(
  app,
  name           = "CodebookMergeApp-screenshots",
  width          = 1200,
  height         = 800,
  load_timeout   = 30000
)
on.exit(driver$stop(), add = TRUE)

# Use an HTML selector rather than shinytest2's symbolic selectors, which are
# not compatible with newer chromote versions.
viewport <- list(selector = "html")
driver$wait_for_idle()
driver$get_screenshot(overview_png, screenshot_args = viewport)

# Select the first conflict so the resolution controls are illustrated.
driver$set_inputs(codebook_merge_tabs = "Harmonization")
driver$wait_for_idle()
driver$click(selector = "#harmonization_table tbody tr:first-child")
driver$wait_for_idle()
driver$get_screenshot(harmonization_png, screenshot_args = viewport)

# Save the default observed resolution and generate reproducible merge rules.
driver$click("save_resolution")
driver$set_inputs(codebook_merge_tabs = "Export")
driver$wait_for_idle()
driver$click("generate_rules")
driver$wait_for_idle()
driver$get_screenshot(export_png, screenshot_args = viewport)

message("Wrote ", overview_png)
message("Wrote ", harmonization_png)
message("Wrote ", export_png)
