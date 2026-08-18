# =============================================================================
# Regenerate interactive dashboard reference screenshots
# =============================================================================
#
# The ExploreDatasetComparison() and ExploreMergeValidation() reference pages
# use static previews because pkgdown cannot host their live htmltools/reactable
# dashboards. This script writes the previews to man/figures/.
#
# Run from the package root with:
#
#   source("data-raw/ExploreDashboards-screenshot.R")
#
# Requirements (development only): chromote and a runnable headless Chrome.
# =============================================================================

if (!requireNamespace("chromote", quietly = TRUE)) {
  stop("chromote is required to regenerate dashboard screenshots. ",
       "Install it with install.packages('chromote').")
}

if (!exists("ExploreDatasetComparison", mode = "function")) {
  library(SciDataReportR)
}

CaptureDashboard <- function(dashboard, filename) {
  html_file <- tempfile("scidr-dashboard-", fileext = ".html")
  on.exit(unlink(html_file), add = TRUE)

  htmltools::save_html(dashboard, html_file, background = "white")

  browser <- chromote::Chromote$new()
  session <- browser$new_session(width = 1200, height = 800)

  tryCatch(
    {
      session$go_to(paste0("file://", normalizePath(html_file, winslash = "/")), delay = 1)
      session$screenshot(
        filename = filename,
        cliprect = c(0, 0, 1200, 800),
        delay = 1
      )
    },
    finally = {
      try(session$close(), silent = TRUE)
      browser$close()
    }
  )
}

out_dir <- "man/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Dataset-comparison example: two removed records, one added record, one added
# variable, and changes to two biomarkers.
data("SampleData", package = "SciDataReportR")

old_data <- cbind(id = seq_len(nrow(SampleData)), SampleData)
new_data <- old_data[-c(1, 2), ]
new_data$MMP7[1:12] <- new_data$MMP7[1:12] * 1.15
new_data$tau[20:35] <- new_data$tau[20:35] + 5
new_data$QualityReview <- ifelse(seq_len(nrow(new_data)) %% 3 == 0, "Review", "Pass")
new_data$Smoker <- NULL
new_data <- rbind(
  new_data,
  transform(new_data[1, ], id = max(old_data$id) + 1, QualityReview = "New record")
)

comparison <- CompareDatasets(old_data, new_data, keys = "id")
CaptureDashboard(
  ExploreDatasetComparison(comparison, TopN = 8),
  file.path(out_dir, "ExploreDatasetComparison.png")
)

# Merge-validation example: a clean one-to-one merge with complete key coverage.
set.seed(1)
left <- data.frame(id = 1:50, x = rnorm(50))
right <- data.frame(id = 1:50, y = rnorm(50))
merged <- merge(left, right, by = "id")
validation <- ValidateMerge(left, right, merged, keys = "id")
CaptureDashboard(
  ExploreMergeValidation(validation),
  file.path(out_dir, "ExploreMergeValidation.png")
)

message("Wrote dashboard screenshots to ", out_dir)
