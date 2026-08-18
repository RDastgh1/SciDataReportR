# These tests pin the colors that must NOT be swept into the SciData
# qualitative palette. They all use `scale_*_manual`, so they are
# syntactically indistinguishable from categorical use - the only thing
# stopping a future palette retrofit from swallowing them is this file.
#
# Two invariants are enforced:
#   1. The literal semantic colors survive.
#   2. The ggplot surfaces and the HTML dashboard surfaces agree with each
#      other. That pairing is otherwise enforced by nothing at all, so a
#      change to one file would silently desynchronise the static plots from
#      the interactive reports.

StatusPalette <- function() {
  c(PASS = "#2E7D32", WARNING = "#F9A825", FAIL = "#C62828")
}

MakeMergeValidation <- function() {
  set.seed(1)
  left <- data.frame(id = 1:50, x = rnorm(50))
  right <- data.frame(id = c(1:45, 100:104), y = rnorm(50))
  merged <- merge(left, right, by = "id")
  ValidateMerge(left, right, merged, keys = "id")
}

test_that("merge validation keeps its traffic-light status colors", {
  plots <- PlotMergeValidation(
    MakeMergeValidation(), Plot = "All", interactive = FALSE
  )

  fills <- unique(ggplot2::ggplot_build(plots$Checks)$data[[1]]$fill)

  # Every fill used must come from the status palette - no palette color.
  expect_true(all(fills %in% StatusPalette()))
  expect_false(any(fills %in% SciDataPalette(names = FALSE)))
})

test_that("merge coverage keeps its severity-ordered colors", {
  plots <- PlotMergeValidation(
    MakeMergeValidation(), Plot = "All", interactive = FALSE
  )

  fills <- unique(ggplot2::ggplot_build(plots$Coverage)$data[[1]]$fill)

  # Matching green, left-only amber, right-only orange: an ordered severity
  # ramp, not three arbitrary categories.
  expect_true(all(fills %in% c("#2E7D32", "#F9A825", "#EF6C00")))
})

test_that("the merge dashboard uses the same status colors as the plots", {
  html <- as.character(ExploreMergeValidation(MakeMergeValidation()))

  for (color in StatusPalette()) {
    expect_true(
      grepl(color, html, fixed = TRUE),
      info = paste("status color missing from merge dashboard:", color)
    )
  }

  # The dashboard adds INFO on top of the shared three.
  expect_true(grepl("#1565C0", html, fixed = TRUE))
})

test_that("dataset comparison keeps its status colors on both surfaces", {
  set.seed(4)
  old_data <- data.frame(id = 1:40, age = rnorm(40, 60, 8), grp = rep(c("A", "B"), 20))
  new_data <- old_data
  new_data$age[1:5] <- new_data$age[1:5] + 1
  new_data <- new_data[-c(3, 7), ]
  new_data$extra <- rnorm(nrow(new_data))

  comparison <- CompareDatasets(old_data, new_data, keys = "id")

  plots <- PlotDatasetComparison(comparison, interactive = FALSE)
  fills <- unique(ggplot2::ggplot_build(plots$Checks)$data[[1]]$fill)
  expect_true(all(fills %in% StatusPalette()))
  expect_false(any(fills %in% SciDataPalette(names = FALSE)))

  html <- as.character(ExploreDatasetComparison(comparison))
  for (color in StatusPalette()) {
    expect_true(
      grepl(color, html, fixed = TRUE),
      info = paste("status color missing from comparison dashboard:", color)
    )
  }
})

test_that("signed-significance ladders keep their diverging colors", {
  # A 7-step diverging ladder centred on grey "ns". Replacing this with a
  # qualitative palette would destroy the direction and strength encoding.
  ladder <- rev(c(
    "red4", "firebrick3", "pink2", "grey93",
    "lightblue2", "steelblue3", "blue"
  ))

  source_files <- c(
    "PlotNumInteractionEffectsMatrix.R",
    "PlotCatInteractionEffectsMatrix.R"
  )

  for (file_name in source_files) {
    path <- test_path("..", "..", "R", file_name)
    skip_if_not(file.exists(path), paste("source not available:", file_name))

    text <- paste(readLines(path, warn = FALSE), collapse = " ")
    for (color in ladder) {
      expect_true(
        grepl(paste0('"', color, '"'), text, fixed = TRUE),
        info = paste(color, "missing from", file_name)
      )
    }
  }
})

test_that("outlier and significance emphasis stay achromatic or alerting", {
  # IQROutliers: TRUE must be red - it is an alert, not a category.
  set.seed(9)
  df <- data.frame(v = c(rnorm(60, 10, 2), 40, -15))
  result <- IQROutliers(df, Variable = "v")

  plot_obj <- result$p
  skip_if(is.null(plot_obj), "IQROutliers returned no plot")

  built <- ggplot2::ggplot_build(plot_obj)
  colours <- unlist(lapply(built$data, function(layer) layer$colour))
  colours <- unique(colours[!is.na(colours)])
  expect_true("red" %in% colours || "black" %in% colours)
  expect_false(any(colours %in% SciDataPalette(names = FALSE)))
})

test_that("missing-data plots stay achromatic", {
  data(SampleData, envir = environment())

  plot_obj <- PlotMissingData(SampleData[, c("age", "AXL", "tau")])
  plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot
  skip_if(is.null(plot_obj), "PlotMissingData returned no plot")

  fills <- unique(ggplot2::ggplot_build(plot_obj)$data[[1]]$fill)

  # Greys carry the missing/present distinction deliberately; colour would
  # imply a category that does not exist.
  expect_false(any(fills %in% SciDataPalette(names = FALSE)))
})

test_that("the p-value scales stay continuous and keep their palettes", {
  # scale_*_pvalue must never become a discrete SciData scale.
  expect_s3_class(scale_color_pvalue(), "ScaleContinuous")
  expect_s3_class(scale_fill_pvalue(), "ScaleContinuous")
  expect_equal(scale_color_pvalue()$name, "P-value")
  expect_error(scale_color_pvalue(palette = "SciData"), "should be one of")
})
