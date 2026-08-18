test_that("PlotMissingData preserves the default observation axis", {
  df_Test <- data.frame(
    id = 1:4,
    MarkerOne = c(1, NA, 3, 4),
    MarkerTwo = c(NA, 2, 3, 4)
  )

  p <- PlotMissingData(df_Test, variables = c("MarkerOne", "MarkerTwo"), Relabel = FALSE)

  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$get_scales("x")$name, "Observations")
  expect_equal(sort(unique(p$data$.x_plot)), 1:4)
})

test_that("PlotMissingData supports labelled x variables", {
  df_Test <- data.frame(
    FACS_date = as.Date("2026-01-01") + 0:3,
    Group = c("A", "B", NA, "A"),
    Marker = c(1, NA, 3, 4)
  )
  sjlabelled::set_label(df_Test$FACS_date) <- "FACS date"

  p_date <- PlotMissingData(df_Test, x_var = "FACS_date", Relabel = FALSE)
  p_group <- PlotMissingData(df_Test, x_var = "Group", Relabel = FALSE)

  expect_equal(p_date$scales$get_scales("x")$name, "FACS date")
  expect_true(inherits(p_date$data$.x_plot, "Date"))
  expect_true(is.factor(p_group$data$.x_plot))
  expect_true("Missing" %in% levels(p_group$data$.x_plot))
  expect_false("FACS_date" %in% as.character(p_date$data$variable))
})

test_that("PlotMissingData validates x_var and facet_by", {
  df_Test <- data.frame(Group = c("A", "B"), Marker = c(1, NA))

  expect_error(PlotMissingData(df_Test, x_var = "Unknown"), "x_var was not found")
  expect_error(PlotMissingData(df_Test, facet_by = c("Group", "Marker")), "facet_by must be a single")
})

test_that("PlotMissingData facets with subgroup-specific percentages", {
  df_Test <- data.frame(
    Visit = 1:4,
    Site = c("A", "A", NA, NA),
    Marker = c(1, NA, NA, NA)
  )

  p <- PlotMissingData(
    df_Test,
    x_var = "Visit",
    facet_by = "Site",
    Relabel = FALSE
  )

  expect_s3_class(p$facet, "FacetWrap")
  expect_true(all(c("A", "Missing") %in% unique(p$data$.facet_value)))
  expect_true(any(grepl("50%", as.character(p$data$.variable_display), fixed = TRUE)))
  expect_true(any(grepl("100%", as.character(p$data$.variable_display), fixed = TRUE)))
  expect_false("Visit" %in% as.character(p$data$variable))
  expect_false("Site" %in% as.character(p$data$variable))
})
