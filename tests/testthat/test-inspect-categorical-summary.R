test_that("InspectCategoricalSummary summarizes common categorical types", {
  df <- data.frame(
    character_var = c("A", "B", "A", NA),
    factor_var = factor(c("low", "high", "low", NA), levels = c("low", "high", "mid")),
    logical_var = c(TRUE, FALSE, TRUE, NA),
    numeric_var = c(1, 2, 3, 4)
  )

  result <- InspectCategoricalSummary(df, Plot = FALSE)

  expect_named(result, c("Summary", "Plot"))
  expect_s3_class(result$Summary, "tbl_df")
  expect_null(result$Plot)
  expect_setequal(
    unique(result$Summary$Variable),
    c("character_var", "factor_var", "logical_var")
  )
})

test_that("missing values can be included or excluded", {
  df <- data.frame(group = c("A", "A", "B", NA))

  included <- InspectCategoricalSummary(df, Variables = "group", Plot = FALSE)$Summary
  excluded <- InspectCategoricalSummary(
    df,
    Variables = "group",
    IncludeMissing = FALSE,
    Plot = FALSE
  )$Summary

  expect_true(any(included$Missing))
  expect_equal(included$TotalN[1], 4)
  expect_equal(included$Percent[included$Level == "A"], 0.5)
  expect_false(any(excluded$Missing))
  expect_equal(excluded$TotalN[1], 3)
  expect_equal(excluded$Percent[excluded$Level == "A"], 2 / 3)
})

test_that("codebook labels are applied with fallback to variable names", {
  df <- data.frame(group = c("A", "B"), status = c("yes", "no"))
  codebook <- data.frame(
    Variable = "group",
    Label = "Treatment group"
  )

  summary <- InspectCategoricalSummary(df, Codebook = codebook, Plot = FALSE)$Summary

  expect_equal(unique(summary$Label[summary$Variable == "group"]), "Treatment group")
  expect_equal(unique(summary$Label[summary$Variable == "status"]), "status")
})

test_that("factor level order can be preserved", {
  df <- data.frame(
    rating = factor(c("medium", "low", "high"), levels = c("low", "medium", "high"))
  )

  summary <- InspectCategoricalSummary(
    df,
    Variables = "rating",
    SortLevelsBy = "None",
    Plot = FALSE
  )$Summary

  expect_equal(summary$Level, c("low", "medium", "high"))
})

test_that("high-cardinality variables collapse to other in plots", {
  df <- data.frame(group = paste0("level_", 1:6))

  p <- InspectCategoricalSummary(
    df,
    Variables = "group",
    MaxLevels = 3,
    Plot = TRUE
  )$Plot

  built <- ggplot2::ggplot_build(p)
  expect_s3_class(p, "ggplot")
  expect_true(any(built$plot$data$LevelLabelPlot == "(Other)"))
})

test_that("Plot TRUE returns ggplot and Plot FALSE returns NULL", {
  df <- data.frame(group = c("A", "B", "A"))

  expect_s3_class(
    InspectCategoricalSummary(df, Variables = "group", Plot = TRUE)$Plot,
    "ggplot"
  )
  expect_null(
    InspectCategoricalSummary(df, Variables = "group", Plot = FALSE)$Plot
  )
})

test_that("validation errors are clear", {
  df <- data.frame(group = c("A", "B"))

  expect_error(
    InspectCategoricalSummary(df, Variables = "missing", Plot = FALSE),
    "Variables not found"
  )
  expect_error(
    InspectCategoricalSummary(df, Codebook = data.frame(Name = "group"), Plot = FALSE),
    "`Codebook` must contain"
  )
})

test_that("labelled value labels are used when available", {
  df <- data.frame(score = c(1, 2, 1, NA))
  df$score <- labelled::labelled(df$score, labels = c(Low = 1, High = 2))

  summary <- InspectCategoricalSummary(df, Variables = "score", Plot = FALSE)$Summary

  expect_equal(summary$LevelLabel[summary$Level == "1"], "Low")
  expect_equal(summary$LevelLabel[summary$Level == "2"], "High")
})
