GetOneSidedMainRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}

test_that("MakeComparisonTable keeps two-sided defaults and respects factor-order direction", {
  set.seed(1601)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 35), levels = c("Control", "Case")),
    outcome = c(stats::rnorm(35), stats::rnorm(35, 1.1))
  )

  tbl_default <- MakeComparisonTable(df_Test, group_var = "group", variables = "outcome")
  tbl_directed <- MakeComparisonTable(
    df_Test,
    group_var = "group",
    variables = "outcome",
    alternative = "greater"
  )

  expected_two_sided <- stats::t.test(outcome ~ group, data = df_Test)$p.value
  expected_greater <- stats::t.test(
    df_Test$outcome[df_Test$group == "Case"],
    df_Test$outcome[df_Test$group == "Control"],
    alternative = "greater"
  )$p.value

  expect_equal(GetOneSidedMainRow(tbl_default, "outcome")$p.value, expected_two_sided, tolerance = 1e-12)
  expect_equal(GetOneSidedMainRow(tbl_directed, "outcome")$p.value, expected_greater, tolerance = 1e-12)
  expect_match(GetOneSidedMainRow(tbl_directed, "outcome")$Notes, "Case > Control", fixed = TRUE)
})

test_that("binary directional testing uses the second outcome level as the event", {
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 30), levels = c("Control", "Case")),
    outcome = factor(
      c(rep("No", 24), rep("Yes", 6), rep("No", 9), rep("Yes", 21)),
      levels = c("No", "Yes")
    )
  )

  tbl <- MakeComparisonTable(
    df_Test,
    group_var = "group",
    variables = "outcome",
    CatMethod = "fisher",
    alternative = "greater"
  )
  expected <- stats::fisher.test(
    matrix(c(6, 21, 24, 9), nrow = 2, byrow = TRUE),
    alternative = "less"
  )$p.value

  row <- GetOneSidedMainRow(tbl, "outcome")
  expect_equal(row$p.value, expected, tolerance = 1e-12)
  expect_match(row$Test, "One-sided Fisher", fixed = TRUE)
  expect_match(row$Notes, "Case > Control", fixed = TRUE)
})

test_that("adjusted continuous and binary tests use directed group coefficients", {
  set.seed(1602)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 100), levels = c("Control", "Case")),
    age = stats::rnorm(200)
  )
  df_Test$outcome <- 0.8 * (df_Test$group == "Case") + 0.5 * df_Test$age + stats::rnorm(200)
  df_Test$binary <- factor(
    ifelse(stats::rbinom(200, 1, stats::plogis(-0.2 + 1.1 * (df_Test$group == "Case") + 0.3 * df_Test$age)) == 1, "Yes", "No"),
    levels = c("No", "Yes")
  )

  tbl <- MakeComparisonTable(
    df_Test,
    group_var = "group",
    variables = c("outcome", "binary"),
    covariates = "age",
    alternative = "greater"
  )
  fit_continuous <- stats::lm(outcome ~ group + age, data = df_Test)
  fit_binary <- stats::glm(binary ~ group + age, data = df_Test, family = stats::binomial())

  expect_equal(
    unname(GetOneSidedMainRow(tbl, "outcome")$p.value),
    stats::pt(summary(fit_continuous)$coefficients["groupCase", "t value"], df = stats::df.residual(fit_continuous), lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    unname(GetOneSidedMainRow(tbl, "binary")$p.value),
    stats::pnorm(summary(fit_binary)$coefficients["groupCase", "z value"], lower.tail = FALSE),
    tolerance = 1e-12
  )
})

test_that("one-sided requests retain global inference for multi-group comparisons", {
  set.seed(1603)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 20)),
    outcome = stats::rnorm(60)
  )

  expect_warning(
    tbl <- MakeComparisonTable(df_Test, group_var = "group", variables = "outcome", alternative = "greater"),
    "only available for two-group"
  )
  expect_equal(
    GetOneSidedMainRow(tbl, "outcome")$p.value,
    summary(stats::aov(outcome ~ group, data = df_Test))[[1]]["group", "Pr(>F)"],
    tolerance = 1e-12
  )
})

test_that("PlotSplitViolin exposes its directional comparison metadata", {
  set.seed(1604)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 25), levels = c("Control", "Case")),
    outcome = c(stats::rnorm(25), stats::rnorm(25, 1))
  )

  p <- PlotSplitViolin(df_Test, outcome, group, alternative = "greater")
  comparison <- attr(p, "comparison")

  expect_identical(comparison$Alternative, "greater")
  expect_identical(comparison$Contrast, "Case > Control")
  expect_equal(
    unname(comparison$PValue),
    stats::t.test(df_Test$outcome[df_Test$group == "Case"], df_Test$outcome[df_Test$group == "Control"], alternative = "greater")$p.value,
    tolerance = 1e-12
  )
})
