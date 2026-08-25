GetAdjustedMainRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}

test_that("adjusted ANCOVA p-values support non-syntactic two-group names", {
  skip_if_not_installed("car")
  skip_if_not_installed("emmeans")

  set.seed(201)
  df_Test <- data.frame(
    "Group Status" = factor(rep(c("Control", "Case"), each = 30)),
    Age = stats::rnorm(60),
    check.names = FALSE
  )
  df_Test$outcome <- 0.8 * (df_Test[["Group Status"]] == "Case") +
    0.4 * df_Test$Age +
    stats::rnorm(60)

  fit <- stats::lm(outcome ~ `Group Status` + Age, data = df_Test)
  expected <- as.numeric(car::Anova(fit, type = 2)["`Group Status`", "Pr(>F)"])

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "Group Status",
    variables = "outcome",
    covariates = "Age",
    AddEffectSize = TRUE
  )
  got <- GetAdjustedMainRow(tbl, "outcome")

  expect_identical(got$Test, "ANCOVA (Type II)")
  expect_true(is.finite(got$p.value))
  expect_equal(got$p.value, expected, tolerance = 1e-12)
  expect_false(is.na(got$p.value_fmt) || !nzchar(got$p.value_fmt))
  expect_true(is.finite(got$effect_size))
})

test_that("adjusted ANCOVA p-values support non-syntactic multi-group names", {
  skip_if_not_installed("car")

  set.seed(202)
  df_Test <- data.frame(
    "Group Status" = factor(rep(c("Control", "Group B", "Group C"), each = 24)),
    Age = stats::rnorm(72),
    check.names = FALSE
  )
  df_Test$outcome <- c(Control = 0, `Group B` = 0.5, `Group C` = 1.0)[df_Test[["Group Status"]]] +
    0.5 * df_Test$Age +
    stats::rnorm(72)

  fit <- stats::lm(outcome ~ `Group Status` + Age, data = df_Test)
  expected <- as.numeric(car::Anova(fit, type = 2)["`Group Status`", "Pr(>F)"])

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "Group Status",
    variables = "outcome",
    covariates = "Age",
    AddEffectSize = TRUE
  )
  got <- GetAdjustedMainRow(tbl, "outcome")

  expect_identical(got$Test, "ANCOVA (Type II)")
  expect_equal(got$p.value, expected, tolerance = 1e-12)
  expect_false(is.na(got$p.value_fmt) || !nzchar(got$p.value_fmt))
})

test_that("adjusted categorical p-values support non-syntactic group names", {
  skip_if_not_installed("nnet")

  set.seed(203)
  n <- 180
  df_Test <- data.frame(
    "Group Status" = factor(rep(c("Control", "Case"), each = n / 2)),
    Age = stats::rnorm(n),
    check.names = FALSE
  )
  prob_binary <- stats::plogis(
    -0.2 + 0.7 * (df_Test[["Group Status"]] == "Case") + 0.3 * df_Test$Age
  )
  df_Test$binary_outcome <- factor(stats::rbinom(n, 1, prob_binary))
  df_Test$multi_outcome <- factor(sample(
    c("Low", "Middle", "High"),
    size = n,
    replace = TRUE,
    prob = c(0.35, 0.40, 0.25)
  ))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "Group Status",
    variables = c("binary_outcome", "multi_outcome"),
    covariates = "Age"
  )
  got_binary <- GetAdjustedMainRow(tbl, "binary_outcome")
  got_multi <- GetAdjustedMainRow(tbl, "multi_outcome")

  expect_identical(got_binary$Test, "Logistic regression (LR)")
  expect_true(is.finite(got_binary$p.value))
  expect_false(is.na(got_binary$p.value_fmt) || !nzchar(got_binary$p.value_fmt))
  expect_identical(got_multi$Test, "Multinomial LR")
  expect_true(is.finite(got_multi$p.value))
  expect_false(is.na(got_multi$p.value_fmt) || !nzchar(got_multi$p.value_fmt))
})
