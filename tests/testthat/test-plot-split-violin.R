test_that("PlotSplitViolin applies parametric covariate adjustment", {
  skip_if_not_installed("emmeans")

  set.seed(42)
  n <- 200
  df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = n / 2))
  )
  df_Test$age <- stats::rnorm(
    n,
    mean = ifelse(df_Test$group == "B", 2, 0)
  )
  df_Test$outcome <- 2 * df_Test$age + stats::rnorm(n)

  plot_unadjusted <- PlotSplitViolin(df_Test, outcome, group)
  plot_adjusted <- PlotSplitViolin(
    df_Test,
    outcome,
    group,
    covariates = "age"
  )

  result_unadjusted <- attr(plot_unadjusted, "comparison")
  result_adjusted <- attr(plot_adjusted, "comparison")
  expected_fit <- stats::lm(outcome ~ group + age, data = df_Test)
  expected_p <- as.data.frame(
    emmeans::contrast(
      emmeans::emmeans(expected_fit, "group"),
      method = "revpairwise"
    )
  )$p.value[1]

  expect_lt(result_unadjusted$PValue, 0.001)
  expect_gt(result_adjusted$PValue, 0.05)
  expect_equal(result_adjusted$PValue, expected_p, tolerance = 1e-12)
  expect_equal(result_adjusted$Covariates, "age")
})

test_that("PlotSplitViolin uses adjusted residuals for nonparametric tests", {
  skip_if_not_installed("emmeans")

  set.seed(84)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 60)),
    age = c(stats::rnorm(60, 0), stats::rnorm(60, 2))
  )
  df_Test$outcome <- 3 * df_Test$age + stats::rt(120, df = 4)

  plot_adjusted <- PlotSplitViolin(
    df_Test,
    outcome,
    group,
    covariates = "age",
    nonparametric = TRUE,
    p_label = "p_value"
  )
  result <- attr(plot_adjusted, "comparison")
  expected_residuals <- stats::residuals(
    stats::lm(outcome ~ age, data = df_Test)
  )
  expected_p <- stats::wilcox.test(
    expected_residuals ~ df_Test$group,
    exact = FALSE
  )$p.value

  expect_equal(result$PValue, expected_p, tolerance = 1e-12)
  expect_match(result$Method, "adjusted residuals")
})

test_that("PlotSplitViolin exposes exact p-values and validates covariates", {
  skip_if_not_installed("emmeans")

  df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 10)),
    value = c(seq_len(10), seq_len(10) + 2)
  )

  plot_p <- PlotSplitViolin(
    df_Test,
    value,
    group,
    p_label = "both",
    show_ns = TRUE
  )

  expect_true(is.finite(attr(plot_p, "comparison")$PValue))
  expect_error(
    PlotSplitViolin(df_Test, value, group, covariates = "missing_age"),
    "Covariate\\(s\\) not found"
  )
})
