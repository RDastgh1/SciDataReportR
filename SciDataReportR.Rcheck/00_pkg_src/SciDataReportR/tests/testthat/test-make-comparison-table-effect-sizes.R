GetMainEffectRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}

GetComparisonCaption <- function(tbl) {
  as.character(tbl$table_styling$caption)
}

test_that("two-group parametric comparisons retain absolute Cohen's d", {
  skip_if_not_installed("effectsize")

  set.seed(101)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 24)),
    outcome = c(stats::rnorm(24), stats::rnorm(24, mean = 0.7))
  )

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    AddEffectSize = TRUE
  )

  expected <- abs(
    effectsize::cohens_d(outcome ~ group, data = df_Test)$Cohens_d
  )
  got <- GetMainEffectRow(tbl, "outcome")

  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "|d|")
})

test_that("two-group ANCOVA reports adjusted d from EMM and residual SD", {
  skip_if_not_installed("emmeans")

  set.seed(102)
  df_Test <- data.frame(
    `group status` = factor(
      c(rep("Control group", 25), rep("Case group", 19)),
      levels = c("Control group", "Case group")
    ),
    covariate = stats::rnorm(44),
    check.names = FALSE
  )
  df_Test$outcome <- 0.8 * (df_Test[["group status"]] == "Case group") +
    0.6 * df_Test$covariate +
    stats::rnorm(44)
  df_Test$covariate[c(3, 31)] <- NA_real_

  df_CC <- df_Test[stats::complete.cases(
    df_Test[, c("outcome", "group status", "covariate")]
  ), , drop = FALSE]
  fit <- stats::lm(
    outcome ~ `group status` + covariate,
    data = df_CC
  )
  emm <- emmeans::emmeans(fit, ~ `group status`)
  emm_pair <- summary(emmeans::contrast(emm, method = "pairwise"))
  expected <- abs(as.numeric(emm_pair$estimate[1]) / stats::sigma(fit))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group status",
    variables = "outcome",
    covariates = "covariate",
    AddEffectSize = TRUE
  )
  got <- GetMainEffectRow(tbl, "outcome")

  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "adjusted |d| (EMM / residual SD)")
})

test_that("multi-group ANOVA reports Cohen's f", {
  skip_if_not_installed("effectsize")

  set.seed(103)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), times = c(18, 25, 21))),
    outcome = stats::rnorm(64)
  )
  df_Test$outcome <- df_Test$outcome +
    c(A = 0, B = 0.5, C = 1)[df_Test$group]

  fit <- stats::aov(outcome ~ group, data = df_Test)
  eta <- effectsize::eta_squared(fit, partial = FALSE)$Eta2[1]
  expected <- sqrt(eta / (1 - eta))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    AddEffectSize = TRUE
  )
  got <- GetMainEffectRow(tbl, "outcome")

  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "Cohen's f")
})

test_that("multi-group Type II ANCOVA reports partial Cohen's f", {
  skip_if_not_installed("car")
  skip_if_not_installed("effectsize")

  set.seed(104)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), times = c(20, 27, 16))),
    covariate = stats::rnorm(63)
  )
  df_Test$outcome <- c(A = 0, B = 0.4, C = 1.0)[df_Test$group] +
    0.7 * df_Test$covariate +
    stats::rnorm(63)
  df_Test$covariate[c(2, 29, 61)] <- NA_real_

  df_CC <- df_Test[stats::complete.cases(df_Test), , drop = FALSE]
  fit <- stats::lm(outcome ~ group + covariate, data = df_CC)
  a2 <- car::Anova(fit, type = 2)
  eta <- effectsize::eta_squared(a2, partial = TRUE)
  eta_group <- eta$Eta2_partial[eta$Parameter == "group"]
  expected <- sqrt(eta_group / (1 - eta_group))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    covariates = "covariate",
    AddEffectSize = TRUE
  )
  got <- GetMainEffectRow(tbl, "outcome")

  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "partial Cohen's f")
})

test_that("nonparametric and categorical effect sizes are preserved", {
  set.seed(105)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 18)),
    outcome = stats::rexp(54),
    category = factor(rep(c("No", "Yes"), 27))
  )

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome", "category"),
    Parametric = FALSE,
    AddEffectSize = TRUE
  )

  got_outcome <- GetMainEffectRow(tbl, "outcome")
  got_category <- GetMainEffectRow(tbl, "category")
  H <- as.numeric(stats::kruskal.test(outcome ~ group, data = df_Test)$statistic)
  expected_epsilon <- (H - 3 + 1) / (nrow(df_Test) - 3)

  expect_equal(got_outcome$effect_size, expected_epsilon, tolerance = 1e-10)
  expect_identical(got_outcome$ES_Method, "epsilon-squared")
  expect_identical(got_category$ES_Method, "Cramer's V")
})

test_that("effect-size caption is dynamic and explains mixed scales", {
  set.seed(106)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 15)),
    outcome_two = stats::rnorm(45),
    outcome_three = stats::rnorm(45),
    category = factor(rep(c("No", "Yes", "No"), 15))
  )
  df_Test$outcome_two[df_Test$group == "C"] <- NA_real_

  tbl_without_es <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome_two", "outcome_three", "category")
  )
  tbl_with_es <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome_two", "outcome_three", "category"),
    AddEffectSize = TRUE
  )

  caption_without_es <- GetComparisonCaption(tbl_without_es)
  caption_with_es <- GetComparisonCaption(tbl_with_es)

  expect_false(grepl("Effect sizes:", caption_without_es, fixed = TRUE))
  expect_true(grepl("|d|: 0.2/0.5/0.8", caption_with_es, fixed = TRUE))
  expect_true(grepl("Cohen's f: 0.1/0.25/0.4", caption_with_es, fixed = TRUE))
  expect_true(grepl("Cramer's V: 0.1/0.3/0.5", caption_with_es, fixed = TRUE))
  expect_true(grepl("d and f are not numerically equivalent", caption_with_es, fixed = TRUE))
  expect_true(grepl("not clinical-importance cutoffs", caption_with_es, fixed = TRUE))
})

test_that("insufficient adjusted data returns a missing effect size", {
  df_Test <- data.frame(
    group = factor(c("A", "A", "B", "B")),
    outcome = c(1, 2, 3, 4),
    covariate = c(NA, NA, 1, NA)
  )

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    covariates = "covariate",
    AddEffectSize = TRUE
  )
  got <- GetMainEffectRow(tbl, "outcome")

  expect_true(is.na(got$effect_size))
})
