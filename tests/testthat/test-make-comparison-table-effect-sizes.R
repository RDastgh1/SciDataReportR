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

test_that("two-group Wilcoxon reports absolute rank-biserial correlation", {
  df_Test <- data.frame(
    group = factor(
      c(rep("A", 4), rep("B", 4), "A", "B", NA),
      levels = c("A", "B")
    ),
    outcome = c(1, 2, 2, 4, 3, 4, 5, 5, NA, 6, 7)
  )

  tbl <- suppressWarnings(MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    Parametric = FALSE,
    AddEffectSize = TRUE
  ))
  got <- GetMainEffectRow(tbl, "outcome")

  df_CC <- df_Test[stats::complete.cases(df_Test), , drop = FALSE]
  ranks <- rank(df_CC$outcome, ties.method = "average")
  n_first <- sum(df_CC$group == "A")
  n_second <- sum(df_CC$group == "B")
  u_first <- sum(ranks[df_CC$group == "A"]) - n_first * (n_first + 1) / 2
  expected <- abs(2 * u_first / (n_first * n_second) - 1)

  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "|rank-biserial r|")

  df_Reordered <- df_Test
  df_Reordered$group <- factor(df_Reordered$group, levels = c("B", "A"))
  got_reordered <- GetMainEffectRow(
    suppressWarnings(MakeComparisonTable(
      data = df_Reordered,
      group_var = "group",
      variables = "outcome",
      Parametric = FALSE,
      AddEffectSize = TRUE
    )),
    "outcome"
  )

  expect_equal(got_reordered$effect_size, expected, tolerance = 1e-10)
})

test_that("multi-group nonparametric reports rank f and retains categorical effect sizes", {
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
  df_Test$ranks <- rank(df_Test$outcome, ties.method = "average")
  expected_f <- effectsize::cohens_f(stats::aov(ranks ~ group, data = df_Test))$Cohens_f[1]

  expect_equal(got_outcome$effect_size, expected_f, tolerance = 1e-10)
  expect_identical(got_outcome$ES_Method, "Cohen's f on ranks")
  expect_identical(got_category$ES_Method, "Cramer's V")
  tab <- table(df_Test$category, df_Test$group)
  expected_v <- if (requireNamespace("DescTools", quietly = TRUE)) {
    as.numeric(DescTools::CramerV(tab, method = "bias.corrected"))
  } else {
    sqrt(as.numeric(stats::chisq.test(tab, correct = FALSE)$statistic) /
      (sum(tab) * (min(dim(tab)) - 1)))
  }
  expect_equal(got_category$effect_size, expected_v)
})

test_that("identical group rank distributions give zero rank f", {
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 3)),
    outcome = rep(c(1, 2, 3), 3)
  )

  got <- GetMainEffectRow(
    MakeComparisonTable(
      data = df_Test,
      group_var = "group",
      variables = "outcome",
      Parametric = FALSE,
      AddEffectSize = TRUE
    ),
    "outcome"
  )

  expect_equal(got$effect_size, 0)
  expect_identical(got$ES_Method, "Cohen's f on ranks")
})

test_that("captions stay focused when effect sizes are included", {
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

  expect_identical(
    GetComparisonCaption(tbl_without_es),
    GetComparisonCaption(tbl_with_es)
  )
  expect_false(grepl("Effect sizes:", GetComparisonCaption(tbl_with_es), fixed = TRUE))
  expect_false(grepl("Cohen's f:", GetComparisonCaption(tbl_with_es), fixed = TRUE))
  expect_false(grepl("Cramer's V:", GetComparisonCaption(tbl_with_es), fixed = TRUE))
})

test_that("caption and methods footnote describe the selected analysis", {
  set.seed(107)
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 20)),
    age = stats::rnorm(40),
    outcome = stats::rnorm(40),
    category = factor(rep(c("No", "Yes"), 20))
  )
  labelled::var_label(df_Test$group) <- "Diagnostic group"
  labelled::var_label(df_Test$age) <- "Age (years)"

  tbl_unadjusted <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome", "category"),
    CatMethod = "auto"
  )
  tbl_adjusted <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome", "category"),
    covariates = "age",
    AddPairwise = TRUE,
    PairwiseMethod = "holm",
    ShowNotes = "always"
  )

  caption_unadjusted <- GetComparisonCaption(tbl_unadjusted)
  caption_adjusted <- GetComparisonCaption(tbl_adjusted)
  note_unadjusted <- tbl_unadjusted$table_styling$footnote_header %>%
    dplyr::filter(.data$column == "p.value_fmt") %>%
    dplyr::pull("footnote")
  note_adjusted <- tbl_adjusted$table_styling$footnote_header %>%
    dplyr::filter(.data$column == "p.value_fmt") %>%
    dplyr::pull("footnote")

  expect_match(caption_unadjusted, "Comparison by Diagnostic group", fixed = TRUE)
  expect_false(grepl("adjusted for", caption_unadjusted, fixed = TRUE))
  expect_false(grepl("Pairwise", caption_unadjusted, fixed = TRUE))
  expect_false(grepl("auto", caption_unadjusted, fixed = TRUE))
  expect_match(note_unadjusted, "chi-square or Fisher's exact", fixed = TRUE)

  expect_match(caption_adjusted, "p-values adjusted for Age (years)", fixed = TRUE)
  expect_match(caption_adjusted, "holm-adjusted within outcome", fixed = TRUE)
  expect_false(grepl("multinomial_LR", caption_adjusted, fixed = TRUE))
  expect_match(note_adjusted, "Type II ANCOVA", fixed = TRUE)
  expect_match(note_adjusted, "multinomial likelihood-ratio", fixed = TRUE)
  expect_match(note_adjusted, "holm-adjusted within outcome", fixed = TRUE)

  html <- gt::as_raw_html(gtsummary::as_gt(tbl_adjusted))
  expect_match(html, "Analysis notes", fixed = TRUE)
  expect_match(html, "multinomial likelihood-ratio", fixed = TRUE)
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
  expect_identical(got$ES_Method, "Adjusted effect unavailable")
})


test_that("rank f uses complete cases, ties, and observed group levels", {
  df_Test <- data.frame(
    group = factor(c(rep("A", 4), rep("B", 5), rep("C", 6), "D", NA),
      levels = c("A", "B", "C", "D", "Unused")),
    outcome = c(1, 1, 2, 3, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, NA, 8)
  )
  df_CC <- droplevels(df_Test[stats::complete.cases(df_Test), ])
  df_CC$ranks <- rank(df_CC$outcome, ties.method = "average")
  expected <- effectsize::cohens_f(stats::aov(ranks ~ group, data = df_CC))$Cohens_f[1]
  tbl <- suppressWarnings(MakeComparisonTable(
    df_Test, "group", "outcome", Parametric = FALSE, AddEffectSize = TRUE,
    AddPairwise = TRUE
  ))
  tbl_without_es <- suppressWarnings(MakeComparisonTable(
    df_Test, "group", "outcome", Parametric = FALSE, AddPairwise = TRUE
  ))
  got <- GetMainEffectRow(tbl, "outcome")
  expect_equal(got$effect_size, expected, tolerance = 1e-10)
  expect_identical(got$ES_Method, "Cohen's f on ranks")
  expect_equal(got$p.value, stats::kruskal.test(outcome ~ group, data = df_CC)$p.value)
  vars_Unchanged <- intersect(names(tbl_without_es$table_body), names(tbl$table_body))
  expect_equal(tbl$table_body[, vars_Unchanged], tbl_without_es$table_body[, vars_Unchanged])
  expect_match(gt::as_raw_html(gtsummary::as_gt(tbl)), "Cohen&#39;s f on ranks|Cohen's f on ranks")
})

test_that("rank f preserves infinite perfect separation", {
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 4)),
    outcome = rep(c(1, 2, 3), each = 4)
  )
  got <- GetMainEffectRow(suppressWarnings(MakeComparisonTable(
    df_Test, "group", "outcome", Parametric = FALSE, AddEffectSize = TRUE
  )), "outcome")
  expect_identical(got$effect_size, Inf)
  expect_identical(got$ES_Method, "Cohen's f on ranks")
})

test_that("undefined rank effects are missing for constant or insufficient data", {
  df_Test <- data.frame(
    group = factor(c(rep(c("A", "B", "C"), each = 4), NA, NA)),
    constant = c(rep(1, 12), 2, 3),
    single_group = c(1:4, rep(NA_real_, 10))
  )
  tbl <- suppressWarnings(MakeComparisonTable(
    df_Test, "group", c("constant", "single_group"),
    Parametric = FALSE, AddEffectSize = TRUE
  ))
  expect_true(is.na(GetMainEffectRow(tbl, "constant")$effect_size))
  expect_true(is.na(GetMainEffectRow(tbl, "single_group")$effect_size))
  expect_identical(GetMainEffectRow(tbl, "single_group")$ES_Method, "Insufficient groups")
})
