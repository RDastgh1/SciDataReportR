GetComparisonTableMainRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}

GetComparisonTablePairwiseP <- function(tbl, variable, contrast) {
  headers <- tbl$table_styling$header
  column <- headers$column[headers$label == paste0("**", contrast, "**")]
  row <- GetComparisonTableMainRow(tbl, variable)
  as.numeric(row[[column]])
}

PairwiseLogisticLR <- function(data, outcome, group, covariates, pairs) {
  raw_p <- vapply(pairs, function(pair) {
    df_pair <- data[data[[group]] %in% pair, , drop = FALSE]
    df_pair[[group]] <- droplevels(df_pair[[group]])
    full <- stats::glm(
      stats::reformulate(c(group, covariates), response = outcome),
      data = df_pair,
      family = stats::binomial()
    )
    reduced <- stats::glm(
      stats::reformulate(covariates, response = outcome),
      data = df_pair,
      family = stats::binomial()
    )
    stats::anova(reduced, full, test = "Chisq")$`Pr(>Chi)`[2]
  }, numeric(1))
  stats::p.adjust(raw_p, method = "holm")
}

PairwiseMultinomialLR <- function(data, outcome, group, covariates, pairs) {
  skip_if_not_installed("nnet")

  raw_p <- vapply(pairs, function(pair) {
    df_pair <- data[data[[group]] %in% pair, , drop = FALSE]
    df_pair[[group]] <- droplevels(df_pair[[group]])
    full <- nnet::multinom(
      stats::reformulate(c(group, covariates), response = outcome),
      data = df_pair,
      trace = FALSE
    )
    reduced <- nnet::multinom(
      stats::reformulate(covariates, response = outcome),
      data = df_pair,
      trace = FALSE
    )
    lr <- 2 * (as.numeric(stats::logLik(full)) - as.numeric(stats::logLik(reduced)))
    df_lr <- attr(stats::logLik(full), "df") - attr(stats::logLik(reduced), "df")
    stats::pchisq(lr, df = df_lr, lower.tail = FALSE)
  }, numeric(1))
  stats::p.adjust(raw_p, method = "holm")
}

test_that("adjusted continuous pairwise contrasts match emmeans Bonferroni results", {
  skip_if_not_installed("emmeans")

  set.seed(601)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 32)),
    age = stats::rnorm(96)
  )
  df_Test$outcome <- c(A = 0, B = 0.6, C = 1.2)[df_Test$group] +
    0.5 * df_Test$age + stats::rnorm(96)

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    covariates = "age",
    AddPairwise = TRUE,
    PairwiseMethod = "bonferroni"
  )

  fit <- stats::lm(outcome ~ group + age, data = df_Test)
  expected <- as.data.frame(summary(
    emmeans::contrast(emmeans::emmeans(fit, ~group), method = "pairwise"),
    adjust = "bonferroni"
  ))

  expect_equal(
    GetComparisonTablePairwiseP(tbl, "outcome", "A - B"),
    expected$p.value[expected$contrast == "A - B"],
    tolerance = 1e-10
  )
  expect_equal(
    GetComparisonTablePairwiseP(tbl, "outcome", "B - C"),
    expected$p.value[expected$contrast == "B - C"],
    tolerance = 1e-10
  )
})

test_that("referent continuous contrasts correct only the displayed family", {
  skip_if_not_installed("emmeans")

  set.seed(602)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 26)),
    outcome = stats::rnorm(78)
  )
  df_Test$outcome <- df_Test$outcome + c(A = 0, B = 0.5, C = 1)[df_Test$group]

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    AddPairwise = TRUE,
    PairwiseMethod = "bonferroni",
    Referent = "A"
  )

  expected_raw <- c(
    stats::t.test(outcome ~ group, data = subset(df_Test, group %in% c("A", "B")))$p.value,
    stats::t.test(outcome ~ group, data = subset(df_Test, group %in% c("A", "C")))$p.value
  )
  expected <- stats::p.adjust(expected_raw, method = "bonferroni")

  expect_equal(GetComparisonTablePairwiseP(tbl, "outcome", "A - B"), expected[1], tolerance = 1e-10)
  expect_equal(GetComparisonTablePairwiseP(tbl, "outcome", "A - C"), expected[2], tolerance = 1e-10)
  expect_false(any(grepl("B - C", tbl$table_styling$header$label, fixed = TRUE)))
})

test_that("unadjusted binary pairwise tests use the requested Pearson or Fisher method", {
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 20)),
    outcome = factor(c(rep(c("No", "Yes"), c(15, 5)), rep(c("No", "Yes"), c(10, 10)), rep(c("No", "Yes"), c(4, 16))))
  )
  pairs <- list(c("A", "B"), c("A", "C"), c("B", "C"))

  for (method in c("chisq", "fisher")) {
    tbl <- MakeComparisonTable(
      data = df_Test,
      group_var = "group",
      variables = "outcome",
      AddPairwise = TRUE,
      PairwiseMethod = "none",
      CatMethod = method
    )

    expected <- vapply(pairs, function(pair) {
      sub <- df_Test[df_Test$group %in% pair, , drop = FALSE]
      sub$group <- droplevels(sub$group)
      tab <- table(sub$outcome, sub$group)
      if (method == "chisq") {
        stats::chisq.test(tab, correct = FALSE)$p.value
      } else {
        stats::fisher.test(tab)$p.value
      }
    }, numeric(1))

    expect_equal(GetComparisonTablePairwiseP(tbl, "outcome", "A - B"), expected[1], tolerance = 1e-10)
    expect_equal(GetComparisonTablePairwiseP(tbl, "outcome", "B - C"), expected[3], tolerance = 1e-10)
  }
})

test_that("adjusted binary and multicategory pairwise tests use likelihood-ratio tests", {
  skip_if_not_installed("nnet")

  set.seed(603)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 80)),
    age = stats::rnorm(240)
  )
  group_effect <- c(A = -0.7, B = 0, C = 0.7)[df_Test$group]
  df_Test$binary <- factor(ifelse(stats::rbinom(240, 1, stats::plogis(group_effect + 0.4 * df_Test$age)) == 1, "Yes", "No"))
  df_Test$multi <- factor(sample(c("Low", "Middle", "High"), 240, replace = TRUE, prob = c(0.3, 0.4, 0.3)))
  df_Test$multi[df_Test$group == "C" & df_Test$age > 0] <- "High"
  pairs <- list(c("A", "B"), c("A", "C"), c("B", "C"))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("binary", "multi"),
    covariates = "age",
    AddPairwise = TRUE,
    PairwiseMethod = "holm"
  )

  binary_expected <- PairwiseLogisticLR(df_Test, "binary", "group", "age", pairs)
  multi_expected <- PairwiseMultinomialLR(df_Test, "multi", "group", "age", pairs)

  expect_equal(GetComparisonTablePairwiseP(tbl, "binary", "A - C"), binary_expected[2], tolerance = 1e-10)
  expect_equal(GetComparisonTablePairwiseP(tbl, "multi", "B - C"), multi_expected[3], tolerance = 1e-10)
  expect_identical(GetComparisonTableMainRow(tbl, "binary")$Test, "Logistic regression (LR)")
  expect_identical(GetComparisonTableMainRow(tbl, "multi")$Test, "Multinomial LR")
})

test_that("robust covariate-adjusted continuous tests use HC3 global and pairwise inference", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("sandwich")
  skip_if_not_installed("car")

  set.seed(604)
  df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 35)),
    age = stats::rnorm(105)
  )
  df_Test$outcome <- c(A = 0, B = 0.5, C = 1.1)[df_Test$group] +
    0.7 * df_Test$age + stats::rnorm(105, sd = exp(df_Test$age / 3))

  tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    covariates = "age",
    Parametric = FALSE,
    AddPairwise = TRUE,
    PairwiseMethod = "holm"
  )

  fit <- stats::lm(outcome ~ group + age, data = df_Test)
  V <- sandwich::vcovHC(fit, type = "HC3")
  expected_global <- car::linearHypothesis(
    fit,
    hypothesis.matrix = c("groupB = 0", "groupC = 0"),
    vcov. = V,
    test = "F"
  )$`Pr(>F)`[2]
  emm <- emmeans::emmeans(fit, ~group, vcov. = V)
  expected_pairs <- as.data.frame(summary(emmeans::contrast(emm, method = "pairwise"), adjust = "holm"))

  expect_equal(GetComparisonTableMainRow(tbl, "outcome")$p.value, expected_global, tolerance = 1e-10)
  expect_equal(
    GetComparisonTablePairwiseP(tbl, "outcome", "A - C"),
    expected_pairs$p.value[expected_pairs$contrast == "A - C"],
    tolerance = 1e-10
  )
  expect_identical(GetComparisonTableMainRow(tbl, "outcome")$Test, "Robust ANCOVA (HC3 Wald)")
})

test_that("pairwise styling survives tbl_merge and the merged table renders", {
  skip_if_not_installed("gt")

  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  Labelled$ApoeGroup <- factor(
    dplyr::case_when(
      Labelled$Genotype %in% c("E2E2", "E2E3") ~ "E2 carrier",
      Labelled$Genotype == "E3E3" ~ "E3/E3",
      TRUE ~ "E4 carrier"
    ),
    levels = c("E2 carrier", "E3/E3", "E4 carrier")
  )

  merged <- MakeFacetCatComparisonTable(
    data = Labelled,
    FacetVariables = c("Diagnosis", "ApoeGroup"),
    variables = c("age", "tau", "p_tau"),
    AddPairwise = TRUE
  )

  # tbl_merge() suffixes the pairwise columns. The styling rules must follow
  # them: a rule left pointing at the pre-merge name makes the table
  # unprintable with "Column `pw_1` not found".
  text_format <- merged$table_styling$text_format
  pairwise_rules <- text_format[grepl("^pw_", text_format$column), ]
  expect_gt(nrow(pairwise_rules), 0)

  for (i in seq_len(nrow(pairwise_rules))) {
    referenced <- all.vars(rlang::quo_get_expr(pairwise_rules$rows[[i]]))
    referenced <- setdiff(referenced, ".data")
    expect_true(
      all(referenced %in% names(merged$table_body)),
      info = paste("rule for", pairwise_rules$column[i])
    )
  }

  # Rendering is the check that matters: the failure only surfaces on print.
  expect_no_error(rendered <- gt::as_raw_html(gtsummary::as_gt(merged)))
  expect_true(grepl("font-weight: bold", as.character(rendered), fixed = TRUE))
})

test_that("an unmerged pairwise table still renders and bolds significant cells", {
  skip_if_not_installed("gt")

  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  Labelled$ApoeGroup <- factor(
    dplyr::case_when(
      Labelled$Genotype %in% c("E2E2", "E2E3") ~ "E2 carrier",
      Labelled$Genotype == "E3E3" ~ "E3/E3",
      TRUE ~ "E4 carrier"
    ),
    levels = c("E2 carrier", "E3/E3", "E4 carrier")
  )

  tbl <- MakeComparisonTable(
    data = Labelled,
    group_var = "ApoeGroup",
    variables = c("age", "tau", "p_tau"),
    AddPairwise = TRUE
  )

  expect_no_error(rendered <- gt::as_raw_html(gtsummary::as_gt(tbl)))
  expect_true(grepl("font-weight: bold", as.character(rendered), fixed = TRUE))
})
