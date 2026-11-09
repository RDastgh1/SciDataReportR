MiningTestData <- function() {
  set.seed(1)
  n <- 60
  data.frame(
    sex = factor(sample(c("F", "M"), n, TRUE)),
    race = factor(sample(c("A", "B", "C"), n, TRUE)),
    emp = sample(c(TRUE, FALSE), n, TRUE),
    age = rnorm(n, 60, 8), edu = rnorm(n, 13, 2),
    y1 = rnorm(n), y2 = rnorm(n)
  )
}

test_that("all mixed-type requested pairs are tested under both settings", {
  df_Test <- MiningTestData()
  attr(df_Test$emp, "label") <- "Employed"
  for (parametric in c(TRUE, FALSE)) {
    res <- PlotMiningMatrix(df_Test, c("sex", "y1", "y2"),
                            c("race", "emp", "age", "edu"), Parametric = parametric)
    df_Table <- res$Unadjusted$PvalTable
    expect_equal(nrow(df_Table), 12L)
    expect_equal(nrow(unique(df_Table[c("XVar", "YVar")])), 12L)
    expect_true(all(is.finite(df_Table$p)))
    expect_true(all(is.finite(df_Table$EffectSize)))
    expect_equal(nrow(res$Unadjusted$Missing), 0L)
    expect_identical(unname(res$Metadata$DisplayLabels["emp"]), "Employed")
    expect_identical(res$p, res$Unadjusted)
    for (vars in list(c("age", "sex"), c("sex", "age"))) {
      single <- PlotMiningMatrix(df_Test, vars[1], vars[2], Parametric = parametric)
      expect_equal(nrow(single$Unadjusted$PvalTable), 1L)
      expect_true(is.finite(single$Unadjusted$PvalTable$p))
    }
  }
})

test_that("Kruskal eta squared uses pairwise complete n and observed levels", {
  df_Test <- MiningTestData()
  df_Test$age[1:5] <- NA_real_
  res <- PlotAnovaRelationshipsMatrix(df_Test, "race", "age", Parametric = FALSE)
  df_Complete <- df_Test[!is.na(df_Test$age), ]
  h <- unname(stats::kruskal.test(age ~ race, df_Complete)$statistic)
  expect_equal(res$Unadjusted$PvalTable$eta2_H, (h - 3 + 1) / (55 - 3))
  expect_true("ges" %in% names(PlotAnovaRelationshipsMatrix(
    df_Test, "race", "age")$Unadjusted$PvalTable))
})

test_that("missing pairs are reported even when every test is skipped", {
  df_Test <- MiningTestData()
  df_Test$age <- 1
  expect_warning(res <- PlotMiningMatrix(df_Test, "age", "sex"),
                 "Missing requested pairs: age vs sex")
  expect_equal(res$Unadjusted$Missing, data.frame(XVar = "age", YVar = "sex"))
  expect_equal(nrow(res$Unadjusted$PvalTable), 0L)
})

test_that("branch errors name variables and preserve other branches", {
  testthat::local_mocked_bindings(
    PlotAnovaRelationshipsMatrix = function(...) stop("forced failure"),
    .package = "SciDataReportR"
  )
  warnings <- character()
  res <- withCallingHandlers(
    PlotMiningMatrix(MiningTestData(), "y1", c("sex", "age")),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("ANOVA branch failed for sex, y1: forced failure", warnings)))
  expect_true(any(grepl("Missing requested pairs: y1 vs sex", warnings)))
  expect_equal(res$Unadjusted$PvalTable$YVar, "age")
})

test_that("missing effect sizes retain p-values and neutral significance markers", {
  testthat::local_mocked_bindings(
    PlotAnovaRelationshipsMatrix = function(...) list(Unadjusted = list(
      PvalTable = data.frame(ContinuousVariable = "age", CategoricalVariable = "sex",
                            p = 0.01, ges = NA_real_))),
    .package = "SciDataReportR"
  )
  res <- PlotMiningMatrix(MiningTestData(), "age", "sex")
  expect_equal(res$Unadjusted$PvalTable$p, 0.01)
  expect_true(is.na(res$Unadjusted$PvalTable$EffectSize))
  expect_equal(nrow(res$Unadjusted$Missing), 0L)
  df_Layer <- ggplot2::ggplot_build(res$Unadjusted$plot)$data[[1]]
  expect_equal(df_Layer$colour, "grey70")
  expect_equal(df_Layer$shape, 15)
})

test_that("per-predictor FDR groups the symmetrized table by YVar", {
  df_Test <- MiningTestData()
  vars_All <- c("age", "edu", "sex", "race")
  res <- PlotMiningMatrix(df_Test, vars_All, fdr_scope = "per_predictor")
  df_Table <- res$Unadjusted$PvalTable
  expected <- ave(df_Table$p, df_Table$YVar, FUN = function(p) stats::p.adjust(p, "fdr"))
  expect_equal(df_Table$p_adj, expected)
})

test_that("SampleData example retains historical statistics and structure", {
  df_Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
  res <- PlotMiningMatrix(df_Labelled,
    outcome_vars = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
    predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                       "Apolipoprotein_B", "C_Reactive_Protein", "Cortisol", "Insulin"))
  df_Baseline <- dget(testthat::test_path("fixtures", "mining-sample-baseline.R"))
  df_Existing <- dplyr::inner_join(
    df_Baseline[c("XVar", "YVar")], res$Unadjusted$PvalTable,
    by = c("XVar", "YVar")
  )
  # Newly covered reverse mixed pairs expand the example; old cells stay identical.
  # FDR uses the expanded family, so its adjusted values intentionally change.
  vars_Compare <- setdiff(names(df_Baseline), "p_adj")
  expect_equal(as.data.frame(df_Existing[vars_Compare]), df_Baseline[vars_Compare])
  expect_equal(nrow(res$Unadjusted$PvalTable), 30L)
  expect_s3_class(res$Unadjusted$plot, "ggplot")
})

test_that("two-valued numeric responses remain testable and covariates only exclude", {
  df_Test <- MiningTestData()
  df_Test$y1 <- rep(c(0, 1), 30)
  for (parametric in c(TRUE, FALSE)) {
    res <- PlotMiningMatrix(df_Test, "y1", c("sex", "race", "age"),
                            covariates = "age", Parametric = parametric)
    expect_equal(nrow(res$Unadjusted$PvalTable), 2L)
    expect_true(all(is.finite(res$Unadjusted$PvalTable$p)))
    expect_equal(nrow(res$Unadjusted$Missing), 0L)
  }
})
