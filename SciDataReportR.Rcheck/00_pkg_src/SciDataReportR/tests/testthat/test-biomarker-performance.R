test_that("EvaluateBiomarkerPerformance handles binary outcomes and labels", {
  skip_if_not_installed("pROC")
  set.seed(901)
  df_Test <- data.frame(
    `Disease Cohort` = factor(sample(c("Control", "Disease"), 90, TRUE)),
    `NfL value` = rnorm(90),
    Age = rnorm(90),
    Sex = factor(sample(c("Female", "Male"), 90, TRUE)),
    check.names = FALSE
  )
  result <- EvaluateBiomarkerPerformance(df_Test, "Disease Cohort", "NfL value", c("Age", "Sex"), PositiveLevel = "Disease", CIBootstrapR = 10)
  expect_equal(result$Metadata$PositiveLevel, "Disease")
  expect_true(all(c("Models", "PerformanceTable", "ValidationPredictions", "Plots") %in% names(result)))
  expect_equal(result$PerformanceTable$N[1], 90)
  expect_true(nrow(result$ThresholdTable) >= 1)
})

test_that("categorical biomarkers, continuous outcomes, and pairwise screening are stable", {
  skip_if_not_installed("pROC")
  set.seed(902)
  df_Test <- data.frame(
    Outcome = factor(sample(c("No", "Yes"), 80, TRUE)),
    ContinuousOutcome = rnorm(80),
    Acrocyanosis = factor(sample(c("Absent", "Present"), 80, TRUE)),
    Pattern = factor(sample(c("None", "Mild", "Severe"), 80, TRUE)),
    NfL = rnorm(80),
    Age = rnorm(80)
  )
  df_Test$NfL[1:5] <- NA_real_
  categorical <- EvaluateBiomarkerPerformance(df_Test, "Outcome", "Acrocyanosis", "Age", CIBootstrapR = 10)
  multi <- EvaluateBiomarkerPerformance(df_Test, "Outcome", "Pattern", "Age", CIBootstrapR = 10)
  continuous <- EvaluateBiomarkerPerformance(df_Test, "ContinuousOutcome", "NfL", "Age", CIBootstrapR = 10)
  expect_true(nrow(categorical$RawBiomarkerThreshold) == 1)
  expect_null(multi$RawBiomarkerThreshold)
  expect_true(all(is.na(continuous$PerformanceTable$AUC)))
  screen <- ScreenBiomarkerPerformance(df_Test, c("Outcome", "ContinuousOutcome"), c("NfL", "Acrocyanosis"), "Age", CIBootstrapR = 10)
  expect_true(all(c("FailureTable", "HeatmapValue") %in% c(names(screen), names(screen$Plots$Heatmap$data))))
  expect_s3_class(add_biomarker_values(screen$Plots$Heatmap), "ggplot")
})

test_that("cross-validation returns one audit row per analysis row and model", {
  skip_if_not_installed("pROC")
  set.seed(903)
  df_Test <- data.frame(Outcome = factor(sample(c("No", "Yes"), 60, TRUE)), Biomarker = rnorm(60), Age = rnorm(60))
  result <- EvaluateBiomarkerPerformance(df_Test, "Outcome", "Biomarker", "Age", Validation = "cross_validation", CVFolds = 5, CIBootstrapR = 10)
  expect_true(all(result$ValidationPredictions$Fold %in% 1:5))
  expect_equal(nrow(result$ValidationPredictions), nrow(df_Test) * 3)
})

test_that("screening returns labelled continuous-biomarker panels and ROC facets", {
  skip_if_not_installed("pROC")
  set.seed(904)
  df_Test <- data.frame(
    Outcome = factor(sample(c("Control", "Impaired"), 100, TRUE)),
    Strong = rnorm(100),
    Weak = rnorm(100),
    Genotype = factor(sample(c("E3E3", "E3E4", "E4E4"), 100, TRUE)),
    Age = rnorm(100)
  )
  codebook <- data.frame(
    Variable = c("Outcome", "Strong", "Weak", "Genotype", "Age"),
    Label = c("Diagnosis", "Strong protein", "Weak protein", "Genotype", "Age")
  )

  screen <- ScreenBiomarkerPerformance(
    df_Test,
    outcome_vars = "Outcome",
    biomarker_vars = c("Strong", "Weak", "Genotype"),
    covariates = "Age",
    PositiveLevel = "Impaired",
    codebook = codebook,
    CIBootstrapR = 10
  )

  expect_s3_class(screen$Plots$BiomarkerPanels, "ggplot")
  expect_s3_class(screen$Plots$ROCFacets, "ggplot")
  expect_setequal(unique(screen$Plots$BiomarkerPanels$data$BiomarkerLabel), c("Strong protein", "Weak protein"))
  expect_true(all(c("AUC", "AUC_Lower", "AUC_Upper", "N") %in% names(screen$Plots$BiomarkerPanels$data)))
})
