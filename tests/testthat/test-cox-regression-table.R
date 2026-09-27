test_that("MakeCoxRegressionTable matches survival::coxph", {
  skip_if_not_installed("survival")
  set.seed(1201)
  n <- 180
  x <- stats::rnorm(n)
  age <- stats::rnorm(n, 55, 8)
  event_time <- stats::rexp(n, rate = exp(0.55 * x + 0.015 * age) / 20)
  censor_time <- stats::rexp(n, rate = 1 / 25)
  df <- data.frame(
    time = pmin(event_time, censor_time),
    event = as.integer(event_time <= censor_time),
    x = x,
    age = age
  )
  attr(df$event, "label") <- "Clinical event"
  attr(df$x, "label") <- "Marker X"

  expected <- survival::coxph(
    survival::Surv(time, event) ~ x + age,
    data = df
  )
  result <- MakeCoxRegressionTable(
    data = df,
    time_var = "time",
    event_vars = "event",
    predictor_vars = "x",
    covariates = "age",
    FDR = FALSE,
    ReturnModels = TRUE
  )
  row <- result$Results[1, ]
  expected_summary <- summary(expected)

  expect_named(result, c("FormattedTable", "LargeTable", "Results", "ModelSummaries", "Metadata"))
  expect_s3_class(result$FormattedTable, "gt_tbl")
  expect_s3_class(result$LargeTable, "gt_tbl")
  expect_s3_class(result$ModelSummaries$event$x, "coxph")
  expect_equal(row$Estimate, exp(stats::coef(expected)[["x"]]), tolerance = 1e-10)
  expect_equal(row$StdError, expected_summary$coefficients["x", "se(coef)"], tolerance = 1e-10)
  expect_equal(row$PValue, expected_summary$coefficients["x", "Pr(>|z|)"], tolerance = 1e-10)
  expect_equal(row$ConfLow, exp(stats::coef(expected)[["x"]] - stats::qnorm(0.975) * row$StdError), tolerance = 1e-10)
  expect_equal(row$ConfHigh, exp(stats::coef(expected)[["x"]] + stats::qnorm(0.975) * row$StdError), tolerance = 1e-10)
  expect_equal(row$N, expected$n)
  expect_equal(row$Events, sum(df$event))
  expect_equal(row$Concordance, unname(expected_summary$concordance[["C"]]), tolerance = 1e-10)
  expect_equal(row$PH_PValue, unname(survival::cox.zph(expected)$table["GLOBAL", "p"]), tolerance = 1e-10)
  expect_equal(row$OutcomeFamily, "survival")
  expect_equal(row$EffectType, "HR")
  expect_equal(row$ReferenceValue, 1)
  expect_equal(row$OutcomeLabel, "Clinical event")
  expect_equal(row$TermLabel, "Marker X")
  expect_false("FDR" %in% names(result$Results))
})

test_that("MakeCoxRegressionTable handles categorical references and labels", {
  skip_if_not_installed("survival")
  set.seed(1202)
  n <- 150
  group <- factor(rep(c("Control", "Dose A", "Dose B"), length.out = n))
  df <- data.frame(
    time = stats::rexp(n, exp(c(Control = 0, `Dose A` = 0.3, `Dose B` = -0.2)[group])),
    event = stats::rbinom(n, 1, 0.75),
    "Treatment Arm" = group,
    check.names = FALSE
  )
  attr(df[["Treatment Arm"]], "label") <- "Treatment group"

  default_result <- MakeCoxRegressionTable(
    df, "time", "event", "Treatment Arm",
    FDR = FALSE
  )
  list_result <- MakeCoxRegressionTable(
    df, "time", "event", "Treatment Arm",
    reference_levels = list("Treatment Arm" = "Dose B"),
    FDR = FALSE
  )
  vector_result <- MakeCoxRegressionTable(
    df, "time", "event", "Treatment Arm",
    reference_levels = c("Treatment Arm" = "Dose B"),
    FDR = FALSE
  )

  expect_setequal(default_result$Results$Level, c("Dose A", "Dose B"))
  expect_equal(list_result$Results, vector_result$Results)
  expect_setequal(list_result$Results$Level, c("Control", "Dose A"))
  expect_setequal(
    list_result$Results$TermLabel,
    c("Treatment group : Control", "Treatment group : Dose A")
  )
  expect_error(
    MakeCoxRegressionTable(df, "time", "event", "Treatment Arm",
                           reference_levels = list("Treatment Arm" = "Absent")),
    "was not found"
  )
})

test_that("MakeCoxRegressionTable keeps predictor terms separate from similarly named covariates", {
  skip_if_not_installed("survival")
  set.seed(1207)
  n <- 100
  df <- data.frame(
    time = stats::rexp(n),
    event = stats::rbinom(n, 1, 0.7),
    age = stats::rnorm(n),
    age2 = stats::rnorm(n)
  )

  result <- MakeCoxRegressionTable(
    df, "time", "event", "age",
    covariates = "age2",
    FDR = FALSE
  )

  expect_equal(nrow(result$Results), 1)
  expect_equal(result$Results$Term, "age")
})

test_that("MakeCoxRegressionTable survey path matches survey::svycoxph", {
  skip_if_not_installed("survival")
  skip_if_not_installed("survey")
  set.seed(1203)
  n <- 120
  df <- data.frame(
    time = stats::rexp(n),
    event = stats::rbinom(n, 1, 0.7),
    x = stats::rnorm(n),
    weight = stats::runif(n, 0.5, 2)
  )
  design <- survey::svydesign(ids = ~1, weights = ~weight, data = df)
  expected <- survey::svycoxph(
    survival::Surv(time, event) ~ x,
    design = design
  )
  result <- MakeCoxRegressionTable(
    df, "time", "event", "x",
    design = design,
    CheckPH = TRUE,
    ReturnModels = TRUE,
    FDR = FALSE
  )

  expect_equal(result$Results$Estimate, exp(stats::coef(expected)[["x"]]), tolerance = 1e-10)
  expect_equal(result$Results$StdError, summary(expected)$coefficients["x", "se(coef)"], tolerance = 1e-10)
  expect_equal(result$Results$PValue, summary(expected)$coefficients["x", "Pr(>|z|)"], tolerance = 1e-10)
  expect_equal(result$Results$Concordance, unname(summary(expected)$concordance[["C"]]), tolerance = 1e-10)
  expect_true(is.na(result$Results$PH_PValue))
  expect_equal(result$Metadata$Outcomes$Method, "svycoxph")
  expect_s3_class(result$ModelSummaries$event$x, "svycoxph")
})

test_that("MakeCoxRegressionTable supports event coding variants and missingness", {
  skip_if_not_installed("survival")
  set.seed(1204)
  n <- 100
  numeric_event <- rep(c(0, 1), length.out = n)
  df <- data.frame(
    time = stats::rexp(n),
    numeric_event = numeric_event,
    logical_event = as.logical(numeric_event),
    factor_event = factor(ifelse(numeric_event == 1, "Event", "Censored"),
                          levels = c("Censored", "Event")),
    x = stats::rnorm(n),
    age = stats::rnorm(n)
  )
  df$x[c(2, 9)] <- NA
  df$age[7] <- NA
  result <- MakeCoxRegressionTable(
    df,
    "time",
    c("numeric_event", "logical_event", "factor_event"),
    "x",
    covariates = "age",
    CheckPH = FALSE,
    FDR = FALSE
  )

  expect_equal(result$Results$N, rep(n - 3, 3))
  expect_true(all(is.na(result$Results$PH_PValue)))
  expect_equal(result$Metadata$Outcomes$EventLevel, c("1", "TRUE", "Event"))
  expect_error(
    MakeCoxRegressionTable(transform(df, bad = rep(c(0, 2), 50)), "time", "bad", "x"),
    "coded only 0/1"
  )
  expect_error(
    MakeCoxRegressionTable(transform(df, bad = factor(rep(c("A", "B", "C", "A"), 25))),
                           "time", "bad", "x"),
    "two-level factor"
  )

  ordered_event <- df
  ordered_event$factor_event <- ordered(
    ordered_event$factor_event,
    levels = c("Censored", "Event")
  )
  ordered_result <- MakeCoxRegressionTable(
    ordered_event,
    "time",
    "factor_event",
    "x",
    TreatOrdinalAs = "Continuous",
    FDR = FALSE
  )
  expect_equal(ordered_result$Metadata$Outcomes$EventLevel, "Event")
})

test_that("MakeCoxRegressionTable applies FDR and standardization", {
  skip_if_not_installed("survival")
  set.seed(1205)
  n <- 140
  df <- data.frame(
    time = stats::rexp(n),
    event = stats::rbinom(n, 1, 0.7),
    x1 = stats::rnorm(n, 10, 3),
    x2 = stats::rnorm(n),
    age = stats::rnorm(n, 50, 10)
  )
  result <- MakeCoxRegressionTable(
    df, "time", "event", c("x1", "x2"),
    covariates = "age",
    Standardize = TRUE,
    FDR = TRUE,
    FDRAlpha = 0.2,
    ReturnModels = TRUE
  )

  expect_equal(result$Results$FDR, ApplyFDRCorrection(result$Results$PValue))
  expect_equal(result$Results$Significant, result$Results$FDR < 0.2)
  model_frame <- result$ModelSummaries$event$x1$model
  expect_equal(as.numeric(scale(model_frame$x1)), model_frame$x1, tolerance = 1e-10)
  expect_equal(as.numeric(scale(model_frame$age)), model_frame$age, tolerance = 1e-10)
})

test_that("Cox results plot alone and bind with logistic results", {
  skip_if_not_installed("survival")
  set.seed(1206)
  n <- 130
  df <- data.frame(
    time = stats::rexp(n),
    event = stats::rbinom(n, 1, 0.7),
    diagnosis = factor(sample(c("Control", "Case"), n, replace = TRUE),
                       levels = c("Control", "Case")),
    x = stats::rnorm(n)
  )
  cox_result <- MakeCoxRegressionTable(df, "time", "event", "x")
  logistic_result <- MakeUnivariateRegressionTable(
    df,
    outcome_vars = "diagnosis",
    predictor_vars = "x"
  )
  combined <- dplyr::bind_rows(logistic_result$Results, cox_result$Results)

  expect_s3_class(PlotForestFromTable(cox_result), "ggplot")
  expect_s3_class(PlotForestFromTable(cox_result$Results), "ggplot")
  expect_s3_class(PlotForestFromTable(combined), "ggplot")
  expect_equal(unique(combined$ReferenceValue), 1)
})
