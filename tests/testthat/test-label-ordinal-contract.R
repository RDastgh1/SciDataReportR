test_that("RevalueData retains ordinal score metadata from numeric codebook codes", {
  df_Test <- data.frame(severity = c(1, 2, 3, NA_real_))
  df_Codebook <- data.frame(
    Variable = "severity",
    Recode = "Yes",
    Code = "1=Mild;2=Moderate;3=Severe",
    Type = "Ordinal",
    Label = "Clinical severity"
  )

  out <- RevalueData(df_Test, df_Codebook)$RevaluedData

  expect_true(is.ordered(out$severity))
  expect_identical(attr(out$severity, "scidr_type"), "ordinal")
  expect_equal(unname(attr(out$severity, "scidr_ordinal_scores")), c(1, 2, 3))
  expect_equal(as.numeric(ConvertOrdinalToNumeric(out)$severity), c(1, 2, 3, NA_real_))
})

test_that("ordinal rank scoring works when no numeric code mapping is available", {
  df_Test <- data.frame(severity = ordered(
    c("Mild", "Severe", NA), levels = c("Mild", "Moderate", "Severe")
  ))
  sjlabelled::set_label(df_Test$severity) <- "Clinical severity"

  out <- ConvertOrdinalToNumeric(df_Test)

  expect_equal(as.numeric(out$severity), c(1, 3, NA_real_))
  expect_identical(attr(out$severity, "scidr_ordinal_score_source"), "rank")
  expect_identical(sjlabelled::get_label(out$severity), "Clinical severity")
})

test_that("MakeTable1 includes ordinal variables twice when requested", {
  skip_if_not_installed("gtsummary")
  df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 4)),
    severity = ordered(
      c("Mild", "Moderate", "Severe", "Moderate", "Mild", "Severe", "Moderate", "Mild"),
      levels = c("Mild", "Moderate", "Severe")
    )
  )
  sjlabelled::set_label(df_Test$severity) <- "Clinical severity"

  tbl <- MakeTable1(df_Test, variables = "severity", TreatOrdinalAs = "Both")
  labels <- tbl$table_body$label[tbl$table_body$row_type == "label"]

  expect_true(any(labels == "Clinical severity (categorical)"))
  expect_true(any(labels == "Clinical severity (continuous)"))
})

test_that("Relabel controls table display labels", {
  skip_if_not_installed("gtsummary")
  df_Test <- data.frame(score = c(1, 2, 3))
  sjlabelled::set_label(df_Test$score) <- "Clinical score"

  labelled <- MakeTable1(df_Test, variables = "score", Relabel = TRUE)
  raw <- MakeTable1(df_Test, variables = "score", Relabel = FALSE)

  expect_identical(labelled$table_body$label[1], "Clinical score")
  expect_identical(raw$table_body$label[1], "score")
})

test_that("PlotMiningMatrix labels both ordinal representations", {
  df_Test <- data.frame(
    group = factor(rep(c("A", "B"), each = 8)),
    severity = rep(0:3, 4),
    marker = stats::rnorm(16)
  )
  df_Codebook <- data.frame(
    Variable = "severity",
    Recode = "Yes",
    Code = "0=None;1=Mild;2=Moderate;3=Severe",
    Type = "Ordinal",
    Label = "Clinical severity"
  )
  df_Revalued <- RevalueData(df_Test, df_Codebook)$RevaluedData

  out <- PlotMiningMatrix(
    df_Revalued,
    outcome_vars = c("group", "severity", "marker"),
    TreatOrdinalAs = "Both"
  )
  labels <- unname(out$Metadata$DisplayLabels)
  result_labels <- unique(c(
    as.character(out$Unadjusted$PvalTable$XLabel),
    as.character(out$Unadjusted$PvalTable$YLabel)
  ))

  expect_true("Clinical severity (categorical)" %in% labels)
  expect_true("Clinical severity (continuous)" %in% labels)
  expect_true("Clinical severity (categorical)" %in% result_labels)
  expect_true("Clinical severity (continuous)" %in% result_labels)
})

test_that("matrix helpers accept ordered and unordered categorical variables together", {
  df_Test <- data.frame(
    group = factor(rep(c("Control", "Disease"), each = 8)),
    severity = ordered(
      rep(c("None", "Mild", "Moderate", "Severe"), 4),
      levels = c("None", "Mild", "Moderate", "Severe")
    ),
    marker = seq_len(16)
  )

  anova_out <- PlotAnovaRelationshipsMatrix(
    df_Test,
    CatVars = c("group", "severity"),
    ContVars = "marker",
    Relabel = FALSE
  )
  chisq_out <- PlotChiSqCovar(
    df_Test,
    predictor_vars = c("group", "severity"),
    outcome_vars = c("group", "severity"),
    Relabel = FALSE
  )

  expect_gt(nrow(anova_out$Unadjusted$PvalTable), 0)
  expect_gt(nrow(chisq_out$p$data), 0)
})
