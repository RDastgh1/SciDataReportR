MakeFeatureDifferentialData <- function() {
  set.seed(31)
  n <- 24
  features <- matrix(rnorm(n * 40), nrow = n, dimnames = list(NULL, paste0("F", 1:40)))
  df_wide <- data.frame(
    `Sample ID` = paste0("S", seq_len(n)),
    `Study group` = rep(c("Control", "Case A", "Case B"), each = 8),
    Age = round(runif(n, 30, 70)),
    Sex = rep(c("F", "M"), n / 2),
    features,
    check.names = FALSE
  )
  df_wide$F1 <- df_wide$F1 + 5 * (df_wide$`Study group` == "Case A")
  df_wide$F2 <- df_wide$F2 - 5 * (df_wide$`Study group` == "Case B")
  names(df_wide)[names(df_wide) == "F3"] <- "HLA-E"
  df_wide
}

test_that("feature differential abundance reproduces limma and orients contrasts", {
  skip_if_not_installed("limma")
  df_wide <- MakeFeatureDifferentialData()
  vars_features <- setdiff(names(df_wide), c("Sample ID", "Study group", "Age", "Sex"))
  da <- AnalyzeFeatureDifferentialAbundance(
    df_wide, vars_features, "Study group",
    covariates = c("Age", "Sex"), reference = "Control", feature_label = "protein"
  )

  expect_s3_class(da, "feature_differential_abundance")
  expect_identical(da$Metadata$Contrasts, c("Case A vs Control", "Case B vs Control"))
  expect_equal(nrow(da$ResultsTable), 80)
  expect_match(da$Metadata$Methods, "Differential protein abundance")

  group <- factor(df_wide$`Study group`, levels = c("Control", "Case A", "Case B"))
  design <- stats::model.matrix(~ group + Age + Sex, data = df_wide)
  fit <- limma::eBayes(limma::lmFit(t(as.matrix(df_wide[vars_features])), design))
  expected <- limma::topTable(fit, coef = "groupCase A", number = Inf, sort.by = "none", confint = TRUE)
  df_a <- da$ResultsTable[da$ResultsTable$Contrast == "Case A vs Control", ]
  expect_identical(df_a$Variable, vars_features)
  expect_equal(df_a$Effect, expected$logFC)
  expect_equal(df_a$TStatistic, expected$t)
  expect_equal(df_a$PValue, expected$P.Value)
  expect_equal(df_a$FDR, expected$adj.P.Val)
  expect_equal(df_a$ConfLow, expected$CI.L)
  expect_equal(df_a$ConfHigh, expected$CI.R)
  expect_true(df_a$FDRSignificant[df_a$Variable == "F1"])
  df_b <- da$ResultsTable[da$ResultsTable$Contrast == "Case B vs Control", ]
  expect_identical(df_b$Direction[df_b$Variable == "F2"], "Lower")

  custom <- AnalyzeFeatureDifferentialAbundance(
    df_wide, vars_features, "Study group",
    contrasts = list(c("Control", "Case A"), c("Case B", "Case A"))
  )
  expect_identical(custom$Metadata$Contrasts, c("Control vs Case A", "Case B vs Case A"))
  expect_lt(custom$ResultsTable$Effect[[1]], -2)
  expect_output(print(da), "Case A vs Control: ")
  expect_error(AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Study group", contrasts = list(c("Case A", "Nope"))), "Unknown group")
  expect_error(AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Study group", reference = "Nope"), "reference")
})

test_that("feature differential abundance handles numeric predictors and exclusions", {
  skip_if_not_installed("limma")
  df_wide <- MakeFeatureDifferentialData()
  vars_features <- setdiff(names(df_wide), c("Sample ID", "Study group", "Age", "Sex"))
  df_wide$F4 <- 0.1 * df_wide$Age + rnorm(nrow(df_wide), sd = 0.3)
  df_wide$Age[1:2] <- NA
  df_wide$F5[3:6] <- NA
  df_wide$F6 <- 1
  df_wide$F7[-(1:4)] <- NA

  da <- AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Age", covariates = "Sex")
  expect_identical(da$Metadata$Contrasts, "Age (per unit)")
  expect_equal(da$SampleAudit$SamplesDroppedIncompleteDesign, 2)
  expect_equal(da$Metadata$NSamples, 22)
  expect_setequal(da$ExcludedFeatures$Variable, c("F6", "F7"))
  df_results <- da$ResultsTable
  expect_equal(df_results$Effect[df_results$Variable == "F4"], 0.1, tolerance = 0.2)
  expect_identical(unique(df_results$EffectType), "Slope")
  expect_equal(df_results$N[df_results$Variable == "F5"], 18)

  df_confounded <- df_wide
  df_confounded$Copy <- df_confounded$`Study group`
  expect_error(AnalyzeFeatureDifferentialAbundance(df_confounded, vars_features, "Study group", covariates = "Copy"), "full rank")
  expect_error(AnalyzeFeatureDifferentialAbundance(df_wide, "F1", "Age"), "at least two")
  expect_error(AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Missing"), "Missing columns")
})

test_that("blocked feature designs gain power from paired samples", {
  skip_if_not_installed("limma")
  skip_if_not_installed("statmod")
  set.seed(32)
  n_participants <- 12
  participant_effect <- matrix(rnorm(n_participants * 30, sd = 2), nrow = n_participants)
  make_visit <- function(shift) participant_effect + matrix(rnorm(n_participants * 30, sd = 0.4), nrow = n_participants) + shift
  features <- rbind(make_visit(0), make_visit(matrix(rep(c(0.6, rep(0, 29)), each = n_participants), nrow = n_participants)))
  colnames(features) <- paste0("F", 1:30)
  df_wide <- data.frame(
    Participant = rep(paste0("ID", seq_len(n_participants)), 2),
    Visit = rep(c("Baseline", "Followup"), each = n_participants),
    features
  )
  vars_features <- paste0("F", 1:30)

  unpaired <- AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Visit")
  paired <- AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Visit", block_var = "Participant")
  expect_gt(paired$Metadata$BlockCorrelation, 0.8)
  expect_lt(paired$ResultsTable$PValue[[1]], unpaired$ResultsTable$PValue[[1]])
  expect_lt(paired$ResultsTable$PValue[[1]], 0.01)
  expect_match(paired$Metadata$Methods, "duplicateCorrelation")
  expect_error(
    AnalyzeFeatureDifferentialAbundance(df_wide[c(1:6, 19:24), ], vars_features, "Visit", block_var = "Participant"),
    "no repeated samples"
  )
})

test_that("PlotFeatureDifferentialAbundance draws a volcano per contrast", {
  skip_if_not_installed("limma")
  df_wide <- MakeFeatureDifferentialData()
  vars_features <- setdiff(names(df_wide), c("Sample ID", "Study group", "Age", "Sex"))
  da <- AnalyzeFeatureDifferentialAbundance(df_wide, vars_features, "Study group", reference = "Control")

  volcano <- PlotFeatureDifferentialAbundance(da)
  expect_s3_class(volcano, "ggplot")
  expect_identical(volcano$labels$title, "Case A vs Control")
  expect_identical(as.character(volcano$data$Status[volcano$data$Variable == "F1"]), "Higher")
  second <- PlotFeatureDifferentialAbundance(da, contrast = "Case B vs Control", min_effect = 1)
  expect_identical(as.character(second$data$Status[second$data$Variable == "F2"]), "Lower")
  # Drawing evaluates the aesthetics, which building the object does not.
  built <- ggplot2::ggplot_build(volcano)
  expect_equal(nrow(built$data[[2]]), 40)
  expect_silent(ggplot2::ggplot_build(second))
  expect_silent(ggplot2::ggplot_build(PlotFeatureDifferentialAbundance(da, label_top_n = 0)))
  expect_error(PlotFeatureDifferentialAbundance(da, contrast = "Nope"), "must be one of")
  expect_error(PlotFeatureDifferentialAbundance(da, label_by = "Nope"), "label column")
})
