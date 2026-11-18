test_that("split violin captures bare columns before inference inspection", {
  skip_if_not_installed("emmeans")
  set.seed(328)
  df_Test <- data.frame(outcome = rnorm(60), group = factor(rep(c("A", "B"), 30)), age = rnorm(60))
  for (nonparametric in c(FALSE, TRUE)) {
    expect_no_warning(result <- PlotSplitViolin(df_Test, outcome, group,
      covariates = "age", nonparametric = nonparametric))
    expect_identical(result$Results$Outcome, "outcome")
    expect_identical(result$Results$Predictor, "group")
    expect_true(all(result$Results$Status == "ok"))
    expect_no_warning(quoted <- PlotSplitViolin(df_Test, "outcome", "group",
      covariates = "age", nonparametric = nonparametric))
    expect_equal(attr(result, "comparison")$PValue, attr(quoted, "comparison")$PValue)
  }
})

test_that("data-frame pathway results do not probe absent list components", {
  set.seed(329)
  df_Test <- data.frame(marker = rnorm(60, 10), group = factor(rep(c("A", "B"), 30)))
  expect_no_warning(result <- calculate_pathway_results(df_Test, "group", metabolites = "marker"))
  expect_true(all(attr(result, "Results")$Status == "ok"))
})

test_that("ANOVA significance sizes use an explicit ordinal scale", {
  set.seed(330)
  df_Test <- data.frame(marker = rnorm(60), group = factor(rep(c("A", "B"), 30)))
  result <- PlotAnovaRelationshipsMatrix(df_Test, "group", "marker")
  expect_no_warning(ggplot2::ggplot_build(result$Unadjusted$plot))
  expect_no_warning(ggplot2::ggplot_build(result$FDRCorrected$plot))
})

LatentWarningFixture <- function() {
  set.seed(810)
  as.data.frame(replicate(4, factor(sample(c("No", "Yes"), 100, TRUE))))
}

test_that("latent class estimates survive the upstream standard-error matrix defect", {
  skip_if_not_installed("poLCA")
  df_Test <- LatentWarningFixture()
  df_Training <- as.data.frame(lapply(df_Test, function(x) as.integer(factor(x))))
  PoLCA <- poLCA::poLCA
  set.seed(439)
  reference <- PoLCA(cbind(V1, V2, V3, V4) ~ 1, df_Training,
    nclass = 2, nrep = 1, verbose = FALSE, calc.se = FALSE)
  for (condition_type in c("warning", "error")) local({
    # The same upstream defect is a warning normally and an error during checks.
    local_mocked_bindings(poLCA = function(..., calc.se = TRUE) {
      if (calc.se) {
        call <- quote(matrix(sqrt(diag(VCE.beta)), nrow = 1, ncol = 1))
        message <- "data length differs from size of matrix: [9 != 1 x 1]"
        if (condition_type == "error") stop(simpleError(message, call = call))
        warning(simpleWarning(message, call = call))
      }
      PoLCA(..., calc.se = calc.se)
    }, .package = "poLCA")
    expect_no_warning(result <- CreateClusterModel_LatentClass(df_Test,
      method = "finalize", final_k = 2, nrep = 1, seed = 437))
    expect_equal(result$ModelInfo$lca_model$llik, reference$llik)
    expect_equal(result$ModelInfo$lca_model$posterior, reference$posterior)
    expect_false(result$ModelInfo$StandardErrors$Available)
    expect_match(result$ModelInfo$StandardErrors$Reason, "data length differs", fixed = TRUE)
    expect_true(all(is.na(result$ModelInfo$lca_model$probs.se)))
  })
})

test_that("latent class failures retain unrelated engine diagnostics", {
  skip_if_not_installed("poLCA")
  local_mocked_bindings(poLCA = function(...) stop("planted input failure"), .package = "poLCA")
  expect_error(CreateClusterModel_LatentClass(LatentWarningFixture(),
    method = "finalize", final_k = 2, nrep = 1), "planted input failure", fixed = TRUE)
})
