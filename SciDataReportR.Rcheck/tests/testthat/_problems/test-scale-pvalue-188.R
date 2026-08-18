# Extracted from test-scale-pvalue.R:188

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(42)
df_Anova <- data.frame(
    Group = factor(rep(c("Control", "Case"), each = 20)),
    Outcome = c(stats::rnorm(20), stats::rnorm(20, mean = 0.8))
  )
effect_result <- PlotAnovaRelationshipsMatrix(
    df_Anova,
    CatVars = "Group",
    ContVars = "Outcome",
    Relabel = FALSE,
    Parametric = TRUE
  )
pvalue_result <- PlotAnovaRelationshipsMatrix(
    df_Anova,
    CatVars = "Group",
    ContVars = "Outcome",
    Relabel = FALSE,
    Parametric = FALSE
  )
effect_scale <- effect_result$Unadjusted$plot$scales$get_scales("colour")
raw_scale <- pvalue_result$Unadjusted$plot$scales$get_scales("colour")
adjusted_scale <- pvalue_result$FDRCorrected$plot$scales$get_scales("colour")
raw_mapping <- rlang::get_expr(
    pvalue_result$Unadjusted$plot$layers[[1]]$mapping$colour
  )
adjusted_mapping <- rlang::get_expr(
    pvalue_result$FDRCorrected$plot$layers[[1]]$mapping$colour
  )
expect_identical(effect_scale$name, "Effect Size")
expect_identical(raw_scale$name, "P-value")
expect_identical(adjusted_scale$name, "FDR-adjusted P-value")
