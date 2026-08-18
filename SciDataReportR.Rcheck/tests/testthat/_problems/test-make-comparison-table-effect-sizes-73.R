# Extracted from test-make-comparison-table-effect-sizes.R:73

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
GetMainEffectRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}
GetComparisonCaption <- function(tbl) {
  as.character(tbl$table_styling$caption)
}

# test -------------------------------------------------------------------------
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
