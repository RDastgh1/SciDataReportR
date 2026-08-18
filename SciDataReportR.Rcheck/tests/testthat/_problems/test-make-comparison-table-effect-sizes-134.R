# Extracted from test-make-comparison-table-effect-sizes.R:134

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
