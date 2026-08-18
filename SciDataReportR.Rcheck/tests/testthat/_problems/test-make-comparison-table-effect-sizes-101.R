# Extracted from test-make-comparison-table-effect-sizes.R:101

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
skip_if_not_installed("effectsize")
set.seed(103)
df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), times = c(18, 25, 21))),
    outcome = stats::rnorm(64)
  )
df_Test$outcome <- df_Test$outcome +
    c(A = 0, B = 0.5, C = 1)[df_Test$group]
fit <- stats::aov(outcome ~ group, data = df_Test)
eta <- effectsize::eta_squared(fit, partial = FALSE)$Eta2[1]
expected <- sqrt(eta / (1 - eta))
tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = "outcome",
    AddEffectSize = TRUE
  )
got <- GetMainEffectRow(tbl, "outcome")
expect_equal(got$effect_size, expected, tolerance = 1e-10)
expect_identical(got$ES_Method, "Cohen's f")
