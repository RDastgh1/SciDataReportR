# Extracted from test-make-comparison-table-adjusted-pvalues.R:67

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
GetAdjustedMainRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}

# test -------------------------------------------------------------------------
skip_if_not_installed("car")
set.seed(202)
df_Test <- data.frame(
    "Group Status" = factor(rep(c("Control", "Group B", "Group C"), each = 24)),
    Age = stats::rnorm(72),
    check.names = FALSE
  )
df_Test$outcome <- c(Control = 0, `Group B` = 0.5, `Group C` = 1.0)[df_Test[["Group Status"]]] +
    0.5 * df_Test$Age +
    stats::rnorm(72)
fit <- stats::lm(outcome ~ `Group Status` + Age, data = df_Test)
expected <- as.numeric(car::Anova(fit, type = 2)["`Group Status`", "Pr(>F)"])
tbl <- MakeComparisonTable(
    data = df_Test,
    group_var = "Group Status",
    variables = "outcome",
    covariates = "Age",
    AddEffectSize = TRUE
  )
got <- GetAdjustedMainRow(tbl, "outcome")
expect_identical(got$Test, "ANCOVA (Type II)")
expect_equal(got$p.value, expected, tolerance = 1e-12)
expect_false(is.na(got$p.value_fmt) || !nzchar(got$p.value_fmt))
