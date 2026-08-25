# Extracted from test-make-comparison-table-adjusted-pvalues.R:34

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
skip_if_not_installed("emmeans")
set.seed(201)
df_Test <- data.frame(
    "Group Status" = factor(rep(c("Control", "Case"), each = 30)),
    Age = stats::rnorm(60),
    check.names = FALSE
  )
df_Test$outcome <- 0.8 * (df_Test[["Group Status"]] == "Case") +
    0.4 * df_Test$Age +
    stats::rnorm(60)
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
expect_true(is.finite(got$p.value))
