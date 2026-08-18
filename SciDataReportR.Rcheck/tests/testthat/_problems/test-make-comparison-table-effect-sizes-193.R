# Extracted from test-make-comparison-table-effect-sizes.R:193

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
set.seed(106)
df_Test <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = 15)),
    outcome_two = stats::rnorm(45),
    outcome_three = stats::rnorm(45),
    category = factor(rep(c("No", "Yes", "No"), 15))
  )
df_Test$outcome_two[df_Test$group == "C"] <- NA_real_
tbl_without_es <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome_two", "outcome_three", "category")
  )
tbl_with_es <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome_two", "outcome_three", "category"),
    AddEffectSize = TRUE
  )
caption_without_es <- GetComparisonCaption(tbl_without_es)
caption_with_es <- GetComparisonCaption(tbl_with_es)
expect_false(grepl("Effect sizes:", caption_without_es, fixed = TRUE))
expect_true(grepl("|d|: 0.2/0.5/0.8", caption_with_es, fixed = TRUE))
expect_true(grepl("Cohen's f: 0.1/0.25/0.4", caption_with_es, fixed = TRUE))
expect_true(grepl("Cramer's V: 0.1/0.3/0.5", caption_with_es, fixed = TRUE))
expect_true(grepl("d and f are not numerically equivalent", caption_with_es, fixed = TRUE))
