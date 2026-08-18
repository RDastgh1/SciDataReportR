# Extracted from test-make-comparison-table-effect-sizes.R:233

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
set.seed(107)
df_Test <- data.frame(
    group = factor(rep(c("Control", "Case"), each = 20)),
    age = stats::rnorm(40),
    outcome = stats::rnorm(40),
    category = factor(rep(c("No", "Yes"), 20))
  )
labelled::var_label(df_Test$group) <- "Diagnostic group"
labelled::var_label(df_Test$age) <- "Age (years)"
tbl_unadjusted <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome", "category"),
    CatMethod = "auto"
  )
tbl_adjusted <- MakeComparisonTable(
    data = df_Test,
    group_var = "group",
    variables = c("outcome", "category"),
    covariates = "age",
    AddPairwise = TRUE,
    PairwiseMethod = "holm",
    ShowNotes = "always"
  )
caption_unadjusted <- GetComparisonCaption(tbl_unadjusted)
caption_adjusted <- GetComparisonCaption(tbl_adjusted)
note_unadjusted <- tbl_unadjusted$table_styling$footnote_header %>%
    dplyr::filter(.data$column == "p.value_fmt") %>%
    dplyr::pull("footnote")
note_adjusted <- tbl_adjusted$table_styling$footnote_header %>%
    dplyr::filter(.data$column == "p.value_fmt") %>%
    dplyr::pull("footnote")
expect_match(caption_unadjusted, "Comparison by Diagnostic group", fixed = TRUE)
expect_false(grepl("adjusted for", caption_unadjusted, fixed = TRUE))
expect_false(grepl("Pairwise", caption_unadjusted, fixed = TRUE))
