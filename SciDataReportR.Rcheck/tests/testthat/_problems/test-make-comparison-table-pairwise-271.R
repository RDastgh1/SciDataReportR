# Extracted from test-make-comparison-table-pairwise.R:271

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
GetComparisonTableMainRow <- function(tbl, variable) {
  tbl$table_body %>%
    dplyr::filter(.data$variable == .env$variable, .data$row_type == "label") %>%
    dplyr::slice(1)
}
GetComparisonTablePairwiseP <- function(tbl, variable, contrast) {
  headers <- tbl$table_styling$header
  column <- headers$column[headers$label == paste0("**", contrast, "**")]
  row <- GetComparisonTableMainRow(tbl, variable)
  as.numeric(row[[column]])
}
PairwiseLogisticLR <- function(data, outcome, group, covariates, pairs) {
  raw_p <- vapply(pairs, function(pair) {
    df_pair <- data[data[[group]] %in% pair, , drop = FALSE]
    df_pair[[group]] <- droplevels(df_pair[[group]])
    full <- stats::glm(
      stats::reformulate(c(group, covariates), response = outcome),
      data = df_pair,
      family = stats::binomial()
    )
    reduced <- stats::glm(
      stats::reformulate(covariates, response = outcome),
      data = df_pair,
      family = stats::binomial()
    )
    stats::anova(reduced, full, test = "Chisq")$`Pr(>Chi)`[2]
  }, numeric(1))
  stats::p.adjust(raw_p, method = "holm")
}
PairwiseMultinomialLR <- function(data, outcome, group, covariates, pairs) {
  skip_if_not_installed("nnet")

  raw_p <- vapply(pairs, function(pair) {
    df_pair <- data[data[[group]] %in% pair, , drop = FALSE]
    df_pair[[group]] <- droplevels(df_pair[[group]])
    full <- nnet::multinom(
      stats::reformulate(c(group, covariates), response = outcome),
      data = df_pair,
      trace = FALSE
    )
    reduced <- nnet::multinom(
      stats::reformulate(covariates, response = outcome),
      data = df_pair,
      trace = FALSE
    )
    lr <- 2 * (as.numeric(stats::logLik(full)) - as.numeric(stats::logLik(reduced)))
    df_lr <- attr(stats::logLik(full), "df") - attr(stats::logLik(reduced), "df")
    stats::pchisq(lr, df = df_lr, lower.tail = FALSE)
  }, numeric(1))
  stats::p.adjust(raw_p, method = "holm")
}

# test -------------------------------------------------------------------------
skip_if_not_installed("gt")
data(SampleData, envir = environment())
data(SampleVariableTypes, envir = environment())
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
Labelled$ApoeGroup <- factor(
    dplyr::case_when(
      Labelled$Genotype %in% c("E2E2", "E2E3") ~ "E2 carrier",
      Labelled$Genotype == "E3E3" ~ "E3/E3",
      TRUE ~ "E4 carrier"
    ),
    levels = c("E2 carrier", "E3/E3", "E4 carrier")
  )
merged <- MakeFacetCatComparisonTable(
    data = Labelled,
    FacetVariables = c("Diagnosis", "ApoeGroup"),
    variables = c("age", "tau", "p_tau"),
    AddPairwise = TRUE
  )
text_format <- merged$table_styling$text_format
pairwise_rules <- text_format[grepl("^pw_", text_format$column), ]
expect_gt(nrow(pairwise_rules), 0)
for (i in seq_len(nrow(pairwise_rules))) {
    referenced <- all.vars(rlang::quo_get_expr(pairwise_rules$rows[[i]]))
    referenced <- setdiff(referenced, ".data")
    expect_true(
      all(referenced %in% names(merged$table_body)),
      info = paste("rule for", pairwise_rules$column[i])
    )
  }
expect_no_error(rendered <- gt::as_raw_html(gtsummary::as_gt(merged)))
