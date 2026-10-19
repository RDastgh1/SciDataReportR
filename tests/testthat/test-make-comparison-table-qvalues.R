MakeQComparisonData <- function() {
  set.seed(831)
  df_Test <- data.frame(
    Group = factor(rep(c("Control", "Case", "Other"), each = 40)),
    Facet = factor(rep(c("A", "B"), 60)),
    Age = stats::rnorm(120)
  )
  df_Test$Continuous <- 2 * as.numeric(df_Test$Group) +
    0.5 * df_Test$Age + stats::rnorm(120)
  df_Test$Binary <- factor(stats::rbinom(120, 1, 0.4))
  df_Test$Category <- factor(sample(c("Low", "Middle", "High"), 120, replace = TRUE))
  # Varies overall, but only one group has observations: no valid group test.
  df_Test$Untestable <- c(stats::rnorm(40), rep(NA_real_, 80))
  df_Test$Continuous[c(2, 9)] <- NA_real_
  df_Test$Category[c(3, 10)] <- NA
  df_Test
}

test_that("addQ counts each main outcome once and excludes missing tests", {
  df_Test <- MakeQComparisonData()
  tbl <- MakeComparisonTable(
    df_Test, "Group", c("Continuous", "Binary", "Category", "Untestable"),
    IncludeMissing = TRUE, addQ = TRUE
  )
  df_Main <- tbl$table_body %>% dplyr::filter(.data$.is_main_row)
  valid <- is.finite(df_Main$p.value)

  expect_equal(sum(valid), 3L)
  expect_equal(df_Main$q.value[valid], stats::p.adjust(df_Main$p.value[valid], "fdr"))
  expect_true(all(is.na(df_Main$q.value[!valid])))
  expect_true(all(is.na(tbl$table_body$q.value[!tbl$table_body$.is_main_row])))
  expect_true(any(tbl$table_body$row_type == "missing"))
  expect_silent(gt::as_raw_html(gtsummary::as_gt(tbl)))
})

test_that("addQ uses reported covariate-adjusted main p-values", {
  df_Test <- MakeQComparisonData()
  tbl <- MakeComparisonTable(
    df_Test, "Group", c("Continuous", "Binary", "Category"),
    covariates = "Age", addQ = TRUE
  )
  df_Main <- tbl$table_body %>% dplyr::filter(.data$.is_main_row)

  expect_true(all(is.finite(df_Main$p_adj)))
  expect_equal(df_Main$p.value, df_Main$p_adj)
  expect_equal(df_Main$q.value, stats::p.adjust(df_Main$p_adj, "fdr"))
})

test_that("addQ preserves default output, pairwise tests, and effect sizes", {
  df_Test <- MakeQComparisonData()
  args <- list(
    data = df_Test, group_var = "Group",
    variables = c("Continuous", "Binary", "Category"),
    AddPairwise = TRUE, PairwiseMethod = "holm", AddEffectSize = TRUE
  )
  tbl_Default <- do.call(MakeComparisonTable, args)
  tbl_False <- do.call(MakeComparisonTable, c(args, list(addQ = FALSE)))
  tbl_Q <- do.call(MakeComparisonTable, c(args, list(addQ = TRUE)))

  expect_identical(tbl_Default$table_body, tbl_False$table_body)
  expect_identical(tbl_Default$table_styling, tbl_False$table_styling)
  expect_false("q.value" %in% names(tbl_Default$table_body))
  expect_identical(
    tbl_Default$table_body,
    tbl_Q$table_body[, names(tbl_Default$table_body)]
  )
  expect_true(any(grepl("^pw_", names(tbl_Q$table_body))))
  expect_true("effect_size" %in% names(tbl_Q$table_body))
  expect_silent(gt::as_raw_html(gtsummary::as_gt(tbl_Q)))
})

test_that("q-values have independent highlighting, precision, and placement", {
  df_Test <- MakeQComparisonData()
  tbl <- MakeComparisonTable(
    df_Test, "Group", c("Continuous", "Binary", "Category"),
    p_digits = 4, addQ = TRUE
  )
  headers <- tbl$table_styling$header
  visible <- headers$column[!headers$hide]
  expect_identical(visible[match("p.value_fmt", visible) + 1L], "q.value")
  expect_identical(headers$label[headers$column == "q.value"], "**q-value**")

  df_Main <- tbl$table_body %>% dplyr::filter(.data$.is_main_row)
  expect_true(any(df_Main$q.value <= 0.05))
  expect_true(any(df_Main$q.value > 0.05))
  rules <- tbl$table_styling$text_format
  q_rule <- rules$rows[[which(rules$column == "q.value" & rules$format_type == "bold")]]
  expect_equal(
    rlang::eval_tidy(q_rule, data = tbl$table_body),
    tbl$table_body$.is_main_row & !is.na(tbl$table_body$q.value) &
      tbl$table_body$q.value <= 0.05
  )
  formatter <- tbl$table_styling$fmt_fun$fmt_fun[[
    which(tbl$table_styling$fmt_fun$column == "q.value")
  ]]
  expect_identical(formatter(c(0.01234, NA_real_)),
                   c(format.pval(0.01234, digits = 4, eps = 1e-4), NA_character_))
  expect_match(paste(tbl$table_styling$footnote_header$footnote, collapse = " "),
               "Benjamini-Hochberg")
})

test_that("addQ handles one or zero valid tests and summary-only tables", {
  df_Test <- MakeQComparisonData()
  tbl_One <- MakeComparisonTable(df_Test, "Group", c("Continuous", "Untestable"), addQ = TRUE)
  df_One <- tbl_One$table_body %>% dplyr::filter(.data$.is_main_row)
  expect_equal(df_One$q.value, df_One$p.value)

  tbl_None <- MakeComparisonTable(df_Test, "Group", "Untestable", addQ = TRUE)
  expect_true(all(is.na(tbl_None$table_body$q.value)))
  expect_silent(gt::as_raw_html(gtsummary::as_gt(tbl_None)))

  tbl_Overall <- MakeComparisonTable(df_Test, variables = "Continuous", addQ = TRUE)
  expect_false("q.value" %in% names(tbl_Overall$table_body))
  df_Test$Group <- factor(rep("Control", nrow(df_Test)))
  tbl_Single <- MakeComparisonTable(df_Test, "Group", "Continuous", addQ = TRUE)
  expect_false("q.value" %in% names(tbl_Single$table_body))
})

test_that("faceted tables correct separately and render with pairwise columns", {
  df_Test <- MakeQComparisonData()
  tbl <- MakeFacetCatComparisonTable(
    df_Test, c("Group", "Facet"), c("Continuous", "Binary", "Category"),
    AddPairwise = TRUE, AddEffectSize = TRUE, addQ = TRUE
  )
  for (facet_index in seq_along(c("Group", "Facet"))) {
    p_values <- tbl$table_body[[paste0("p.value_", facet_index)]]
    q_values <- tbl$table_body[[paste0("q.value_", facet_index)]]
    main_rows <- tbl$table_body[[paste0(".is_main_row_", facet_index)]]
    expect_equal(q_values[main_rows], stats::p.adjust(p_values[main_rows], "fdr"))
    expect_true(all(is.na(q_values[!main_rows])))
  }
  expect_silent(gt::as_raw_html(gtsummary::as_gt(tbl)))
})

test_that("addQ requires a single nonmissing logical value", {
  df_Test <- MakeQComparisonData()
  for (invalid in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(
      MakeComparisonTable(df_Test, "Group", "Continuous", addQ = invalid),
      "addQ must be a single nonmissing logical value", fixed = TRUE
    )
  }
})
