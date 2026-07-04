make_baseline <- function() {
  data.frame(
    id = 1:5,
    age = c(50, 61, 45, 58, 63),
    stringsAsFactors = FALSE
  )
}

test_that("safe_merge performs a clean exact merge", {
  baseline <- make_baseline()
  labs <- data.frame(
    id = 1:5,
    glucose = c(90, 110, 100, 95, 105)
  )

  m <- safe_merge(baseline, labs, by = "id", name = "Baseline + labs")

  expect_named(m, c("data", "validation", "log", "summary"))

  expect_s3_class(m$data, "data.frame")
  expect_equal(nrow(m$data), nrow(baseline))
  expect_true("glucose" %in% names(m$data))

  expect_type(m$validation, "list")
  expect_true("Checks" %in% names(m$validation))

  expect_s3_class(m$log, "tbl_df")
  expect_equal(nrow(m$log), 1)
  expect_equal(m$log$Merge, "Baseline + labs")
  expect_equal(m$log$Status, "PASS")
  expect_true(m$log$ReadyForAnalysis)
  expect_equal(m$log$RowsBefore, 5)
  expect_equal(m$log$RowsAfter, 5)
  expect_equal(m$log$ColsBefore, 2)
  expect_equal(m$log$ColsAfter, 3)
  expect_equal(m$log$MatchedKeys, 5)
  expect_equal(m$log$DuplicateKeyGroups, 0)
  expect_equal(m$log$UnresolvedDupVars, 0)

  expect_s3_class(m$summary, "knitr_kable")
  expect_true(any(grepl("Baseline \\+ labs", m$summary)))
})

test_that("safe_merge flags a duplicate key as FAIL", {
  baseline <- make_baseline()
  labs <- data.frame(
    id = c(1, 2, 2, 3),
    glucose = c(90, 110, 115, 100)
  )

  m <- safe_merge(baseline, labs, by = "id", name = "Duplicated labs")

  expect_equal(m$log$Status, "FAIL")
  expect_false(m$log$ReadyForAnalysis)
  expect_gt(m$log$DuplicateKeyGroups, 0)
  expect_gt(m$log$RowsAfter, m$log$RowsBefore)

  fail_checks <- m$validation$Checks$Check[m$validation$Checks$Status == "FAIL"]
  expect_true("Duplicate Keys" %in% fail_checks)
})

test_that("safe_merge supports closest_time merges with a date column", {
  visits <- data.frame(
    id = c(1, 1, 2, 3),
    visit_date = as.Date(c("2024-01-10", "2024-06-15", "2024-02-01", "2024-03-20")),
    score = c(10, 12, 8, 15)
  )
  labs <- data.frame(
    id = c(1, 1, 2, 3),
    lab_date = as.Date(c("2024-01-12", "2024-06-10", "2024-02-05", "2024-03-25")),
    glucose = c(90, 95, 110, 100)
  )

  m <- safe_merge(
    visits,
    labs,
    by = "id",
    name = "Visits + closest labs",
    method = "closest_time",
    time_var_before = "visit_date",
    time_var_add = "lab_date",
    is_date = TRUE
  )

  expect_named(m, c("data", "validation", "log", "summary"))
  expect_equal(nrow(m$data), nrow(visits))
  expect_true("glucose" %in% names(m$data))

  # Each visit should pick up the nearest lab draw for that id
  merged_id1 <- m$data[m$data$id == "1", ]
  expect_equal(sort(merged_id1$glucose), c(90, 95))

  # Repeated ids in the visit data are legitimate repeated visits, but
  # ValidateMerge audits duplicate keys on `by` alone, so this is expected
  # to surface as a duplicate-key FAIL (documented caveat).
  expect_equal(m$log$Status, "FAIL")
  expect_gt(m$log$DuplicateKeyGroups, 0)
})

test_that("safe_merge validates its arguments", {
  baseline <- make_baseline()
  labs <- data.frame(id = 1:5, glucose = 1:5)

  expect_error(
    safe_merge(baseline, labs, by = "id"),
    "name is required"
  )

  expect_error(
    safe_merge(baseline, labs, by = "id", name = c("a", "b")),
    "single non-empty character"
  )

  expect_error(
    safe_merge(baseline, labs, by = 1, name = "bad by"),
    "character vector"
  )

  expect_error(
    safe_merge(
      baseline, labs,
      by = "id", name = "no time vars",
      method = "closest_time"
    ),
    "time_var_before and time_var_add"
  )
})
