test_that("merge_detail prints kable sections for a problematic merge", {
  baseline <- data.frame(
    id = 1:5,
    age = c(50, 61, 45, 58, 63),
    sex = c("F", "M", "F", "F", "M")
  )
  labs <- data.frame(
    id = c(1, 2, 2, 6),
    sex = c("F", "M", "M", "F"),
    glucose = c(90, 110, 115, 100)
  )

  m <- safe_merge(baseline, labs, by = "id", name = "Messy merge")

  output <- paste(capture.output(merge_detail(m, TopN = 3)), collapse = "\n")

  expect_match(output, "Messy merge: validation checks")
  expect_match(output, "keys only in left data")
  expect_match(output, "keys only in right data")
  expect_match(output, "overlapping non-key variables")

  expect_invisible(merge_detail(m))
})

test_that("merge_detail skips empty sections", {
  baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
  labs <- data.frame(id = 1:4, glucose = c(90, 110, 100, 95))

  m <- safe_merge(baseline, labs, by = "id", name = "Clean merge")

  output <- paste(capture.output(merge_detail(m)), collapse = "\n")

  expect_match(output, "Clean merge: validation checks")
  expect_no_match(output, "keys only in left data")
  expect_no_match(output, "keys only in right data")
  expect_no_match(output, "overlapping non-key variables")
  expect_no_match(output, "suspicious duplicate-variable conflicts")
})

test_that("merge_detail rejects invalid input", {
  expect_error(merge_detail(list(data = 1)), "safe_merge")

  baseline <- data.frame(id = 1:4, age = c(50, 61, 45, 58))
  labs <- data.frame(id = 1:4, glucose = c(90, 110, 100, 95))
  m <- safe_merge(baseline, labs, by = "id", name = "Clean merge")

  expect_error(merge_detail(m, TopN = 0), "positive")
})
