test_that("ReplaceMissingCode handles several markers for one variable", {
  df <- data.frame(
    id = 1:6,
    age = c(34, 999, 52, -7, 41, -8),
    score = c(10, -9, -9, 15, 20, 12)
  )

  codebook <- data.frame(
    Variable = c("age", "score"),
    MissingCode = c("999, -7, -8", "-9")
  )

  out <- ReplaceMissingCode(df, codebook)

  expect_equal(out$age, c(34, NA, 52, NA, 41, NA))
  expect_equal(out$score, c(10, NA, NA, 15, 20, 12))
  expect_equal(out$id, df$id)
})

test_that("ReplaceMissingCode accepts commas, semicolons, and spacing variants", {
  df <- data.frame(age = c(34, 999, -7, -8, 41))
  expected <- c(34, NA, NA, NA, 41)

  for (spec in c("999, -7, -8", "999,-7,-8", "999; -7; -8", "999;-7;-8", " 999 , -7 ,-8 ")) {
    out <- ReplaceMissingCode(df, data.frame(Variable = "age", MissingCode = spec))
    expect_equal(out$age, expected, info = spec)
  }
})

test_that("ReplaceMissingCode accepts a list-column of codes", {
  df <- data.frame(age = c(34, 999, -7, 41))
  codebook <- data.frame(Variable = "age", MissingCode = I(list(c(999, -7))))

  out <- ReplaceMissingCode(df, codebook)

  expect_equal(out$age, c(34, NA, NA, 41))
})

test_that("ReplaceMissingCode is a no-op when no variable carries a code", {
  df <- data.frame(id = 1:3, age = c(34, 999, 52))
  codebook <- data.frame(
    Variable = c("id", "age"),
    MissingCode = c(NA_character_, NA_character_)
  )

  expect_silent(out <- ReplaceMissingCode(df, codebook))
  expect_equal(out, df)

  # Blank strings are skipped the same way as NA.
  blank <- data.frame(Variable = c("id", "age"), MissingCode = c("", ""))
  expect_equal(ReplaceMissingCode(df, blank), df)
})

test_that("ReplaceMissingCode warns instead of failing on variables absent from data", {
  df <- data.frame(id = 1:3, age = c(34, 999, 52))
  codebook <- data.frame(
    Variable = c("age", "not_selected"),
    MissingCode = c("999", "-9")
  )

  expect_warning(out <- ReplaceMissingCode(df, codebook), "not_selected")
  # The variables that are present are still cleaned.
  expect_equal(out$age, c(34, NA, 52))
})

test_that("ReplaceMissingCode matches text markers on character and factor columns", {
  df_chr <- data.frame(
    grp = c("A", "Unknown", "B", "Unknown"),
    stringsAsFactors = FALSE
  )
  out_chr <- ReplaceMissingCode(
    df_chr, data.frame(Variable = "grp", MissingCode = "Unknown")
  )
  expect_equal(out_chr$grp, c("A", NA, "B", NA))

  df_fct <- data.frame(grp = factor(c("A", "Refused", "B", "Refused", "C")))
  out_fct <- ReplaceMissingCode(
    df_fct, data.frame(Variable = "grp", MissingCode = "Refused")
  )
  expect_equal(as.character(out_fct$grp), c("A", NA, "B", NA, "C"))
  # The marker stops being a level, so it cannot show up as an empty category.
  expect_false("Refused" %in% levels(out_fct$grp))
})

test_that("ReplaceMissingCode leaves unrelated empty factor levels alone", {
  df <- data.frame(
    grp = factor(c("A", "Refused", "B"), levels = c("A", "B", "C", "Refused"))
  )

  out <- ReplaceMissingCode(df, data.frame(Variable = "grp", MissingCode = "Refused"))

  expect_equal(levels(out$grp), c("A", "B", "C"))
})

test_that("ReplaceMissingCode does not coerce a numeric column to match text", {
  df <- data.frame(age = c(34, 999, 52))

  # A text marker cannot match a numeric column; nothing should change and no
  # coercion warning should be raised.
  expect_silent(
    out <- ReplaceMissingCode(df, data.frame(Variable = "age", MissingCode = "Unknown"))
  )
  expect_equal(out$age, df$age)
})

test_that("ReplaceMissingCode validates the codebook", {
  df <- data.frame(age = c(34, 999))

  expect_error(ReplaceMissingCode(df, "not a data frame"), "data frame")
  expect_error(
    ReplaceMissingCode(df, data.frame(Variable = "age")),
    "MissingCode"
  )
  expect_error(
    ReplaceMissingCode(df, data.frame(MissingCode = "999")),
    "Variable"
  )
})
