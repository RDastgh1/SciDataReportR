test_that("RevalueData flags unrecognized Recode values", {
  df_Test <- data.frame(Status = c(0, 1))
  df_Codebook <- data.frame(
    Variable = "Status",
    Recode = "yess",
    Code = "0=No;1=Yes",
    Type = "Categorical"
  )

  out <- RevalueData(df_Test, df_Codebook)

  expect_true(any(grepl("Status: Unrecognized Recode value 'yess'", out$warninglist, fixed = TRUE)))
  expect_equal(out$errors, data.frame(Variable = character(), Error = character()))
})

test_that("RevalueData identifies the variable that fails", {
  df_Test <- data.frame(Bad = c(0, 1), Good = c(1, 2))
  df_Codebook <- data.frame(
    Variable = c("Bad", "Good"),
    Recode = c("Yes", "No"),
    Code = c("0=No;1=Yes", NA_character_),
    Type = c("Categorical", "Double")
  )

  testthat::local_mocked_bindings(
    set_labels = function(...) stop("Invalid codebook mapping"),
    .package = "sjlabelled"
  )

  expect_error(RevalueData(df_Test, df_Codebook), "Error revaluing variable 'Bad'")

  out <- RevalueData(
    df_Test,
    df_Codebook,
    on_error = "warn"
  )

  expect_equal(out$errors$Variable, "Bad")
  expect_true(grepl("Invalid codebook mapping", out$errors$Error, fixed = TRUE))
  expect_true(any(grepl("Error revaluing variable 'Bad'", out$warninglist, fixed = TRUE)))
  expect_equal(out$RevaluedData$Good, c(1, 2))
})
test_that("RevalueData needs only a Variable column in the codebook", {
  df <- data.frame(year = c(2005L, 2007L), hours = c(7, 99))
  codebook <- data.frame(
    Variable = c("year", "hours"),
    Label = c("Cycle", "Sleep hours"),
    MissingCode = c(NA, "77,99")
  )

  result <- RevalueData(df, codebook)$RevaluedData

  expect_identical(result$year, c(2005L, 2007L), ignore_attr = TRUE)
  expect_identical(attr(result$year, "label"), "Cycle")
  expect_equal(result$hours, c(7, NA), ignore_attr = TRUE)

  minimal <- RevalueData(df, data.frame(Variable = "year"))$RevaluedData
  expect_identical(minimal$year, c(2005L, 2007L))
})
