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
