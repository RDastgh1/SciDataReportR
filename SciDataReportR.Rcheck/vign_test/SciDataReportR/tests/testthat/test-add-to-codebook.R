test_that("AddToCodebook populates existing user-defined columns", {
  df_codebook <- data.frame(
    Variable = "sex",
    Label = "Sex assigned at birth",
    Type = "Double",
    Domain = "Clinical",
    stringsAsFactors = FALSE
  )

  actual <- AddToCodebook(
    df_codebook,
    "age",
    "Age at enrollment",
    "Double",
    Domain = "Clinical"
  )

  expect_identical(names(actual), names(df_codebook))
  expect_identical(actual$Domain[[2]], "Clinical")
})

test_that("AddToCodebook warns and creates new user-defined columns", {
  df_codebook <- data.frame(Variable = "sex", stringsAsFactors = FALSE)

  expect_warning(
    actual <- AddToCodebook(df_codebook, "age", Domain = "Clinical"),
    "not an existing codebook column"
  )

  expect_identical(names(actual), c("Variable", "Label", "Domain"))
  expect_true(is.na(actual$Domain[[1]]))
  expect_identical(actual$Domain[[2]], "Clinical")
})

test_that("AddToCodebook validates the codebook and VariableName", {
  expect_error(AddToCodebook(list(), "age"), "data frame")
  expect_error(AddToCodebook(data.frame(x = 1), "age"), "Variable")

  df_codebook <- data.frame(Variable = "age", stringsAsFactors = FALSE)
  expect_error(AddToCodebook(df_codebook, NA_character_), "non-missing")
  expect_error(AddToCodebook(df_codebook, " "), "non-empty")
  expect_error(AddToCodebook(df_codebook, "age"), "already exists")
})

test_that("AddToCodebook validates custom field inputs and standard-field collisions", {
  df_codebook <- data.frame(Variable = "age", stringsAsFactors = FALSE)

  expect_error(AddToCodebook(df_codebook, "sex", "a", Domain = c("A", "B")), "single atomic")
  expect_error(AddToCodebook(df_codebook, "sex", "a", Domain = list("A")), "single atomic")
  expect_error(AddToCodebook(df_codebook, "sex", "a", "Double", Type = "Double"), "formal")
})

test_that("AddToCodebook warns for unseen values and storage type mismatches", {
  df_codebook <- data.frame(
    Variable = c("sex", "site"),
    Type = c("Categorical", "Categorical"),
    Recode = c(1, 1),
    Exclude = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  expect_warning(
    AddToCodebook(
      df_codebook,
      "age",
      VariableType = "Double"
    ),
    "`Type` has not previously appeared"
  )
  expect_warning(
    AddToCodebook(df_codebook, "age_recode", VariableRecode = 4),
    "`Recode` has not previously appeared"
  )
  expect_warning(
    actual <- AddToCodebook(df_codebook, "age_exclude", VariableExclude = 3),
    "`Exclude` has not previously appeared"
  )
  expect_warning(
    AddToCodebook(
      df_codebook,
      "age_2",
      VariableRecode = "1"
    ),
    "different storage type"
  )
  expect_identical(actual$Exclude[[3]], 3)
  expect_type(actual$Exclude, "double")
})

test_that("AddToCodebook accepts values when the old reference column is all missing", {
  df_codebook <- data.frame(
    Variable = "age",
    Recode = NA_real_,
    stringsAsFactors = FALSE
  )

  expect_no_warning(AddToCodebook(df_codebook, "sex", VariableRecode = 4))
})

test_that("AddToCodebook creates supplied legacy fields and supports CB", {
  df_codebook <- data.frame(Variable = "age", stringsAsFactors = FALSE)

  actual <- AddToCodebook(df_codebook, "sex", VariableRecode = 1)
  expect_true("Recode" %in% names(actual))
  expect_identical(actual$Recode[[2]], 1)

  expect_warning(
    actual_deprecated <- AddToCodebook(CB = df_codebook, VariableName = "sex"),
    "deprecated"
  )
  expect_identical(actual_deprecated$Variable[[2]], "sex")
})
