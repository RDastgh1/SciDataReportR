test_that("ReadSciData reads delimited files through the fast path", {
  skip_if_not_installed("data.table")

  path <- tempfile(fileext = ".csv")
  writeLines(
    c(
      "id,value,label",
      "1,2.5,A",
      "2,3.5,B"
    ),
    path
  )

  df <- ReadSciData(path, inspect = FALSE)

  expect_s3_class(df, "tbl_df")
  expect_equal(names(df), c("id", "value", "label"))
  expect_equal(nrow(df), 2)
  expect_equal(df$value, c(2.5, 3.5))
  expect_true(!is.null(attr(df, "scidata_source")))
  expect_null(attr(df, "scidata_inspection"))
})

test_that("ReadSciData can use the readr delimited path when requested", {
  path <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "id\tvalue",
      "1\t2",
      "2\t3"
    ),
    path
  )

  df <- ReadSciData(path, inspect = FALSE, fast_delimited = FALSE)

  expect_s3_class(df, "tbl_df")
  expect_equal(names(df), c("id", "value"))
  expect_equal(nrow(df), 2)
})
