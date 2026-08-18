test_that("FormattedDataDictionary reports a missing codebook dependency", {
  skip_if_not_installed("gt")

  testthat::local_mocked_bindings(
    ScidrPackageAvailable = function(package) {
      if (identical(package, "codebook")) return(FALSE)
      base::requireNamespace(package, quietly = TRUE)
    },
    .package = "SciDataReportR"
  )

  expect_error(
    FormattedDataDictionary(data.frame(Value = 1)),
    "'codebook' package is required"
  )
})

test_that("FormattedDataDictionary returns a gt table when dependencies are installed", {
  skip_if_not_installed("gt")
  skip_if_not_installed("codebook")

  out <- FormattedDataDictionary(data.frame(Value = c(1, 2, 3)))

  expect_s3_class(out, "gt_tbl")
})
