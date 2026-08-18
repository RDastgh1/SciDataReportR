# Extracted from test-formatted-data-dictionary.R:10

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("gt")
testthat::local_mocked_bindings(
    ScidrPackageAvailable = function(package) {
      if (identical(package, "codebook")) return(FALSE)
      base::requireNamespace(package, quietly = TRUE)
    },
    .package = "SciDataReportR"
  )
