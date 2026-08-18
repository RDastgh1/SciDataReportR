# Extracted from test-revalue-data-diagnostics.R:26

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
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
