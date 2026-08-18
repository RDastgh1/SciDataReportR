# Extracted from test-revalue-data-diagnostics.R:10

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
df_Test <- data.frame(Status = c(0, 1))
df_Codebook <- data.frame(
    Variable = "Status",
    Recode = "yess",
    Code = "0=No;1=Yes",
    Type = "Categorical"
  )
out <- RevalueData(df_Test, df_Codebook)
