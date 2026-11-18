test_that("partial regression treats unusual column names literally", {
  set.seed(427)
  df_Test <- data.frame(x = rnorm(60), y = rnorm(60), z = rnorm(60), w = rnorm(60))
  df_Test$y <- df_Test$y + 0.7 * df_Test$x + 0.4 * df_Test$z
  df_Test$z[3] <- NA_real_
  vars_Unusual <- c("% marker", "marker value", "marker+value", "marker`value")
  vars_Roles <- c("x", "y", "z", "w")

  for (adjusted in c(FALSE, TRUE)) {
    vars_Covariates <- if (adjusted) c("z", "w") else NULL
    reference <- PlotPartialRegressionScatter(df_Test, "x", "y",
      covariates = vars_Covariates, Relabel = FALSE)

    # Exercise each name in every role, then all unusual names together.
    name_sets <- lapply(vars_Roles, function(role) {
      lapply(vars_Unusual, function(unusual) {
        vars_Names <- vars_Roles
        vars_Names[vars_Roles == role] <- unusual
        vars_Names
      })
    })
    name_sets <- c(unlist(name_sets, recursive = FALSE), list(vars_Unusual))

    for (vars_Names in name_sets) {
      df_Renamed <- df_Test
      names(df_Renamed) <- vars_Names
      result <- PlotPartialRegressionScatter(df_Renamed, vars_Names[1], vars_Names[2],
        covariates = if (adjusted) vars_Names[3:4] else NULL, Relabel = FALSE)

      expect_equal(result$r2, reference$r2)
      expect_equal(result$p_value, reference$p_value)
      expect_equal(result$n, reference$n)
      expect_gt(nrow(result$Results), 0)
      expect_true(all(result$Results$Status == "ok"))
      for (field in c("Estimate", "ConfLow", "ConfHigh", "PValue", "N")) {
        expect_equal(result$Results[[field]], reference$Results[[field]])
      }
      expect_equal(result$Results$Predictor, vars_Names[1])
      expect_equal(result$Results$Outcome, vars_Names[2])
      expect_s3_class(suppressMessages(ggplot2::ggplot_build(result$plot)), "ggplot_built")
    }
  }
})
