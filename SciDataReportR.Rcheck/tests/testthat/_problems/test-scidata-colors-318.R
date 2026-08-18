# Extracted from test-scidata-colors.R:318

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_not_installed("plotly")
data(SampleData, envir = environment())
data(SampleVariableTypes, envir = environment())
labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
pca <- suppressWarnings(CreatePCAObject(
    data = labelled,
    VarsToReduce = c("AXL", "Adiponectin", "Ferritin", "MMP7", "tau", "p_tau")
  ))
plot_obj <- suppressWarnings(plotPCA(
    pca, Var = "Genotype", Mode = "2D", Components = c("RC1", "RC2")
  ))
built <- suppressWarnings(plotly::plotly_build(plot_obj))
to_hex <- function(x) {
    parts <- regmatches(x, regexec("rgba[(]([0-9]+),([0-9]+),([0-9]+)", x))[[1]]
    if (length(parts) == 4) {
      toupper(grDevices::rgb(
        as.numeric(parts[2]), as.numeric(parts[3]), as.numeric(parts[4]),
        maxColorValue = 255
      ))
    } else {
      NA_character_
    }
  }
trace_colors <- vapply(
    built$x$data,
    function(trace) {
      if (is.null(trace$marker$color)) NA_character_
      else to_hex(as.character(trace$marker$color)[[1]])
    },
    character(1)
  )
trace_colors <- trace_colors[!is.na(trace_colors)]
expect_setequal(trace_colors, .SciDataColorValues(length(trace_colors)))
