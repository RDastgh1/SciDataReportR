MakeSplitCorrelationData <- function() {
  set.seed(42)
  df_Test <- data.frame(
    Compartment = factor(rep(c("Plasma", "CSF"), each = 30),
      levels = c("CSF", "Plasma", "Unused")),
    x = rnorm(60), z = rnorm(60), y = rnorm(60), age = rnorm(60)
  )
  df_Test$y <- df_Test$x + df_Test$age + df_Test$y
  df_Test$x <- sjlabelled::set_label(df_Test$x, "Shared label")
  df_Test$z <- sjlabelled::set_label(df_Test$z, "Shared label")
  df_Test$y <- sjlabelled::set_label(df_Test$y, "Cognitive score")
  df_Test$Compartment <- sjlabelled::set_label(df_Test$Compartment, "Specimen type")
  df_Test
}

test_that("split results agree with independent stratum calls and FDR scopes", {
  df_Test <- MakeSplitCorrelationData()
  for (scope in c("matrix", "per_outcome", "per_predictor")) {
    for (covariates in list(NULL, "age")) {
      out <- PlotSplitCorrelationHeatmap(
        df_Test, c("z", "x"), c("y", "age"), "Compartment",
        covariates = covariates, fdr_scope = scope, method = "spearman"
      )
      expect_identical(out$Metadata$SplitLevels, c("CSF", "Plasma"))
      expect_identical(out$Metadata$AxisOrder$rows, c("z", "x"))
      expect_identical(out$Metadata$AxisOrder$columns, c("y", "age"))
      for (level in out$Metadata$SplitLevels) {
        expected <- PlotCorrelationsHeatmap(
          df_Test[df_Test$Compartment == level, ], c("z", "x"), c("y", "age"),
          covariates = covariates, fdr_scope = scope, method = "spearman"
        )
        for (branch in c("Unadjusted", "FDRCorrected")) {
          for (stat in c("r", "p", "npairs")) {
            expect_equal(out[[branch]][[stat]][[level]], expected[[branch]][[stat]])
          }
        }
        df_Level <- out$Results[out$Results$Split == level, ]
        expect_equal(df_Level$R, as.vector(expected$Unadjusted$r))
        expect_equal(df_Level$P_adj, as.vector(expected$FDRCorrected$p))
      }
      expect_identical(out$p, out$Unadjusted)
      expect_identical(out$p_fdr, out$FDRCorrected)
      expect_identical(out$Metadata$FDRFamily, "within_stratum")
    }
  }
})

test_that("split validation, selection and exclusion counts are explicit", {
  df_Test <- MakeSplitCorrelationData()
  expect_error(PlotSplitCorrelationHeatmap(list(), split_var = "Compartment"), "data frame")
  expect_error(PlotSplitCorrelationHeatmap(df_Test, split_var = "missing"), "split_var")
  for (levels in list(c("CSF", "CSF"), c("CSF", "Unused"), NA_character_, character())) {
    expect_error(PlotSplitCorrelationHeatmap(df_Test, "x", "y", "Compartment",
      split_levels = levels), "two distinct observed")
  }
  df_Test$Compartment <- as.character(df_Test$Compartment)
  expect_identical(PlotSplitCorrelationHeatmap(df_Test, "x", "y", "Compartment")$
    Metadata$SplitLevels, c("Plasma", "CSF"))
  df_Test$Compartment[1] <- "Other"
  df_Test$Compartment[2] <- NA_character_
  expect_error(PlotSplitCorrelationHeatmap(df_Test, "x", "y", "Compartment"), "Specify `split_levels`")
  expect_warning(out <- PlotSplitCorrelationHeatmap(
    df_Test, "x", "y", "Compartment", split_levels = c("CSF", "Plasma")
  ), "1 rows with missing split")
  expect_equal(out$Metadata$MissingSplitRows, 1)
  expect_equal(out$Metadata$UnselectedSplitRows, 1)
  expect_equal(out$Metadata$StratumRows, c(CSF = 30L, Plasma = 28L))
  df_Test <- df_Test[!is.na(df_Test$Compartment), ]
  expect_error(PlotSplitCorrelationHeatmap(df_Test, "typo", "y", "Compartment",
    split_levels = c("CSF", "Plasma")), "typo")
})

test_that("unavailable correlations retain their triangles and counts", {
  df_Test <- MakeSplitCorrelationData()
  df_Test$x[df_Test$Compartment == "CSF"] <- NA_real_
  df_Test$z[df_Test$Compartment == "Plasma"] <- 1
  out <- PlotSplitCorrelationHeatmap(df_Test, c("x", "z"), "y", "Compartment")
  expect_equal(nrow(out$Results), 4L)
  expect_equal(nrow(out$Metadata$Unavailable), 2L)
  expect_equal(out$Results$nPairs[out$Results$Split == "CSF" & out$Results$XVar == "x"], 0)
  expect_equal(out$Results$nPairs[out$Results$Split == "Plasma" & out$Results$XVar == "z"], 30)
  df_Built <- ggplot2::ggplot_build(out$Unadjusted$plot)$data[[1]]
  expect_equal(sum(df_Built$fill == "grey95"), 6L)
  expect_equal(nrow(df_Built), 12L)
})

test_that("duplicate labels retain distinct cells, correct vertices and stars", {
  df_Test <- MakeSplitCorrelationData()
  out <- PlotSplitCorrelationHeatmap(df_Test, c("z", "x"), "y", "Compartment",
    IncludeTrendSymbol = TRUE)
  expect_equal(length(unique(out$Results$CellId)), 2L)
  expect_true(all(out$Results$XLabel == "Shared label"))
  expect_true(all(out$Results$YLabel == "Cognitive score"))
  expect_true(all(grepl("Specimen type", out$Results$Tooltip)))
  expect_identical(out$Results$stars, SciDataReportR:::.FormatPValueStars(
    out$Results$P, IncludeTrendSymbol = TRUE))
  expect_identical(out$Results$stars_FDR, SciDataReportR:::.FormatPValueStars(
    out$Results$P_adj, IncludeTrendSymbol = TRUE))
  df_Poly <- out$Unadjusted$plot$data
  expect_equal(df_Poly$x[1:3], c(0.52, 0.52, 1.48))
  expect_equal(df_Poly$y[1:3], c(1.52, 2.48, 2.48))
  expect_equal(df_Poly$x[7:9], c(0.52, 1.48, 1.48))
  expect_equal(df_Poly$y[7:9], c(1.52, 1.52, 2.48))
  expect_identical(out$Unadjusted$plot$scales$get_scales("shape")$name, "Specimen type")
  expect_identical(out$Unadjusted$plot$scales$get_scales("fill")$limits, c(-1, 1))
  for (branch in c("Unadjusted", "FDRCorrected")) {
    df_Stars <- out[[branch]]$plot$layers[[2]]$data
    expect_true(all(df_Stars$y[df_Stars$Triangle == "top-left"] >
      df_Stars$y[df_Stars$Triangle == "bottom-right"]))
    expect_silent(ggplot2::ggplotGrob(out[[branch]]$plot))
  }
  no_labels <- PlotSplitCorrelationHeatmap(df_Test, "x", "y", "Compartment", Relabel = FALSE)
  expect_true(all(no_labels$Results$XLabel == "x"))
})

test_that("numeric splits are excluded from default predictor selection", {
  df_Test <- MakeSplitCorrelationData()
  df_Test$Compartment <- as.integer(df_Test$Compartment)
  expect_warning(
    out <- PlotSplitCorrelationHeatmap(df_Test, split_var = "Compartment"),
    "converting locally to a factor"
  )
  expect_false("Compartment" %in% out$Metadata$AxisOrder$rows)
  expect_equal(out$Metadata$AxisOrder$rows, out$Metadata$AxisOrder$columns)
  expect_true(all(is.na(out$Results$R[out$Results$XVar == out$Results$YVar])))
})

test_that("optional widgets retain tooltip content and leave static layers intact", {
  skip_if_not_installed("plotly")
  df_Test <- MakeSplitCorrelationData()
  out <- PlotSplitCorrelationHeatmap(df_Test, c("z", "x"), "y", "Compartment",
    interactive = "plotly")
  for (branch in c("Unadjusted", "FDRCorrected")) {
    expect_s3_class(out$Interactive$Plotly[[branch]], "plotly")
    expect_s3_class(out[[branch]]$plot, "ggplot")
    expect_null(out[[branch]]$plot$layers[[1]]$mapping$text)
    widget <- plotly::plotly_build(out$Interactive$Plotly[[branch]])
    tooltip <- unlist(lapply(widget$x$data, function(trace) trace$text))
    for (term in c("CSF", "Plasma", "Shared label", "Cognitive score",
      "FDR p-value", "Complete pairs")) {
      expect_true(any(grepl(term, tooltip, fixed = TRUE)))
    }
  }
})
