test_that("SciDataPalette keeps its public contract", {
  full <- SciDataPalette()

  expect_length(full, 34)
  expect_named(full)
  expect_equal(unname(full[1:2]), c("#0B1F5E", "#E87400"))
  expect_null(names(SciDataPalette(names = FALSE)))
  expect_length(SciDataPalette(3), 3)

  # Asking an explicitly-sized qualitative palette for more than it has is a
  # user error and must stay one. The internals bypass this deliberately.
  expect_error(SciDataPalette(35), "cannot exceed")
  expect_error(SciDataPalette(0), "positive integer")
  expect_error(SciDataPalette(2.5), "positive integer")
})

test_that(".SciDataColorValues returns the palette verbatim up to its length", {
  base <- SciDataPalette(names = FALSE)

  expect_identical(.SciDataColorValues(34), base)
  expect_identical(.SciDataColorValues(1), base[1])
  expect_identical(.SciDataColorValues(10), base[1:10])
})

test_that(".SciDataColorValues extends past the palette without losing the anchors", {
  base <- SciDataPalette(names = FALSE)

  # Extending must not disturb the branded colors that came before it.
  expect_identical(.SciDataColorValues(40)[1:34], base)
  expect_identical(.SciDataColorValues(170)[1:34], base)
})

test_that(".SciDataColorValues returns exactly n colors for any n", {
  for (n in c(1, 2, 33, 34, 35, 68, 100, 170)) {
    expect_length(.SciDataColorValues(n), n)
  }
})

test_that(".SciDataColorValues stays distinct through its designed range", {
  for (n in c(34, 35, 68, 102, 136, 170)) {
    expect_equal(anyDuplicated(.SciDataColorValues(n)), 0, info = paste("n =", n))
  }
})

test_that(".SciDataColorValues warns rather than failing past its distinct range", {
  expect_warning(values <- .SciDataColorValues(200), "being reused")
  expect_length(values, 200)
})

test_that(".SciDataColorValues is deterministic", {
  expect_identical(.SciDataColorValues(50), .SciDataColorValues(50))
})

test_that(".SciDataColorValues validates n", {
  expect_error(.SciDataColorValues(0), "positive integer")
  expect_error(.SciDataColorValues(-1), "positive integer")
  expect_error(.SciDataColorValues(2.5), "positive integer")
  expect_error(.SciDataColorValues(c(1, 2)), "positive integer")
})

test_that(".SciDataNamedValues holds reserved levels out of the palette", {
  values <- .SciDataNamedValues(
    c("1", "2", "3", "Noise"),
    fixed = c(Noise = "grey60")
  )

  expect_named(values, c("1", "2", "3", "Noise"))
  expect_equal(unname(values[["Noise"]]), "grey60")

  # The reserved level must not consume a palette slot: the real categories
  # take positions 1..3, not 1, 2, 4.
  expect_equal(
    unname(values[c("1", "2", "3")]),
    .SciDataColorValues(3)
  )
})

test_that(".SciDataNamedValues works with no reserved levels and with only reserved levels", {
  plain <- .SciDataNamedValues(c("a", "b"))
  expect_equal(unname(plain), .SciDataColorValues(2))

  all_fixed <- .SciDataNamedValues("Noise", fixed = c(Noise = "grey60"))
  expect_equal(unname(all_fixed), "grey60")
})

test_that(".SciDataContrastText picks a legible color for each background", {
  # The palette spans 1.4:1 to 14:1 against black, so a fixed label color is
  # illegible on part of it.
  expect_equal(
    unname(.SciDataContrastText(c("#0B1F5E", "#641B68"))),
    c("white", "white")
  )
  expect_equal(
    unname(.SciDataContrastText(c("#F2B134", "#8FD8B8", "#E7C6B2"))),
    c("black", "black", "black")
  )
  expect_equal(unname(.SciDataContrastText("white")), "black")
  expect_equal(unname(.SciDataContrastText("black")), "white")
})

test_that(".SciDataContrastText covers the whole palette and returns both values", {
  labels <- .SciDataContrastText(SciDataPalette(names = FALSE))

  expect_length(labels, 34)
  expect_setequal(unique(labels), c("black", "white"))
})

test_that("the internal scales survive more categories than the palette holds", {
  # This is the direct test of the decision that package plots must never fail
  # merely because a variable has many levels.
  df <- data.frame(g = factor(paste0("L", seq_len(40))), y = seq_len(40))

  fill_plot <- ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, fill = g)) +
    ggplot2::geom_col() +
    .SciDataFillScale()

  colour_plot <- ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, colour = g)) +
    ggplot2::geom_point() +
    .SciDataColourScale()

  expect_no_error(ggplot2::ggplot_build(fill_plot))
  expect_no_error(ggplot2::ggplot_build(colour_plot))
})

test_that("the internal scales map the first categories to the palette", {
  df <- data.frame(g = factor(c("a", "b", "c")), y = 1:3)

  built <- ggplot2::ggplot_build(
    ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, fill = g)) +
      ggplot2::geom_col() +
      .SciDataFillScale()
  )

  expect_equal(unique(built$data[[1]]$fill), .SciDataColorValues(3))
})

test_that("the exported scales map to the palette", {
  df <- data.frame(g = factor(c("a", "b", "c")), y = 1:3)

  built_fill <- ggplot2::ggplot_build(
    ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, fill = g)) +
      ggplot2::geom_col() +
      scale_fill_SciData()
  )
  expect_equal(unique(built_fill$data[[1]]$fill), SciDataPalette(3, names = FALSE))

  built_colour <- ggplot2::ggplot_build(
    ggplot2::ggplot(df, ggplot2::aes(x = g, y = y, colour = g)) +
      ggplot2::geom_point() +
      scale_color_SciData()
  )
  expect_equal(unique(built_colour$data[[1]]$colour), SciDataPalette(3, names = FALSE))
})

test_that(".ClusterPalette draws clusters from the SciData palette", {
  values <- .ClusterPalette(factor(c("1", "2", "3")))

  expect_equal(unname(values), .SciDataColorValues(3))
  expect_named(values, c("1", "2", "3"))
})

test_that(".ClusterPalette holds the noise group in grey", {
  values <- .ClusterPalette(factor(c("1", "2", "3", "Noise")))

  expect_equal(unname(values[["Noise"]]), "grey60")
  # Noise must not consume a palette slot.
  expect_equal(unname(values[c("1", "2", "3")]), .SciDataColorValues(3))
})

test_that(".ClusterPalette still finds noise after PlotClusterMap relabels levels", {
  # PlotClusterMap rewrites levels to "Noise (n = 12, 10%)" before the scale
  # sees them. An equality test on "Noise" would silently hand noise a real
  # cluster colour and shift every other cluster by one.
  relabelled <- factor(c(
    "1 (n = 40, 33%)", "2 (n = 40, 33%)", "Noise (n = 12, 10%)"
  ))

  values <- .ClusterPalette(relabelled)

  expect_equal(unname(values[["Noise (n = 12, 10%)"]]), "grey60")
  expect_equal(
    unname(values[c("1 (n = 40, 33%)", "2 (n = 40, 33%)")]),
    .SciDataColorValues(2)
  )
})

test_that(".ClusterPalette survives more clusters than the palette holds", {
  values <- .ClusterPalette(factor(paste0("C", seq_len(40))))

  expect_length(values, 40)
  expect_equal(anyDuplicated(values), 0)
})

test_that("PlotCategoricalDistributions fills from the SciData palette", {
  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  plot_obj <- PlotCategoricalDistributions(
    labelled, variables = c("Diagnosis", "sex")
  )
  plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot

  fills <- unique(ggplot2::ggplot_build(plot_obj)$data[[1]]$fill)

  expect_true(all(fills %in% .SciDataColorValues(170)))
})

test_that("PlotCategoricalDistributions labels stay legible on every fill", {
  # Labels are drawn inside the bars. The palette opens on Navy at 1.4:1
  # against black, so a fixed label colour would be unreadable on the first
  # category - and nothing else in the suite would catch it.
  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  plot_obj <- PlotCategoricalDistributions(
    labelled, variables = c("Diagnosis", "sex", "Genotype")
  )
  plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot

  built <- ggplot2::ggplot_build(plot_obj)
  text_layer <- built$data[[2]]

  # Both label colours are in use, i.e. the choice is actually adapting.
  expect_setequal(unique(text_layer$colour), c("black", "white"))

  # The text layer carries the resolved fill it sits on, so every label can be
  # checked against what the contrast helper would pick for that background.
  expect_gt(nrow(text_layer), 0)
  expect_equal(
    text_layer$colour,
    unname(.SciDataContrastText(text_layer$fill))
  )
})

test_that("PlotCategoricalDistributions keeps missing levels grey", {
  df <- data.frame(
    grp = factor(c("A", "B", "A", NA, "B", NA)),
    other = factor(c("X", "Y", "X", "Y", NA, "X"))
  )

  plot_obj <- PlotCategoricalDistributions(df, variables = c("grp", "other"))
  plot_obj <- if (inherits(plot_obj, "ggplot")) plot_obj else plot_obj$plot

  fills <- unique(ggplot2::ggplot_build(plot_obj)$data[[1]]$fill)

  expect_true("grey70" %in% fills)
})

test_that("CreateZScorePlot survives more variable categories than the old palette held", {
  # The previous `classcolors` vector was exactly 20 colours and went straight
  # into scale_color_manual(), so 21+ variable categories hard-failed with
  # "Insufficient values in manual scale". This is that regression test.
  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  variables <- getNumVars(labelled)[1:25]
  categories <- paste0("Domain", seq_along(variables))

  plot_obj <- suppressWarnings(CreateZScorePlot(
    data = labelled,
    TargetVar = "Diagnosis",
    variables = variables,
    VariableCategories = categories
  ))

  built <- suppressWarnings(ggplot2::ggplot_build(plot_obj))
  mapped <- unique(unlist(lapply(built$data, function(layer) layer$colour)))
  mapped <- mapped[!is.na(mapped)]

  expect_true(all(.SciDataColorValues(25) %in% mapped))
})

test_that("plotPCA hands plotly exactly-sized colours so it cannot interpolate", {
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

  # plotly silently interpolates a `colors` vector shorter than the level
  # count, which would yield a ramp *through* the palette rather than the
  # palette itself. Exact equality is what proves that did not happen.
  expect_setequal(trace_colors, .SciDataColorValues(length(trace_colors)))
})

test_that("plotPCA leaves a continuous colour variable to plotly's own scale", {
  skip_if_not_installed("plotly")

  data(SampleData, envir = environment())
  data(SampleVariableTypes, envir = environment())
  labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

  pca <- suppressWarnings(CreatePCAObject(
    data = labelled,
    VarsToReduce = c("AXL", "Adiponectin", "Ferritin", "MMP7", "tau", "p_tau")
  ))

  plot_obj <- suppressWarnings(plotPCA(
    pca, Var = "age", Mode = "2D", Components = c("RC1", "RC2")
  ))

  expect_no_error(suppressWarnings(plotly::plotly_build(plot_obj)))
})

test_that("every figure in every clustering pipeline actually builds", {
  # A palette scale can be constructed from the wrong column and still look
  # fine until ggplot resolves it, at which point it fails with
  # "Insufficient values in manual scale". Building each figure is the only
  # way to catch that, and it has to cover every pipeline: this exact bug
  # reached a vignette because only the K-means pipeline had been exercised.
  data(SimulatedPhenotypeData, envir = environment())

  training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
  numeric_vars <- paste0("Var", 1:12)
  categorical_vars <- paste0("CatVar", 1:3)

  BuildEveryPlot <- function(x, path = "") {
    if (inherits(x, "ggplot")) {
      err <- tryCatch(
        {
          invisible(ggplot2::ggplot_build(x))
          NULL
        },
        error = function(e) paste0(path, ": ", conditionMessage(e))
      )
      return(err)
    }
    if (is.list(x)) {
      return(unlist(lapply(
        names(x), function(n) BuildEveryPlot(x[[n]], paste0(path, "$", n))
      )))
    }
    NULL
  }

  models <- list(
    Mclust = function() {
      CreateClusterModel_MClust(
        training, numeric_vars, method = "finalize",
        final_k = 4, final_model = 1)
    },
    KMeans = function() {
      CreateClusterModel_KMeans(
        training, numeric_vars, method = "finalize", final_k = 4, nstart = 10)
    },
    HDBSCAN = function() {
      CreateClusterModel_HDBSCAN(
        training, c("DensityX", "DensityY"), method = "finalize",
        final_minPts = 10, final_cluster_selection_epsilon = 0)
    },
    LatentClass = function() {
      CreateClusterModel_LatentClass(
        training, categorical_vars, method = "finalize", final_k = 3)
    }
  )

  for (name in names(models)) {
    model <- suppressWarnings(models[[name]]())
    failures <- BuildEveryPlot(model)
    expect_null(failures, info = paste(name, "->", paste(failures, collapse = "; ")))
  }
})

test_that("the cluster scales reject an empty cluster vector with a clear message", {
  # NULL is what you get from `df$ColumnThatDoesNotExist`, which is how the
  # uncertainty-map bug happened.
  expect_error(.ClusterFillScale(NULL), "does not exist")
  expect_error(.ClusterColourScale(NULL), "does not exist")
  expect_error(.ClusterColourScale(factor(character(0))), "does not exist")
})
