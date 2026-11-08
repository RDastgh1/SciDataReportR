# A hand-built screen so the selection rule can be tested without fitting maps
FakeSOMGridScreen <- function(n, p, sides, qe = NULL, one_person = NULL, empty = NULL,
    te = NULL, trust = NULL, jitter_sd = 0.25) {
  k <- length(sides)
  x <- list(
    table = data.frame(Side = sides, Nodes = sides^2, PeoplePerNode = n / sides^2,
      QuantizationError = if (is.null(qe)) rep(5, k) else qe,
      OnePersonSmooth = if (is.null(one_person)) rep(0, k) else one_person,
      EmptySmooth = if (is.null(empty)) rep(0, k) else empty,
      TopographicError = if (is.null(te)) rep(0.3, k) else te,
      Trustworthiness = if (is.null(trust)) rep(0.95, k) else trust),
    settings = list(n = n, p = p, jitter_sd = rep_len(jitter_sd, p), refits = 2L))
  class(x) <- "SOMGridScreen"
  x
}

test_that("the 5 * sqrt(n) target is chosen when it lies inside the allowed range", {
  # n = 405: target round(sqrt(5 * sqrt(405))) = 10; p = 5: floor 8, cap 12
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12))
  expect_identical(x$selection$target, 10L)
  expect_identical(x$selection$floor, 8L)
  expect_identical(x$selection$side, 10L)
  expect_identical(x$selection$ceiling_binding, "cap")
  expect_true(x$table$Chosen[x$table$Side == 10])
})

test_that("the sparsity ceiling moves the target down and is reported", {
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12,
    one_person = c(0, 0, 0, .02, .05, .08, .15)))
  expect_identical(x$selection$ceiling, 11L)
  expect_identical(x$selection$side, 10L)
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12,
    one_person = c(0, 0, 0, .02, .12, .2, .3)))
  expect_identical(x$selection$ceiling_binding, "sparsity")
  expect_identical(x$selection$side, 9L)
  expect_match(x$selection$reason, "moved down to the sparsity ceiling")
  # Empty nodes have their own cap
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12, empty = c(0, 0, 0, 0, .06, .1, .2)))
  expect_identical(x$selection$side, 9L)
})

test_that("the noise ceiling stops at the first map finer than measurement noise", {
  # noise distance = 0.25 * sqrt(5) = 0.559
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12,
    qe = c(.83, .76, .69, .63, .58, .53, .49)))
  expect_equal(x$settings$noise_distance, 0.25 * sqrt(5))
  expect_identical(x$selection$ceilings[["noise"]], 11)
  expect_identical(x$selection$side, 10L)
  # Larger noise lowers the ceiling below the target
  x <- SelectSOMGrid(x, noise_sd = 0.3)
  expect_identical(x$selection$ceiling_binding, "noise")
  expect_identical(x$selection$side, 9L)
})

test_that("large samples are capped at 30 nodes per variable, or max_nodes", {
  # NHANES-like: n = 18196, p = 15 -> target 26, cap floor(sqrt(450)) = 21
  x <- SelectSOMGrid(FakeSOMGridScreen(18196, 15, 11:22))
  expect_identical(x$selection$target, 26L)
  expect_identical(x$selection$side, 21L)
  expect_identical(x$selection$ceiling_binding, "cap")
  expect_match(x$selection$ceiling_label, "30-nodes-per-variable")
  x <- SelectSOMGrid(x, max_nodes = 300)
  expect_identical(x$selection$side, 17L)
  expect_match(x$selection$ceiling_label, "300-node cap")
})

test_that("small samples move up to the floor, and a floor above the ceiling warns", {
  # n = 120, p = 10: target round(sqrt(5 * sqrt(120))) = 7, floor 10
  x <- SelectSOMGrid(FakeSOMGridScreen(120, 10, 5:12))
  expect_identical(x$selection$side, 10L)
  expect_match(x$selection$reason, "moved up to the variable floor")
  expect_warning(x <- SelectSOMGrid(FakeSOMGridScreen(120, 10, 5:12,
    one_person = c(0, 0, 0, .2, .3, .4, .5, .6))), "variable floor")
  expect_identical(x$selection$ceiling, 7L)
  expect_identical(x$selection$side, 7L)
  expect_identical(x$selection$allowed, c(7L, 10L))
})

test_that("largest_feasible takes the top of the allowed range", {
  x <- SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12, one_person = c(0, 0, 0, .02, .05, .08, .15)),
    select = "largest_feasible")
  expect_identical(x$selection$side, 11L)
})

test_that("map checks warn for a folded map, low trustworthiness, or a too-small screen", {
  expect_warning(SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12, te = c(.3, .3, .3, .3, .6, .3, .3))),
    "folded")
  expect_warning(SelectSOMGrid(FakeSOMGridScreen(405, 5, 6:12, trust = rep(.7, 7))),
    "Trustworthiness")
  expect_warning(SelectSOMGrid(FakeSOMGridScreen(18196, 15, 11:15)), "largest side screened")
})

test_that("selection settings are validated", {
  x <- FakeSOMGridScreen(405, 5, 6:12)
  expect_error(SelectSOMGrid(list()), "SOMGridScreen")
  expect_error(SelectSOMGrid(x, min_nodes_per_variable = 0), "min_nodes_per_variable")
  expect_error(SelectSOMGrid(x, max_nodes_per_variable = 5), "at least min_nodes_per_variable")
  expect_error(SelectSOMGrid(x, max_nodes = 2), "max_nodes")
})

test_that("candidate k with too few nodes per cluster give one warning", {
  expect_warning(low <- SciDataReportR:::.WarnSOMNodesPerCluster(100, c(2, 5, 11, 12), 10),
    "k = 11, 12")
  expect_identical(low, c(11, 12))
  expect_silent(SciDataReportR:::.WarnSOMNodesPerCluster(100, 2:10, 10))
  expect_silent(SciDataReportR:::.WarnSOMNodesPerCluster(100, 2:20, NULL))
})

test_that("ScreenSOMGrid measures real maps and trains the pipeline's map", {
  skip_if_not_installed("kohonen")
  set.seed(4411)
  df <- EmbeddedTestData()
  names(df) <- c("A", "B", "C")
  screen <- suppressWarnings(ScreenSOMGrid(df, c("A", "B", "C"), sides = 4:7, refits = 2))
  expect_s3_class(screen, "SOMGridScreen")
  expect_identical(screen$table$Side, 4:7)
  expect_true(all(c("QuantizationError", "VarianceCaptured", "TopographicError", "Trustworthiness",
    "Continuity", "EmptyNodes", "OnePersonNodes", "LowDensityNodes", "KappaSeed", "KappaJitter",
    "OnePersonSmooth", "EmptySmooth", "PassFloor", "PassSparsity", "FinerThanNoise", "PassCap",
    "Chosen") %in% names(screen$table)))
  expect_true(all(diff(screen$table$QuantizationError) < 0))
  expect_equal(nrow(screen$refits), 4L * 2L)
  expect_equal(nrow(screen$pairs), 4L * 1L)
  expect_identical(sum(screen$table$Chosen), 1L)
  expect_s3_class(plot(screen), "ggplot")
  expect_equal(nrow(SOMGridDecisionTable(screen)), 4L)
  expect_output(print(screen), "Chosen")

  # Deterministic for a given seed
  again <- suppressWarnings(ScreenSOMGrid(df, c("A", "B", "C"), sides = 4:7, refits = 2))
  expect_equal(again$table, screen$table)

  # The first seed refit is exactly the map the pipeline trains with the same seed
  zmat <- as.matrix(CreateZScoreObject(df, c("A", "B", "C"), names_prefix = "Z_")$ZScores)
  m1 <- SciDataReportR:::.FitSOMGridScreen(zmat, 5L, 934521L, "hexagonal", "gaussian", 100L)
  q1 <- SciDataReportR:::.SOMGridMapQuality(m1, zmat, 3L)
  expect_equal(unname(q1$metrics[["QuantizationError"]]),
    screen$refits$QuantizationError[screen$refits$Side == 5 & screen$refits$Refit == 1])

  expect_error(ScreenSOMGrid(df, c("A", "B", "C"), jitter_sd = c(.1, .2)), "jitter_sd")
  expect_error(ScreenSOMGrid(df, "Z"), "variables")
})

test_that("the pipeline validates grid_screen_args", {
  df <- EmbeddedAB()
  expect_error(CreateClusterModel_SOM_MClust(df, c("A", "B"), k_range = 2, models = 1,
    grid_screen_args = list(seed_som = 1)), "grid_screen_args")
  expect_error(CreateClusterModel_SOM_MClust(df, c("A", "B"), k_range = 2, models = 1,
    grid_screen_args = list(1)), "grid_screen_args")
})
