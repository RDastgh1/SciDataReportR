test_that("default method benchmark candidates match the planned comparison", {
  candidates <- DefaultClusterBenchmarkCandidates()
  expect_true(all(c("som_mclust", "mclust", "kmeans", "hdbscan", "ward_rf") %in%
    candidates$Method))
  expect_setequal(sort(unique(candidates$Classes[candidates$Method == "kmeans"])), 3:5)
  expect_setequal(sort(unique(candidates$NodeWeighting[candidates$Method == "som_mclust"])),
    c("occupancy", "uniform"))
  expect_setequal(sort(unique(candidates$MinPts[candidates$Method == "hdbscan"])),
    c(20L, 50L, 100L))
})

test_that("common benchmark projects every complete participant from shared draws", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("tidyLPA")
  data("SampleData", package = "SciDataReportR")
  vars_Clustering <- c("age", "AXL", "Adiponectin")
  candidates <- dplyr::tibble(
    CandidateID = c("mclust_2", "kmeans_2"),
    Method = c("mclust", "kmeans"), Classes = c(2L, 2L),
    Model = c(1L, NA_integer_)
  )
  benchmark <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = vars_Clustering, candidates = candidates,
    stability_resamples = 2L, stability_seed = 6001L, stability_cores = 1L
  ))
  n_complete <- sum(stats::complete.cases(SampleData[, vars_Clustering, drop = FALSE]))
  expect_true(all(benchmark$replicates$ParticipantsCompared == n_complete))
  expect_true(all(benchmark$subsamples %>% dplyr::count(Replicate) %>% dplyr::pull(n) ==
    floor(n_complete * .90)))
  expect_identical(benchmark$settings$scoring_cohort, "all_reference_participants")
  expect_identical(benchmark$settings$perturbation, "90% subsample without replacement")

  repeated <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = vars_Clustering, candidates = candidates,
    stability_resamples = 2L, stability_seed = 6001L, stability_cores = 1L
  ))
  expect_identical(benchmark$subsamples, repeated$subsamples)
  expect_identical(benchmark$replicates, repeated$replicates)
})

test_that("ARI ignores cluster label permutations", {
  reference <- c(1L, 1L, 2L, 2L, 3L, 3L)
  permuted <- c(3L, 3L, 1L, 1L, 2L, 2L)
  expect_equal(mclust::adjustedRandIndex(reference, permuted), 1)
})

test_that("all-noise HDBSCAN candidates are not phenotype candidates", {
  skip_if_not_installed("dbscan")
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "all_noise", Method = "hdbscan",
    MinPts = 300L, Epsilon = 0)
  benchmark <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 1L, stability_seed = 6002L, stability_cores = 1L
  ))
  expect_identical(benchmark$decision_table$Eligibility, "not_applicable")
  expect_true(is.na(benchmark$decision_table$StabilityRank))
  expect_true(all(benchmark$replicates$Status == "not_applicable"))
})

test_that("Ward plus RF records its surrogate projection provenance", {
  skip_if_not_installed("ranger")
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "ward_2", Method = "ward_rf",
    Classes = 2L, WardN = 50L, RF_Trees = 20L)
  benchmark <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 1L, stability_seed = 6004L, stability_cores = 1L
  ))
  native <- benchmark$decision_table$NativeMetrics[[1]]
  expect_true(all(c("WardN", "WardSeed", "RFSeed", "ProjectionPolicy") %in% names(native)))
  expect_identical(native$ProjectionPolicy,
    "Ward.D2 labels projected by random-forest surrogate")
})

test_that("benchmark confirmation retains exactly selected settings", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("tidyLPA")
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "mclust_2", Method = "mclust",
    Classes = 2L, Model = 1L)
  benchmark <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 1L, stability_seed = 6003L, stability_cores = 1L
  ))
  confirmed <- suppressWarnings(ConfirmClusterMethodStability(benchmark,
    candidates = "mclust_2", confirm_resamples = 2L))
  expect_identical(confirmed$decision_table$CandidateID, "mclust_2")
  expect_identical(confirmed$decision_table$Classes, 2L)
  expect_identical(confirmed$decision_table$Model, 1L)
  expect_identical(confirmed$settings$refits_requested, 2L)
  expect_true(isTRUE(confirmed$settings$keep_models))
  expect_identical(names(confirmed$models), "mclust_2")
})

test_that("benchmark defaults SOM candidates to standard soft-3 projection", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("tidyLPA")
  skip_if_not_installed("aweSOM")
  skip_if_not_installed("kohonen")
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "som_standard", Method = "som_mclust",
    Classes = 2L, Model = 1L, NodeWeighting = "uniform", SOM_XDim = 3L,
    SOM_YDim = 3L, SOM_Topology = "hexagonal")
  benchmark <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 1L, stability_seed = 6006L, stability_cores = 1L,
    keep_models = TRUE
  ))
  expect_identical(benchmark$models$som_standard$Specification$projection$assignment,
    "soft")
  expect_identical(benchmark$models$som_standard$Specification$projection$soft_n_nodes,
    3L)
})

test_that("benchmark screen drops full models but retains worker audit", {
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "km", Method = "kmeans", Classes = 2L)
  benchmark <- CompareClusterMethodStability(SampleData,
    variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 1L, stability_seed = 6010L, stability_cores = 1L)
  expect_length(benchmark$models, 0L)
  expect_true(all(c("ResolvedWorkers", "WorkerPolicy") %in% names(benchmark$decision_table)))
  expect_true(all(c("bytes_per_worker", "safe_memory_budget_bytes") %in%
    names(benchmark$settings$worker_estimates$km)))
})

test_that("benchmark stability is deterministic in serial and parallel", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("tidyLPA")
  skip_if(future::availableCores() < 2L,
    "This test environment permits only one parallel worker.")
  skip_if(isTRUE(pkgload::is_dev_package("SciDataReportR")),
    "Multisession workers load the installed package; this assertion runs after package installation.")
  data("SampleData", package = "SciDataReportR")
  candidates <- dplyr::tibble(CandidateID = "mclust_2", Method = "mclust",
    Classes = 2L, Model = 1L)
  serial <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 2L, stability_seed = 6005L, stability_cores = 1L
  ))
  parallel <- suppressWarnings(CompareClusterMethodStability(
    SampleData, variables = c("age", "AXL", "Adiponectin"), candidates = candidates,
    stability_resamples = 2L, stability_seed = 6005L, stability_cores = 2L
  ))
  expect_identical(serial$subsamples, parallel$subsamples)
  expect_identical(serial$replicates, parallel$replicates)
})
