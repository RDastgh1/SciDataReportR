test_that("BuildFeatureWGCNA validates selected wide features", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  set.seed(1)
  df_data <- data.frame(SampleID = paste0("S", seq_len(30)), matrix(rnorm(30 * 9), ncol = 9), check.names = FALSE)
  names(df_data)[-1] <- paste0("Protein", seq_len(9))
  obj <- BuildFeatureWGCNA(df_data, variables = names(df_data)[-1], sample_id = "SampleID", soft_power = 4, min_module_size = 3, deep_split = 3, seed = 1)
  expect_s3_class(obj, "feature_wgcna_obj")
  expect_equal(nrow(obj$modules), 9)
  expect_equal(attr(obj$datExpr, "sample_ids"), df_data$SampleID)
  expect_equal(obj$params$deep_split, 3)
  expect_error(BuildFeatureWGCNA(df_data, c("Protein1", "Protein1", "Protein2")), "unique")
  expect_error(BuildFeatureWGCNA(df_data, names(df_data)[-1], soft_power = 4, min_module_size = 3, deep_split = 5), "deep_split")
  df_data$Protein1 <- 1
  expect_error(BuildFeatureWGCNA(df_data, names(df_data)[-1], soft_power = 4, min_module_size = 3), "constant")
})

test_that("AnalyzeFeatureModuleTraits requires an exact sample-ID join", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  set.seed(2)
  df_data <- data.frame(SampleID = paste0("S", seq_len(30)), matrix(rnorm(30 * 9), ncol = 9), check.names = FALSE)
  names(df_data)[-1] <- paste0("Protein", seq_len(9))
  obj <- BuildFeatureWGCNA(df_data, names(df_data)[-1], sample_id = "SampleID", soft_power = 4, min_module_size = 3)
  df_traits <- data.frame(SampleID = rev(df_data$SampleID), Age = seq_len(30), Group = rep(c("A", "B"), 15))
  result <- AnalyzeFeatureModuleTraits(obj, df_traits, sample_id = "SampleID", outcome = c("Age", "Group"))
  expect_s3_class(result, "feature_module_trait_obj")
  expect_true(all(result$join_audit$Matched))
  df_traits$SampleID[1] <- "other"
  expect_error(AnalyzeFeatureModuleTraits(obj, df_traits, sample_id = "SampleID", outcome = "Age"), "incomplete")
})

CreateWGCNAFixture <- function() {
  set.seed(104)
  signals <- matrix(rnorm(60 * 4), ncol = 4)
  signals[, 2] <- 0.8 * signals[, 1] + 0.6 * signals[, 2]
  signals[, 4] <- -signals[, 3] + 0.3 * signals[, 4]
  dynamic <- c(rep(1:4, c(6, 5, 4, 3)), 0L, 0L)
  df_data <- as.data.frame(vapply(seq_along(dynamic), function(i) {
    if (dynamic[i] == 0L) rnorm(60) else signals[, dynamic[i]] + rnorm(60, sd = 0.1)
  }, numeric(60)))
  names(df_data) <- paste0("Protein", seq_along(dynamic))
  list(data = df_data, dynamic = dynamic)
}

FitWGCNAFixture <- function(fixture, ...) {
  BuildFeatureWGCNA(fixture$data, names(fixture$data), soft_power = 4,
                    min_module_size = 3, deep_split = 3, seed = 104, ...)
}

test_that("NULL preserves height-based merging and records non-grey counts", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  implicit <- FitWGCNAFixture(fixture, keep_tom = TRUE)
  explicit <- FitWGCNAFixture(fixture, keep_tom = TRUE, n_modules = NULL)
  expect_identical(implicit, explicit)
  dynamic <- dynamicTreeCut::cutreeDynamic(implicit$gene_tree, distM = 1 - implicit$tom,
    deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 3, verbose = 0)
  legacy <- WGCNA::mergeCloseModules(implicit$datExpr, WGCNA::labels2colors(dynamic),
    cutHeight = 0.25, verbose = 0)
  expect_identical(implicit$module_colors, legacy$colors)
  expect_null(implicit$params$n_modules)
  expect_equal(implicit$n_modules, length(setdiff(unique(legacy$colors), "grey")))
})

test_that("exact counts use average-linkage eigengenes and preserve grey", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  testthat::local_mocked_bindings(
    cutreeDynamic = function(dendro, distM, deepSplit, pamRespectsDendro, minClusterSize) {
      expect_equal(deepSplit, 3)
      expect_equal(minClusterSize, 3)
      fixture$dynamic
    }, .package = "dynamicTreeCut")
  testthat::local_mocked_bindings(
    mergeCloseModules = function(...) stop("Height-based merging must not run"), .package = "WGCNA")
  original_colors <- WGCNA::labels2colors(fixture$dynamic)
  initial_mes <- WGCNA::moduleEigengenes(fixture$data, original_colors, excludeGrey = TRUE)$eigengenes
  expected_groups <- stats::cutree(stats::hclust(stats::as.dist(1 - stats::cor(initial_mes)),
    method = "average"), k = 2)
  expected_assignments <- unname(expected_groups[paste0("ME", original_colors[original_colors != "grey"])])
  for (k in c(1L, 2L, 4L)) {
    obj <- FitWGCNAFixture(fixture, n_modules = k)
    expect_equal(obj$params$n_modules, k)
    expect_equal(obj$n_modules, k)
    expect_equal(length(setdiff(unique(obj$module_colors), "grey")), k)
    expect_identical(obj$module_colors == "grey", original_colors == "grey")
    non_grey_sizes <- obj$module_sizes[WGCNA::labels2colors(seq_len(k))]
    expect_true(all(diff(as.numeric(non_grey_sizes)) <= 0))
    recomputed <- WGCNA::moduleEigengenes(obj$datExpr, obj$module_colors)$eigengenes
    recomputed <- recomputed[, order(names(recomputed)), drop = FALSE]
    expect_equal(obj$eigengenes, recomputed)
    correlations <- stats::cor(obj$datExpr, recomputed, use = "pairwise.complete.obs")
    expected_kme <- correlations[cbind(seq_along(obj$module_colors),
      match(paste0("ME", obj$module_colors), colnames(correlations)))]
    expect_equal(obj$modules$ModuleMembership, as.numeric(expected_kme))
    expect_identical(obj$modules$Module, obj$module_colors)
    if (k == 2L) {
      actual <- obj$module_colors[original_colors != "grey"]
      expect_identical(outer(actual, actual, "=="), outer(expected_assignments, expected_assignments, "=="))
      expect_equal(obj$module_colors[1], "turquoise")
      expect_equal(obj$module_colors[7], "turquoise")
      expect_identical(obj$module_colors,
        FitWGCNAFixture(fixture, n_modules = k, merge_cut_height = 0.9)$module_colors)
    }
  }
})

test_that("exact-count colors resolve size ties by original color", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  fixture$dynamic <- c(rep(1:3, each = 6), 0L, 0L)
  testthat::local_mocked_bindings(cutreeDynamic = function(...) fixture$dynamic, .package = "dynamicTreeCut")
  obj <- FitWGCNAFixture(fixture, n_modules = 3)
  expect_equal(obj$module_colors[c(1, 7, 13)], c("brown", "turquoise", "blue"))
  expect_identical(obj$module_colors, FitWGCNAFixture(fixture, n_modules = 3)$module_colors)
})

test_that("exact counts handle one module and reject insufficient modules", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  testthat::local_mocked_bindings(cutreeDynamic = function(...) fixture$dynamic, .package = "dynamicTreeCut")
  expect_error(FitWGCNAFixture(fixture, n_modules = 5),
    "4 non-grey modules.*n_modules = 5.*Raise `deep_split` or lower `min_module_size`")
  fixture$dynamic <- c(rep(2L, 18), 0L, 0L)
  obj <- FitWGCNAFixture(fixture, n_modules = 1)
  expect_identical(obj$module_colors, c(rep("turquoise", 18), "grey", "grey"))
  expect_equal(obj$n_modules, 1)
  fixture$dynamic[] <- 0L
  expect_error(FitWGCNAFixture(fixture, n_modules = 1),
    "0 non-grey modules.*n_modules = 1.*Raise `deep_split` or lower `min_module_size`")
})

test_that("n_modules requires a finite positive scalar whole number", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  invalid <- list(0, -1, 1.5, NA_real_, NaN, Inf, -Inf, TRUE, "2", c(1, 2), numeric(0), 2 + 1i)
  for (value in invalid) {
    expect_error(FitWGCNAFixture(fixture, n_modules = value), "`n_modules` must be NULL or a single finite positive whole number", fixed = TRUE)
  }
})

test_that("exact counts work with the real dynamic tree cut", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  fixture <- CreateWGCNAFixture()
  obj <- FitWGCNAFixture(fixture, n_modules = 2, keep_tom = TRUE)
  dynamic <- dynamicTreeCut::cutreeDynamic(obj$gene_tree, distM = 1 - obj$tom,
    deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 3, verbose = 0)
  expect_gte(length(setdiff(unique(dynamic), 0L)), 2)
  expect_equal(obj$n_modules, 2)
  expect_identical(obj$module_colors == "grey", unname(dynamic == 0L))
})
