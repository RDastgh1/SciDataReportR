test_that("DeriveFreesurferVolumes preserves automatic ICV volume derivations", {
  df_volumes <- data.frame(
    EstimatedTotalIntraCranialVol = c(1000, 2000),
    eTIV = c(1000, 2000),
    BrainSegVol = c(800, 1600),
    "Left-Hippocampus" = c(10, 20),
    "Right-Hippocampus" = c(20, 40),
    lh_fusiform_volume = c(4, 8),
    rh_fusiform_volume = c(6, 12),
    check.names = FALSE
  )

  actual <- DeriveFreesurferVolumes(df_volumes, verbose = FALSE)

  expect_equal(actual$Hippocampus_total, c(30, 60))
  expect_equal(actual$fusiform_volume_total, c(10, 20))
  expect_equal(actual$BrainSegVol_icv, c(0.8, 0.8))
  expect_equal(actual$Left_Hippocampus_icv, c(0.01, 0.01))
  expect_equal(actual$Hippocampus_total_icv, c(0.03, 0.03))
  expect_identical(
    attr(actual, "Freesurfer_derivation_log")$icv_source,
    "EstimatedTotalIntraCranialVol"
  )
  expect_identical(
    attr(actual, "Freesurfer_derivation_log")$bilateral_method,
    "sum"
  )
})

test_that("DeriveFreesurferVolumes accepts and validates an explicit ICV variable", {
  df_volumes <- data.frame(
    EstimatedTotalIntraCranialVol = c(1000, 2000),
    eTIV = c(1200, 2400),
    CustomICV = c(500, 1000),
    "Left-Hippocampus" = c(10, 20),
    "Right-Hippocampus" = c(20, 40),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  actual <- DeriveFreesurferVolumes(
    df_volumes,
    icv_var = "CustomICV",
    verbose = FALSE
  )

  expect_equal(actual$Hippocampus_total_icv, c(0.06, 0.06))
  expect_identical(
    attr(actual, "Freesurfer_derivation_log")$icv_source,
    "CustomICV"
  )
  expect_error(
    DeriveFreesurferVolumes(df_volumes, icv_var = "MissingICV"),
    "not found"
  )

  df_volumes$CharacterICV <- c("500", "1000")
  expect_error(
    DeriveFreesurferVolumes(df_volumes, icv_var = "CharacterICV"),
    "must be numeric"
  )
})

test_that("DeriveFreesurferVolumes derives cortical thickness means without ICV", {
  df_thickness <- data.frame(
    SubjectID = c("A", "B", "C"),
    lh_fusiform_thickness = c(2, NA, 4),
    rh_fusiform_thickness = c(4, 6, 8),
    lh_insula_thickness = c(3, 4, 5),
    rh_insula_thickness = c(5, 6, 7),
    lh_unmatched_thickness = c(1, 2, 3),
    stringsAsFactors = FALSE
  )

  actual <- DeriveFreesurferVolumes(
    df_thickness,
    derive_icv_ratios = FALSE,
    bilateral_method = "mean",
    verbose = FALSE
  )

  expect_named(
    actual,
    c("fusiform_thickness_total", "insula_thickness_total")
  )
  expect_equal(actual$fusiform_thickness_total, c(3, NA, 6))
  expect_equal(actual$insula_thickness_total, c(4, 5, 6))
  expect_false(any(grepl("_icv$", names(actual))))
  expect_identical(
    attr(actual, "Freesurfer_derivation_log")$derive_icv_ratios,
    FALSE
  )
  expect_true(
    "lh_unmatched_thickness" %in%
      attr(actual, "Freesurfer_derivation_log")$unmatched_dkt_left
  )
})

test_that("DeriveFreesurferVolumes skips ICV checks when ratios are disabled", {
  df_thickness <- data.frame(
    lh_fusiform_thickness = c(2, 3),
    rh_fusiform_thickness = c(4, 5)
  )

  expect_no_error(
    actual <- DeriveFreesurferVolumes(
      df_thickness,
      icv_var = "NotPresent",
      derive_icv_ratios = FALSE,
      bilateral_method = "mean",
      verbose = FALSE
    )
  )
  expect_equal(actual$fusiform_thickness_total, c(3, 4))
})
