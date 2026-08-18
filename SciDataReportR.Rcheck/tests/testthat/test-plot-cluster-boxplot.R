test_that("PlotClusterBoxplot preserves ordinary numeric plotting data", {
  df_Test <- data.frame(
    Cluster = rep(c("A", "B"), each = 4),
    MarkerOne = seq_len(8),
    MarkerTwo = seq(10, 80, by = 10)
  )

  p <- PlotClusterBoxplot(
    data = df_Test,
    ClusterVar = "Cluster",
    variables = c("MarkerOne", "MarkerTwo"),
    Relabel = FALSE
  )

  expect_equal(
    p$data$Score,
    c(rbind(df_Test$MarkerOne, df_Test$MarkerTwo))
  )
  expect_type(p$data$Score, "double")
  expect_equal(levels(p$data$Label), c("MarkerOne", "MarkerTwo"))
})

test_that("PlotClusterBoxplot handles labelled numeric variables", {
  skip_if_not_installed("haven")

  df_Test <- data.frame(
    Cluster = rep(c("A", "B"), each = 4),
    MarkerOne = haven::labelled(
      seq_len(8),
      labels = c(Low = 1, High = 8),
      label = "First marker label"
    ),
    MarkerTwo = haven::labelled(
      seq(10, 80, by = 10),
      labels = c(Low = 10, High = 80),
      label = "Second marker label"
    )
  )

  expect_no_error(
    p <- PlotClusterBoxplot(
      data = df_Test,
      ClusterVar = "Cluster",
      variables = c("MarkerOne", "MarkerTwo")
    )
  )

  expect_type(p$data$Score, "double")
  expect_equal(
    levels(p$data$Label),
    c("First marker label", "Second marker label")
  )
})

test_that("PlotClusterBoxplot handles sjlabelled numeric variables", {
  df_Test <- data.frame(
    Cluster = rep(c("A", "B"), each = 4),
    MarkerOne = sjlabelled::set_label(seq_len(8), "First marker label"),
    MarkerTwo = sjlabelled::set_label(
      seq(10, 80, by = 10),
      "Second marker label"
    )
  )

  expect_no_error(
    p <- PlotClusterBoxplot(
      data = df_Test,
      ClusterVar = "Cluster",
      variables = c("MarkerOne", "MarkerTwo")
    )
  )

  expect_type(p$data$Score, "double")
  expect_equal(
    levels(p$data$Label),
    c("First marker label", "Second marker label")
  )
})

test_that("PlotClusterBoxplot codebook labels override variable labels", {
  skip_if_not_installed("haven")

  df_Test <- data.frame(
    Cluster = rep(c("A", "B"), each = 4),
    MarkerOne = haven::labelled(
      seq_len(8),
      label = "Variable attribute label"
    ),
    MarkerTwo = haven::labelled(
      seq(10, 80, by = 10),
      label = "Another attribute label"
    )
  )
  df_Codebook <- data.frame(
    Variable = c("MarkerOne", "MarkerTwo"),
    Label = c("Codebook marker one", "Codebook marker two")
  )

  p <- PlotClusterBoxplot(
    data = df_Test,
    ClusterVar = "Cluster",
    variables = c("MarkerOne", "MarkerTwo"),
    codebook = df_Codebook,
    Relabel = TRUE
  )

  expect_equal(
    levels(p$data$Label),
    c("Codebook marker one", "Codebook marker two")
  )
})
