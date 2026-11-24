# A CHARSXP cannot be built by ordinary R code, but it survives serialization,
# so one is made here by removing the one-element character-vector wrapper from
# a serialized string and unserializing the bare element that remains. This
# needs no compiler and no R version pinning: the wrapper's bytes are located in
# the stream this R wrote.
MakeCharsxp <- function(text = "Female") {
  bytes <- serialize(text, NULL)
  wrapper <- as.raw(c(0, 0, 0, 16, 0, 0, 0, 1))
  at <- which(vapply(seq_len(length(bytes) - 7L), function(i)
    identical(bytes[i:(i + 7L)], wrapper), logical(1)))[1]
  expect_false(is.na(at))
  unserialize(bytes[-(at:(at + 7L))])
}

test_that("the CHARSXP fixture is what the failing backend value was", {
  x <- MakeCharsxp()
  expect_identical(typeof(x), "char")
  # This is the exact abort the conditional backend's values produced.
  expect_error(inherits(x, "Date"), "CHARSXP")
  expect_error(as.character(x), "cannot coerce")
})

test_that("ScidrRecoverValues restores a CHARSXP wherever it is nested", {
  finding <- list(suspicous_value = list(column = "Target"),
    conditions = list(list(column = "Sex", comparison = "=", value_comp = MakeCharsxp("Female")),
      list(column = "Age", comparison = "<=", value_comp = 60)),
    group_statistics = list(n_obs = 40L, categ_maj = MakeCharsxp('a "quoted" \\ value')))
  recovered <- SciDataReportR:::ScidrRecoverValues(finding)
  expect_identical(recovered$conditions[[1]]$value_comp, "Female")
  # deparse() reads a CHARSXP verbatim, so quotes and backslashes survive unescaped.
  expect_identical(recovered$group_statistics$categ_maj, 'a "quoted" \\ value')
  # Everything else is left exactly as it was, structure and names included.
  expect_identical(recovered$conditions[[2]], finding$conditions[[2]])
  expect_identical(recovered$suspicous_value, finding$suspicous_value)
  expect_identical(names(recovered), names(finding))
  expect_identical(recovered$group_statistics$n_obs, 40L)
})

test_that("ScidrReviewFormat renders a CHARSXP instead of aborting", {
  expect_identical(SciDataReportR:::ScidrReviewFormat(MakeCharsxp("Female")), "Female")
})

test_that("a backend condition value R cannot touch still renders a queue and a page", {
  skip_if_not_installed("outliertree")
  set.seed(31)
  n <- 160
  df_Data <- data.frame(ID = sprintf("P%03d", seq_len(n)),
    Sex = rep(c("Female", "Male"), each = n / 2),
    Measure = c(stats::rnorm(n / 2, 10, 1), stats::rnorm(n / 2, 30, 1)))
  df_Data$Measure[1] <- 30
  result <- DetectDataAnomalies(df_Data, c("Sex", "Measure"), id_var = "ID", methods = "conditional")
  rows <- which(vapply(result$CellResults$Conditions, length, integer(1)) > 0L)
  skip_if(!length(rows), "the backend reported no conditions for this fixture")
  # Poison one stored condition the way the backend did, then confirm the queue
  # and the page are built rather than aborting on the whole run.
  result$CellResults$Conditions[[rows[1]]][[1]]$value_comp <- MakeCharsxp("Female")
  queue <- CreateAnomalyReviewQueue(result)
  poisoned <- queue$Evidence$Conditions[vapply(queue$Evidence$Conditions, length, integer(1)) > 0L]
  expect_true(any(vapply(poisoned, function(conditions)
    any(vapply(conditions, function(condition) identical(condition$value_comp, "Female"), logical(1))), logical(1))))
  expect_true(all(vapply(queue$Evidence$Conditions, function(conditions)
    all(vapply(conditions, function(condition) typeof(condition$value_comp) != "char", logical(1))), logical(1))))
  html <- as.character(ExploreDataAnomalies(result))
  expect_true(any(grepl("Female", html, fixed = TRUE)))
})
