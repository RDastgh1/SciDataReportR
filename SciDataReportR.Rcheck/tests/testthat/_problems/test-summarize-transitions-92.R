# Extracted from test-summarize-transitions.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "SciDataReportR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
MakeTransitionData <- function() {
  tibble::tibble(
    ParticipantID = rep(paste0("P", 1:5), each = 4),
    VisitOrder = rep(1:4, times = 5),
    # P1 develops, P2 resolves, P3 never positive, P4 positive throughout,
    # P5 develops then resolves.
    Status = c(
      0, 0, 1, 1,
      1, 1, 0, 0,
      0, 0, 0, 0,
      1, 1, 1, 1,
      0, 1, 1, 0
    )
  )
}

# test -------------------------------------------------------------------------
none_positive <- tibble::tibble(
    ParticipantID = rep(c("P1", "P2"), each = 3),
    VisitOrder = rep(1:3, times = 2),
    Status = rep(0, 6)
  )
out_none <- SummarizeTransitions(
    data = none_positive,
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )
expect_s3_class(out_none$Plots$ConditionCascade, "ggplot")
