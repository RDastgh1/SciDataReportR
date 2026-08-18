# Extracted from test-summarize-transitions.R:30

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
out <- SummarizeTransitions(
    data = MakeTransitionData(),
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )
expect_named(
    out,
    c("participant_summary", "condition_summary", "Plots")
  )
expect_named(out$Plots, c("ConditionCascade", "TransitionPatterns"))
expect_s3_class(out$Plots$ConditionCascade, "ggplot")
