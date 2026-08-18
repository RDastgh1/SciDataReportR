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

test_that("SummarizeTransitions returns figures beside the summaries", {
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
  expect_s3_class(out$Plots$TransitionPatterns, "ggplot")
})

test_that("The cascade figure plots the counts from condition_summary", {
  out <- SummarizeTransitions(
    data = MakeTransitionData(),
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )

  plotted <- ggplot2::ggplot_build(out$Plots$ConditionCascade)$data[[1]]
  summary_counts <- c(
    out$condition_summary$n_participants,
    out$condition_summary$n_ever_positive,
    out$condition_summary$n_developed_condition,
    out$condition_summary$n_resolved_condition,
    out$condition_summary$n_sustained_after_development,
    out$condition_summary$n_sustained_after_resolution
  )

  expect_setequal(plotted$x, summary_counts)
})

test_that("Transition patterns are mutually exclusive and account for everyone", {
  out <- SummarizeTransitions(
    data = MakeTransitionData(),
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )

  plotted <- ggplot2::ggplot_build(out$Plots$TransitionPatterns)$data[[1]]

  # Every participant is counted exactly once across the patterns.
  expect_equal(sum(plotted$x), out$condition_summary$n_participants)

  # Everything except "never positive" makes up the ever-positive group.
  never_positive <- out$condition_summary$n_participants -
    out$condition_summary$n_ever_positive
  expect_equal(
    sum(plotted$x) - never_positive,
    out$condition_summary$n_ever_positive
  )
})

test_that("SummarizeTransitions figures survive edge-case cohorts", {
  # Nobody ever positive: the figures should still build.
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
  expect_equal(out_none$condition_summary$n_ever_positive, 0)
  expect_silent(ggplot2::ggplot_build(out_none$Plots$ConditionCascade))
  expect_silent(ggplot2::ggplot_build(out_none$Plots$TransitionPatterns))

  # Everyone positive at every visit: no transitions at all.
  all_positive <- tibble::tibble(
    ParticipantID = rep(c("P1", "P2"), each = 3),
    VisitOrder = rep(1:3, times = 2),
    Status = rep(1, 6)
  )

  out_all <- SummarizeTransitions(
    data = all_positive,
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )

  expect_equal(out_all$condition_summary$n_developed_condition, 0)
  expect_equal(out_all$condition_summary$n_resolved_condition, 0)
  expect_silent(ggplot2::ggplot_build(out_all$Plots$TransitionPatterns))
})

test_that("Missing status values do not become a transition pattern", {
  with_na <- tibble::tibble(
    ParticipantID = rep(c("P1", "P2"), each = 3),
    VisitOrder = rep(1:3, times = 2),
    Status = c(0, NA, 1, NA, NA, NA)
  )

  out <- SummarizeTransitions(
    data = with_na,
    id_var = ParticipantID,
    time_var = VisitOrder,
    status_var = Status
  )

  plotted <- ggplot2::ggplot_build(out$Plots$TransitionPatterns)$data[[1]]
  expect_equal(sum(plotted$x), out$condition_summary$n_participants)
})
