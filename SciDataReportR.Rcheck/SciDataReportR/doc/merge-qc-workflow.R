## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
has_survival <- requireNamespace("survival", quietly = TRUE)
has_widgets <- requireNamespace("reactable", quietly = TRUE) &&
  requireNamespace("htmltools", quietly = TRUE)
knitr::opts_chunk$set(eval = has_survival)

## ----setup, message = FALSE---------------------------------------------------
library(SciDataReportR)
library(survival)

## -----------------------------------------------------------------------------
pbc_baseline <- pbc[, c("id", "age", "sex", "trt", "stage")]
pbc_labs <- pbcseq[, c("id", "day", "bili", "albumin", "platelet")]

m_pbc <- safe_merge(
  pbc_baseline,
  pbc_labs,
  by = "id",
  name = "pbc baseline + follow-up labs"
)

## -----------------------------------------------------------------------------
m_pbc$log

## -----------------------------------------------------------------------------
m_pbc$summary

## ----results = "asis"---------------------------------------------------------
merge_detail(m_pbc, TopN = 5)

## -----------------------------------------------------------------------------
pbc_enroll <- pbc_baseline
pbc_enroll$enroll_day <- 0

m_closest <- safe_merge(
  pbc_enroll,
  pbc_labs,
  by = "id",
  name = "closest lab to enrollment",
  method = "closest_time",
  time_var_before = "enroll_day",
  time_var_add = "day"
)

m_closest$log

## -----------------------------------------------------------------------------
demographics <- data.frame(
  id = 1:6, # integer
  sex = c("F", "M", "F", "F", "M", "F"),
  age = c(54, 61, 47, 66, 58, 50)
)

device_data <- data.frame(
  id = c(1, 2, 2, 5, 7), # double, with id 2 duplicated
  sex = c("F", "F", "M", "M", "F"), # disagrees with demographics for id 2/5
  score = c(0.82, 0.75, 0.71, 0.64, 0.90)
)

m_synth <- safe_merge(
  demographics,
  device_data,
  by = "id",
  name = "demographics + device data"
)

m_synth$log

## -----------------------------------------------------------------------------
m_synth$summary

## ----results = "asis"---------------------------------------------------------
merge_detail(m_synth)

## ----eval = has_survival && has_widgets---------------------------------------
ExploreMergeValidation(
  m_synth$validation,
  Title = "Demographics + device data",
  Detail = "Compact"
)

## -----------------------------------------------------------------------------
merge_summary_table(
  list(m_pbc$log, m_closest$log, m_synth$log),
  flagged_only = TRUE
)

