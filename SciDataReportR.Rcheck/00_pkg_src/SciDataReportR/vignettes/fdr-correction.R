## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)


## ----setup, message = FALSE---------------------------------------------------
library(SciDataReportR)


## -----------------------------------------------------------------------------
pm <- matrix(
  c(0.005, 0.04, 0.03, 0.01, NA, 0.002),
  nrow = 2,
  dimnames = list(c("pred1", "pred2"), c("out1", "out2", "out3"))
)
pm


## -----------------------------------------------------------------------------
ApplyFDRCorrection(pm, fdr_scope = "matrix")


## -----------------------------------------------------------------------------
ApplyFDRCorrection(pm, fdr_scope = "per_outcome", outcome_margin = 2)


## -----------------------------------------------------------------------------
res_matrix <- PlotCorrelationsHeatmap(
  mtcars,
  predictor_vars = c("disp", "hp", "drat", "wt"),
  outcome_vars   = c("mpg", "qsec"),
  fdr_scope      = "matrix"
)

res_per_outcome <- PlotCorrelationsHeatmap(
  mtcars,
  predictor_vars = c("disp", "hp", "drat", "wt"),
  outcome_vars   = c("mpg", "qsec"),
  fdr_scope      = "per_outcome"
)


## -----------------------------------------------------------------------------
round(res_matrix$p$p, 5)          # res$p is an alias for res$Unadjusted
all.equal(res_matrix$p$p, res_per_outcome$p$p)


## -----------------------------------------------------------------------------
round(res_matrix$p_fdr$p, 5)       # one family of 8 tests
round(res_per_outcome$p_fdr$p, 5)  # two families of 4 tests each


## ----fig.show = "hold"--------------------------------------------------------
res_matrix$p_fdr$plot
res_per_outcome$p_fdr$plot

