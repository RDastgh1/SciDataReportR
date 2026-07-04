# Integration tests: every FDR-using heatmap/matrix function must
#   (a) reproduce historical behavior with fdr_scope = "matrix"
#       (checked against an independent stats::p.adjust() recomputation of
#        the returned raw p-values - the definition of the old behavior), and
#   (b) produce per-outcome-grouped, generally different results with
#       fdr_scope = "per_outcome".

make_fdr_data <- function(n = 70) {
  set.seed(99)
  d <- data.frame(
    a1 = rnorm(n), a2 = rnorm(n), a3 = rnorm(n),
    b1 = rnorm(n), b2 = rnorm(n),
    g1 = factor(sample(c("x", "y"), n, TRUE)),
    g2 = factor(sample(c("lo", "hi"), n, TRUE)),
    g3 = factor(sample(c("A", "B"), n, TRUE))
  )
  d$b1 <- d$b1 + 0.8 * d$a1                       # strong signal for one outcome
  d$a2[d$g1 == "x"] <- d$a2[d$g1 == "x"] + 0.9    # group effect
  d
}

manual_matrix_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  out[ok] <- stats::p.adjust(p[ok], method = "fdr")
  out
}

manual_grouped_adjust <- function(p, ids) {
  out <- rep(NA_real_, length(p))
  for (idx in split(seq_along(p), as.character(ids))) {
    out[idx] <- manual_matrix_adjust(p[idx])
  }
  out
}

test_that("PlotCorrelationsHeatmap fdr_scope works (outcomes = columns)", {
  d <- make_fdr_data()
  rm_ <- PlotCorrelationsHeatmap(d, predictor_vars = c("a1", "a2", "a3"),
                                 outcome_vars = c("b1", "b2"))
  # matrix scope == whole-matrix p.adjust of the returned raw p matrix
  expect_equal(as.vector(rm_$FDRCorrected$p),
               manual_matrix_adjust(as.vector(rm_$Unadjusted$p)))

  rp <- PlotCorrelationsHeatmap(d, predictor_vars = c("a1", "a2", "a3"),
                                outcome_vars = c("b1", "b2"),
                                fdr_scope = "per_outcome")
  expected_cols <- apply(rm_$Unadjusted$p, 2, manual_matrix_adjust)
  dimnames(expected_cols) <- dimnames(rm_$Unadjusted$p)
  expect_equal(rp$FDRCorrected$p, expected_cols)
  expect_false(isTRUE(all.equal(rp$FDRCorrected$p, rm_$FDRCorrected$p)))
  # standardized aliases present and pointing at the same objects
  expect_identical(rm_$p, rm_$Unadjusted)
  expect_identical(rm_$p_fdr, rm_$FDRCorrected)
})

test_that("PlotPointCorrelationsHeatmap fdr_scope works (outcomes = continuous vars)", {
  d <- make_fdr_data()
  rm_ <- PlotPointCorrelationsHeatmap(d, CatVars = c("g1", "g2", "g3"),
                                      ContVars = c("a1", "a2", "b1"))
  tab <- rm_$Unadjusted$PvalTable
  expect_equal(tab$p.adj, manual_matrix_adjust(tab$p_value))

  rp <- PlotPointCorrelationsHeatmap(d, CatVars = c("g1", "g2", "g3"),
                                     ContVars = c("a1", "a2", "b1"),
                                     fdr_scope = "per_outcome")
  tabp <- rp$Unadjusted$PvalTable
  expect_equal(tabp$p.adj,
               manual_grouped_adjust(tabp$p_value, tabp$ContinuousVariable))
  expect_false(isTRUE(all.equal(tabp$p.adj, tab$p.adj)))
  expect_identical(rm_$p, rm_$Unadjusted)
  expect_identical(rm_$p_fdr, rm_$FDRCorrected)
})

test_that("PlotPhiHeatmap fdr_scope works (outcomes = y-axis variables)", {
  d <- make_fdr_data()
  rm_ <- PlotPhiHeatmap(d, CatVars = c("g1", "g2", "g3"))
  tab <- rm_$Unadjusted$PvalTable
  expect_equal(tab$p.adj, manual_matrix_adjust(tab$p_value))

  rp <- PlotPhiHeatmap(d, CatVars = c("g1", "g2", "g3"), fdr_scope = "per_outcome")
  tabp <- rp$Unadjusted$PvalTable
  expect_equal(tabp$p.adj, manual_grouped_adjust(tabp$p_value, tabp$YVar))
  expect_identical(rm_$p, rm_$Unadjusted)
  expect_identical(rm_$p_fdr, rm_$FDRCorrected)
})

test_that("PlotChiSqCovar fdr_scope works (outcomes = y-axis variables)", {
  d <- make_fdr_data()
  rm_ <- PlotChiSqCovar(d, predictor_vars = c("g1", "g2"),
                        outcome_vars = c("g3", "g1"))
  tab <- rm_$p$data
  expect_equal(tab$pval.adj, manual_matrix_adjust(tab$pval))

  rp <- PlotChiSqCovar(d, predictor_vars = c("g1", "g2"),
                       outcome_vars = c("g3", "g1"), fdr_scope = "per_outcome")
  tabp <- rp$p$data
  expect_equal(tabp$pval.adj, manual_grouped_adjust(tabp$pval, tabp$YVar))
  # standardized aliases (flat family: p is the unadjusted plot)
  expect_identical(rm_$p_fdr, rm_$p_FDR)
  expect_identical(rm_$pvaltable_fdr, rm_$pvaltable_FDR)
})

test_that("PlotAnovaRelationshipsMatrix fdr_scope works (outcomes = ContVars)", {
  d <- make_fdr_data()
  rm_ <- PlotAnovaRelationshipsMatrix(d, CatVars = c("g1", "g2"),
                                      ContVars = c("a1", "a2", "b1"))
  tab <- rm_$Unadjusted$PvalTable
  expect_equal(tab$p.adj, manual_matrix_adjust(tab$p))

  rp <- PlotAnovaRelationshipsMatrix(d, CatVars = c("g1", "g2"),
                                     ContVars = c("a1", "a2", "b1"),
                                     fdr_scope = "per_outcome")
  tabp <- rp$Unadjusted$PvalTable
  expect_equal(tabp$p.adj, manual_grouped_adjust(tabp$p, tabp$ContinuousVariable))
  expect_false(isTRUE(all.equal(tabp$p.adj, tab$p.adj)))
  expect_identical(rm_$p, rm_$Unadjusted)
  expect_identical(rm_$p_fdr, rm_$FDRCorrected)
})

test_that("PlotMiningMatrix fdr_scope works (outcomes = x-axis / XVar)", {
  d <- make_fdr_data()
  rm_ <- PlotMiningMatrix(d, outcome_vars = c("a1", "a2", "g3"),
                          predictor_vars = c("b1", "b2", "g1"))
  tab <- rm_$Unadjusted$PvalTable
  expect_equal(tab$p_adj, manual_matrix_adjust(tab$p))

  rp <- PlotMiningMatrix(d, outcome_vars = c("a1", "a2", "g3"),
                         predictor_vars = c("b1", "b2", "g1"),
                         fdr_scope = "per_outcome")
  tabp <- rp$Unadjusted$PvalTable
  expect_equal(tabp$p_adj, manual_grouped_adjust(tabp$p, tabp$XVar))
  expect_false(isTRUE(all.equal(tabp$p_adj, tab$p_adj)))
  expect_identical(rm_$p, rm_$Unadjusted)
})

test_that("PlotInteractionEffectsMatrix fdr_scope works (outcomes = Outcome)", {
  d <- make_fdr_data()
  rm_ <- PlotInteractionEffectsMatrix(d, interVar = "g1",
                                      outcome_vars = c("a1", "a2"),
                                      predictor_vars = c("b1", "b2"))
  tab <- rm_$raw_data
  expect_equal(tab$P_FDR, manual_matrix_adjust(tab$P))

  rp <- PlotInteractionEffectsMatrix(d, interVar = "g1",
                                     outcome_vars = c("a1", "a2"),
                                     predictor_vars = c("b1", "b2"),
                                     fdr_scope = "per_outcome")
  tabp <- rp$raw_data
  expect_equal(tabp$P_FDR, manual_grouped_adjust(tabp$P, tabp$Outcome))
  expect_false(isTRUE(all.equal(tabp$P_FDR, tab$P_FDR)))
  expect_identical(rm_$p, rm_$Unadjusted)
  expect_identical(rm_$p_fdr, rm_$FDRCorrected)
})

test_that("PlotCatInteractionEffectsMatrix fdr_scope works (outcomes = Y)", {
  d <- make_fdr_data()
  rm_ <- PlotCatInteractionEffectsMatrix(d, predictor_vars = c("g2", "g3"),
                                         outcome_vars = c("g3", "g2"),
                                         interVar = "g1")
  tab <- rm_$p$data
  expect_equal(tab$P_FDR, manual_matrix_adjust(tab$P))

  rp <- PlotCatInteractionEffectsMatrix(d, predictor_vars = c("g2", "g3"),
                                        outcome_vars = c("g3", "g2"),
                                        interVar = "g1", fdr_scope = "per_outcome")
  tabp <- rp$p$data
  expect_equal(tabp$P_FDR, manual_grouped_adjust(tabp$P, tabp$Y))
  expect_identical(rm_$p_fdr, rm_$p_FDR)
  expect_identical(rm_$pvals_fdr, rm_$pvals_FDR)
})

test_that("PlotNumInteractionEffectsMatrix fdr_scope works (outcomes = Y)", {
  d <- make_fdr_data()
  # Covariates are optional: without them the plain interaction model is
  # fit (silent all-NA behavior was fixed in 19.15.0).
  rm_ <- PlotNumInteractionEffectsMatrix(d, predictor_vars = c("a1", "a2"),
                                         outcome_vars = c("b1", "b2"),
                                         interVar = "g1")
  tab <- rm_$p$data
  expect_false(anyNA(tab$P))
  expect_equal(tab$P_FDR, manual_matrix_adjust(tab$P))
  # one cell cross-checked against a direct lm fit
  m_direct <- lm(b1 ~ a1 * g1, data = d)
  p_direct <- utils::tail(as.data.frame(summary(m_direct)$coefficients)$`Pr(>|t|)`, 1)
  expect_equal(tab$P[tab$X == "a1" & tab$Y == "b1"], p_direct)

  # single and multiple covariates also fit real models
  rc <- PlotNumInteractionEffectsMatrix(d, predictor_vars = c("a1", "a2"),
                                        outcome_vars = c("b1", "b2"),
                                        interVar = "g1", covariates = c("a3", "b2"))
  expect_false(anyNA(rc$p$data$P))

  # missing interVar errors clearly instead of returning all-NA
  expect_error(
    PlotNumInteractionEffectsMatrix(d, predictor_vars = "a1", outcome_vars = "b1"),
    "interVar must be provided"
  )

  rp <- PlotNumInteractionEffectsMatrix(d, predictor_vars = c("a1", "a2"),
                                        outcome_vars = c("b1", "b2"),
                                        interVar = "g1", fdr_scope = "per_outcome")
  tabp <- rp$p$data
  expect_equal(tabp$P_FDR, manual_grouped_adjust(tabp$P, tabp$Y))
  expect_false(isTRUE(all.equal(tabp$P_FDR, tab$P_FDR)))
  expect_identical(rm_$p_fdr, rm_$p_FDR)
})

test_that("PlotDirectionalHeatmaps threads fdr_scope to its sub-analyses", {
  d <- make_fdr_data()
  cont <- c("a1", "a2", "a3", "b1")

  key_cols <- function(df) {
    df <- as.data.frame(df)[, c("XVar", "YVar", "stars", "stars_FDR")]
    df[order(df$XVar, df$YVar), , drop = FALSE]
  }

  for (scope in c("matrix", "per_outcome")) {
    dir_res <- PlotDirectionalHeatmaps(d, variables = cont, fdr_scope = scope)
    cor_res <- PlotCorrelationsHeatmap(d, predictor_vars = cont, fdr_scope = scope)
    expect_equal(
      key_cols(dir_res$Unadjusted$data),
      key_cols(cor_res$Unadjusted$plot$data),
      ignore_attr = TRUE
    )
  }

  res <- PlotDirectionalHeatmaps(d, variables = cont)
  expect_identical(res$p, res$Unadjusted)
  expect_identical(res$p_fdr, res$FDRCorrected)
})
