.ScidrDirectionLabel <- function(group, alternative) {
  levels_group <- levels(droplevels(as.factor(group)))

  if (length(levels_group) != 2 || identical(alternative, "two.sided")) {
    return("two-sided")
  }

  if (identical(alternative, "greater")) {
    paste0(levels_group[2], " > ", levels_group[1])
  } else {
    paste0(levels_group[2], " < ", levels_group[1])
  }
}

.ScidrDirectedContinuousP <- function(x, group, Parametric, alternative) {
  group <- droplevels(as.factor(group))
  levels_group <- levels(group)

  if (length(levels_group) != 2 || identical(alternative, "two.sided")) {
    return(NA_real_)
  }

  x_first <- x[group == levels_group[1]]
  x_second <- x[group == levels_group[2]]

  tryCatch(
    if (isTRUE(Parametric)) {
      stats::t.test(x_second, x_first, alternative = alternative, var.equal = FALSE)$p.value
    } else {
      stats::wilcox.test(x_second, x_first, alternative = alternative, exact = FALSE)$p.value
    },
    error = function(e) NA_real_
  )
}

.ScidrDirectedBinaryP <- function(outcome, group, alternative, method = c("auto", "chisq", "fisher")) {
  method <- match.arg(method)
  outcome <- droplevels(as.factor(outcome))
  group <- droplevels(as.factor(group))

  if (nlevels(outcome) != 2 || nlevels(group) != 2 || identical(alternative, "two.sided")) {
    return(list(p = NA_real_, label = NA_character_))
  }

  levels_group <- levels(group)
  event <- levels(outcome)[2]
  success <- outcome == event
  n_group <- as.numeric(table(group))
  x_group <- c(
    sum(success[group == levels_group[1]]),
    sum(success[group == levels_group[2]])
  )

  # prop.test() and fisher.test() orient their alternatives as level 1 minus
  # level 2. The public API deliberately defines direction as level 2 vs level 1.
  base_alternative <- if (identical(alternative, "greater")) "less" else "greater"
  tab <- matrix(
    c(x_group[1], x_group[2], n_group[1] - x_group[1], n_group[2] - x_group[2]),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c(event, "Other"), levels_group)
  )

  chi <- tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)), error = function(e) NULL)
  use_fisher <- identical(method, "fisher") || is.null(chi) || any(chi$expected < 5)

  if (use_fisher) {
    p <- tryCatch(stats::fisher.test(tab, alternative = base_alternative)$p.value, error = function(e) NA_real_)
    return(list(p = p, label = "One-sided Fisher exact"))
  }

  p <- tryCatch(
    stats::prop.test(x = x_group, n = n_group, alternative = base_alternative, correct = FALSE)$p.value,
    error = function(e) NA_real_
  )
  list(p = p, label = "One-sided two-sample proportion test")
}

.ScidrDirectedCoefficientP <- function(fit, group_term, alternative, robust = FALSE) {
  if (identical(alternative, "two.sided")) return(NA_real_)

  mm <- tryCatch(stats::model.matrix(fit), error = function(e) NULL)
  if (is.null(mm)) return(NA_real_)

  group_pos <- match(group_term, attr(stats::terms(fit), "term.labels"))
  coefficient_index <- which(attr(mm, "assign") == group_pos)
  coefficient_index <- coefficient_index[coefficient_index != 1]

  if (length(coefficient_index) != 1) return(NA_real_)

  estimate <- stats::coef(fit)[coefficient_index]
  covariance <- if (isTRUE(robust)) {
    tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  } else {
    tryCatch(stats::vcov(fit), error = function(e) NULL)
  }
  if (is.null(covariance)) return(NA_real_)

  se <- sqrt(covariance[coefficient_index, coefficient_index])
  if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(NA_real_)

  statistic <- estimate / se
  if (isTRUE(robust) || inherits(fit, "glm")) {
    if (identical(alternative, "greater")) {
      return(stats::pnorm(statistic, lower.tail = FALSE))
    }
    return(stats::pnorm(statistic, lower.tail = TRUE))
  }

  if (identical(alternative, "greater")) {
    return(stats::pt(statistic, df = stats::df.residual(fit), lower.tail = FALSE))
  }
  stats::pt(statistic, df = stats::df.residual(fit), lower.tail = TRUE)
}
