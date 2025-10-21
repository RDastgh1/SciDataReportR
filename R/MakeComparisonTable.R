#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, Pairwise Contrasts, and Optional Overall Summary
#'
#' @description
#' Creates a comparison table (via **gtsummary**) that
#'   • summarises variables by a grouping factor,
#'   • optionally adjusts continuous outcomes for covariates (ANCOVA),
#'   • adds p-values, effect sizes, and pair-wise contrasts,
#'   • automatically drops constant variables,
#'   • can suppress the “Unknown” row or add an overall *N* column, and
#'   • can return an overall-only summary when requested or when no valid group exists.
#'
#' @param DataFrame      Data frame with the raw data.
#' @param Variables      Character vector of columns to compare.
#' @param CompVariable   Grouping / comparison variable (character scalar). Optional if `IncludeOverallStats = TRUE`.
#' @param Covariates     Optional covariates (default `NULL`).
#' @param ValueDigits    Digits in summary stats (default `2`).
#' @param pDigits        Digits in p-values (default `3`).
#' @param AddEffectSize  Logical; add effect sizes? (default `FALSE`).
#' @param EffectSizeDigits Digits in effect sizes (default `2`).
#' @param AddPairwise    Logical; add pair-wise contrasts? (default `FALSE`).
#' @param PairwiseMethod P-adjustment method (default `"bonferroni"`).
#' @param Parametric     `TRUE` = parametric, `FALSE` = non-parametric tests.
#' @param IncludeOverallN Logical; add overall *N*? (default `FALSE`).
#' @param IncludeMissing Logical; include “Unknown” row? (default `FALSE`).
#' @param suppress_warnings Suppress gtsummary warnings? (default `FALSE`).
#' @param Referent       Optional reference level for contrasts.
#' @param IncludeOverallStats Logical; if `TRUE` (default `FALSE`) or if `CompVariable`
#'                            is not present in `DataFrame`, return overall-only stats.
#' @param ShowPositiveBinaryOnLabel Logical; if `TRUE` (default) show only the
#'                            “positive” level (TRUE/1/YES/Yes) of binary
#'                            categorical/factor/logical variables on the label row.
#'
#' @return A **gtsummary::tbl_summary** object.
#' @export
MakeComparisonTable <- function(
    DataFrame,
    Variables,
    CompVariable         = NULL,
    Covariates           = NULL,
    ValueDigits          = 2,
    pDigits              = 3,
    AddEffectSize        = FALSE,
    EffectSizeDigits     = 2,
    AddPairwise          = FALSE,
    PairwiseMethod       = "bonferroni",
    Parametric           = TRUE,
    IncludeOverallN      = FALSE,
    IncludeMissing       = FALSE,
    suppress_warnings    = FALSE,
    Referent             = NULL,
    IncludeOverallStats  = FALSE,
    ShowPositiveBinaryOnLabel = TRUE
) {
  ## ── required packages ───────────────────────────────────────────────────
  req_pkgs <- c(
    "gtsummary", "dplyr", "car", "emmeans", "broom", "effectsize",
    "purrr", "tidyr", "rlang", "tibble", "tidyselect"
  )
  ok <- vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok))
    stop("Please install: ", paste(req_pkgs[!ok], collapse = ", "))

  ## ── mode detection ──────────────────────────────────────────────────────
  comp_present <- !is.null(CompVariable) && is.character(CompVariable) &&
    length(CompVariable) == 1 && CompVariable %in% names(DataFrame)
  overall_mode <- isTRUE(IncludeOverallStats) || !comp_present

  ## ── sanity checks ───────────────────────────────────────────────────────
  if (!overall_mode && !CompVariable %in% names(DataFrame))
    stop("Grouping variable not found: ", CompVariable)
  if (!all(Variables %in% names(DataFrame)))
    stop("Variable(s) not found: ",
         paste(setdiff(Variables, names(DataFrame)), collapse = ", "))
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame)))
    stop("Covariate(s) not found: ",
         paste(setdiff(Covariates, names(DataFrame)), collapse = ", "))

  ## ── data prep ───────────────────────────────────────────────────────────
  cols <- c(Variables, Covariates)
  if (!overall_mode) cols <- c(CompVariable, cols)
  df <- dplyr::select(DataFrame, tidyselect::all_of(cols))
  if (!overall_mode) df[[CompVariable]] <- factor(df[[CompVariable]])

  ## drop constants ---------------------------------------------------------
  keep <- Variables[sapply(Variables, function(v)
    length(unique(df[[v]][!is.na(df[[v]])])) >= 2)]
  drop <- setdiff(Variables, keep)
  if (length(drop))
    warning("Dropping constant variable(s): ", paste(drop, collapse = ", "))
  Variables <- keep
  if (!length(Variables))
    stop("No variables left to summarize after dropping constants.")

  ## types for continuous vars ---------------------------------------------
  numeric_vars <- Variables[sapply(df[Variables], is.numeric)]
  type_list <- if (length(numeric_vars) > 0) {
    rlang::set_names(as.list(rep("continuous", length(numeric_vars))), numeric_vars)
  } else NULL

  ## ── binary-positive display map for tbl_summary(value=) -----------------
  value_list <- NULL
  if (isTRUE(ShowPositiveBinaryOnLabel)) {
    pos_tokens <- c("TRUE", "1", "YES", "Yes")
    vmap <- list()
    for (v in Variables) {
      x <- df[[v]]
      if (is.numeric(x)) next
      ux <- unique(x[!is.na(x)])
      if (length(ux) != 2) next

      if (is.logical(x)) {
        vmap[[v]] <- TRUE
      } else if (is.factor(x)) {
        levs <- levels(droplevels(x))
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      } else {
        levs <- sort(unique(as.character(ux)))
        hit  <- levs[levs %in% pos_tokens]
        if (length(hit) >= 1) vmap[[v]] <- hit[1]
      }
    }
    if (length(vmap)) value_list <- vmap
  }

  ## ── statistic templates ─────────────────────────────────────────────────
  stat_cont <- if (Parametric) "{mean} ({sd})" else "{median} [{p25}, {p75}]"
  stat_cat  <- "{n} ({p}%)"

  ## ── OVERALL-ONLY MODE ───────────────────────────────────────────────────
  if (overall_mode) {
    tbl <- gtsummary::tbl_summary(
      df,
      include   = tidyselect::all_of(Variables),
      missing   = if (IncludeMissing) "ifany" else "no",
      statistic = list(
        gtsummary::all_continuous()  ~ stat_cont,
        gtsummary::all_categorical() ~ stat_cat
      ),
      digits    = list(gtsummary::all_continuous() ~ ValueDigits),
      type      = type_list,
      value     = value_list
    )
    if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
    if (suppress_warnings) tbl <- suppressWarnings(tbl)

    cap <- sprintf("Overall Summary (%s)", if (Parametric) "parametric" else "non-parametric")
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  ## ── GROUPED SUMMARY ─────────────────────────────────────────────────────
  tbl <- gtsummary::tbl_summary(
    df,
    by        = CompVariable,
    include   = tidyselect::all_of(Variables),
    missing   = if (IncludeMissing) "ifany" else "no",
    statistic = list(
      gtsummary::all_continuous()  ~ stat_cont,
      gtsummary::all_categorical() ~ stat_cat
    ),
    digits    = list(gtsummary::all_continuous() ~ ValueDigits),
    type      = type_list,
    value     = value_list
  )
  if (IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if (suppress_warnings) tbl <- suppressWarnings(tbl)

  if (nlevels(df[[CompVariable]]) < 2) {
    cap <- sprintf("Overall Summary (only one level in '%s')", CompVariable)
    return(tbl %>% gtsummary::modify_caption(cap))
  }

  ## ── helpers for downstream steps ────────────────────────────────────────
  .cramers_v_rc <- function(tab) {
    tab <- tab[rowSums(tab) > 0, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
    if (requireNamespace("DescTools", quietly = TRUE))
      return(as.numeric(DescTools::CramerV(tab, method = "bias.corrected")))
    chi <- suppressWarnings(stats::chisq.test(tab)$statistic)
    nt  <- sum(tab); m <- min(nrow(tab), ncol(tab)) - 1
    sqrt(as.numeric(chi) / (nt * m))
  }
  .pair_label      <- function(a, b) paste(a, b, sep = " - ")
  .canonical_label <- function(lab) {
    p <- trimws(strsplit(lab, " - ", fixed = TRUE)[[1]])
    .pair_label(sort(p)[1], sort(p)[2])
  }

  ## ── p-values & test labels ──────────────────────────────────────────────
  pdat <- purrr::map_dfr(Variables, function(var) {
    df_var <- df[!is.na(df[[var]]), , drop = FALSE]
    df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
    k <- nlevels(df_var[[CompVariable]])

    if (k < 2)
      return(tibble::tibble(variable = var,
                            p_unadj = NA_real_, p_adj = NA_real_,
                            test_label = NA_character_))

    if (is.numeric(df_var[[var]])) {                     # continuous
      fmla_unadj <- reformulate(CompVariable, response = as.name(var))

      if (Parametric) {
        if (!is.null(Covariates)) {                      # ANCOVA (Type II)
          fmla_adj <- reformulate(c(CompVariable, Covariates), response = as.name(var))
          m1   <- stats::lm(fmla_adj, data = df_var)
          pu   <- summary(stats::aov(fmla_unadj, data = df_var))[[1]][CompVariable, "Pr(>F)"]
          pad  <- car::Anova(m1, type = 2)[CompVariable, "Pr(>F)"]
          tl   <- "ANCOVA (Type II)"
        } else if (k == 2) {                             # Welch t
          grp_n <- table(df_var[[CompVariable]])
          if (min(grp_n) < 2) {
            pu <- pad <- NA_real_; tl <- "Insufficient data"
          } else {
            pu <- stats::t.test(fmla_unadj, data = df_var, var.equal = FALSE)$p.value
            pad <- NA_real_; tl <- "Welch t-test"
          }
        } else {                                         # ANOVA
          pu  <- summary(stats::aov(fmla_unadj, data = df_var))[[1]][CompVariable, "Pr(>F)"]
          pad <- NA_real_; tl <- "ANOVA"
        }
      } else {                                           # non-parametric
        pu <- if (k == 2) {
          stats::wilcox.test(fmla_unadj, data = df_var)$p.value
        } else {
          stats::kruskal.test(fmla_unadj, data = df_var)$p.value
        }
        pad <- NA_real_; tl <- if (k == 2) "Wilcoxon rank-sum" else "Kruskal–Wallis"
      }

    } else {                                             # categorical
      tbl0 <- table(df_var[[var]], df_var[[CompVariable]])
      tbl0 <- tbl0[rowSums(tbl0) > 0, , drop = FALSE]
      tbl0 <- tbl0[, colSums(tbl0) > 0, drop = FALSE]

      if (nrow(tbl0) < 2 || ncol(tbl0) < 2) {
        pu <- pad <- NA_real_; tl <- "Insufficient categories"

      } else if (Parametric && !is.null(Covariates) &&
                 nlevels(droplevels(as.factor(df_var[[var]]))) == 2) {
        fmla_adj <- reformulate(c(CompVariable, Covariates), response = as.name(var))
        gm   <- stats::glm(fmla_adj, df_var, family = stats::binomial())
        pad  <- stats::drop1(gm, test = "Chisq")[CompVariable, "Pr(>Chi)"]
        pu   <- suppressWarnings(stats::chisq.test(tbl0)$p.value)
        tl   <- "Logistic regression (LR)"

      } else {
        expct <- suppressWarnings(stats::chisq.test(tbl0)$expected)
        if (any(expct < 5)) {
          pu <- stats::fisher.test(tbl0, simulate.p.value = TRUE, B = 1e4)$p.value
          tl <- "Fisher (sim.)"
        } else {
          pu <- stats::chisq.test(tbl0)$p.value
          tl <- "Chi-squared"
        }
        pad <- NA_real_
      }
    }

    tibble::tibble(variable = var, p_unadj = pu, p_adj = pad, test_label = tl)
  })

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   dplyr::left_join(pdat, by = "variable") %>%
                                   dplyr::mutate(
                                     p.value     = dplyr::coalesce(.data$p_adj, .data$p_unadj),
                                     p.value_fmt = gtsummary::style_pvalue(.data$p.value, digits = pDigits),
                                     p.value_fmt = ifelse(.data$row_type == "label", .data$p.value_fmt, NA_character_),
                                     Test        = ifelse(.data$row_type == "label", .data$test_label, NA_character_)
                                   )) %>%
    gtsummary::modify_table_styling(
      columns = "p.value_fmt",
      rows    = .data$p.value <= 0.05 & .data$row_type == "label",
      text_format = "bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt ~ "**p-value**") %>%
    gtsummary::modify_header(Test ~ "**Test**")

  ## ── effect sizes ────────────────────────────────────────────────────────
  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {
      df_var <- df[!is.na(df[[var]]), , drop = FALSE]
      df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
      k  <- nlevels(df_var[[CompVariable]])
      n  <- nrow(df_var)

      if (n < 2 || k < 2)
        return(tibble::tibble(variable = var,
                              effect_size = NA_real_, es_method = NA_character_))

      if (is.numeric(df_var[[var]])) {                   # continuous
        if (Parametric && k == 2 && is.null(Covariates)) {
          es_val <- tryCatch({
            abs(effectsize::cohens_d(
              reformulate(CompVariable, response = as.name(var)),
              data = df_var
            )$Cohens_d)
          }, error = function(e) NA_real_)
          method <- "|d|"

        } else if (Parametric && !is.null(Covariates)) {
          fmla_es <- reformulate(c(CompVariable, Covariates), response = as.name(var))
          lm_es   <- stats::lm(fmla_es, df_var)
          a2      <- car::Anova(lm_es, type = 2)
          et      <- effectsize::eta_squared(a2, partial = TRUE)
          idx <- if ("Parameter" %in% names(et)) et$Parameter == CompVariable else rownames(et) == CompVariable
          es_val  <- suppressWarnings(et$Eta2_partial[idx][1])
          method  <- "partial η²"

        } else if (Parametric) {
          fmla_es <- reformulate(CompVariable, response = as.name(var))
          es_val  <- effectsize::eta_squared(stats::aov(fmla_es, df_var), partial = FALSE)$Eta2[1]
          method  <- "η²"

        } else {
          H <- stats::kruskal.test(reformulate(CompVariable, response = as.name(var)), data = df_var)$statistic
          es_val <- as.numeric((H - k + 1) / (n - k))
          method <- "ε²"
        }

      } else {                                          # categorical
        tab    <- table(df_var[[var]], df_var[[CompVariable]])
        es_val <- .cramers_v_rc(tab)
        method <- "Cramer's V"
      }

      tibble::tibble(
        variable    = var,
        effect_size = ifelse(is.finite(es_val),
                             round(es_val, EffectSizeDigits), NA_real_),
        es_method   = method
      )
    })

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     dplyr::left_join(es_df, by = "variable") %>%
                                     dplyr::mutate(
                                       effect_size = dplyr::if_else(.data$row_type == "label",
                                                                    .data$effect_size, NA_real_),
                                       ES_Method   = dplyr::if_else(.data$row_type == "label",
                                                                    .data$es_method,   NA_character_)
                                     )) %>%
      gtsummary::modify_fmt_fun(effect_size ~ function(x)
        ifelse(is.na(x), NA_character_,
               formatC(x, digits = EffectSizeDigits, format = "f"))) %>%
      gtsummary::modify_header(effect_size ~ "**Effect Size**") %>%
      gtsummary::modify_header(ES_Method   ~ "**ES Method**")
  }

  ## ── pair-wise contrasts ────────────────────────────────────────────────
  if (AddPairwise && nlevels(df[[CompVariable]]) > 1) {

    lvls   <- levels(df[[CompVariable]])
    combos <- if (!is.null(Referent)) {
      if (!Referent %in% lvls)
        stop("Referent level not found: ", Referent)
      lapply(setdiff(lvls, Referent), function(x) c(Referent, x))
    } else utils::combn(lvls, 2, simplify = FALSE)

    pw <- purrr::map_dfr(Variables, function(var) {

      empty_pw <- tibble::tibble(
        variable = character(), contrast = character(), p_val = numeric()
      )

      ## continuous ---------------------------------------------------------
      if (is.numeric(df[[var]])) {
        if (Parametric) {
          res <- tryCatch(
            stats::pairwise.t.test(
              x = df[[var]],
              g = df[[CompVariable]],
              p.adjust.method = PairwiseMethod,
              pool.sd = FALSE
            ),
            error = function(e) NULL
          )
          if (is.null(res) || is.null(res$p.value)) return(empty_pw)
          r <- as.data.frame(as.table(res$p.value))
          if (!all(c("Var1","Var2","Freq") %in% names(r))) return(empty_pw)
          r <- r %>%
            dplyr::filter(!is.na(Freq)) %>%
            dplyr::transmute(
              variable = var,
              g1 = as.character(Var1),
              g2 = as.character(Var2),
              contrast = .canonical_label(.pair_label(g1, g2)),
              p_val    = Freq
            )
          if (!is.null(Referent))
            r <- dplyr::filter(r, g1 == Referent | g2 == Referent)
          dplyr::select(r, variable, contrast, p_val)

        } else {                                        # Wilcoxon
          if (length(combos) == 0) return(empty_pw)
          purrr::map_dfr(combos, function(cp) {
            idx <- df[[CompVariable]] %in% cp
            res <- tryCatch(
              stats::pairwise.wilcox.test(
                x = df[[var]][idx],
                g = droplevels(df[[CompVariable]][idx]),
                p.adjust.method = PairwiseMethod
              )$p.value,
              error = function(e) NULL
            )
            if (is.null(res)) return(empty_pw)
            rr <- as.data.frame(as.table(res))
            if (!all(c("Var1","Var2","Freq") %in% names(rr))) return(empty_pw)
            rr %>%
              dplyr::transmute(
                variable = var,
                contrast = .canonical_label(.pair_label(Var1, Var2)),
                p_val    = Freq
              )
          })
        }

        ## categorical --------------------------------------------------------
      } else {
        if (length(combos) == 0) return(empty_pw)
        purrr::map_dfr(combos, function(cp) {
          sub   <- df[df[[CompVariable]] %in% cp, ]
          tbl0  <- table(droplevels(as.factor(sub[[var]])),
                         droplevels(sub[[CompVariable]]))
          pv <- if (nrow(tbl0) < 2 || ncol(tbl0) < 2) NA_real_
          else if (any(
            suppressWarnings(stats::chisq.test(tbl0)$expected) < 5))
            stats::fisher.test(tbl0, simulate.p.value = TRUE,
                               B = 1e4)$p.value
          else stats::chisq.test(tbl0)$p.value
          tibble::tibble(
            variable = var,
            contrast = .canonical_label(.pair_label(cp[1], cp[2])),
            p_val    = pv
          )
        }) %>%
          dplyr::mutate(
            p_val = if (PairwiseMethod != "none")
              p.adjust(p_val, PairwiseMethod) else p_val)
      }
    })

    if (nrow(pw) > 0) {
      pw_wide <- tidyr::pivot_wider(
        pw, id_cols = "variable",
        names_from = "contrast",
        values_from = "p_val",
        values_fn = list(p_val = ~ mean(.x, na.rm = TRUE))
      )

      tbl <- tbl %>%
        gtsummary::modify_table_body(~ .x %>%
                                       dplyr::left_join(pw_wide, by = "variable") %>%
                                       dplyr::mutate(
                                         dplyr::across(tidyselect::all_of(
                                           setdiff(names(pw_wide), "variable")),
                                           ~ ifelse(.data$row_type == "label", ., NA_real_))
                                       ))

      for (col in setdiff(names(pw_wide), "variable")) {
        tbl <- tbl %>%
          gtsummary::modify_fmt_fun(
            !!rlang::sym(col) ~
              function(x) gtsummary::style_pvalue(x, digits = pDigits)
          ) %>%
          gtsummary::modify_table_styling(
            columns     = col,
            rows        = .data[[col]] <= 0.05 & .data$row_type == "label",
            text_format = "bold"
          ) %>%
          gtsummary::modify_header(
            !!rlang::sym(col) := paste0("**", col, "**")
          ) %>%
          gtsummary::modify_footnote(
            !!rlang::sym(col) ~ paste0("Pair-wise p-value (", PairwiseMethod,
                                       if (!is.null(Referent))
                                         paste0("; vs ", Referent), ")")
          )
      }
    }
  }

  ## ── caption & return ────────────────────────────────────────────────────
  cap <- sprintf(
    "Comparison Table (%s analysis%s) — see **Test** column for per-variable method",
    if (Parametric) "parametric" else "non-parametric",
    if (!is.null(Covariates) && Parametric)
      paste0("; adjusted for ",
             paste(Covariates, collapse = ", "))
    else ""
  )
  tbl %>% gtsummary::modify_caption(cap)
}
