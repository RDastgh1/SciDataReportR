#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, and Pairwise Contrasts
#'
#' @description
#' Creates a comparison table using {gtsummary} summarizing continuous and categorical variables
#' by a grouping factor, with optional covariate adjustment for continuous outcomes only (ANCOVA),
#' effect-size calculations, pairwise contrasts with p-value adjustments (or logistic regression for binary outcomes),
#' and optional suppression of the "Unknown" (missing) row.  Variables with fewer than
#' two unique non-NA values are dropped automatically before any tests.
#'
#' @param DataFrame A `data.frame` containing the raw data.
#' @param Variables Character vector of column names to include in the comparison.
#' @param CompVariable Character string specifying the grouping (comparison) variable.
#' @param Covariates Optional character vector of covariate column names for adjustment (default `NULL`).
#' @param ValueDigits Integer; number of digits for summary statistics (default `2`).
#' @param pDigits Integer; number of digits for formatted p-values (default `3`).
#' @param AddEffectSize Logical; include an effect-size column? (default `FALSE`).
#' @param EffectSizeDigits Integer; number of digits for effect sizes (default `2`).
#' @param AddPairwise Logical; include pairwise contrast columns? (default `FALSE`).
#' @param PairwiseMethod Character; p-value adjustment method for contrasts (default `"bonferroni"`; `"none"` for unadjusted).
#' @param Parametric Logical; use parametric tests (ANOVA/ANCOVA) if `TRUE`, or nonparametric (Kruskal–Wallis) if `FALSE` (default `TRUE`).
#' @param IncludeOverallN Logical; include a column with the overall N in the table? (default `FALSE`).
#' @param IncludeMissing Logical; include a row summarizing missing data ("Unknown")? (default `FALSE`).
#' @param suppress_warnings Logical; suppress intermediate warnings from gtsummary (default `FALSE`).
#'
#' @return A `gtsummary::tbl_summary` object augmented with formatted p-values, effect sizes, and contrasts.
#' @export
MakeComparisonTable <- function(
    DataFrame,
    Variables,
    CompVariable,
    Covariates       = NULL,
    ValueDigits      = 2,
    pDigits          = 3,
    AddEffectSize    = FALSE,
    EffectSizeDigits = 2,
    AddPairwise      = FALSE,
    PairwiseMethod   = "bonferroni",
    Parametric       = TRUE,
    IncludeOverallN  = FALSE,
    IncludeMissing   = FALSE,
    suppress_warnings= FALSE
) {
  req_pkgs <- c("gtsummary","dplyr","car","emmeans","broom","effectsize","purrr","tidyr","rlang")
  ok <- vapply(req_pkgs, requireNamespace, FUN.VALUE=logical(1), quietly=TRUE)
  if(any(!ok)) stop("Please install: ", paste(req_pkgs[!ok], collapse=", "))

  if(!CompVariable %in% names(DataFrame)) stop("Grouping var not found: ", CompVariable)
  if(!all(Variables %in% names(DataFrame))) stop("Variable(s) not found: ",
                                                 paste(setdiff(Variables, names(DataFrame)), collapse=", "))
  if(!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) stop("Covariate(s) not found: ",
                                                                          paste(setdiff(Covariates, names(DataFrame)), collapse=", "))

  df <- DataFrame %>% dplyr::select(all_of(c(CompVariable, Variables, Covariates)))
  df[[CompVariable]] <- factor(df[[CompVariable]])

  keep <- Variables[sapply(Variables, function(v)
    length(unique(df[[v]][!is.na(df[[v]])])) >= 2
  )]
  dropped <- setdiff(Variables, keep)
  if(length(dropped)) warning("Dropping: ", paste(dropped, collapse=", "))
  Variables <- keep
  if(length(Variables)==0) stop("No variables to compare after dropping constants.")

  # explicitly tag numeric variables as continuous
  numeric_vars <- Variables[sapply(df[Variables], is.numeric)]
  type_list <- rlang::set_names(as.list(rep("continuous", length(numeric_vars))), numeric_vars)

  tbl <- gtsummary::tbl_summary(
    df,
    by      = CompVariable,
    missing = if(IncludeMissing) "ifany" else "no",
    statistic = list(
      all_continuous()  ~ if(Parametric) "{mean} ({sd})" else "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(all_continuous() ~ ValueDigits),
    type   = type_list
  )
  if(IncludeOverallN) tbl <- tbl %>% gtsummary::add_n()
  if(suppress_warnings) tbl <- suppressWarnings(tbl)

  if(nlevels(df[[CompVariable]])<2) {
    warning("Only one level of ",CompVariable,"; returning summary.")
    return(tbl)
  }

  # unadjusted and adjusted p-values
  pvals <- purrr::map_dfr(Variables, function(var) {
    df_var <- df[!is.na(df[[var]]), , drop=FALSE]
    df_var[[CompVariable]] <- droplevels(df_var[[CompVariable]])
    if(nlevels(df_var[[CompVariable]])<2) {
      return(tibble::tibble(variable=var, p_unadj=NA_real_, p_adj=NA_real_))
    }
    if(is.numeric(df_var[[var]])) {
      pu <- if(Parametric) {
        summary(aov(reformulate(CompVariable, response = var), data=df_var))[[1]][CompVariable,"Pr(>F)"]
      } else {
        kruskal.test(reformulate(CompVariable, response = var), data=df_var)$p.value
      }
    } else {
      tbl0 <- table(df_var[[var]], df_var[[CompVariable]])
      tbl0 <- tbl0[rowSums(tbl0)>0, colSums(tbl0)>0, drop=FALSE]
      pu <- if(nrow(tbl0)<2||ncol(tbl0)<2) {
        NA_real_
      } else if(any(suppressWarnings(chisq.test(tbl0)$expected)<5)) {
        fisher.test(tbl0, simulate.p.value=TRUE, B=1e5)$p.value
      } else {
        chisq.test(tbl0)$p.value
      }
    }

    pad <- NA_real_
    if(Parametric && !is.null(Covariates)) {
      if(is.numeric(df_var[[var]])) {
        fmla <- reformulate(c(CompVariable, Covariates), response = var)
        m1   <- lm(fmla, data=df_var)
        pad  <- car::Anova(m1,type=2)[CompVariable,"Pr(>F)"]
      } else if(nlevels(df_var[[var]])==2) {
        fmla <- reformulate(c(CompVariable, Covariates), response = var)
        gm   <- glm(fmla, data=df_var, family=binomial)
        pad  <- drop1(gm, test="Chisq")[CompVariable,"Pr(>Chi)"]
      }
    }

    tibble::tibble(variable=var, p_unadj=pu, p_adj=pad)
  })

  test_note <- if(Parametric)
    "Continuous: ANOVA/ANCOVA; Categorical: chi-square/Fisher or logistic regression adj"
  else
    "Continuous: Kruskal-Wallis; Categorical: chi-square/Fisher"

  tbl <- tbl %>%
    gtsummary::modify_table_body(~ .x %>%
                                   left_join(pvals, by="variable") %>%
                                   mutate(
                                     p.value     = coalesce(p_adj, p_unadj),
                                     p.value_fmt = gtsummary::style_pvalue(p.value, digits=pDigits),
                                     p.value_fmt = ifelse(row_type=="label", p.value_fmt, NA_character_)
                                   )
    ) %>%
    gtsummary::modify_table_styling(
      columns="p.value_fmt",
      rows   = p.value<=0.05 & row_type=="label",
      text_format="bold"
    ) %>%
    gtsummary::modify_header(p.value_fmt~"**p-value**") %>%
    gtsummary::modify_footnote(p.value_fmt~test_note)

  if(AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {
      k <- nlevels(df[[CompVariable]]); n <- nrow(df)
      if(is.numeric(df[[var]])) {
        if(Parametric && k==2 && is.null(Covariates)) {
          es <- abs(effectsize::cohens_d(reformulate(CompVariable, response = var), data=df)$Cohens_d)
          method <- "|Cohen's d|"
        } else if(Parametric && !is.null(Covariates)) {
          m1 <- lm(reformulate(c(CompVariable, Covariates), response = var), data=df)
          et <- effectsize::eta_squared(car::Anova(m1,type=2), partial=TRUE)
          es <- et$Eta2_partial[1]; method <- "Partial η²"
        } else if(Parametric) {
          m0 <- aov(reformulate(CompVariable, response = var), data=df)
          et <- effectsize::eta_squared(m0, partial=FALSE)
          es <- et$Eta2[1]; method <- "η²"
        } else {
          H  <- kruskal.test(reformulate(CompVariable, response = var), data=df)$statistic
          es <- (H - k + 1)/(n - k); method <- "ε²"
        }
      } else {
        tbl0 <- table(df[[var]], df[[CompVariable]])
        chi  <- suppressWarnings(chisq.test(tbl0)$statistic)
        es   <- sqrt(chi/(sum(tbl0)*(min(nrow(tbl0),ncol(tbl0))-1)))
        method <- "Cramer's V"
      }
      tibble::tibble(variable=var, effect_size=round(es,EffectSizeDigits), es_method=method)
    })

    cont_methods <- unique(es_df$es_method[es_df$es_method!="Cramer's V"])
    footnote_es  <- sprintf("Effect size: Continuous—%s; Categorical—Cramer's V",
                            paste(cont_methods, collapse=", "))

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     left_join(es_df, by="variable") %>%
                                     mutate(effect_size = ifelse(row_type=="label", effect_size, NA_real_))
      ) %>%
      gtsummary::modify_fmt_fun(
        effect_size~function(x) formatC(x, digits=EffectSizeDigits, format="f")
      ) %>%
      gtsummary::modify_header(effect_size~"**Effect Size**") %>%
      gtsummary::modify_footnote(effect_size~footnote_es)
  }

  if(AddPairwise && nlevels(df[[CompVariable]])>1) {
    pw <- purrr::map_dfr(Variables, function(var) {
      if(is.numeric(df[[var]])) {
        if(Parametric) {
          fit <- aov(reformulate(CompVariable, response = var), data=df)
          em  <- emmeans::emmeans(fit, CompVariable)
          ct  <- emmeans::contrast(em, method="pairwise", adjust=PairwiseMethod)
          r   <- broom::tidy(ct)
          pcol<- intersect(c("adj.p.value","p.value"), names(r))[1]
          r$contrast <- sapply(r$contrast, function(x)
            paste(sort(trimws(strsplit(x," - ")[[1]])), collapse="-"))
          tibble::tibble(variable=var, contrast=r$contrast, p_val=r[[pcol]])
        } else {
          combos <- combn(levels(df[[CompVariable]]),2, simplify=FALSE)
          purrr::map_dfr(combos, function(cp) {
            res <- pairwise.wilcox.test(
              df[[var]][df[[CompVariable]]%in%cp],
              df[[CompVariable]][df[[CompVariable]]%in%cp],
              p.adjust.method=PairwiseMethod
            )$p.value
            as.data.frame(as.table(res)) %>%
              transmute(variable=var,
                        contrast=paste(Var1,Var2,sep="-"),
                        p_val=Freq)
          })
        }
      } else {
        combos <- combn(levels(df[[CompVariable]]),2, simplify=FALSE)
        purrr::map_dfr(combos, function(cp) {
          sub  <- df[df[[CompVariable]]%in%cp, ]
          x_sub<- droplevels(as.factor(sub[[var]]))
          g_sub<- droplevels(sub[[CompVariable]])
          tbl0 <- table(x_sub, g_sub)
          tbl0 <- tbl0[rowSums(tbl0)>0, colSums(tbl0)>0, drop=FALSE]
          if(nrow(tbl0)<2||ncol(tbl0)<2) {
            pv <- NA_real_
          } else if(any(suppressWarnings(chisq.test(tbl0)$expected)<5)) {
            pv <- fisher.test(tbl0, simulate.p.value=TRUE, B=1e5)$p.value
          } else {
            pv <- chisq.test(tbl0)$p.value
          }
          tibble::tibble(variable=var,
                         contrast=paste(sort(cp),collapse="-"),
                         p_val=pv)
        }) %>% dplyr::mutate(
          p_val = if(PairwiseMethod!="none") p.adjust(p_val, method=PairwiseMethod) else p_val
        )
      }
    })

    pw_wide <- tidyr::pivot_wider(
      pw,
      id_cols    = "variable",
      names_from  = "contrast",
      values_from = "p_val",
      values_fn   = list(p_val=mean)
    )

    tbl <- tbl %>%
      gtsummary::modify_table_body(~ .x %>%
                                     left_join(pw_wide, by="variable") %>%
                                     dplyr::mutate(across(
                                       all_of(setdiff(names(pw_wide),"variable")),
                                       ~ ifelse(row_type=="label", ., NA_real_)
                                     ))
      )

    for(col in setdiff(names(pw_wide),"variable")) {
      tbl <- tbl %>%
        gtsummary::modify_fmt_fun(
          !!rlang::sym(col)~function(x) gtsummary::style_pvalue(x,digits=pDigits)
        ) %>%
        gtsummary::modify_table_styling(
          columns     = col,
          rows        = .data[[col]] <= 0.05 & row_type=="label",
          text_format = "bold"
        ) %>%
        gtsummary::modify_header(!!rlang::sym(col):=paste0("**",col,"**")) %>%
        gtsummary::modify_footnote(
          !!rlang::sym(col) ~ paste0("Pairwise p-value (",PairwiseMethod,")")
        )
    }
  }

  cap <- sprintf(
    "Comparison Table (%s analysis%s)",
    if(Parametric) "parametric" else "non-parametric",
    if(!is.null(Covariates)&&Parametric)
      paste0("; adjusted for ",paste(Covariates,collapse=",")) else ""
  )
  tbl %>% gtsummary::modify_caption(cap)
}
