#' @title Make Comparison Table with Covariate Adjustment, Effect Sizes, and Pairwise Contrasts
#'
#' @description
#' Creates a comparison table using {gtsummary} summarizing continuous and categorical variables
#' by a grouping factor, with optional covariate adjustment for continuous outcomes only (ANCOVA),
#' effect-size calculations, and pairwise contrasts with p-value adjustments (or logistic regression for binary outcomes).
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
#' @param PairwiseMethod Character; p-value adjustment method for contrasts (e.g. "bonferroni", default). Use "none" for unadjusted.
#' @param Parametric Logical; use parametric tests (ANOVA/ANCOVA) if `TRUE`, or nonparametric (Kruskal–Wallis) if `FALSE` (default `TRUE`).
#' @param suppress_warnings Logical; suppress intermediate warnings from gtsummary (default `FALSE`).
#'
#' @return A `gtsummary::tbl_summary` object augmented with formatted p-values, effect sizes, and contrasts.
#' @export
MakeComparisonTable <- function(
    DataFrame,
    Variables,
    CompVariable,
    Covariates = NULL,
    ValueDigits = 2,
    pDigits = 3,
    AddEffectSize = FALSE,
    EffectSizeDigits = 2,
    AddPairwise = FALSE,
    PairwiseMethod = "bonferroni",
    Parametric = TRUE,
    suppress_warnings = FALSE
) {
  # ensure required packages
  required_pkgs <- c("gtsummary","dplyr","car","emmeans","broom","effectsize","purrr","tidyr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing_pkgs)) stop("Install packages: ", paste(missing_pkgs, collapse=", "))

  # validate inputs
  if (!CompVariable %in% names(DataFrame)) stop("Grouping variable not found: ", CompVariable)
  if (!all(Variables %in% names(DataFrame))) stop("Variable(s) not found: ", paste(setdiff(Variables, names(DataFrame)), collapse=", "))
  if (!is.null(Covariates) && !all(Covariates %in% names(DataFrame))) stop("Covariate(s) not found: ", paste(setdiff(Covariates, names(DataFrame)), collapse=", "))

  # prepare data
  df <- DataFrame %>%
    dplyr::select(all_of(c(CompVariable, Variables, Covariates))) %>%
    dplyr::select_if(~ !(is.factor(.) && nlevels(.) > 15))
  df[[CompVariable]] <- factor(df[[CompVariable]])

  # base summary
  tbl <- gtsummary::tbl_summary(
    df,
    by = CompVariable,
    statistic = list(
      all_continuous() ~ if (Parametric) "{mean} ({sd})" else "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p})"
    ),
    digits = list(all_continuous() ~ ValueDigits)
  ) %>% gtsummary::add_n()
  if (suppress_warnings) tbl <- suppressWarnings(tbl)

  # compute p-values
  pvals <- purrr::map_dfr(Variables, function(var) {
    # unadjusted
    if (is.numeric(df[[var]])) {
      pu <- if (Parametric) summary(aov(as.formula(paste(var, '~', CompVariable)), df))[[1]][CompVariable,'Pr(>F)']
      else kruskal.test(as.formula(paste(var, '~', CompVariable)), df)$p.value
    } else {
      tbl0 <- table(df[[var]], df[[CompVariable]])
      pu <- if (any(suppressWarnings(chisq.test(tbl0)$expected) < 5)) fisher.test(tbl0)$p.value else chisq.test(tbl0)$p.value
    }
    # adjusted
    pad <- NA_real_
    if (Parametric && !is.null(Covariates)) {
      if (is.numeric(df[[var]])) {
        pad <- car::Anova(lm(as.formula(paste(var, '~', CompVariable, '+', paste(Covariates, collapse='+'))), df), type=2)[CompVariable,'Pr(>F)']
      } else if (nlevels(df[[var]]) == 2) {
        glm_mod <- glm(as.formula(paste(var, '~', CompVariable, '+', paste(Covariates, collapse='+'))), df, family=binomial)
        pad <- drop1(glm_mod, test='Chisq')[CompVariable,'Pr(>Chi)']
      }
    }
    tibble::tibble(variable=var, p_unadj=pu, p_adj=pad)
  })

  # add p-values
  test_note <- if (Parametric) 'Continuous: ANOVA/ANCOVA; Categorical: chi-square/Fisher or logistic adj' else 'Continuous: Kruskal-Wallis; Categorical: chi-square/Fisher'
  tbl <- tbl %>%
    modify_table_body(~ .x %>%
                        left_join(pvals, by='variable') %>%
                        mutate(
                          p.value = coalesce(p_adj, p_unadj),
                          p.value_fmt = gtsummary::style_pvalue(p.value, digits=pDigits),
                          p.value_fmt = ifelse(row_type=='label', p.value_fmt, NA_character_)
                        )
    ) %>%
    modify_table_styling(columns='p.value_fmt', rows = p.value <= 0.05 & row_type=='label', text_format='bold') %>%
    modify_header(p.value_fmt ~ '**p-value**') %>%
    modify_footnote(p.value_fmt ~ test_note)

  # effect sizes
  if (AddEffectSize) {
    es_df <- purrr::map_dfr(Variables, function(var) {
      if (is.numeric(df[[var]])) {
        k <- nlevels(df[[CompVariable]]); n <- nrow(df)
        if (Parametric && k==2 && is.null(Covariates)) {
          d <- effectsize::cohens_d(as.formula(paste(var,'~',CompVariable)), df); method <- "Cohen's d"; es <- d$Cohens_d
        } else if (Parametric && !is.null(Covariates)) {
          et <- effectsize::eta_squared(car::Anova(lm(as.formula(paste(var,'~',CompVariable,'+',paste(Covariates,collapse='+'))),df),type=2), partial=TRUE)
          method <- 'Partial η²'; es <- et$Eta2_partial[1]
        } else if (Parametric) {
          et <- effectsize::eta_squared(aov(as.formula(paste(var,'~',CompVariable)), df), partial=FALSE); method <- 'η²'; es <- et$Eta2[1]
        } else {
          H <- kruskal.test(as.formula(paste(var,'~',CompVariable)), df)$statistic; method <- 'ε²'; es <- (H - k +1)/(n-k)
        }
      } else {
        tbl0 <- table(df[[var]], df[[CompVariable]]); chi <- suppressWarnings(chisq.test(tbl0)$statistic);
        method <- "Cramer's V"; es <- sqrt(chi/(sum(tbl0)*(min(nrow(tbl0),ncol(tbl0))-1)))
      }
      tibble::tibble(variable=var, effect_size=round(es,EffectSizeDigits), es_method=method)
    })
    tbl <- tbl %>%
      modify_table_body(~ .x %>% left_join(es_df, by='variable') %>%
                          mutate(effect_size = ifelse(row_type=='label', effect_size, NA_real_))
      ) %>%
      modify_fmt_fun(effect_size ~ function(x) formatC(x, digits=EffectSizeDigits, format='f')) %>%
      modify_header(effect_size ~ '**Effect Size**') %>%
      modify_footnote(effect_size ~ sprintf('Effect size: Continuous—%s; Categorical—Cramer\'s V', unique(es_df$es_method)))
  }

  # pairwise contrasts
  if (AddPairwise && nlevels(df[[CompVariable]])>1) {
    pw <- purrr::map_dfr(Variables, function(var) {
      if (is.numeric(df[[var]])) {
        if (Parametric) {
          fit<-aov(as.formula(paste(var,'~',CompVariable)),df)
          r<-broom::tidy(emmeans::contrast(emmeans::emmeans(fit,CompVariable),'pairwise',adjust=PairwiseMethod))
          pcol<-intersect(c('adj.p.value','p.value'),names(r))[1]
          r$contrast<-sapply(r$contrast,function(x)paste(sort(trimws(strsplit(x,' - ')[[1]])),collapse='-'))
          tibble::tibble(variable=var, contrast=r$contrast, p_val=r[[pcol]])
        } else {
          combos<-combn(levels(df[[CompVariable]]),2, simplify=FALSE)
          purrr::map_dfr(combos,function(cp){
            m<-pairwise.wilcox.test(df[[var]][df[[CompVariable]]%in%cp], df[[CompVariable]][df[[CompVariable]]%in%cp], p.adjust.method=PairwiseMethod)$p.value
            dfm<-as.data.frame(as.table(m))
            tibble::tibble(variable=var, contrast=paste(dfm$Var1,dfm$Var2,sep='-'), p_val=dfm$Freq)
          })
        }
      } else {
        combos<-combn(levels(df[[CompVariable]]),2, simplify=FALSE)
        purrr::map_dfr(combos,function(cp){
          tbl0<-table(df[[var]][df[[CompVariable]]%in%cp], df[[CompVariable]][df[[CompVariable]]%in%cp])
          pv<-if(any(suppressWarnings(chisq.test(tbl0)$expected)<5)) fisher.test(tbl0)$p.value else chisq.test(tbl0)$p.value
          tibble::tibble(variable=var, contrast=paste(sort(cp),collapse='-'), p_val=pv)
        }) %>% dplyr::mutate(p_val = if(PairwiseMethod!='none') p.adjust(p_val, method=PairwiseMethod) else p_val)
      }
    })
    pw_wide<-tidyr::pivot_wider(pw, id_cols='variable', names_from='contrast', values_from='p_val', values_fn=list(p_val=mean))
    tbl <- tbl %>%
      modify_table_body(~ .x %>% left_join(pw_wide, by='variable') %>%
                          dplyr::mutate(across(all_of(names(pw_wide)[-1]), ~ ifelse(row_type=='label', ., NA_real_)))
      )
    for(col in names(pw_wide)[-1]) tbl<-tbl %>%
      modify_fmt_fun(!!sym(col) ~ function(x) style_pvalue(x, digits=pDigits)) %>%
      modify_table_styling(columns=col, rows=.data[[col]]<=0.05 & row_type=='label', text_format='bold') %>%
      modify_header(!!sym(col) := paste0('**',col,'**')) %>%
      modify_footnote(!!sym(col) ~ paste0('Pairwise p-value (', PairwiseMethod, ')'))
  }

  # final caption
  cap <- sprintf('Comparison Table (%s analysis%s)', if(Parametric) 'parametric' else 'non-parametric',
                 if(!is.null(Covariates) && Parametric) sprintf('; adjusted for %s', paste(Covariates, collapse=',')) else '')
  tbl %>% modify_caption(cap)
}
