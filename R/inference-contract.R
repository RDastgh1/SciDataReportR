# Shared provenance and tidy results for inferential public functions.
# Collect numerical objects where they are calculated; never parse display text.

ScidrInferenceValue <- function(frame, candidates, default = NULL) {
  for (name in candidates) {
    if (!exists(name, envir = frame, inherits = FALSE)) next
    value <- tryCatch(get(name, envir = frame, inherits = FALSE), error = function(e) NULL)
    if (rlang::is_missing(value)) next
    if (!is.null(value) && !inherits(value, 'deprecated')) return(value)
  }
  default
}

ScidrStartInference <- function(function_name, call, frame) {
  context <- new.env(parent = emptyenv())
  context$Function <- function_name
  context$Call <- call
  context$Frame <- frame
  context$Observations <- list()
  context$ObservationStore <- new.env(hash = TRUE, parent = emptyenv())
  context$ObservationCount <- 0L
  context$Failures <- list()
  context$Estimates <- list()
  context$PairwisePValues <- list()
  context$Intervals <- list()
  context$Children <- list()
  context$Parent <- ScidrCurrentInference()
  context$PVar <- ScidrInferenceValue(frame, 'PVar', 'PValue')
  context$Pthresh <- ScidrInferenceValue(frame, 'Pthresh', 0.05)
  if (!is.character(context$PVar) || length(context$PVar) != 1L || is.na(context$PVar)) stop('PVar must name one p-value column.', call. = FALSE)
  if (!is.null(context$Parent)) {
    if (!'PVar' %in% names(call)) context$PVar <- context$Parent$PVar
    if (!'Pthresh' %in% names(call)) context$Pthresh <- context$Parent$Pthresh
  }
  if (!'Pthresh' %in% names(call)) {
    for (name in intersect(c('Alpha', 'FDRAlpha', 'label_q', 'adjusted_significance_threshold'), names(call))) {
      context$Pthresh <- ScidrInferenceValue(frame, name, context$Pthresh)
      break
    }
  }
  if (context$PVar %in% c('P', 'p', 'p.value', 'pval')) context$PVar <- 'PValue'
  if (context$PVar %in% c('FDR', 'p_adj', 'P_adj', 'p.adj', 'pval.adj', 'PAdjusted')) context$PVar <- 'PValueFDR'
  if (!is.character(context$PVar) || length(context$PVar) != 1L ||
      !context$PVar %in% c('PValue', 'PValueFDR', 'P', 'p', 'p.value', 'p_adj', 'P_adj', 'FDR', 'AdjustedPValue')) {
    stop('PVar must identify a raw or adjusted p-value column.', call. = FALSE)
  }
  if (!is.numeric(context$Pthresh) || length(context$Pthresh) != 1L ||
      is.na(context$Pthresh) || context$Pthresh <= 0 || context$Pthresh >= 1) {
    stop('Pthresh must be a single number strictly between 0 and 1.', call. = FALSE)
  }
  context$InputAxes <- ScidrInferenceAxes(context, frame)
  assign('PVar', context$PVar, envir = frame)
  assign('Pthresh', context$Pthresh, envir = frame)
  context
}

ScidrInferenceP <- function(raw, adjusted = NULL) {
  context <- ScidrCurrentInference()
  if (!is.null(context) && context$PVar %in% c('PValueFDR', 'FDR', 'P_adj', 'p_adj', 'AdjustedPValue')) {
    if (is.null(adjusted) || (any(is.finite(raw)) && all(is.na(adjusted[is.finite(raw)])))) return(stats::p.adjust(raw, method = 'BH'))
    return(adjusted)
  }
  raw
}

ScidrPairwiseDisplay <- function(df_Results, scope, star_p, IncludeTrendSymbol) {
  if (!nrow(df_Results)) return(df_Results)
  # BH is independent of any legacy Holm/Bonferroni/uncorrected column.
  df_Results$PValueFDR <- .ApplyTidyPAdjustment(df_Results, p_col = 'PValue',
    adjust_scope = if (scope == 'none') 'matrix' else scope, method = 'BH')
  context <- ScidrCurrentInference()
  selected <- if (star_p == 'none') rep(NA_real_, nrow(df_Results)) else if (star_p == 'raw') df_Results$PValue else if (context$PVar == 'AdjustedPValue') df_Results$AdjustedPValue else df_Results$PValueFDR
  df_Results$SignificanceLabel <- .FormatPValueStars(selected, IncludeTrendSymbol = IncludeTrendSymbol)
  # Conventional star tiers are retained; the requested cutoff controls inclusion.
  df_Results$SignificanceLabel[is.na(selected) | selected >= context$Pthresh] <- ''
  df_Results
}

ScidrStyleComparisonInference <- function(output, results, context, frame) {
  if (!inherits(output, 'gtsummary') || !'.is_main_row' %in% names(output$table_body)) return(output)
  if (!'TestRole' %in% names(results)) results$TestRole <- rep('Omnibus', nrow(results))
  df_Main <- results[results$TestRole %in% 'Omnibus', , drop = FALSE]
  columns <- 'p.value_fmt'
  map <- ScidrInferenceValue(frame, 'map_cols')
  if (is.data.frame(map)) columns <- c(columns, map$col_safe)
  for (column in intersect(columns, names(output$table_body))) {
    output <- gtsummary::modify_table_styling(output, columns = tidyselect::all_of(column), undo_text_format = 'bold')
    df_Selected <- if (column == 'p.value_fmt') df_Main else results[results$TestRole %in% 'Pairwise' & results$Term %in% map$contrast_label[map$col_safe == column], , drop = FALSE]
    vars_Significant <- df_Selected$Outcome[df_Selected$Significant]
    rows <- rlang::new_quosure(rlang::expr(.data$.is_main_row & .data$variable %in% !!vars_Significant), baseenv())
    output <- gtsummary::modify_table_styling(output, columns = tidyselect::all_of(column), rows = !!rows, text_format = 'bold')
  }
  output
}

ScidrSignificantPairs <- function(object, p_var, threshold, context, frame) {
  df_Results <- object$Results
  if (!is.data.frame(df_Results)) {
    df_Results <- object$Unadjusted$plot$data
    df_Results <- ScidrNormalizeInference(df_Results, context, frame)
  }
  p_column <- if (p_var %in% c('p', 'P', 'p.value')) 'PValue' else if (p_var %in% c('p_adj', 'P_adj', 'FDR', 'p.adj')) 'PValueFDR' else p_var
  if (!p_column %in% names(df_Results)) stop('Requested PVar is unavailable: ', p_var, call. = FALSE)
  df_Results <- df_Results[!is.na(df_Results[[p_column]]) & df_Results[[p_column]] < threshold, , drop = FALSE]
  df_Results$XVar <- df_Results$Predictor
  df_Results$YVar <- df_Results$Outcome
  df_Results$CategoricalVariable <- df_Results$Predictor
  df_Results$ContinuousVariable <- df_Results$Outcome
  df_Results$R <- df_Results$Estimate
  df_Results$P <- df_Results$PValue
  df_Results$P_adj <- df_Results$PValueFDR
  df_Results$nPairs <- df_Results$N
  df_Results
}

ScidrRecordInferenceFailure <- function(outcome, predictor, reason, n = NA_real_, stage = 'estimation') {
  context <- ScidrCurrentInference()
  if (is.null(context)) return(invisible(NULL))
  context$Failures[[length(context$Failures) + 1L]] <- tibble::tibble(
    Outcome = outcome, Predictor = predictor, Status = 'failed', Reason = reason, N = n, Stage = stage)
  invisible(NULL)
}

ScidrDropUnavailableOutcomes <- function(data, outcome_vars, predictor_vars, covariates = NULL, excluded_outcomes = character()) {
  vars_Unavailable <- outcome_vars[vapply(outcome_vars, function(variable) {
    x <- data[[variable]]
    supported <- is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x)
    supported && !variable %in% excluded_outcomes && length(unique(stats::na.omit(x))) < 2L
  }, logical(1))]
  for (outcome in vars_Unavailable) for (predictor in predictor_vars) {
    n <- sum(stats::complete.cases(data[, unique(c(outcome, predictor, covariates)), drop = FALSE]))
    ScidrRecordInferenceFailure(outcome, predictor, 'Outcome has fewer than two usable values.', n)
  }
  setdiff(outcome_vars, vars_Unavailable)
}

ScidrTryInference <- function(expr, fallback = NULL) {
  context <- ScidrCurrentInference()
  n_failures <- if (is.null(context)) 0L else length(context$Failures)
  value <- tryCatch(force(expr), error = function(e) {
    if (is.null(context)) stop(e)
    pair <- ScidrInferencePair(context, parent.frame())
    axes <- ScidrInferenceAxes(context, context$Frame)
    if (is.na(pair$Outcome) && length(axes$Outcome) == 1L) pair$Outcome <- axes$Outcome
    if (is.na(pair$Predictor) && length(axes$Predictor) == 1L) pair$Predictor <- axes$Predictor
    ScidrRecordInferenceFailure(pair$Outcome, pair$Predictor, conditionMessage(e))
    fallback
  })
  if (is.null(value) && !is.null(context) && length(context$Failures) == n_failures) {
    pair <- ScidrInferencePair(context, parent.frame())
    ScidrRecordInferenceFailure(pair$Outcome, pair$Predictor, 'The calculation returned no usable result.')
  }
  value
}

ScidrFailRequestedInference <- function(function_name, reason, group_var = NULL) {
  context <- ScidrCurrentInference()
  frame <- context$Frame
  axes <- ScidrInferenceAxes(context, frame)
  df_Data <- ScidrInferenceData(frame)
  ScidrValidateVariables(df_Data, axes$Outcome, 'outcome_vars')
  ScidrValidateVariables(df_Data, axes$Predictor, 'predictor_vars')
  covariates <- ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars', 'covars'))
  ScidrValidateVariables(df_Data, covariates, 'covariates')
  for (outcome in axes$Outcome) for (predictor in axes$Predictor) {
    if (outcome == predictor) next
    variables <- unique(c(outcome, predictor, group_var, covariates))
    n <- sum(stats::complete.cases(df_Data[, variables, drop = FALSE]))
    ScidrRecordInferenceFailure(outcome, predictor, reason, n)
  }
  ScidrEmptyLegacyInference(function_name)
}

ScidrRecordInferenceIntervals <- function(model, intervals, standardized_coefficients = NULL) {
  context <- ScidrCurrentInference()
  if (!is.null(context)) context$Intervals[[length(context$Intervals) + 1L]] <- list(Model = model,
    Intervals = intervals, StandardizedCoefficients = standardized_coefficients)
  invisible(NULL)
}

ScidrStoreInferenceObservation <- function(context, observation) {
  context$ObservationCount <- context$ObservationCount + 1L
  assign(as.character(context$ObservationCount), observation, envir = context$ObservationStore)
  invisible(NULL)
}

ScidrIndexInference <- function(context) {
  keys <- as.character(seq_len(context$ObservationCount))
  context$Observations <- if (length(keys)) unname(mget(keys, envir = context$ObservationStore, inherits = FALSE)) else list()
  context$ByPair <- new.env(hash = TRUE, parent = emptyenv())
  context$ByOutcome <- new.env(hash = TRUE, parent = emptyenv())
  if (!length(keys)) return(invisible(NULL))
  outcomes <- vapply(context$Observations, `[[`, character(1), 'Outcome')
  predictors <- vapply(context$Observations, `[[`, character(1), 'Predictor')
  for (index in split(seq_along(keys), paste(outcomes, predictors, sep = '\r'))) {
    observation <- context$Observations[[index[1]]]
    assign(paste(observation$Outcome, observation$Predictor, sep = '\r'), context$Observations[index], envir = context$ByPair)
  }
  for (index in split(seq_along(keys), outcomes)) assign(outcomes[index[1]], context$Observations[index], envir = context$ByOutcome)
  invisible(NULL)
}

ScidrCurrentInference <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (exists('.scidr_inference', envir = frames[[i]], inherits = FALSE)) {
      return(get('.scidr_inference', envir = frames[[i]], inherits = FALSE))
    }
  }
  NULL
}

ScidrInferenceLocal <- function(context, candidates) {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    value <- ScidrInferenceValue(frames[[i]], candidates)
    if (!is.null(value)) return(value)
    if (identical(frames[[i]], context$Frame)) break
  }
  NULL
}

ScidrInferencePair <- function(context, frame) {
  Find <- function(candidates) {
    value <- ScidrInferenceValue(frame, candidates)
    if (is.character(value) && length(value) == 1L) return(value)
    NA_character_
  }
  outcome <- Find(c('YVar', 'yVar', 'outcomeVar', 'outcome', 'outcome_var', 'yname', 'vy', 'trait_name', 'metab'))
  predictor <- Find(c('xVar', 'predictorVar', 'predictor', 'this_var', 'xname', 'vx', 'module_name', 'var', 'v', 'comp'))
  keys <- ScidrInferenceValue(frame, c('keys', '.y'))
  if (is.data.frame(keys) && nrow(keys) == 1L) {
    for (name in c('ContinuousVariable', 'YVar')) if (name %in% names(keys)) outcome <- as.character(keys[[name]][1])
    for (name in c('CategoricalVariable', 'XVar')) if (name %in% names(keys)) predictor <- as.character(keys[[name]][1])
  }
  if (context$Function %in% c('MakeComparisonTable', 'MakePairwiseHeatmap', 'PlotPairwiseMiningMatrix', 'PlotSplitViolin')) {
    value <- ScidrInferenceValue(frame, c('var', 'v', 'variable', 'var_nm'))
    if (is.character(value) && length(value) == 1L) outcome <- value
    value <- ScidrInferenceValue(context$Frame, c('group_var', 'CompVariable'))
    if (is.character(value) && length(value) == 1L) predictor <- value
  }
  if (is.na(outcome) || is.na(predictor)) {
    frames <- sys.frames()
    for (i in rev(seq_along(frames))) {
      if (identical(frames[[i]], context$Frame)) break
      if (is.na(outcome)) {
        value <- ScidrInferenceValue(frames[[i]], c('YVar', 'yVar', 'outcomeVar', 'outcome', 'outcome_var', 'yname', 'vy', 'trait_name'))
        if (is.character(value) && length(value) == 1L) outcome <- value
      }
      if (is.na(predictor)) {
        value <- ScidrInferenceValue(frames[[i]], c('xVar', 'predictorVar', 'predictor', 'this_var', 'xname', 'vx', 'module_name', 'var', 'v'))
        if (is.character(value) && length(value) == 1L) predictor <- value
      }
    }
  }
  if (context$Function %in% c('MakeComparisonTable', 'MakePairwiseHeatmap', 'PlotPairwiseMiningMatrix', 'PlotZScore')) {
    value <- ScidrInferenceLocal(context, c('var', 'v', 'variable'))
    if (is.character(value) && length(value) == 1L) outcome <- value
    value <- ScidrInferenceValue(context$Frame, c('group_var', 'CompVariable', 'TargetVar'))
    if (is.character(value) && length(value) == 1L) predictor <- value
  }
  if (context$Function %in% c('plotSigAssociations', 'plotSigCorrelations')) {
    df_Pairs <- ScidrInferenceValue(context$Frame, 'SigVarCombos')
    index <- ScidrInferenceValue(context$Frame, 'i')
    if (is.data.frame(df_Pairs) && length(index) == 1L && index <= nrow(df_Pairs)) {
      outcome <- df_Pairs$Outcome[index]
      predictor <- df_Pairs$Predictor[index]
    }
  }
  list(Outcome = outcome, Predictor = predictor)
}

ScidrObserveInference <- function(expr) {
  context <- ScidrCurrentInference()
  if (is.null(context)) return(force(expr))
  frame <- parent.frame()
  pair <- ScidrInferencePair(context, frame)
  axes <- context$ResolvedAxes
  if (is.null(axes)) {
    axes <- ScidrInferenceAxes(context, context$Frame)
    context$ResolvedAxes <- axes
  }
  if (is.na(pair$Outcome) && length(axes$Outcome) == 1L) pair$Outcome <- axes$Outcome
  if (is.na(pair$Predictor) && length(axes$Predictor) == 1L) pair$Predictor <- axes$Predictor
  expression <- substitute(expr)
  observation <- c(pair, list(Expression = paste(deparse(expression), collapse = ' '),
    Engine = paste(deparse(expression[[1]]), collapse = ''), Object = NULL,
    Reason = NA_character_, Status = 'ok'))
  value <- tryCatch(force(expr), error = function(e) {
    observation$Status <- 'failed'
    observation$Reason <- conditionMessage(e)
    ScidrStoreInferenceObservation(context, observation)
    stop(e)
  })
  observation$Object <- value
  observation$Formula <- tryCatch(paste(deparse(stats::formula(value)), collapse = ' '), error = function(e) NA_character_)
  if (inherits(value, c('lm', 'glm', 'coxph'))) {
    formula_vars <- all.vars(stats::formula(value))
    outcomes <- intersect(formula_vars[1], axes$Outcome)
    predictors <- intersect(formula_vars[-1], axes$Predictor)
    if (length(outcomes) == 1L) observation$Outcome <- outcomes
    if (length(predictors) == 1L) observation$Predictor <- predictors
  }
  if (is.na(observation$Outcome) && length(axes$Outcome) == 1L) observation$Outcome <- axes$Outcome
  if (is.na(observation$Predictor) && length(axes$Predictor) == 1L) observation$Predictor <- axes$Predictor
  if (inherits(value, 'lm') || inherits(value, 'glm') || inherits(value, 'coxph')) {
    observation$N <- tryCatch(stats::nobs(value), error = function(e) NA_real_)
    observation$ReferenceLevels <- value$xlevels
  }
  if (inherits(value, 'htest')) {
    x <- tryCatch(eval(expression[[2]], frame), error = function(e) NULL)
    y <- if (length(expression) >= 3L) tryCatch(eval(expression[[3]], frame), error = function(e) NULL) else NULL
    if (is.atomic(x) && !is.null(x) && is.atomic(y) && length(x) == length(y)) observation$N <- sum(stats::complete.cases(x, y))
    if (is.table(x)) observation$N <- sum(x)
    pair_levels <- ScidrInferenceLocal(context, 'cp')
    if (context$Function == 'MakeComparisonTable' && is.character(pair_levels) && length(pair_levels) == 2L) {
      estimates <- value$estimate
      estimate <- if (length(estimates) == 2L) estimates[1] - estimates[2] else if (length(estimates) == 1L) estimates[1] else NA_real_
      ci <- if (length(value$conf.int) == 2L) as.numeric(value$conf.int) else c(NA_real_, NA_real_)
      df_Subset <- ScidrInferenceLocal(context, 'sub')
      context$Estimates[[length(context$Estimates) + 1L]] <- tibble::tibble(
        Outcome = observation$Outcome, Predictor = observation$Predictor,
        Term = paste(pair_levels, collapse = ' - '), ReferenceLevel = pair_levels[2], ComparisonLevel = pair_levels[1],
        Estimate = as.numeric(estimate), ConfLow = ci[1], ConfHigh = ci[2], PValue = value$p.value,
        N = if (is.data.frame(df_Subset)) nrow(df_Subset) else NA_real_, TestRole = 'Pairwise',
        EffectType = if (grepl('Fisher', value$method)) 'Odds ratio' else 'Difference', Test = value$method)
    }
  }
  if (inherits(value, 'emmGrid')) {
    df_Contrasts <- tryCatch(as.data.frame(summary(value, infer = c(TRUE, TRUE), adjust = 'none')), error = function(e) NULL)
    if (is.data.frame(df_Contrasts) && all(c('contrast', 'estimate', 'p.value') %in% names(df_Contrasts))) {
      fit <- ScidrInferenceValue(frame, c('fit', 'model'))
      n <- if (inherits(fit, c('lm', 'glm'))) stats::nobs(fit) else NA_real_
      context$Estimates[[length(context$Estimates) + 1L]] <- tibble::tibble(
        Outcome = observation$Outcome, Predictor = observation$Predictor, Term = df_Contrasts$contrast,
        Estimate = df_Contrasts$estimate, ConfLow = ScidrInferenceColumn(df_Contrasts, c('lower.CL', 'asymp.LCL')),
        ConfHigh = ScidrInferenceColumn(df_Contrasts, c('upper.CL', 'asymp.UCL')), PValue = df_Contrasts$p.value,
        N = n, TestRole = 'Pairwise', EffectType = 'Mean difference')
    }
  }
  # Scaling attributes survive in the analysis frame before model.frame strips them.
  observation$Standardization <- list()
  for (name in c('ModelData', 'model_data', 'df_model', 'df_Model', 'DataSub')) {
    candidate <- ScidrInferenceValue(frame, name)
    if (!is.data.frame(candidate)) next
    for (variable in names(candidate)) {
      sd_used <- attr(candidate[[variable]], 'scaled:scale', exact = TRUE)
      if (!is.null(sd_used)) observation$Standardization[[variable]] <- list(
        SD = as.numeric(sd_used), Center = as.numeric(attr(candidate[[variable]], 'scaled:center', exact = TRUE)), N = nrow(candidate))
    }
  }
  scale_used <- attr(value, 'scaled:scale', exact = TRUE)
  if (!is.null(scale_used)) {
    centers <- attr(value, 'scaled:center', exact = TRUE)
    variables <- names(scale_used)
    if (is.null(variables)) {
      variable <- ScidrInferenceValue(frame, c('nm', 'variable', 'col'))
      if (is.null(variable)) variable <- if (grepl('outcome|YVar|scidr_y', paste(deparse(expression), collapse = ' '))) observation$Outcome else observation$Predictor
      variables <- rep(variable, length(scale_used))
    }
    for (i in seq_along(scale_used)) if (!is.na(variables[i])) observation$Standardization[[variables[i]]] <- list(
      SD = as.numeric(scale_used[i]), Center = as.numeric(centers[i]), N = if (is.matrix(value)) sum(!is.na(value[, i])) else sum(!is.na(value)))
  }
  # These engines explicitly convert scale() output to a vector.
  for (item in if (is.null(scale_used)) list(c('predictor_sd', 'this_var'), c('outcome_sd', 'outcome_var'), c('ref_sd', 'v')) else list()) {
    sd_used <- ScidrInferenceValue(frame, item[1])
    variable <- ScidrInferenceValue(frame, item[2])
    if (is.numeric(sd_used) && length(sd_used) == 1L && is.character(variable) && length(variable) == 1L) {
      if (item[1] == 'ref_sd' && !inherits(value, 'lm')) next
      ref_values <- ScidrInferenceValue(frame, 'ref_values')
      observation$Standardization[[variable]] <- list(SD = sd_used,
        Center = if (item[1] == 'ref_sd') ScidrInferenceValue(frame, 'ref_mean', NA_real_) else NA_real_,
        N = if (item[1] == 'ref_sd' && is.numeric(ref_values)) sum(!is.na(ref_values)) else ScidrInferenceValue(frame, c('final_n'), NA_real_))
    }
  }
  if (inherits(value, 'htest')) observation$Object <- structure(value[intersect(names(value),
    c('estimate', 'p.value', 'conf.int', 'method', 'parameter', 'statistic'))], class = 'htest')
  ScidrStoreInferenceObservation(context, observation)
  value
}

ScidrEmptyInferenceResults <- function() {
  tibble::tibble(Outcome = character(), Predictor = character(), OutcomeLabel = character(),
    PredictorLabel = character(), Group = character(), Term = character(), Estimate = double(),
    ConfLow = double(), ConfHigh = double(), PValue = double(), PValueFDR = double(), N = double(),
    TermLabel = character(), EffectType = character(), ReferenceValue = double(), Significant = logical(),
    Status = character(), Reason = character())
}

ScidrEmptyLegacyInference <- function(function_name) {
  df_Results <- ScidrEmptyInferenceResults()
  blank <- ggplot2::ggplot() + ggplot2::theme_void()
  table <- if (requireNamespace('gt', quietly = TRUE)) gt::gt(tibble::tibble(Outcome = character())) else NULL
  context <- ScidrCurrentInference()
  frame <- if (is.null(context)) parent.frame() else context$Frame
  relabel <- ScidrInferenceValue(frame, 'Relabel', TRUE)
  covariates <- ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars'))
  switch(function_name,
    MultivariableRegressionTable = list(Models = list(), FormattedTable = table, LargeTable = data.frame(),
      RegressionMatrix = data.frame(), VariableImportanceMatrix = data.frame(), Predictions = data.frame(),
      Diagnostics = data.frame(), ModelSummary = data.frame(), Multicollinearity = list(), Plots = list(),
      Results = df_Results, Metadata = list()),
    DiagnosticLikelihoodRatioTable = list(FormattedTable = table, LargeTable = table, BinaryFormattedTable = NULL,
      BinaryLargeTable = NULL, Results = df_Results, BinarySummary = tibble::tibble(), Metadata = list()),
    EvaluateBiomarkerPerformance = list(Models = list(), PerformanceTable = tibble::tibble(), RegressionTable = tibble::tibble(),
      ThresholdTable = tibble::tibble(), RawBiomarkerThreshold = NULL, Predictions = tibble::tibble(), ROCData = tibble::tibble(),
      CalibrationData = tibble::tibble(), Validation = tibble::tibble(), ValidationPredictions = tibble::tibble(),
      Plots = list(ROC = NULL, Probability = NULL, Calibration = NULL, CalibrationGrouped = NULL,
        BiomarkerDistribution = NULL, PredictedObserved = NULL, BiomarkerRelationship = NULL, Residuals = NULL),
      Results = df_Results, Metadata = list()),
    PlotBetaProfile = list(Plot = blank, ResultsTable = tibble::tibble(), Results = df_Results, Metadata = list()),
    ScreenBiomarkerPerformance = list(PerformanceTable = tibble::tibble(), RegressionTable = tibble::tibble(),
      ThresholdTable = tibble::tibble(), FailureTable = tibble::tibble(), AxisOrder = list(), Evaluations = list(),
      Plots = list(Heatmap = blank, AUC = blank, IncrementalAUC = blank, R2 = blank,
        IncrementalR2 = blank, BiomarkerPanels = list(), ROCFacets = blank), Metadata = list()),
    PlotPointCorrelationsHeatmap = {
      branch <- list(PvalTable = data.frame(), plot = blank)
      list(Unadjusted = branch, FDRCorrected = branch, method = 'R_pb', Relabel = relabel,
        AxisOrder = list(), Covariates = covariates, BinaryMapping = ScidrInferenceValue(frame, 'binary_map'), p = branch, p_fdr = branch)
    },
    PlotPhiHeatmap = {
      branch <- list(PvalTable = data.frame(), plot = blank)
      list(Unadjusted = branch, FDRCorrected = branch, method = 'Phi', Relabel = relabel,
        AxisOrder = list(), Triangle = ScidrInferenceValue(frame, 'triangle', 'full')[1], BinaryMapping = ScidrInferenceValue(frame, 'binary_map'), p = branch, p_fdr = branch)
    },
    Plot2GroupStats = list(plot = blank, table = NULL, pvaltable = data.frame(), data_used = tibble::tibble()),
    MakePairwiseHeatmap = list(Plot = blank, Results = df_Results, Models = list(), Settings = list(),
      ScalingParameters = data.frame(), Warnings = character(), AxisOrder = list()),
    PlotPairwiseMiningMatrix = list(Plots = list(Continuous = blank, Categorical = blank), Results = df_Results,
      ContinuousAudit = df_Results, CategoryLevelResults = tibble::tibble(), Settings = list(), Models = list(), Warnings = character()),
    PlotSplitCorrelationHeatmap = {
      branch <- list(r = list(), p = list(), npairs = list(), plot = blank)
      list(Unadjusted = branch, FDRCorrected = branch, Correlations = list(), Results = df_Results,
        Metadata = list(), p = branch, p_fdr = branch)
    },
    PlotCorrelationComparisons = {
      branch <- list(delta_r = matrix(numeric(), 0L, 0L), cohens_q = matrix(numeric(), 0L, 0L),
        z = matrix(numeric(), 0L, 0L), p = matrix(numeric(), 0L, 0L),
        r_comparison = matrix(numeric(), 0L, 0L), r_reference = matrix(numeric(), 0L, 0L),
        n_comparison = matrix(numeric(), 0L, 0L), n_reference = matrix(numeric(), 0L, 0L), plot = blank)
      list(Correlations = list(Comparison = NULL, Reference = NULL), Unadjusted = branch, FDRCorrected = branch,
        Results = df_Results, AxisOrder = list(), DirectionReversal = matrix(logical(), 0L, 0L), Metadata = list(), p = branch, p_fdr = branch)
    },
    PlotBlandAltman = list(plot = blank, stats = NULL),
    PlotPartialRegressionScatter = list(plot = blank, Relabel = relabel, Covariates = covariates,
      r2 = NA_real_, p_value = NA_real_, n = 0L, equation = NA_character_),
    blank)
}

ScidrInferenceData <- function(frame) {
  ScidrInferenceValue(frame, c('data', 'Data', 'DataFrame', 'trait_data'), data.frame())
}

ScidrInferenceAxes <- function(context, frame) {
  df_Data <- ScidrInferenceData(frame)
  vars_Outcomes <- ScidrInferenceValue(frame, c('outcome_vars', 'OutcomeVars', 'outcomeVars', 'ContVars', 'yVars', 'event_vars', 'outcome_var', 'OutcomeVar', 'outcomeVar', 'DepVar', 'Variable2', 'var_nm'))
  vars_Predictors <- ScidrInferenceValue(frame, c('predictor_vars', 'PredictorVars', 'predictorVars', 'CatVars', 'xVars', 'biomarker_vars', 'biomarker_var', 'predictor_var', 'predictorVar', 'IndepVar', 'Variable1'))
  group_var <- ScidrInferenceValue(frame, c('group_var', 'CompVariable', 'TargetVar', 'GroupVariable'))
  vars_Comparisons <- ScidrInferenceValue(frame, c('variables', 'Variables', 'var_nm', 'Var'))
  if (is.symbol(vars_Comparisons)) vars_Comparisons <- as.character(vars_Comparisons)
  if (is.null(group_var)) {
    group <- ScidrInferenceValue(frame, 'Group')
    if (is.symbol(group)) group_var <- as.character(group)
  }
  if (is.null(vars_Outcomes) && !is.null(vars_Comparisons)) vars_Outcomes <- vars_Comparisons
  if (is.null(vars_Predictors) && !is.null(group_var)) vars_Predictors <- group_var
  if (context$Function == 'CreateStatisticsTable') vars_Outcomes <- setdiff(names(df_Data), group_var)
  if (context$Function == 'PlotPValueComparisons' && is.null(vars_Comparisons)) vars_Outcomes <- setdiff(names(df_Data), group_var)
  if (context$Function %in% c('PlotSplitCorrelationHeatmap', 'PlotCorrelationComparisons') &&
      is.null(ScidrInferenceValue(frame, 'outcome_vars')) && is.null(ScidrInferenceValue(frame, 'predictor_vars')) && is.data.frame(df_Data)) {
    vars_Outcomes <- names(df_Data)[vapply(df_Data, is.numeric, logical(1))]
    vars_Predictors <- vars_Outcomes
  }
  if (context$Function == 'PlotAssociations') {
    vars_Outcomes <- ScidrInferenceValue(frame, 'Var2')
    vars_Predictors <- ScidrInferenceValue(frame, 'Var1')
  }
  if (context$Function == 'calculate_pathway_results') {
    vars_Outcomes <- ScidrInferenceValue(frame, c('metabolites', 'metabolite_vars'))
    vars_Predictors <- ScidrInferenceValue(frame, 'comparison_var')
  }
  if (is.null(vars_Predictors) && context$Function %in% c('PlotMiningMatrix', 'PlotNumInteractionEffectsMatrix')) vars_Predictors <- vars_Outcomes
  if (is.null(vars_Outcomes) && !is.null(vars_Predictors)) vars_Outcomes <- vars_Predictors
  if (is.null(vars_Outcomes) && context$Function == 'CreateStatisticsTable') vars_Outcomes <- setdiff(names(df_Data), group_var)
  list(Outcome = as.character(vars_Outcomes), Predictor = as.character(vars_Predictors))
}

ScidrInferenceColumn <- function(data, candidates, default = NA_real_) {
  out <- rep(if (is.character(default)) NA_character_ else NA_real_, nrow(data))
  for (name in candidates) {
    if (!name %in% names(data) || is.list(data[[name]])) next
    value <- if (is.character(default)) as.character(data[[name]]) else suppressWarnings(as.numeric(data[[name]]))
    use <- is.na(out) & !is.na(value)
    out[use] <- value[use]
  }
  out[is.na(out)] <- default
  out
}

ScidrNormalizeInference <- function(data, context, frame) {
  if (!is.data.frame(data) || !nrow(data)) return(ScidrEmptyInferenceResults())
  df_Results <- tibble::as_tibble(data)
  axes <- ScidrInferenceAxes(context, frame)
  single_outcome <- if (length(axes$Outcome) == 1L) axes$Outcome else NA_character_
  single_predictor <- if (length(axes$Predictor) == 1L) axes$Predictor else NA_character_
  outcome_names <- c('Outcome', 'ContinuousVariable', 'YVar', 'Y', 'Trait')
  predictor_names <- c('Predictor', 'CategoricalVariable', 'XVar', 'X', 'Variable', 'variable', 'Module', 'metabolite')
  if (context$Function == 'PlotMiningMatrix') {
    outcome_names <- c('Outcome', 'XVar')
    predictor_names <- c('Predictor', 'YVar')
  }
  if (context$Function == 'calculate_pathway_results') {
    outcome_names <- c('Outcome', 'Metabolite')
    predictor_names <- 'Predictor'
  }
  if (context$Function %in% c('MakeComparisonTable', 'MakeFacetCatComparisonTable', 'CreateStatisticsTable',
      'MakePairwiseHeatmap', 'PlotPairwiseMiningMatrix', 'Plot2GroupStats', 'PlotPValueComparisons', 'PlotZScore')) {
    outcome_names <- c('Outcome', 'Variable', 'variable')
    predictor_names <- 'Predictor'
  }
  Set <- function(name, candidates, default = NA_real_) {
    if (!name %in% names(df_Results)) df_Results[[name]] <<- ScidrInferenceColumn(df_Results, candidates, default)
  }
  Set('Outcome', outcome_names, NA_character_)
  Set('Predictor', predictor_names, NA_character_)
  df_Results$Outcome[is.na(df_Results$Outcome)] <- single_outcome
  df_Results$Predictor[is.na(df_Results$Predictor)] <- single_predictor
  Set('Group', c('Stratum', 'Split'), 'Pooled')
  Set('Term', c('term', 'ResultLevel', 'contrast', 'Contrast'), NA_character_)
  df_Results$Term[is.na(df_Results$Term)] <- df_Results$Predictor[is.na(df_Results$Term)]
  if (context$Function == 'MultivariableRegressionTable') {
    placeholder <- df_Results$Term %in% c('Continuous', 'continuous', '')
    df_Results$Term[placeholder] <- df_Results$Predictor[placeholder]
  }
  Set('Estimate', c('Beta', 'estimate', 'LikelihoodRatio', 'DeltaR', 'correlation', 'R', 'r', 'Phi', 'EffectSize', 'effect_size', 'Effect', 'C', 'EstimatedMeanDifference', 'EffectMagnitude', 'ges', 'eta2_H', 'cramers_v'))
  if (context$Function == 'MultivariableRegressionTable' && 'Effect' %in% names(df_Results)) df_Results$Estimate <- df_Results$Effect
  Set('ConfLow', c('CILow', 'LowerCI', 'conf.low', 'lower.CL', 'LRLowerCI', 'ConfidenceIntervalLower'))
  Set('ConfHigh', c('CIHigh', 'UpperCI', 'conf.high', 'upper.CL', 'LRUpperCI', 'ConfidenceIntervalUpper'))
  Set('PValue', c('p.value', 'p_value', 'pval', 'p', 'P', 'p_val'))
  Set('N', c('nPairs', 'n_pairs', 'npairs', 'SampleSize', 'n', 'n.obs'))
  if (all(c('n.obs.1', 'n.obs.2') %in% names(df_Results))) df_Results$N <- df_Results$n.obs.1 + df_Results$n.obs.2
  if (context$Function == 'PlotCorrelationComparisons' && all(c('NComparison', 'NReference') %in% names(df_Results))) df_Results$N <- df_Results$NComparison + df_Results$NReference
  Set('PValueFDR', c('FDR', 'p.adj', 'p_adj', 'P_adj', 'PAdjusted', 'P_FDR', 'q.value'))
  method <- ScidrInferenceValue(frame, c('p_adjust_method', 'AdjustMethod', 'adjust_method'), 'fdr')
  if (!method %in% c('fdr', 'BH')) df_Results$PValueFDR <- NA_real_
  df_Data <- ScidrInferenceData(frame)
  codebook <- ScidrInferenceValue(frame, c('codebook', 'Codebook'))
  display_labels <- if (ncol(df_Data)) ScidrDisplayLabels(df_Data, names(df_Data)) else character()
  if (is.data.frame(codebook) && all(c('Variable', 'Label') %in% names(codebook))) {
    valid <- !is.na(codebook$Variable) & !is.na(codebook$Label)
    display_labels[as.character(codebook$Variable[valid])] <- as.character(codebook$Label[valid])
  }
  Label <- function(variables) {
    labels <- unname(display_labels[variables])
    labels[is.na(labels)] <- variables[is.na(labels)]
    labels
  }
  Set('OutcomeLabel', character(), NA_character_)
  Set('PredictorLabel', character(), NA_character_)
  for (name in c('Outcome', 'Predictor', 'OutcomeLabel', 'PredictorLabel', 'Term')) {
    if (is.list(df_Results[[name]])) df_Results[[name]] <- vapply(df_Results[[name]], function(x) if (length(x)) as.character(x[1]) else NA_character_, character(1))
    df_Results[[name]] <- as.character(df_Results[[name]])
  }
  missing_labels <- is.na(df_Results$OutcomeLabel)
  df_Results$OutcomeLabel[missing_labels] <- Label(df_Results$Outcome[missing_labels])
  missing_labels <- is.na(df_Results$PredictorLabel)
  df_Results$PredictorLabel[missing_labels] <- Label(df_Results$Predictor[missing_labels])
  Set('TermLabel', character(), NA_character_)
  missing_labels <- is.na(df_Results$TermLabel)
  vars_SortedNames <- names(df_Data)[order(nchar(names(df_Data)), decreasing = TRUE)]
  RenderTerm <- function(term, predictor_label) {
    if (is.na(term)) return(predictor_label)
    parts <- strsplit(term, ':', fixed = TRUE)[[1]]
    vars_Names <- vars_SortedNames
    parts <- vapply(parts, function(part) {
      for (variable in vars_Names) {
        for (prefix in c(paste0('`', variable, '`'), variable)) {
          if (startsWith(part, prefix)) {
            suffix <- substring(part, nchar(prefix) + 1L)
            return(paste0(Label(variable), if (nzchar(suffix)) paste0(' : ', suffix) else ''))
          }
        }
      }
      part
    }, character(1))
    paste(parts, collapse = ' \u00d7 ')
  }
  indices <- which(missing_labels)
  keys <- paste(df_Results$Term[indices], df_Results$PredictorLabel[indices], sep = '\r')
  unique_indices <- indices[!duplicated(keys)]
  rendered <- vapply(unique_indices, function(i) RenderTerm(df_Results$Term[i], df_Results$PredictorLabel[i]), character(1))
  df_Results$TermLabel[missing_labels] <- rendered[match(keys, unique(keys))]
  Set('EffectType', c('Test', 'test'), NA_character_)
  if (context$Function == 'DiagnosticLikelihoodRatioTable') {
    df_Results$EffectType <- rep('Likelihood ratio', nrow(df_Results))
    if (is.null(ScidrInferenceValue(frame, 'stratify_by'))) df_Results$Group <- rep('Pooled', nrow(df_Results))
  }
  Set('ReferenceValue', character())
  ratio <- grepl('odds|hazard|likelihood', df_Results$EffectType, ignore.case = TRUE)
  df_Results$ReferenceValue[is.na(df_Results$ReferenceValue)] <- ifelse(ratio[is.na(df_Results$ReferenceValue)], 1, 0)
  Set('Status', character(), NA_character_)
  df_Results$Status[df_Results$Status %in% c('available', 'success')] <- 'ok'
  df_Results$Status[df_Results$Status %in% c('unavailable', 'error')] <- 'failed'
  # Missing p-values are legitimate for descriptive effects and penalized models.
  available <- !is.na(df_Results$Estimate) | is.finite(df_Results$PValue)
  df_Results$Status[is.na(df_Results$Status)] <- ifelse(available[is.na(df_Results$Status)], 'ok', 'failed')
  Set('Reason', c('Note', 'Error'), NA_character_)
  df_Results$Reason[df_Results$Status == 'failed' & is.na(df_Results$Reason)] <- 'Relationship could not be estimated from the available observations.'
  df_Results
}

ScidrInferenceModelRows <- function(observation) {
  model <- observation$Object
  if (!inherits(model, c('lm', 'glm', 'coxph', 'svycoxph'))) return(NULL)
  coefficients <- tryCatch(as.data.frame(summary(model)$coefficients), error = function(e) NULL)
  if (is.null(coefficients) || !nrow(coefficients)) return(NULL)
  terms <- rownames(coefficients)
  keep <- terms != '(Intercept)'
  coefficients <- coefficients[keep, , drop = FALSE]
  terms <- terms[keep]
  if (!length(terms)) return(NULL)
  estimate <- coefficients[[1]]
  raw_estimate <- estimate
  se_name <- intersect(c('Std. Error', 'se(coef)', 'robust se'), names(coefficients))
  se <- if (length(se_name)) coefficients[[tail(se_name, 1)]] else rep(NA_real_, length(estimate))
  p_name <- grep('^Pr\\(', names(coefficients), value = TRUE)
  p <- if (length(p_name)) coefficients[[p_name[1]]] else rep(NA_real_, length(estimate))
  critical <- if (inherits(model, 'glm') || inherits(model, 'coxph')) stats::qnorm(.975) else stats::qt(.975, stats::df.residual(model))
  ratio <- inherits(model, 'coxph') || (inherits(model, 'glm') && model$family$family == 'binomial')
  lower <- estimate - critical * se
  upper <- estimate + critical * se
  if (!is.null(observation$Intervals)) {
    index <- match(terms, rownames(observation$Intervals))
    lower <- observation$Intervals[index, 1]
    upper <- observation$Intervals[index, 2]
  }
  if (ratio) { estimate <- exp(estimate); lower <- exp(lower); upper <- exp(upper) }
  tibble::tibble(Outcome = observation$Outcome, Predictor = observation$Predictor,
    Term = terms, Estimate = estimate, RawEstimate = raw_estimate, StdError = se, ConfLow = lower, ConfHigh = upper,
    PValue = p, N = if (is.null(observation$N)) NA_real_ else observation$N,
    EffectType = if (inherits(model, 'coxph')) 'Hazard ratio' else if (ratio) 'Odds ratio' else 'Estimate',
    ReferenceValue = if (ratio) 1 else 0)
}

ScidrFullRegressionRows <- function(output, context, frame) {
  if (!is.list(output$Models) || !length(output$Models) || !is.data.frame(output$LargeTable)) return(NULL)
  df_Legacy <- output$LargeTable
  if (!all(names(output$Models) %in% unique(df_Legacy$Outcome))) return(NULL)
  vars_Predictors <- unique(c(ScidrInferenceValue(frame, c('predictor_vars', 'PredictorVars')),
    ScidrInferenceValue(frame, c('covariates', 'Covars'))))
  vars_Predictors <- vars_Predictors[order(nchar(vars_Predictors), decreasing = TRUE)]
  CoefficientRows <- function(outcome, coefficients, standard_errors = NULL, outcome_level = NA_character_, reference_level = NA_character_, ratio = TRUE, n = NA_real_) {
    keep <- !startsWith(names(coefficients), '(Intercept)')
    coefficients <- coefficients[keep]
    se <- if (is.null(standard_errors)) rep(NA_real_, length(coefficients)) else standard_errors[names(coefficients)]
    p <- if (is.null(standard_errors)) rep(NA_real_, length(coefficients)) else 2 * stats::pnorm(abs(coefficients / se), lower.tail = FALSE)
    lower <- coefficients - stats::qnorm(.975) * se
    upper <- coefficients + stats::qnorm(.975) * se
    tibble::tibble(Outcome = outcome, Predictor = NA_character_, Term = names(coefficients),
      RawEstimate = unname(coefficients), Estimate = if (ratio) exp(unname(coefficients)) else unname(coefficients),
      ConfLow = if (ratio) exp(unname(lower)) else unname(lower), ConfHigh = if (ratio) exp(unname(upper)) else unname(upper),
      StdError = unname(se), PValue = unname(p), N = n,
      OutcomeLevel = outcome_level, ReferenceLevel = reference_level,
      EffectType = if (ratio) 'Odds ratio' else 'Estimate', ReferenceValue = if (ratio) 1 else 0)
  }
  rows <- lapply(names(output$Models), function(outcome) {
    model <- output$Models[[outcome]]
    intervals <- Filter(function(item) identical(item$Model, model), context$Intervals)
    df_Outcome <- df_Legacy[df_Legacy$Outcome == outcome, , drop = FALSE]
    n <- if ('SampleSize' %in% names(df_Outcome)) df_Outcome$SampleSize[1] else NA_real_
    if (inherits(model, c('lm', 'glm'))) {
      observation <- list(Object = model, Outcome = outcome, Predictor = NA_character_, N = stats::nobs(model),
        Intervals = if (length(intervals)) intervals[[1]]$Intervals else NULL)
      df_Coefficients <- ScidrInferenceModelRows(observation)
    } else if (inherits(model, 'multinom')) {
      sm <- summary(model)
      df_Coefficients <- dplyr::bind_rows(lapply(rownames(sm$coefficients), function(level) {
        CoefficientRows(outcome, sm$coefficients[level, ], sm$standard.errors[level, ], level, model$lev[1], n = n)
      }))
    } else if (inherits(model, 'polr')) {
      coefficients <- stats::coef(model)
      se <- sqrt(diag(stats::vcov(model)))[names(coefficients)]
      df_Coefficients <- CoefficientRows(outcome, coefficients, se, 'Higher category', 'Lower category', n = n)
    } else if (inherits(model, 'cv.glmnet')) {
      lambda <- if ('Lambda' %in% names(df_Outcome)) df_Outcome$Lambda[1] else model$lambda.min
      coefficients <- stats::coef(model, s = lambda)
      if (is.list(coefficients)) {
        reference <- unique(stats::na.omit(df_Outcome$ReferenceLevel))[1]
        if (is.na(reference)) reference <- names(coefficients)[1]
        df_Coefficients <- dplyr::bind_rows(lapply(setdiff(names(coefficients), reference), function(level) {
          values <- as.matrix(coefficients[[level]])[, 1] - as.matrix(coefficients[[reference]])[, 1]
          CoefficientRows(outcome, values, outcome_level = level, reference_level = reference, n = n)
        }))
      } else {
        ratio <- any(grepl('odds', df_Outcome$EffectType, ignore.case = TRUE))
        df_Coefficients <- CoefficientRows(outcome, as.matrix(coefficients)[, 1], ratio = ratio, n = n)
      }
    } else if (inherits(model, 'ordinalNet')) {
      tuning <- output$Metadata$AnalysisSettings$Tuning[[outcome]]
      values <- stats::coef(model, whichLambda = tuning$BestLambdaIndex)
      df_Coefficients <- CoefficientRows(outcome, values, outcome_level = 'Higher category', reference_level = 'Lower category', n = n)
    } else return(NULL)
    if (is.null(df_Coefficients)) return(NULL)
    df_Coefficients$Predictor <- vapply(df_Coefficients$Term, function(term) {
      for (variable in vars_Predictors) if (startsWith(term, variable) || startsWith(term, paste0('`', variable, '`'))) return(variable)
      term
    }, character(1))
    if ('OutcomeLevel' %in% names(df_Coefficients) && any(!is.na(df_Coefficients$OutcomeLevel)) && 'OutcomeLevel' %in% names(df_Legacy)) {
      index <- match(paste(df_Coefficients$Outcome, df_Coefficients$Predictor, df_Coefficients$OutcomeLevel), paste(df_Legacy$Outcome, df_Legacy$Predictor, df_Legacy$OutcomeLevel))
    } else index <- match(paste(df_Coefficients$Outcome, df_Coefficients$Predictor), paste(df_Legacy$Outcome, df_Legacy$Predictor))
    df_Result <- df_Legacy[index, setdiff(names(df_Legacy), c('FDR', 'Stars', 'TermLabel')), drop = FALSE]
    for (name in names(df_Coefficients)) df_Result[[name]] <- df_Coefficients[[name]]
    df_Result$Effect <- df_Coefficients$Estimate
    df_Result$LowerCI <- df_Coefficients$ConfLow
    df_Result$UpperCI <- df_Coefficients$ConfHigh
    df_Result$StandardError <- df_Coefficients$StdError
    if (length(intervals) && length(intervals[[1]]$StandardizedCoefficients)) {
      standardized <- intervals[[1]]$StandardizedCoefficients
      if (is.matrix(standardized) && 'OutcomeLevel' %in% names(df_Coefficients)) {
        df_Result$StandardizedBeta <- standardized[cbind(match(df_Coefficients$OutcomeLevel, rownames(standardized)), match(df_Coefficients$Term, colnames(standardized)))]
      } else if (!is.list(standardized)) df_Result$StandardizedBeta <- unname(standardized[df_Coefficients$Term])
    } else if ('StandardizedBeta' %in% names(df_Result)) df_Result$StandardizedBeta <- NA_real_
    df_Result
  })
  if (any(vapply(rows, is.null, logical(1)))) return(NULL)
  dplyr::bind_rows(rows)
}

ScidrInferenceSource <- function(output, context, frame) {
  # This matrix computes explicit linear-combination rows, including simple
  # slopes; generic model coefficient extraction would discard those rows.
  if (context$Function == 'PlotInteractionEffectsMatrix' && is.data.frame(output$Results)) return(output$Results)
  if (context$Function == 'MultivariableRegressionTable') {
    if (is.data.frame(output$Results) && nrow(output$Results)) return(output$Results)
    df_Coefficients <- ScidrFullRegressionRows(output, context, frame)
    if (!is.null(df_Coefficients)) return(df_Coefficients)
  }
  if (context$Function == 'MakeComparisonTable' && isTRUE(ScidrInferenceValue(frame, 'overall_mode', FALSE))) return(ScidrEmptyInferenceResults())
  if (context$Function %in% c('plotSigAssociations', 'plotSigCorrelations')) return(ScidrInferenceValue(frame, 'SigVarCombos', ScidrEmptyInferenceResults()))
  if (context$Function == 'EvaluateBiomarkerPerformance') {
    fits <- ScidrInferenceValue(frame, 'fits', list())
    axes <- ScidrInferenceAxes(context, frame)
    df_Rows <- dplyr::bind_rows(lapply(names(fits), function(name) {
      fit <- fits[[name]]
      if (fit$Status == 'not_applicable') return(NULL)
      if (name == 'Adjusted' && !length(ScidrInferenceValue(frame, 'covariates'))) return(NULL)
      observation <- list(Object = fit$Model, Outcome = axes$Outcome[1], Predictor = axes$Predictor[1],
        N = ScidrInferenceValue(frame, 'df_analysis', data.frame()) %>% nrow())
      df_Row <- ScidrInferenceModelRows(observation)
      if (is.null(df_Row)) df_Row <- tibble::tibble(Outcome = axes$Outcome[1], Predictor = axes$Predictor[1], N = observation$N)
      df_Row$Model <- name
      df_Row$Status <- if (fit$Status == 'ok') 'ok' else 'failed'
      df_Row$Reason <- if (fit$Status == 'ok') NA_character_ else fit$Note
      if (fit$Status != 'ok') {
        for (column in intersect(c('Estimate', 'ConfLow', 'ConfHigh', 'PValue', 'StdError'), names(df_Row))) df_Row[[column]] <- NA_real_
      }
      df_Row
    }))
    if (nrow(df_Rows)) return(df_Rows)
  }
  if (context$Function == 'MakeFacetCatComparisonTable') {
    if (length(context$Children)) return(dplyr::bind_rows(lapply(context$Children, `[[`, 'Results')))
    subtables <- ScidrInferenceValue(frame, 'subtables', list())
    df_Results <- dplyr::bind_rows(lapply(subtables, function(table) if (is.data.frame(table$Results)) table$Results else NULL))
    if (nrow(df_Results)) return(df_Results)
    return(ScidrEmptyInferenceResults())
  }
  if (context$Function %in% c('PlotInteractionEffectsMatrix', 'PlotNumInteractionEffectsMatrix',
      'PlotCatInteractionEffectsMatrix', 'PlotInteractionEffectsContinuous', 'EvaluateBiomarkerPerformance', 'PlotPartialRegressionScatter')) {
    observations <- context$Observations
    if (context$Function == 'PlotPartialRegressionScatter') observations <- tail(observations, 1)
    rows <- dplyr::bind_rows(lapply(observations, ScidrInferenceModelRows))
    if (context$Function %in% c('PlotInteractionEffectsMatrix', 'PlotNumInteractionEffectsMatrix', 'PlotCatInteractionEffectsMatrix')) rows <- rows[grepl(':', rows$Term, fixed = TRUE), , drop = FALSE]
    if (nrow(rows)) return(rows)
  }
  if (is.list(output) && !inherits(output, 'ggplot')) {
    for (name in c('Results', 'ResultsTable', 'results')) if (is.data.frame(output[[name]])) return(output[[name]])
    if (is.data.frame(output$LargeTable)) return(output$LargeTable)
  }
  if (is.data.frame(output)) return(output)
  if (inherits(output, 'gtsummary')) {
    df_Table <- output$table_body
    if ('variable' %in% names(df_Table)) df_Table <- df_Table[df_Table$variable %in% names(ScidrInferenceData(frame)), , drop = FALSE]
    if ('.is_main_row' %in% names(df_Table)) df_Table <- df_Table[df_Table$.is_main_row %in% TRUE, , drop = FALSE]
    else if ('row_type' %in% names(df_Table)) df_Table <- df_Table[df_Table$row_type == 'label', , drop = FALSE]
    return(df_Table)
  }
  # Objects retained in the calculation frame precede presentation reshaping.
  for (name in c('stat.test', 'df_Results', 'results', 'results_table', 'pvaltable', 'pvaldata', 'plot.data', 'df_Combined')) {
    value <- ScidrInferenceValue(frame, name)
    if (is.data.frame(value) && any(c('PValue', 'p.value', 'p', 'P', 'p_value', 'pval', 'p_value.x') %in% names(value))) return(value)
  }
  if (context$Function %in% c('PlotInteractionEffectsMatrix', 'PlotNumInteractionEffectsMatrix',
      'PlotCatInteractionEffectsMatrix', 'PlotInteractionEffectsContinuous', 'EvaluateBiomarkerPerformance', 'PlotPartialRegressionScatter')) {
    rows <- lapply(context$Observations, ScidrInferenceModelRows)
    rows <- dplyr::bind_rows(rows)
    if (context$Function %in% c('PlotInteractionEffectsMatrix', 'PlotNumInteractionEffectsMatrix', 'PlotCatInteractionEffectsMatrix')) rows <- rows[grepl(':', rows$Term, fixed = TRUE), , drop = FALSE]
    if (nrow(rows)) return(rows)
  }
  if (context$Function == 'PlotAssociations' && inherits(output, 'ggplot')) {
    statistics <- tryCatch(ggstatsplot::extract_stats(output), error = function(e) NULL)
    if (is.data.frame(statistics$subtitle_data)) return(statistics$subtitle_data)
    if (is.data.frame(statistics$pairwise_comparisons_data)) {
      df_Pairs <- statistics$pairwise_comparisons_data
      df_Pairs$Term <- paste(df_Pairs$group1, df_Pairs$group2, sep = ' - ')
      df_Pairs$ReferenceLevel <- df_Pairs$group2
      df_Pairs$ComparisonLevel <- df_Pairs$group1
      df_Pairs$EffectType <- df_Pairs$effectsize
      return(df_Pairs)
    }
  }
  comparison <- attr(output, 'comparison', exact = TRUE)
  if (is.list(comparison)) {
    if (is.data.frame(comparison$PValueTable)) return(comparison$PValueTable)
    if (is.numeric(comparison$PValue)) return(tibble::tibble(PValue = comparison$PValue, N = comparison$N,
      EffectType = comparison$Method, Reason = if (is.null(comparison$Error)) NA_character_ else comparison$Error))
  }
  if (context$Function == 'PlotCorrelationsHeatmap') {
    mc <- ScidrInferenceValue(frame, 'MC')
    mp <- ScidrInferenceValue(frame, 'MP')
    mn <- ScidrInferenceValue(frame, 'MN')
    if (is.matrix(mc)) return(tibble::tibble(XVar = rep(rownames(mc), ncol(mc)),
      YVar = rep(colnames(mc), each = nrow(mc)), Estimate = as.vector(mc), PValue = as.vector(mp), N = as.vector(mn)))
  }
  if (context$Function == 'PlotBlandAltman' && is.list(output$stats)) {
    ba <- output$stats
    estimates <- c(ba$lower.limit, ba$lines[2], ba$upper.limit)
    ci <- ba$CI.lines
    return(tibble::tibble(Term = c('Lower limit of agreement', 'Bias', 'Upper limit of agreement'),
      Estimate = as.numeric(estimates), ConfLow = as.numeric(ci[c(1, 3, 5)]), ConfHigh = as.numeric(ci[c(2, 4, 6)]), N = length(ba$diffs)))
  }
  if (length(context$Children)) return(dplyr::bind_rows(lapply(context$Children, `[[`, 'Results')))
  ScidrEmptyInferenceResults()
}

ScidrInferenceRequested <- function(context, frame) {
  axes <- context$InputAxes
  if (is.null(axes) || !length(axes$Outcome) || !length(axes$Predictor)) axes <- ScidrInferenceAxes(context, frame)
  if (!length(axes$Outcome) || !length(axes$Predictor)) return(tibble::tibble())
  df_Requested <- tidyr::expand_grid(Outcome = axes$Outcome, Predictor = axes$Predictor)
  df_Requested <- df_Requested[df_Requested$Outcome != df_Requested$Predictor, , drop = FALSE]
  if (context$Function == 'MultivariableRegressionTable') {
    skipped <- ScidrInferenceValue(frame, 'skipped', character())
    df_Requested <- df_Requested[!df_Requested$Outcome %in% skipped, , drop = FALSE]
  }
  strata <- ScidrInferenceValue(frame, 'strata')
  if (is.data.frame(strata) && 'Stratum' %in% names(strata) && length(ScidrInferenceValue(frame, 'stratify_by'))) df_Requested <- tidyr::crossing(df_Requested, Group = strata$Stratum)
  if (context$Function == 'PlotSplitCorrelationHeatmap' && length(ScidrInferenceValue(frame, 'split_levels'))) df_Requested <- tidyr::crossing(df_Requested, Group = ScidrInferenceValue(frame, 'split_levels'))
  covariates <- ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars', 'covars'))
  if (context$Function %in% c('PlotMiningMatrix', 'PlotNumInteractionEffectsMatrix', 'PlotCatInteractionEffectsMatrix')) {
    df_Requested <- df_Requested[!df_Requested$Outcome %in% covariates & !df_Requested$Predictor %in% covariates, , drop = FALSE]
  }
  df_Requested
}

ScidrInferenceFDR <- function(results, frame) {
  if (!nrow(results)) return(results)
  scope <- ScidrInferenceValue(frame, c('fdr_scope', 'adjust_scope'), 'matrix')[1]
  groups <- switch(scope, per_outcome = results$Outcome, per_predictor = results$Predictor,
    per_variable = results$Outcome, per_group = results$Group, rep('matrix', nrow(results)))
  # Preserve existing within-stratum families.
  if (any(c('Split', 'Stratum') %in% names(results))) groups <- paste(results$Group, groups, sep = '\r')
  if ('TestRole' %in% names(results)) {
    pairwise <- results$TestRole %in% 'Pairwise'
    groups[pairwise] <- paste('pairwise', results$Outcome[pairwise], sep = '\r')
  }
  for (index in split(seq_len(nrow(results)), groups)) {
    ok <- index[is.finite(results$PValue[index])]
    if (length(ok) && all(is.na(results$PValueFDR[ok]))) results$PValueFDR[ok] <- stats::p.adjust(results$PValue[ok], method = 'BH')
  }
  results
}

ScidrAttachInference <- function(output, results, metadata) {
  if (inherits(output, 'ggplot')) {
    output$Results <- results
    output$Metadata <- metadata
  } else if (is.data.frame(output) || is.atomic(output) || (is.list(output) && is.null(names(output)))) {
    attr(output, 'Results') <- results
    attr(output, 'Metadata') <- metadata
    class(output) <- unique(c('SciDataReportRInference', class(output)))
  } else {
    output$Results <- results
    output$Metadata <- metadata
  }
  output
}

ScidrEnrichInference <- function(results, context, frame) {
  df_Data <- ScidrInferenceData(frame)
  covariates <- ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars', 'covars'))
  for (i in seq_len(nrow(results))) {
    key <- paste(results$Outcome[i], results$Predictor[i], sep = '\r')
    observations <- if (exists(key, envir = context$ByPair, inherits = FALSE)) get(key, envir = context$ByPair, inherits = FALSE) else list()
    for (observation in observations) {
      if (!identical(observation$Outcome, results$Outcome[i]) || !identical(observation$Predictor, results$Predictor[i])) next
      if (is.na(results$N[i]) && !is.null(observation$N)) results$N[i] <- observation$N
      test <- observation$Object
      if (inherits(test, 'htest') && !is.null(test$conf.int) && length(test$conf.int) == 2L &&
          is.finite(results$PValue[i]) && isTRUE(all.equal(test$p.value, results$PValue[i], tolerance = 1e-6))) {
        estimate <- if (length(test$estimate) == 2L) as.numeric(test$estimate[1] - test$estimate[2]) else if (length(test$estimate) == 1L) as.numeric(test$estimate[1]) else NA_real_
        if (is.na(results$Estimate[i]) && !is.na(estimate)) {
          results$Estimate[i] <- estimate
          results$EffectType[i] <- if (length(test$estimate) == 2L) 'Mean difference' else if (grepl('Fisher', test$method)) 'Odds ratio' else test$method
          if (grepl('Fisher', test$method)) results$ReferenceValue[i] <- 1
        }
        same_scale <- isTRUE(all.equal(results$Estimate[i], estimate, tolerance = 1e-6))
        if (same_scale && is.na(results$ConfLow[i])) results$ConfLow[i] <- test$conf.int[1]
        if (same_scale && is.na(results$ConfHigh[i])) results$ConfHigh[i] <- test$conf.int[2]
      }
    }
    variables <- unique(c(results$Outcome[i], results$Predictor[i], covariates))
    if (is.na(results$N[i]) && all(!is.na(variables)) && all(variables %in% names(df_Data))) {
      results$N[i] <- sum(stats::complete.cases(df_Data[, variables, drop = FALSE]))
    }
    if (results$Status[i] == 'failed' && grepl('Relationship could not|Insufficient usable', results$Reason[i]) &&
        all(c(results$Outcome[i], results$Predictor[i]) %in% names(df_Data))) {
      if (length(unique(stats::na.omit(df_Data[[results$Outcome[i]]]))) < 2L) results$Reason[i] <- 'Outcome has fewer than two usable values.'
      if (length(unique(stats::na.omit(df_Data[[results$Predictor[i]]]))) < 2L) results$Reason[i] <- 'Predictor has fewer than two usable values.'
    }
  }
  results
}

ScidrInferenceRegistry <- function(results, context) {
  vars_Keys <- c('Outcome', 'Predictor', 'Group')
  if (context$Function == 'MultivariableRegressionTable') {
    results$ModelComparison <- if ('OutcomeMode' %in% names(results) && 'OutcomeLevel' %in% names(results)) ifelse(results$OutcomeMode %in% 'one_vs_rest', results$OutcomeLevel, '') else rep('', nrow(results))
    vars_Keys <- c('Outcome', 'ModelComparison', 'Group')
  }
  if (context$Function %in% c('MakeComparisonTable', 'MakeFacetCatComparisonTable') && 'TestRole' %in% names(results)) {
    covariates <- ScidrInferenceValue(context$Frame, c('Covariates', 'covariates'))
    df_Data <- ScidrInferenceData(context$Frame)
    continuous <- ScidrInferenceValue(context$Frame, 'treat_as_continuous')
    if (is.null(continuous)) continuous <- vapply(df_Data, is.numeric, logical(1))
    shared_model <- length(covariates) > 0L & results$Outcome %in% names(continuous)[continuous]
    results$ModelComparison <- ifelse(results$TestRole %in% 'Pairwise' & !shared_model, results$Term, '')
    vars_Keys <- c(vars_Keys, 'ModelComparison')
  }
  if (context$Function %in% c('EvaluateBiomarkerPerformance', 'ScreenBiomarkerPerformance') && 'Model' %in% names(results)) vars_Keys <- c(vars_Keys, 'Model')
  df_Models <- dplyr::distinct(results[, c(vars_Keys, 'Status', 'Reason'), drop = FALSE])
  group_keys <- do.call(paste, c(df_Models[vars_Keys], sep = '\r'))
  groups <- split(seq_len(nrow(df_Models)), group_keys)
  df_Models$Partial <- FALSE
  for (indices in groups) {
    df_Models$Partial[indices[1]] <- any(df_Models$Status[indices] == 'ok') && any(df_Models$Status[indices] == 'failed')
    df_Models$Status[indices[1]] <- if (any(df_Models$Status[indices] == 'ok')) 'ok' else 'failed'
    reasons <- unique(stats::na.omit(df_Models$Reason[indices]))
    df_Models$Reason[indices[1]] <- if (length(reasons)) paste(reasons, collapse = '; ') else NA_character_
  }
  df_Models <- df_Models[!duplicated(group_keys), , drop = FALSE]
  if (!'Predictor' %in% names(df_Models)) df_Models$Predictor <- NA_character_
  if (context$Function %in% c('PlotMiningMatrix', 'PlotCorrelationsHeatmap', 'PlotPhiHeatmap', 'PlotDirectionalHeatmaps', 'PlotSplitCorrelationHeatmap', 'PlotCorrelationComparisons')) {
    keys <- paste(pmin(df_Models$Outcome, df_Models$Predictor), pmax(df_Models$Outcome, df_Models$Predictor), df_Models$Group, sep = '\r')
    df_Models <- df_Models[!duplicated(keys), , drop = FALSE]
  }
  df_Models$ModelID <- paste0('model_', seq_len(nrow(df_Models)))
  df_Models$Formula <- NA_character_
  df_Models$Engine <- context$Function
  for (i in seq_len(nrow(df_Models))) {
    key <- if (is.na(df_Models$Predictor[i])) df_Models$Outcome[i] else paste(df_Models$Outcome[i], df_Models$Predictor[i], sep = '\r')
    index <- if (is.na(df_Models$Predictor[i])) context$ByOutcome else context$ByPair
    observations <- if (!is.na(key) && exists(key, envir = index, inherits = FALSE)) get(key, envir = index, inherits = FALSE) else list()
    matches <- Filter(function(observation) {
      identical(observation$Outcome, df_Models$Outcome[i]) &&
        (is.na(df_Models$Predictor[i]) || identical(observation$Predictor, df_Models$Predictor[i])) &&
        !is.null(observation$Formula) && !is.na(observation$Formula) &&
        !grepl('scale$', observation$Engine)
    }, observations)
    if (length(matches)) {
      observation <- matches[[1]]
      df_Models$Formula[i] <- observation$Formula
      df_Models$Engine[i] <- observation$Engine
    }
  }
  df_Models
}

ScidrInferenceFailureKey <- function(data, function_name) {
  outcome <- data$Outcome
  predictor <- data$Predictor
  if (function_name %in% c('PlotMiningMatrix', 'PlotDirectionalHeatmaps', 'PlotPhiHeatmap')) {
    outcome <- pmin(data$Outcome, data$Predictor)
    predictor <- pmax(data$Outcome, data$Predictor)
  }
  paste(outcome, predictor, data$Group, data$Stage, sep = '\r')
}

#' @export
`$.SciDataReportRInference` <- function(x, name) {
  if (name %in% c('Results', 'Metadata')) return(attr(x, name, exact = TRUE))
  NextMethod('$')
}

ScidrFinishInference <- function(output, context, frame) {
  ScidrIndexInference(context)
  source <- ScidrInferenceSource(output, context, frame)
  df_Results <- ScidrNormalizeInference(source, context, frame)
  child_failures <- dplyr::bind_rows(lapply(seq_along(context$Children), function(i) {
    df_Failed <- context$Children[[i]]$Metadata$Failed
    if (!nrow(df_Failed)) return(NULL)
    if (context$Function == 'PlotSplitCorrelationHeatmap') {
      levels <- ScidrInferenceValue(frame, 'split_levels')
      if (length(levels) >= i) df_Failed$Group <- levels[i]
    }
    df_Failed$ChildIndex <- i
    df_Failed
  }))
  if (length(context$Failures)) {
    df_ExplicitFailures <- ScidrNormalizeInference(dplyr::bind_rows(context$Failures), context, frame)
    for (i in seq_len(nrow(df_ExplicitFailures))) {
      index <- which(df_Results$Outcome %in% df_ExplicitFailures$Outcome[i] & df_Results$Predictor %in% df_ExplicitFailures$Predictor[i] & df_Results$Status == 'failed')
      df_Results$Reason[index] <- df_ExplicitFailures$Reason[i]
    }
    df_ExplicitFailures <- dplyr::anti_join(df_ExplicitFailures, df_Results, by = c('Outcome', 'Predictor'))
    df_Results <- dplyr::bind_rows(df_Results, df_ExplicitFailures)
  }
  if (length(context$Estimates)) {
    df_Estimates <- ScidrNormalizeInference(dplyr::bind_rows(context$Estimates), context, frame)
    if (!'TestRole' %in% names(df_Results)) df_Results$TestRole <- rep('Omnibus', nrow(df_Results))
    if (context$Function %in% c('MakeComparisonTable', 'MakeFacetCatComparisonTable')) df_Results <- dplyr::bind_rows(df_Results, df_Estimates)
    if (context$Function == 'PlotSplitViolin') df_Results <- df_Estimates
  }
  if (length(context$PairwisePValues)) {
    if (!'TestRole' %in% names(df_Results)) df_Results$TestRole <- rep('Omnibus', nrow(df_Results))
    df_Raw <- dplyr::bind_rows(context$PairwisePValues)
    keys <- paste(df_Raw$Outcome, df_Raw$Term, sep = '\r')
    df_Raw <- df_Raw[!duplicated(keys, fromLast = TRUE), , drop = FALSE]
    df_Raw <- ScidrNormalizeInference(df_Raw, context, frame)
    keys <- paste(df_Raw$Outcome, df_Raw$Term, sep = '\r')
    index <- match(paste(df_Results$Outcome, df_Results$Term, sep = '\r'), keys)
    ok <- which(!is.na(index) & df_Results$TestRole %in% 'Pairwise')
    df_Results$PValue[ok] <- df_Raw$PValue[index[ok]]
    df_Results$PValueFDR[ok] <- NA_real_
    df_Missing <- dplyr::anti_join(df_Raw, df_Results[df_Results$TestRole %in% 'Pairwise', , drop = FALSE], by = c('Outcome', 'Term'))
    df_Results <- dplyr::bind_rows(df_Results, df_Missing)
  }
  # Self-pairs represent structure, not fitted relationships.
  if (context$Function %in% c('PlotMiningMatrix', 'PlotCorrelationsHeatmap', 'PlotPhiHeatmap', 'PlotChiSqCovar',
      'PlotDirectionalHeatmaps', 'PlotSplitCorrelationHeatmap', 'PlotCorrelationComparisons')) {
    df_Results <- df_Results[is.na(df_Results$Outcome) | is.na(df_Results$Predictor) | df_Results$Outcome != df_Results$Predictor, , drop = FALSE]
  }
  df_Requested <- ScidrInferenceRequested(context, frame)
  if (nrow(df_Requested)) {
    vars_Keys <- intersect(c('Outcome', 'Predictor', 'Group'), names(df_Requested))
    df_Missing <- dplyr::anti_join(df_Requested, df_Results, by = vars_Keys)
    if (nrow(df_Missing)) {
      df_Missing$Status <- 'failed'
      df_Missing$Reason <- 'Insufficient usable observations, variation, or estimable model terms.'
      df_Data <- ScidrInferenceData(frame)
      df_Missing$N <- vapply(seq_len(nrow(df_Missing)), function(i) {
        variables <- unique(c(df_Missing$Outcome[i], df_Missing$Predictor[i], ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars', 'covars'))))
        if (!all(variables %in% names(df_Data))) return(NA_real_)
        sum(stats::complete.cases(df_Data[, variables, drop = FALSE]))
      }, numeric(1))
      df_Results <- dplyr::bind_rows(df_Results, ScidrNormalizeInference(df_Missing, context, frame))
    }
  }
  for (observation in context$Observations) {
    if (observation$Status != 'failed') next
    matches <- which(df_Results$Outcome %in% observation$Outcome & df_Results$Predictor %in% observation$Predictor)
    if (length(matches)) df_Results$Reason[matches[df_Results$Status[matches] == 'failed']] <- observation$Reason
  }
  df_Results <- ScidrEnrichInference(df_Results, context, frame)
  adjustment_issues <- dplyr::bind_rows(context$Failures)
  if (nrow(child_failures)) {
    adjustment_issues <- dplyr::bind_rows(adjustment_issues,
      child_failures[child_failures$Stage %in% c('adjustment', 'estimation'), , drop = FALSE])
  }
  if (nrow(adjustment_issues) && any(adjustment_issues$Stage %in% 'adjustment')) {
    df_Results$AdjustmentStatus <- rep('adjusted', nrow(df_Results))
    for (i in which(adjustment_issues$Stage %in% 'adjustment')) {
      index <- which(df_Results$Outcome %in% adjustment_issues$Outcome[i] & df_Results$Predictor %in% adjustment_issues$Predictor[i])
      if (context$Function == 'PlotSplitCorrelationHeatmap' && 'Group' %in% names(adjustment_issues)) index <- index[df_Results$Group[index] %in% adjustment_issues$Group[i]]
      if (!'FallbackEstimate' %in% names(df_Results)) df_Results$FallbackEstimate <- rep(NA_real_, nrow(df_Results))
      if (!'FallbackPValue' %in% names(df_Results)) df_Results$FallbackPValue <- rep(NA_real_, nrow(df_Results))
      df_Results$FallbackEstimate[index] <- dplyr::coalesce(df_Results$FallbackEstimate[index], df_Results$Estimate[index])
      df_Results$FallbackPValue[index] <- dplyr::coalesce(df_Results$FallbackPValue[index], df_Results$PValue[index])
      for (name in c('Estimate', 'ConfLow', 'ConfHigh', 'PValue', 'PValueFDR')) df_Results[[name]][index] <- NA_real_
      df_Results$Status[index] <- 'failed'
      df_Results$AdjustmentStatus[index] <- 'unadjusted_fallback'
      df_Results$Reason[index] <- adjustment_issues$Reason[i]
    }
  }
  df_Results <- ScidrInferenceFDR(df_Results, frame)
  p_column <- if (context$PVar %in% c('P', 'p', 'p.value')) 'PValue' else if (context$PVar %in% c('FDR', 'P_adj', 'p_adj')) 'PValueFDR' else context$PVar
  if (!p_column %in% names(df_Results)) stop('Requested PVar is unavailable: ', context$PVar, call. = FALSE)
  df_Results$Significant <- !is.na(df_Results[[p_column]]) & df_Results[[p_column]] < context$Pthresh
  if (context$Function == 'PlotInteractionEffectsMatrix') {
    df_Results$Significant <- !is.na(df_Results[[p_column]]) & df_Results[[p_column]] <= context$Pthresh
  }
  df_Failed <- df_Results[df_Results$Status == 'failed', c('Outcome', 'Predictor', 'Group', 'Term', 'N', 'Reason'), drop = FALSE]
  df_Failed$Stage <- rep('estimation', nrow(df_Failed))
  if (nrow(adjustment_issues)) {
    for (i in seq_len(nrow(df_Failed))) {
      index <- which(adjustment_issues$Outcome %in% df_Failed$Outcome[i] & adjustment_issues$Predictor %in% df_Failed$Predictor[i])
      if (length(index) && any(adjustment_issues$Stage[index] %in% 'adjustment')) df_Failed$Stage[i] <- 'adjustment'
    }
  }
  if (length(context$Failures)) {
    df_Issues <- dplyr::bind_rows(context$Failures)
    df_Issues <- dplyr::anti_join(df_Issues, df_Failed, by = c('Outcome', 'Predictor'))
    if (nrow(df_Issues)) {
      df_Issues$Group <- 'Pooled'
      df_Issues$Term <- NA_character_
      default_stage <- if (context$Function %in% c('plotSigAssociations', 'plotSigCorrelations')) 'plotting' else 'auxiliary'
      df_Issues$Stage[df_Issues$Stage == 'estimation'] <- default_stage
      df_Failed <- dplyr::bind_rows(df_Failed, df_Issues[, names(df_Failed), drop = FALSE])
    }
  }
  for (observation in context$Observations) {
    if (observation$Status != 'failed') next
    if (any(df_Failed$Outcome %in% observation$Outcome & df_Failed$Predictor %in% observation$Predictor)) next
    df_Failed <- dplyr::bind_rows(df_Failed, tibble::tibble(Outcome = observation$Outcome,
      Predictor = observation$Predictor, Group = 'Pooled', Term = NA_character_, N = NA_real_,
      Reason = observation$Reason, Stage = 'auxiliary'))
  }
  if (nrow(child_failures)) {
    keys <- ScidrInferenceFailureKey(child_failures, context$Function)
    additional <- child_failures[!keys %in% ScidrInferenceFailureKey(df_Failed, context$Function), , drop = FALSE]
    df_Failed <- dplyr::bind_rows(df_Failed, additional)
  }
  df_Warnings <- if (nrow(child_failures)) df_Failed[!ScidrInferenceFailureKey(df_Failed, context$Function) %in%
    ScidrInferenceFailureKey(child_failures, context$Function), , drop = FALSE] else df_Failed
  if (nrow(df_Warnings)) {
    messages <- paste0('outcome ', df_Warnings$Outcome, ', predictor ', df_Warnings$Predictor,
      ' (', df_Warnings$Group, '): ', df_Warnings$Reason)
    prefix <- if (context$Function == 'PlotMiningMatrix') 'Missing requested pairs: ' else ''
    warning(context$Function, ': ', prefix, paste(unique(messages), collapse = '\n'), call. = FALSE)
  }
  existing <- if (inherits(output, 'ggplot')) output$Metadata else if (is.list(output) && !is.data.frame(output)) output$Metadata else attr(output, 'Metadata', exact = TRUE)
  if (!is.list(existing)) existing <- list()
  axes <- ScidrInferenceAxes(context, frame)
  if (!length(axes$Outcome)) axes$Outcome <- unique(df_Results$Outcome)
  if (!length(axes$Predictor)) axes$Predictor <- unique(df_Results$Predictor)
  arguments <- list(data = if ('data' %in% names(context$Call)) paste(deparse(context$Call$data), collapse = ' ') else if (length(context$Call) >= 2L) paste(deparse(context$Call[[2]]), collapse = ' ') else NA_character_,
    outcome_vars = axes$Outcome, predictor_vars = axes$Predictor,
    interVar = ScidrInferenceValue(frame, 'interVar'), group_var = ScidrInferenceValue(frame, c('group_var', 'CompVariable', 'TargetVar', 'GroupVariable', 'split_var', 'stratify_by')),
    covariates = ScidrInferenceValue(frame, c('covariates', 'Covariates', 'Covars', 'covars')),
    TreatOrdinalAs = ScidrInferenceValue(frame, 'TreatOrdinalAs', 'not supported'),
    standardize = ScidrInferenceValue(frame, c('Standardize', 'standardize'), FALSE),
    fdr_scope = ScidrInferenceValue(frame, c('fdr_scope', 'adjust_scope'), 'matrix')[1],
    covariate_interactions = ScidrInferenceValue(frame, 'covariate_interactions', 'not supported'),
    PVar = p_column, Pthresh = context$Pthresh)
  if (context$Function %in% c('PlotBetaProfile', 'PlotZScore', 'MakePairwiseHeatmap', 'PlotPairwiseMiningMatrix')) arguments$standardize <- TRUE
  if (context$Function == 'PlotVolcanoEffects') arguments$standardize <- identical(ScidrInferenceValue(frame, 'OutcomeType'), 'continuous')
  options <- list()
  function_object <- get(context$Function, envir = asNamespace('SciDataReportR'))
  for (name in setdiff(names(formals(function_object)), c('data', 'Data', 'DataFrame', 'df', '...', names(arguments)))) {
    value <- ScidrInferenceValue(frame, name)
    if (!is.data.frame(value) && !is.function(value) && !is.environment(value)) options[name] <- list(value)
  }
  arguments$Options <- options
  df_Models <- dplyr::bind_rows(lapply(seq_along(context$Observations), function(i) {
    observation <- context$Observations[[i]]
    tibble::tibble(ModelID = paste0('model_', i), Outcome = observation$Outcome,
      Predictor = observation$Predictor, Formula = if (is.null(observation$Formula)) NA_character_ else observation$Formula,
      Engine = observation$Engine, Status = observation$Status, Reason = observation$Reason)
  }))
  standardization <- dplyr::bind_rows(lapply(seq_along(context$Observations), function(i) {
    observation <- context$Observations[[i]]
    dplyr::bind_rows(lapply(names(observation$Standardization), function(variable) {
      item <- observation$Standardization[[variable]]
      tibble::tibble(ObservationID = i, Outcome = observation$Outcome, Predictor = observation$Predictor,
        Variable = variable, SD = item$SD, Center = item$Center, N = item$N)
    }))
  }))
  reference_levels <- Filter(Negate(is.null), lapply(context$Observations, `[[`, 'ReferenceLevels'))
  if (!nrow(df_Models)) {
    df_Models <- dplyr::distinct(df_Results, Outcome, Predictor, Group, Status, Reason)
    df_Models$ModelID <- paste0('model_', seq_len(nrow(df_Models)))
    df_Models$Formula <- NA_character_
    df_Models$Engine <- context$Function
  }
  df_AuxiliaryModels <- df_Models
  df_Models <- ScidrInferenceRegistry(df_Results, context)
  if (nrow(standardization)) {
    has_scale_observation <- vapply(context$Observations, function(x) grepl('(^|::)scale$', x$Engine), logical(1))
    if (any(has_scale_observation)) standardization <- standardization[standardization$ObservationID %in% which(has_scale_observation), , drop = FALSE]
    standardization$ModelID <- vapply(seq_len(nrow(standardization)), function(i) {
      matches <- which(df_Models$Outcome %in% standardization$Outcome[i] &
        (is.na(df_Models$Predictor) | df_Models$Predictor %in% standardization$Predictor[i]))
      if (length(matches)) df_Models$ModelID[matches[1]] else NA_character_
    }, character(1))
    standardization <- dplyr::distinct(standardization[, setdiff(names(standardization), 'ObservationID'), drop = FALSE])
    standardization$Scope <- rep('Observations used by the original scaling operation', nrow(standardization))
  }
  if (context$Function == 'MakePairwiseHeatmap') {
    scale_object <- ScidrInferenceValue(frame, 'scale_obj')
    parameters <- scale_object$Parameters
    if (is.data.frame(parameters) && 'Variable' %in% names(parameters) &&
        (all(c('Mean', 'SD') %in% names(parameters)) || all(c('Median', 'MAD') %in% names(parameters)))) {
      standardization <- dplyr::bind_rows(lapply(seq_len(nrow(df_Models)), function(i) {
        index <- match(df_Models$Outcome[i], parameters$Variable)
        if (is.na(index)) return(NULL)
        tibble::tibble(ModelID = df_Models$ModelID[i], Outcome = df_Models$Outcome[i], Predictor = df_Models$Predictor[i],
          Variable = parameters$Variable[index], SD = if ('SD' %in% names(parameters)) parameters$SD[index] else NA_real_,
          ScaleFactor = if ('SD' %in% names(parameters)) parameters$SD[index] else parameters$MAD[index],
          Center = if ('Mean' %in% names(parameters)) parameters$Mean[index] else parameters$Median[index],
          ScaleMethod = if ('SD' %in% names(parameters)) 'SD' else 'MAD',
          N = if ('N' %in% names(parameters)) parameters$N[index] else NA_real_, Scope = 'Referent-group scaling parameters')
      }))
    }
  }
  if (context$Function == 'MultivariableRegressionTable' && length(context$Children)) {
    assigned <- character()
    for (child_index in seq_along(context$Children)) {
      child <- context$Children[[child_index]]$Metadata
      if (!identical(child$Function, 'MultivariableRegressionTable')) next
      child_models <- child$Models
      child_scale <- child$Standardization
      for (j in seq_len(nrow(child_models))) {
        index <- which(df_Models$Outcome %in% child_models$Outcome[j] & !df_Models$ModelID %in% assigned)
        if (!length(index)) next
        index <- index[1]
        parent_id <- df_Models$ModelID[index]
        assigned <- c(assigned, parent_id)
        if (is.na(df_Models$Formula[index])) df_Models$Formula[index] <- child_models$Formula[j]
        if (nrow(child_scale)) {
          df_Scale <- child_scale[child_scale$ModelID %in% child_models$ModelID[j], , drop = FALSE]
          df_Scale$ModelID <- rep(parent_id, nrow(df_Scale))
          df_Scale$ChildIndex <- rep(child_index, nrow(df_Scale))
          standardization <- dplyr::bind_rows(standardization, df_Scale)
        }
      }
    }
  }
  if (!ncol(standardization)) standardization <- tibble::tibble(ModelID = character(), Variable = character(),
    SD = double(), Center = double(), N = double(), Scope = character())
  formula_template <- if (length(arguments$interVar)) 'outcome ~ predictor * interVar + covariates' else 'outcome ~ predictor + covariates'
  if (context$Function == 'PlotInteractionEffectsMatrix' && isTRUE(arguments$covariate_interactions)) formula_template <- 'outcome ~ predictor * interVar + covariates * interVar'
  if (context$Function == 'MakeCoxRegressionTable') formula_template <- 'survival::Surv(time_var, event_var) ~ predictor + covariates'
  if (context$Function %in% c('PlotCorrelationsHeatmap', 'PlotSplitCorrelationHeatmap', 'PlotCorrelationComparisons', 'calculate_pathway_results')) formula_template <- 'cor(outcome, predictor); residualize both variables on covariates when supplied'
  if (context$Function == 'PlotBlandAltman') formula_template <- 'difference = Variable1 - Variable2; bias and limits of agreement'
  if (context$Function == 'DiagnosticLikelihoodRatioTable') formula_template <- 'P(test result | outcome positive) / P(test result | outcome negative), within each stratum'
  if (context$Function == 'PlotBetaProfile') formula_template <- 'scale(outcome) ~ scale(predictor) + covariates'
  if (context$Function == 'PlotMiningMatrix') formula_template <- 'Mixed-type screening: correlation, ANOVA/Kruskal-Wallis, or chi-square; supplied covariates are excluded from tested axes'
  if (context$Function == 'PlotChiSqCovar') formula_template <- 'chisq.test(outcome, predictor), using observations complete for the tested variables and supplied covariates'
  if (context$Function == 'PlotPointCorrelationsHeatmap') formula_template <- 'cor(continuous outcome, binary predictor), using observations complete for tested variables and supplied covariates'
  reference_levels <- list(Models = reference_levels,
    Outcomes = if (!is.null(existing$OutcomeLevels)) existing$OutcomeLevels else existing$Outcomes, Predictors = existing$PredictorLevels,
    BinaryPredictors = existing$BinaryPredictorLevels,
    Comparisons = if (nrow(df_Results)) unique(df_Results[, intersect(c('Outcome', 'Predictor', 'Group', 'Level', 'ReferenceLevel', 'ComparisonLevel', 'OutcomePositiveLevel', 'ResultLevel'), names(df_Results)), drop = FALSE]) else tibble::tibble(),
    Referent = ScidrInferenceValue(frame, c('Referent', 'reference_group')))
  if (length(arguments$group_var) == 1L && arguments$group_var %in% names(ScidrInferenceData(frame))) reference_levels$GroupLevels <- levels(factor(ScidrInferenceData(frame)[[arguments$group_var]]))
  if (length(arguments$interVar) == 1L && arguments$interVar %in% names(ScidrInferenceData(frame))) {
    moderator <- ScidrInferenceData(frame)[[arguments$interVar]]
    reference_levels$Moderator <- if (is.numeric(moderator)) list(
      Mean = ScidrInferenceValue(frame, c('interVar_mean', 'mod_mean'), mean(moderator, na.rm = TRUE)),
      SD = ScidrInferenceValue(frame, c('interVar_sd', 'mod_sd'), stats::sd(moderator, na.rm = TRUE))) else levels(factor(moderator))
  }
  family <- paste0('Benjamini-Hochberg adjustment; scope = ', arguments$fdr_scope,
    '. Finite raw p-values only; unavailable requests are excluded. ',
    if (arguments$fdr_scope == 'per_outcome') 'Separate family within each outcome.' else if (arguments$fdr_scope == 'per_predictor') 'Separate family within each predictor.' else 'The returned analysis family is adjusted together.',
    if (context$Function == 'PlotMiningMatrix') ' Legacy adjustment includes the symmetrized table before display subsetting.' else '',
    if (any(c('Split', 'Stratum') %in% names(df_Results))) ' Separate families within each stratum.' else '')
  if (context$Function == 'PlotInteractionEffectsMatrix') family <- paste0(family, ' Interaction tests and simple-slope tests are adjusted as separate families using the same scope.')
  if (context$Function == 'PlotDirectionalHeatmaps') family <- paste0(family,
    ' Correlation, phi, and point-biserial blocks retain separate correction families. Mirrored point-biserial cells reuse their original adjusted values.')
  if ('TestRole' %in% names(df_Results)) family <- paste0(family,
    ' Omnibus tests and pairwise contrasts are separate families; pairwise BH correction is within each outcome. Legacy non-BH contrast corrections retain their original reporting columns.')
  source_metadata <- NULL
  if (context$Function %in% c('plotSigAssociations', 'plotSigCorrelations')) {
    source_object <- ScidrInferenceValue(frame, 'Plot')
    source_metadata <- source_object$Metadata
    if (is.list(source_metadata) && !is.null(source_metadata$Arguments)) {
      for (name in c('covariates', 'interVar', 'group_var', 'TreatOrdinalAs', 'standardize', 'fdr_scope', 'covariate_interactions')) arguments[name] <- source_metadata$Arguments[name]
      formula_template <- source_metadata$FormulaTemplate
      reference_levels <- source_metadata$ReferenceLevels
      standardization <- source_metadata$Standardization
      family <- paste0('Retained from the source analysis before significance selection: ', source_metadata$FDRFamily)
    }
  }
  contract <- list(Function = context$Function,
    PackageVersion = as.character(utils::packageVersion('SciDataReportR')), Call = context$Call,
    Arguments = arguments, FormulaTemplate = formula_template,
    Models = df_Models, ReferenceLevels = reference_levels, Standardization = standardization,
    FDRFamily = family, ModelCounts = list(Attempted = nrow(df_Models), Fit = sum(df_Models$Status == 'ok'), Failed = sum(df_Models$Status == 'failed'), AuxiliaryFailed = sum(df_Failed$Stage == 'auxiliary')),
    Failed = df_Failed, AuxiliaryModels = df_AuxiliaryModels)
  metadata <- existing
  metadata[names(contract)] <- contract
  if (is.list(source_metadata)) metadata$Source <- source_metadata
  if (!is.null(existing$FDRFamily)) metadata$FDRFamilyLegacy <- existing$FDRFamily
  if (length(context$Children)) metadata$Children <- lapply(context$Children, `[[`, 'Metadata')
  if (!is.null(context$Parent)) context$Parent$Children[[length(context$Parent$Children) + 1L]] <- list(Results = df_Results, Metadata = metadata)
  if (context$Function == 'MakeComparisonTable') output <- ScidrStyleComparisonInference(output, df_Results, context, frame)
  output <- ScidrAttachInference(output, df_Results, metadata)
  # Reporting objects can retain closures over this frame. Do not retain the
  # transient collector, its model objects, or the caller's execution frame.
  rm('.scidr_inference', envir = frame)
  output
}
