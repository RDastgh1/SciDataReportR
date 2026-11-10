#' Plot Interaction Effects Matrix
#'
#' Creates a heatmap visualization of interaction effects between continuous variables,
#' showing whether interactions result in slope reversals or maintain the same direction.
#'
#' @param data A data frame containing the variables to be analyzed
#' @param interVar Character string specifying the interaction variable (moderator). Can be categorical or continuous.
#' @param outcome_vars Character vector of outcome variable names (displayed on rows)
#' @param predictor_vars Character vector of predictor variable names (displayed on columns)
#' @param covariates Character vector of covariate names to include in the models
#' @param Relabel Logical indicating whether to use variable labels if available (default: TRUE)
#' @param TreatOrdinalAs How ordinal variables are handled. This interaction
#'   matrix accepts `"Exclude"`, `"Continuous"`, or `"Categorical"`.
#'   Categorical ordinal outcomes are not supported by the linear models.
#' @param Ordinal \strong{Deprecated} (since 20.20.0). Use
#'   \code{TreatOrdinalAs} instead.
#'
#' @return A list containing:
#' \describe{
#'   \item{Unadjusted}{List with unadjusted results including:
#'     \itemize{
#'       \item C: Matrix of interaction coefficients
#'       \item S: Matrix of slope direction indicators (1 = same, -1 = reversed)
#'       \item P: Matrix of p-values
#'       \item D: Matrix of interaction coefficient signs
#'       \item Slope1: Matrix of slopes at low values of interVar (mean - 1SD for continuous, reference group for categorical)
#'       \item Slope2: Matrix of slopes at high values of interVar (mean + 1SD for continuous, comparison group for categorical)
#'       \item Slope1Std, Slope2Std, SlopeDeltaStd: Common-scale standardized slope matrices
#'       \item plot: ggplot object of the heatmap
#'       \item pvaltable: P-value table in wide format
#'     }
#'   }
#'   \item{FDRCorrected}{List with FDR-corrected results (same structure as Unadjusted)}
#'   \item{Relabel}{Logical indicating whether relabeling was applied}
#'   \item{Covariates}{Character vector of covariates used}
#'   \item{interVar}{The interaction variable name}
#'   \item{raw_data}{Processed data frame with all calculated values}
#'   \item{Results}{Long-form interaction and simple-slope inference}
#'   \item{Metadata}{Model provenance, failures, scaling, and analysis choices}
#'   \item{SlopeDifferenceLimits}{Shared symmetric plot limits}
#' }
#'
#' @details
#' Models use `outcome ~ predictor * interVar + covariates` by default.
#' With `covariate_interactions = TRUE`, every covariate also interacts with
#' the moderator. For categorical moderators, point estimates of group slopes
#' then equal separate within-group linear fits on the same complete cases.
#' Their standard errors use the pooled model residual variance, rather than
#' separate group residual variances. The interaction test differs from the
#' common-covariate-effects model.
#'
#' Standardized slopes multiply raw slopes by SD(predictor)/SD(outcome).
#' Both SDs come from all complete cases of the fitted model, pooled across
#' groups. For continuous moderators, simple slopes are evaluated at the
#' fitted-sample moderator mean minus/plus one SD. Standardization of the
#' interaction coefficient uses the same predictor/outcome SD ratio; it does
#' not multiply by the moderator SD. `SlopeDeltaStd` instead describes the
#' difference between the two evaluated standardized slopes.
#'
#' `Results` contains one `Interaction` row and one `Slope: <level>` row per
#' evaluated slope. `Estimate`, `StdError`, `ConfLow`, `ConfHigh`, and `PValue`
#' describe raw-unit inference; `EstimateStd`, `StdErrorStd`, `ConfLowStd`, and
#' `ConfHighStd` use the common standardization. Intervals are 95% t intervals
#' from linear combinations of coefficients and `vcov(model)`. `N` is model N;
#' `NGroup`, `NReference`, and `NComparison` report categorical group counts.
#' `ModeratorValue` records continuous moderator evaluation points.
#'
#' Categorical moderators with more than two levels use a nested-model joint
#' F-test for all predictor-by-moderator terms. All level slopes are returned.
#' The plotted pair is the pair with the largest absolute slope difference;
#' factor-level order sets reference then comparison, with ties resolved by
#' the first pair. `ReferenceLevel`, `ComparisonLevel`, and
#' `LargestLevelDifference` identify this contrast. The interaction row's
#' estimate and CI describe that selected contrast, while `PValue` is the
#' joint F-test; `ContrastPValue` reports the contrast's t-test separately.
#' These contrast intervals are unadjusted for selecting the largest pair.
#' For two levels, the contrast is comparison minus reference. Factor coding
#' does not change the simple slopes or these contrasts.
#'
#' The matrix fill is `SlopeDeltaStd`: comparison minus reference (high minus
#' low for continuous moderators), with a single diverging scale. Positive
#' values are blue, negative red, and zero white. Shared symmetric limits are
#' returned in `SlopeDifferenceLimits`. Only interaction p-values <= `Pthresh`
#' are displayed: raw p-values in `Unadjusted`, BH-adjusted p-values in
#' `FDRCorrected`, regardless of `PVar`. Nonsignificant cells are blank;
#' failed models are grey and recorded with warnings in `Metadata$Failed`.
#' A no-interactions note is displayed when no cell passes the threshold.
#' Hatching marks strict sign reversals (zero is not a reversal); stars
#' describe interaction significance, not significance of either slope.
#'
#' Existing raw-unit matrices, p-value tables, aliases, and `raw_data` columns
#' are retained. `Slope1Std`, `Slope2Std`, and `SlopeDeltaStd` matrices are added
#' to both result branches. For multilevel moderators these matrices describe
#' the selected pair. `raw_data` retains `StrengthChange` as a descriptive
#' magnitude change; it no longer controls fill. `StrengthLimits` is retained
#' for compatibility and is distinct from the new plotted limits.
#' Simple-slope BH corrections form a separate family from interaction tests,
#' using the same `fdr_scope`. `Metadata$PooledScaling` records each model's SDs.
#'
#' @param Pthresh Interaction significance cutoff, default 0.05. Cells are
#'   displayed when their raw (unadjusted plot) or BH-adjusted (FDR plot)
#'   interaction p-value is less than or equal to this cutoff.
#' @param PVar P-value column used for `Results$Significant`, default
#'   `"PValue"`. `"PValueFDR"` uses BH-adjusted values. The named raw and FDR
#'   matrix plots always use their respective p-value columns.
#' @param covariate_interactions Logical, default FALSE. Interact every
#'   covariate with the moderator, allowing moderator-specific covariate effects.
#'
#' @examples
#' \donttest{
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels for readable axes
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' outcomes <- c("age", "ACE_CD143_Angiotensin_Converti",
#'               "ACTH_Adrenocorticotropic_Hormon", "AXL", "Adiponectin",
#'               "Alpha_1_Antichymotrypsin", "Alpha_1_Antitrypsin",
#'               "Alpha_1_Microglobulin")
#' predictors <- c("Alpha_2_Macroglobulin", "Angiopoietin_2_ANG_2",
#'                 "Apolipoprotein_A_IV", "Apolipoprotein_A1",
#'                 "Apolipoprotein_A2", "Apolipoprotein_B", "Apolipoprotein_CI",
#'                 "Apolipoprotein_CIII", "Apolipoprotein_D", "Apolipoprotein_E")
#'
#' # With a categorical interaction variable (Diagnosis)
#' results <- PlotInteractionEffectsMatrix(
#'   data = Labelled,
#'   interVar = "Diagnosis",
#'   outcome_vars = outcomes,
#'   predictor_vars = predictors
#' )
#'
#' # With a continuous interaction variable (age)
#' results_cont <- PlotInteractionEffectsMatrix(
#'   data = Labelled,
#'   interVar = "age",
#'   outcome_vars = setdiff(outcomes, "age"),
#'   predictor_vars = predictors
#' )
#'
#' # Raw p-value interaction matrix
#' results$Unadjusted$plot
#'
#' # FDR-adjusted interaction matrix
#' results$FDRCorrected$plot
#' }
#'
#' @param fdr_scope One of `"matrix"` (default), `"per_outcome"`, or
#'   `"per_predictor"`, passed to [ApplyFDRCorrection()]. `"matrix"` adjusts
#'   all interaction p-values together. `"per_outcome"` adjusts across predictor
#'   columns separately within each outcome row. `"per_predictor"` adjusts
#'   across outcome rows separately within each predictor column.
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param outcomeVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars} instead.
#' @param predictorVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars} instead.
#' @param covars \strong{Deprecated} (since 19.15.0). Use \code{covariates} instead.
#' @export
#' @importFrom dplyr %>% left_join mutate select rename filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_manual theme_minimal theme element_text element_blank labs
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats lm coef p.adjust as.formula sd

#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
PlotInteractionEffectsMatrix <- function(data,
    interVar = NULL,
    outcome_vars = NULL,
    predictor_vars = NULL,
    covariates = NULL,
    Relabel = TRUE,
    Ordinal = lifecycle::deprecated(),
    TreatOrdinalAs = "Exclude",
    Data = lifecycle::deprecated(),
    outcomeVars = lifecycle::deprecated(),
    predictorVars = lifecycle::deprecated(),
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    covars = lifecycle::deprecated(),
    Pthresh = 0.05,
    PVar = "PValue",
    covariate_interactions = FALSE) {
  .scidr_inference <- ScidrStartInference("PlotInteractionEffectsMatrix", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsMatrix(Data)", "PlotInteractionEffectsMatrix(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(outcomeVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsMatrix(outcomeVars)", "PlotInteractionEffectsMatrix(outcome_vars)")
    outcome_vars <- outcomeVars
  }
  outcomeVars <- outcome_vars
  if (lifecycle::is_present(predictorVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsMatrix(predictorVars)", "PlotInteractionEffectsMatrix(predictor_vars)")
    predictor_vars <- predictorVars
  }
  predictorVars <- predictor_vars
  if (lifecycle::is_present(covars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotInteractionEffectsMatrix(covars)", "PlotInteractionEffectsMatrix(covariates)")
    covariates <- covars
  }
  covars <- covariates
  fdr_scope <- match.arg(fdr_scope)
  if (lifecycle::is_present(Ordinal)) {
    lifecycle::deprecate_warn("20.20.0", "PlotInteractionEffectsMatrix(Ordinal)", "PlotInteractionEffectsMatrix(TreatOrdinalAs)")
    TreatOrdinalAs <- if (isTRUE(Ordinal)) "Continuous" else "Exclude"
  }
  TreatOrdinalAs <- match.arg(TreatOrdinalAs, c("Categorical", "Continuous", "Both", "Exclude"))
  if (TreatOrdinalAs == "Both") {
    stop("TreatOrdinalAs = 'Both' is not meaningful for PlotInteractionEffectsMatrix().", call. = FALSE)
  }


  # Check for required packages
  required_packages <- c("dplyr", "tidyr", "ggplot2", "tibble")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(paste("Missing required packages:", paste(missing_packages, collapse = ", ")))
  }

  # Handles p-values with NAs; thresholds come from the package-wide scale.
  get_significance_stars <- function(p_values) {
    ScidrPValueStars(
      p_values,
      tiers = ScidrStarTiersThree,
      ns_label = "",
      na_label = ""
    )
  }

  # Define helper functions if they don't exist
  if (!exists("ReplaceMissingLabels")) {
    ReplaceMissingLabels <- function(data) {
      return(data)
    }
  }

  if (!exists("getNumVars")) {
    getNumVars <- function(data, Ordinal = FALSE) {
      if (Ordinal) {
        vars <- names(data)[sapply(data, function(x) is.numeric(x) || is.ordered(x))]
      } else {
        vars <- names(data)[sapply(data, is.numeric)]
      }
      return(vars)
    }
  }

  # Input validation
  if (is.null(Data) || !is.data.frame(Data)) {
    stop("Data must be a non-null data frame")
  }

  if (is.null(interVar)) {
    stop("interVar must be specified")
  }

  if (!interVar %in% names(Data)) {
    stop(paste("interVar", interVar, "not found in Data"))
  }

  # Check for missing values and replace with appropriate labels
  removediag <- FALSE
  if(is.null(predictorVars)){
    predictorVars <- outcomeVars
    removediag <- TRUE
  }

  # Replace missing labels
  Data <- ReplaceMissingLabels(Data)

  # If outcomeVars is null, get all numeric variables
  if (is.null(outcomeVars)) {
    outcomeVars <- getNumVars(Data, Ordinal = FALSE)
    if (TreatOrdinalAs == "Continuous") {
      outcomeVars <- getNumVars(Data, Ordinal = TRUE)
    }
  }

  if (is.null(predictorVars)) predictorVars <- outcomeVars
  if (interVar %in% covars) {
    warning("Dropping interVar from covariates: the moderator is already in the model.", call. = FALSE)
    covars <- setdiff(covars, interVar)
  }
  # Update outcomeVars and predictorVars to exclude covariates
  if (!is.null(covars)) {
    outcomeVars <- outcomeVars[!outcomeVars %in% covars]
    predictorVars <- predictorVars[!predictorVars %in% covars]
  }

  ScidrValidateVariables(Data, outcomeVars, "outcome_vars")
  ScidrValidateVariables(Data, predictorVars, "predictor_vars")
  ScidrValidateVariables(Data, covars, "covariates")

  ordinal <- ConvertOrdinalToNumeric(
    Data, unique(c(interVar, outcomeVars, predictorVars, covars)),
    TreatOrdinalAs = TreatOrdinalAs, Relabel = Relabel, ReturnMetadata = TRUE
  )
  Data <- ordinal$data
  if (!interVar %in% ordinal$variables) {
    stop("`interVar` is ordinal and was excluded by TreatOrdinalAs = 'Exclude'.", call. = FALSE)
  }
  outcomeVars <- intersect(outcomeVars, ordinal$variables)
  predictorVars <- intersect(predictorVars, ordinal$variables)
  covars <- intersect(covars, ordinal$variables)
  if (TreatOrdinalAs == "Categorical" &&
      any(outcomeVars %in% ordinal$ordinal_variables)) {
    stop(
      "Categorical ordinal outcomes are not supported by PlotInteractionEffectsMatrix(); use TreatOrdinalAs = 'Continuous' or 'Exclude'.",
      call. = FALSE
    )
  }

  if (!is.logical(covariate_interactions) || length(covariate_interactions) != 1L ||
      is.na(covariate_interactions)) {
    stop("covariate_interactions must be TRUE or FALSE.", call. = FALSE)
  }
  covariates <- covars
  if (!length(outcomeVars) || !length(predictorVars)) {
    stop("At least one outcome and predictor are required.", call. = FALSE)
  }
  interVar_is_continuous <- is.numeric(Data[[interVar]])
  if (!interVar_is_continuous) Data[[interVar]] <- factor(Data[[interVar]])

  # These helpers stay local: contrasts must use this fitted model's design and
  # covariance matrix rather than independently fitted or rescaled groups.
  QuoteVariable <- function(variable) if (length(variable)) paste0("`", variable, "`") else character()
  ContrastSummary <- function(weights, model, scale_factor) {
    estimate <- sum(weights * stats::coef(model))
    se <- sqrt(max(0, as.numeric(t(weights) %*% stats::vcov(model) %*% weights)))
    critical <- stats::qt(0.975, stats::df.residual(model))
    p_value <- 2 * stats::pt(abs(estimate / se), stats::df.residual(model), lower.tail = FALSE)
    tibble::tibble(Estimate = estimate, StdError = se,
      ConfLow = estimate - critical * se, ConfHigh = estimate + critical * se,
      EstimateStd = estimate * scale_factor, StdErrorStd = se * scale_factor,
      ConfLowStd = (estimate - critical * se) * scale_factor,
      ConfHighStd = (estimate + critical * se) * scale_factor, PValue = p_value)
  }
  cells <- list()
  result_rows <- list()
  scaling_rows <- list()
  for (outcomeVar in outcomeVars) {
    for (predictorVar in predictorVars) {
      cell <- tibble::tibble(Outcome = outcomeVar, Predictor = predictorVar,
        C = NA_real_, S = NA_real_, P = NA_real_, D = NA_real_,
        Slope1 = NA_real_, Slope2 = NA_real_, Slope1Std = NA_real_,
        Slope2Std = NA_real_, SlopeDeltaStd = NA_real_, N = NA_real_,
        ReferenceLevel = NA_character_, ComparisonLevel = NA_character_,
        NReference = NA_real_, NComparison = NA_real_,
        Status = "failed", Reason = NA_character_)
      model_n <- NA_real_
      if (removediag && outcomeVar == predictorVar) {
        cell$Status <- "skipped"
        cell$Reason <- "Self-comparison was not fitted."
      } else {
        fit <- tryCatch({
          if (outcomeVar == interVar || predictorVar == interVar) {
            stop("The moderator must differ from the outcome and predictor.")
          }
          if (!is.numeric(Data[[outcomeVar]]) || !is.numeric(Data[[predictorVar]])) {
            stop("Interaction slopes require numeric outcomes and predictors.")
          }
          vars_Model <- unique(c(outcomeVar, predictorVar, interVar, covars))
          df_Complete <- Data[stats::complete.cases(Data[, vars_Model, drop = FALSE]),
                              vars_Model, drop = FALSE]
          model_n <- nrow(df_Complete)
          if (model_n < 3L) stop("Insufficient complete observations.")
          if (!interVar_is_continuous) {
            df_Complete[[interVar]] <- droplevels(df_Complete[[interVar]])
            if (nlevels(df_Complete[[interVar]]) < 2L) stop("Moderator has fewer than two observed levels.")
          }
          covar_terms <- QuoteVariable(setdiff(covars, c(outcomeVar, predictorVar)))
          if (covariate_interactions && length(covar_terms)) covar_terms <- paste(covar_terms, "*", QuoteVariable(interVar))
          rhs <- c(paste(QuoteVariable(predictorVar), "*", QuoteVariable(interVar)), covar_terms)
          formula <- stats::as.formula(paste(QuoteVariable(outcomeVar), "~", paste(rhs, collapse = " + ")))
          m <- ScidrObserveInference(stats::lm(formula, data = df_Complete))
          if (anyNA(stats::coef(m)) || stats::df.residual(m) <= 0L) {
            stop("Model is rank deficient or has no residual degrees of freedom.")
          }
          df_Model <- stats::model.frame(m)
          scale_factor <- stats::sd(df_Model[[predictorVar]]) / stats::sd(df_Model[[outcomeVar]])
          if (!is.finite(scale_factor) || scale_factor <= 0) stop("Model-sample SDs must be finite and positive.")

          # Identify the exact predictor-by-moderator term through terms/assign.
          model_terms <- stats::terms(m)
          term_factors <- attr(model_terms, "factors")
          term_variables <- lapply(rownames(term_factors), function(x) all.vars(str2lang(x)))
          interaction_id <- which(vapply(seq_len(ncol(term_factors)), function(k) {
            variables <- unlist(term_variables[term_factors[, k] > 0], use.names = FALSE)
            setequal(variables, c(predictorVar, interVar))
          }, logical(1)))
          design <- stats::model.matrix(m)
          interaction_columns <- which(attr(design, "assign") %in% interaction_id)
          if (!length(interaction_columns)) stop("Predictor-by-moderator interaction is not estimable.")

          if (interVar_is_continuous) {
            mod_mean <- mean(df_Model[[interVar]])
            mod_sd <- stats::sd(df_Model[[interVar]])
            moderator_values <- c(mod_mean - mod_sd, mod_mean + mod_sd)
            slope_labels <- c("Mean - 1 SD", "Mean + 1 SD")
            group_n <- rep(NA_real_, 2)
          } else {
            moderator_values <- levels(df_Model[[interVar]])
            slope_labels <- moderator_values
            group_n <- as.numeric(table(df_Model[[interVar]]))
          }
          slope_weights <- lapply(moderator_values, function(value) {
            # A unit change in predictor yields its derivative for this linear
            # model, independent of covariate values and factor contrast coding.
            df_New <- df_Model[rep(1L, 2L), , drop = FALSE]
            df_New[[predictorVar]] <- c(0, 1)
            if (interVar_is_continuous) df_New[[interVar]] <- rep(value, 2) else {
              df_New[[interVar]] <- factor(rep(value, 2), levels = m$xlevels[[interVar]])
            }
            new_design <- stats::model.matrix(stats::delete.response(model_terms),
              df_New, contrasts.arg = m$contrasts, xlev = m$xlevels)
            unname(new_design[2, ] - new_design[1, ])
          })
          df_Slopes <- dplyr::bind_rows(lapply(slope_weights, ContrastSummary,
            model = m, scale_factor = scale_factor))
          # Preserve factor-level order when orienting a selected pair.
          pair_indices <- utils::combn(seq_along(slope_labels), 2)
          pair_delta <- df_Slopes$Estimate[pair_indices[2, ]] - df_Slopes$Estimate[pair_indices[1, ]]
          selected_pair <- if (length(slope_labels) == 2L) 1L else which.max(abs(pair_delta))
          reference <- pair_indices[1, selected_pair]
          comparison <- pair_indices[2, selected_pair]
          delta_weights <- slope_weights[[comparison]] - slope_weights[[reference]]
          if (interVar_is_continuous) {
            interaction_weights <- numeric(ncol(design))
            interaction_weights[interaction_columns] <- 1
          } else interaction_weights <- delta_weights
          df_Interaction <- ContrastSummary(interaction_weights, m, scale_factor)
          contrast_p <- df_Interaction$PValue
          joint_f <- joint_df <- NA_real_
          if (!interVar_is_continuous && length(slope_labels) > 2L) {
            reduced_terms <- stats::drop.terms(model_terms, interaction_id, keep.response = TRUE)
            reduced <- stats::lm(stats::formula(reduced_terms), data = df_Complete)
            joint <- stats::anova(reduced, m)
            df_Interaction$PValue <- joint$`Pr(>F)`[2]
            joint_f <- joint$F[2]
            joint_df <- joint$Df[2]
          }
          if (!is.finite(df_Interaction$PValue)) stop("Interaction test returned no finite p-value.")
          df_Interaction$Term <- "Interaction"
          df_Interaction$Level <- NA_character_
          df_Interaction$NGroup <- NA_real_
          df_Interaction$ContrastPValue <- contrast_p
          df_Interaction$JointF <- joint_f
          df_Interaction$JointDF <- joint_df
          df_Interaction$Test <- if (is.finite(joint_f)) "Joint interaction F-test" else "Interaction t-test"
          df_Slopes$Term <- paste0("Slope: ", slope_labels)
          df_Slopes$Level <- slope_labels
          df_Slopes$NGroup <- group_n
          df_Slopes$ModeratorValue <- if (interVar_is_continuous) moderator_values else NA_real_
          df_Slopes$Test <- "Simple-slope t-test"
          df_Rows <- dplyr::bind_rows(df_Interaction, df_Slopes)
          df_Rows$Outcome <- outcomeVar
          df_Rows$Predictor <- predictorVar
          df_Rows$Group <- "Pooled"
          df_Rows$N <- nrow(df_Model)
          df_Rows$ReferenceLevel <- slope_labels[reference]
          df_Rows$ComparisonLevel <- slope_labels[comparison]
          df_Rows$NReference <- group_n[reference]
          df_Rows$NComparison <- group_n[comparison]
          df_Rows$LargestLevelDifference <- (df_Slopes$Estimate[comparison] - df_Slopes$Estimate[reference]) * scale_factor
          df_Rows$Status <- "ok"
          df_Rows$Reason <- NA_character_
          df_Rows$EffectType <- ifelse(df_Rows$Term == "Interaction", "Interaction", "Simple slope")
          df_Rows$TermLabel <- df_Rows$Term
          cell$C <- df_Interaction$Estimate
          cell$P <- df_Interaction$PValue
          cell$Slope1 <- df_Slopes$Estimate[reference]
          cell$Slope2 <- df_Slopes$Estimate[comparison]
          cell$Slope1Std <- df_Slopes$EstimateStd[reference]
          cell$Slope2Std <- df_Slopes$EstimateStd[comparison]
          cell$SlopeDeltaStd <- cell$Slope2Std - cell$Slope1Std
          cell$S <- ifelse(sign(cell$Slope1) == sign(cell$Slope2), 1, -1)
          cell$D <- sign(cell$C)
          cell$N <- nrow(df_Model)
          cell$ReferenceLevel <- slope_labels[reference]
          cell$ComparisonLevel <- slope_labels[comparison]
          cell$NReference <- group_n[reference]
          cell$NComparison <- group_n[comparison]
          cell$Status <- "ok"
          list(Rows = df_Rows, Scale = tibble::tibble(Outcome = outcomeVar,
            Predictor = predictorVar, PredictorSD = stats::sd(df_Model[[predictorVar]]),
            OutcomeSD = stats::sd(df_Model[[outcomeVar]]), N = nrow(df_Model)))
        }, error = function(err) err)
        if (inherits(fit, "error")) {
          cell$Reason <- conditionMessage(fit)
          cell$N <- model_n
          ScidrRecordInferenceFailure(outcomeVar, predictorVar, cell$Reason, model_n)
        } else {
          result_rows[[length(result_rows) + 1L]] <- fit$Rows
          scaling_rows[[length(scaling_rows) + 1L]] <- fit$Scale
        }
      }
      cells[[length(cells) + 1L]] <- cell
    }
  }
  m_G <- dplyr::bind_rows(cells)
  m_G$P_FDR <- ApplyFDRCorrection(m_G$P, fdr_scope = fdr_scope,
    outcome_ids = m_G$Outcome, predictor_ids = m_G$Predictor)
  m_G$P_Selected <- ScidrInferenceP(m_G$P, m_G$P_FDR)
  m_G$sig <- get_significance_stars(m_G$P)
  m_G$sig_FDR <- get_significance_stars(m_G$P_FDR)
  # Retain the previous descriptive raw_data columns for downstream code.
  m_G$sign <- ifelse(is.na(m_G$P), "na", ifelse(m_G$P > Pthresh, "ns",
    paste(ifelse(m_G$S == 1, "same", "reversed"), m_G$sig)))
  m_G$sign_FDR <- ifelse(is.na(m_G$P_FDR), "na", ifelse(m_G$P_FDR > Pthresh, "ns",
    paste(ifelse(m_G$S == 1, "same", "reversed"), m_G$sig_FDR)))
  m_G$sigsign <- trimws(paste(m_G$sign, m_G$sig))
  m_G$StrengthChange <- abs(m_G$Slope2Std) - abs(m_G$Slope1Std)
  unchanged <- is.finite(m_G$StrengthChange) & abs(m_G$StrengthChange) <=
    1e-12 * pmax(1, abs(m_G$Slope1Std), abs(m_G$Slope2Std))
  m_G$StrengthChange[unchanged] <- 0
  m_G$DirectionReversal <- m_G$Slope1 * m_G$Slope2 < 0
  m_G$EndingDirection <- ifelse(m_G$Slope2 > 0, "Positive", ifelse(m_G$Slope2 < 0, "Negative", "Zero"))
  DisplayLabels <- function(variables) {
    if (Relabel) unname(ScidrDisplayLabels(Data, variables)) else variables
  }
  m_G$OutcomeLabel <- factor(DisplayLabels(m_G$Outcome), levels = rev(unique(DisplayLabels(outcomeVars))))
  m_G$PredictorLabel <- factor(DisplayLabels(m_G$Predictor), levels = unique(DisplayLabels(predictorVars)))
  m_G$PlotText <- paste0("Outcome: ", m_G$OutcomeLabel, "\nPredictor: ", m_G$PredictorLabel,
    "\n", m_G$ReferenceLevel, ": ", signif(m_G$Slope1Std, 3),
    "\n", m_G$ComparisonLevel, ": ", signif(m_G$Slope2Std, 3),
    "\nStandardized slope difference: ", signif(m_G$SlopeDeltaStd, 3),
    "\nInteraction p: ", signif(m_G$P, 3), "\nFDR p: ", signif(m_G$P_FDR, 3))

  df_Results <- dplyr::bind_rows(result_rows)
  df_Failed <- m_G[m_G$Status != "ok", , drop = FALSE] %>%
    dplyr::transmute(Outcome, Predictor, Term = "Interaction", Group = "Pooled",
      Estimate = NA_real_, StdError = NA_real_, ConfLow = NA_real_, ConfHigh = NA_real_,
      EstimateStd = NA_real_, StdErrorStd = NA_real_, ConfLowStd = NA_real_, ConfHighStd = NA_real_,
      PValue = NA_real_, N, NGroup = NA_real_, ReferenceLevel, ComparisonLevel,
      NReference, NComparison, Status, Reason)
  df_Results <- dplyr::bind_rows(df_Results, df_Failed)
  cell_keys <- paste(m_G$Outcome, m_G$Predictor, sep = "\r")
  row_indices <- match(paste(df_Results$Outcome, df_Results$Predictor, sep = "\r"), cell_keys)
  df_Results$OutcomeLabel <- as.character(m_G$OutcomeLabel[row_indices])
  df_Results$PredictorLabel <- as.character(m_G$PredictorLabel[row_indices])
  df_Results$PValueFDR <- NA_real_
  interaction_rows <- df_Results$Term == "Interaction"
  df_Results$PValueFDR[interaction_rows] <- m_G$P_FDR[row_indices[interaction_rows]]
  # Simple-slope tests form a separate family, never inflating interaction FDR.
  slope_rows <- which(!interaction_rows)
  if (length(slope_rows)) df_Results$PValueFDR[slope_rows] <- ApplyFDRCorrection(
    df_Results$PValue[slope_rows], fdr_scope = fdr_scope,
    outcome_ids = df_Results$Outcome[slope_rows], predictor_ids = df_Results$Predictor[slope_rows])

  slope_limit <- max(abs(m_G$SlopeDeltaStd[is.finite(m_G$SlopeDeltaStd)]), 0)
  if (slope_limit == 0) slope_limit <- 1
  strength_limit <- max(abs(m_G$StrengthChange[is.finite(m_G$StrengthChange)]), 0)
  if (strength_limit == 0) strength_limit <- 1
  if (interVar_is_continuous) {
    direction_label <- "More positive at high moderator (+) / low moderator (-)"
  } else if (nlevels(Data[[interVar]]) == 2L) {
    group_labels <- levels(Data[[interVar]])
    direction_label <- paste0("More positive in ", group_labels[2], " (+) / ", group_labels[1], " (-)")
  } else direction_label <- "Largest level difference: comparison (+) / reference (-); see Results"

  BuildDifferencePlot <- function(p_column, star_column) {
    df_Plot <- m_G
    available <- df_Plot$Status == "ok" & is.finite(df_Plot$SlopeDeltaStd) & is.finite(df_Plot[[p_column]])
    display <- available & df_Plot[[p_column]] <= Pthresh
    df_Plot$FillValue <- NA_real_
    df_Plot$FillValue[display] <- df_Plot$SlopeDeltaStd[display]
    df_Plot$Stars <- ifelse(display, df_Plot[[star_column]], "")
    plot <- ggplot2::ggplot(df_Plot, ggplot2::aes(x = PredictorLabel, y = OutcomeLabel)) +
      ggplot2::geom_tile(fill = "white", colour = "white", linewidth = 0.4) +
      ggplot2::geom_tile(data = df_Plot[df_Plot$Status == "failed", , drop = FALSE],
                        fill = "#BDBDBD", colour = "white", linewidth = 0.4) +
      ggplot2::geom_tile(data = df_Plot[display, , drop = FALSE],
        ggplot2::aes(fill = FillValue), colour = "white", linewidth = 0.4) +
      ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
        midpoint = 0, limits = c(-slope_limit, slope_limit),
        name = paste0("Slope difference (SD per SD)\n", direction_label))
    df_Reversals <- df_Plot[display & !is.na(df_Plot$DirectionReversal) & df_Plot$DirectionReversal, , drop = FALSE]
    if (nrow(df_Reversals)) {
      df_Hatching <- df_Reversals[rep(seq_len(nrow(df_Reversals)), each = 5), ]
      offsets <- rep(seq(-0.64, 0.64, length.out = 5), nrow(df_Reversals))
      cell_x <- as.numeric(df_Hatching$PredictorLabel)
      cell_y <- as.numeric(df_Hatching$OutcomeLabel)
      df_Hatching$XStart <- cell_x + pmax(-0.46, -0.46 - offsets)
      df_Hatching$XEnd <- cell_x + pmin(0.46, 0.46 - offsets)
      df_Hatching$YStart <- cell_y + pmax(-0.46, -0.46 + offsets)
      df_Hatching$YEnd <- cell_y + pmin(0.46, 0.46 + offsets)
      plot <- plot + ggplot2::geom_segment(data = df_Hatching,
        ggplot2::aes(x = XStart, xend = XEnd, y = YStart, yend = YEnd, linetype = "Sign reversal"),
        inherit.aes = FALSE, colour = "black", linewidth = 0.25, alpha = 0.6) +
        ggplot2::scale_linetype_manual(values = c("Sign reversal" = "solid"), name = "Hatching")
    }
    if (!any(display)) {
      note <- paste0("No interactions with ", if (p_column == "P_FDR") "FDR p" else "p", " <= ", Pthresh)
      plot <- plot + ggplot2::annotate("label", x = (length(levels(df_Plot$PredictorLabel)) + 1)/2,
        y = (length(levels(df_Plot$OutcomeLabel)) + 1)/2, label = note)
    }
    plot + ggplot2::geom_text(ggplot2::aes(label = Stars), size = 4) +
      ggplot2::theme_minimal() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = NULL, y = NULL)
  }
  MakeMatrix <- function(column) {
    matrix(m_G[[column]], nrow = length(outcomeVars), byrow = TRUE,
           dimnames = list(outcomeVars, predictorVars))
  }
  M <- lapply(c("C", "S", "P", "D", "Slope1", "Slope2", "Slope1Std", "Slope2Std", "SlopeDeltaStd"), MakeMatrix)
  names(M) <- c("C", "S", "P", "D", "Slope1", "Slope2", "Slope1Std", "Slope2Std", "SlopeDeltaStd")
  M$plot <- BuildDifferencePlot("P", "sig")
  M$pvaltable <- m_G %>% dplyr::select(Outcome, Predictor, P) %>%
    tidyr::pivot_wider(names_from = Outcome, values_from = P)
  M_FDR <- M
  M_FDR$P_adjusted <- MakeMatrix("P_FDR")
  M_FDR$plot <- BuildDifferencePlot("P_FDR", "sig_FDR")
  M_FDR$pvaltable <- m_G %>% dplyr::select(Outcome, Predictor, P_FDR) %>%
    tidyr::pivot_wider(names_from = Outcome, values_from = P_FDR)
  out <- list(Unadjusted = M, FDRCorrected = M_FDR, Relabel = Relabel,
    Covariates = covars, interVar = interVar,
    interVar_type = ifelse(interVar_is_continuous, "continuous", "categorical"),
    StrengthLimits = c(-strength_limit, strength_limit),
    SlopeDifferenceLimits = c(-slope_limit, slope_limit), raw_data = m_G,
    Results = df_Results,
    Metadata = list(PooledScaling = dplyr::bind_rows(scaling_rows),
      ScalingDefinition = "Slope * SD(predictor) / SD(outcome), pooled fitted-model complete cases",
      InteractionTest = if (covariate_interactions) "Moderator-specific covariate effects" else "Common covariate effects",
      PlotPValues = c(Unadjusted = "PValue", FDRCorrected = "PValueFDR"),
      SimpleSlopeFDRFamily = "Separate from interactions; same fdr_scope",
      MultilevelContrast = "Largest absolute pairwise slope difference; factor-level order determines comparison minus reference"))
  out$p <- out$Unadjusted
  out$p_fdr <- out$FDRCorrected
  return(out)
}
