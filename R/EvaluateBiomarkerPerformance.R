.BiomarkerFormula <- function(response, terms) {
  QuoteName <- function(x) paste0("`", gsub("`", "\\\\`", x, fixed = TRUE), "`")
  stats::as.formula(paste(QuoteName(response), "~", if (length(terms)) paste(QuoteName(terms), collapse = " + ") else "1"))
}

.BiomarkerLabel <- function(data, variable, codebook, Relabel) {
  if (!Relabel) return(variable)
  if (!is.null(codebook)) {
    index <- match(variable, codebook$Variable)
    if (!is.na(index) && !is.na(codebook$Label[index]) && nzchar(codebook$Label[index])) return(as.character(codebook$Label[index]))
  }
  value <- tryCatch(sjlabelled::get_label(data[[variable]]), error = function(e) NULL)
  if (!is.null(value) && length(value) && !is.na(value[1]) && nzchar(value[1])) as.character(value[1]) else variable
}

.BiomarkerFit <- function(formula, data, binary) {
  warnings <- character()
  fit <- withCallingHandlers(tryCatch(if (binary) stats::glm(formula, data = data, family = stats::binomial()) else stats::lm(formula, data = data), error = function(e) e), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(fit, "error")) return(list(Model = NULL, Status = "failed", Note = conditionMessage(fit)))
  separation <- binary && (!isTRUE(fit$converged) || any(grepl("0 or 1", warnings, fixed = TRUE)))
  aliased <- any(is.na(stats::coef(fit)))
  list(Model = fit, Status = if (separation || aliased) "unreliable" else "ok", Note = paste(c(if (separation) "Separation or non-convergence detected." else NULL, if (aliased) "Aliased coefficient detected." else NULL, warnings), collapse = " "))
}

.BiomarkerCI <- function(y, pred, metric, R, CILevel, Seed, binary = FALSE) {
  if (R < 10) return(c(NA_real_, NA_real_))
  set.seed(Seed)
  values <- vapply(seq_len(R), function(i) { index <- sample.int(length(y), replace = TRUE); if (binary && length(unique(y[index])) < 2) return(NA_real_); metric(y[index], pred[index]) }, numeric(1))
  values <- values[is.finite(values)]
  if (length(values) < 10) return(c(NA_real_, NA_real_))
  stats::quantile(values, c((1 - CILevel) / 2, 1 - (1 - CILevel) / 2), names = FALSE)
}

.BiomarkerDiagnostic <- function(y, predicted, CILevel) {
  tp <- sum(y == 1 & predicted); fn <- sum(y == 1 & !predicted); fp <- sum(y == 0 & predicted); tn <- sum(y == 0 & !predicted)
  RatioCI <- function(a, b) if (!b) c(NA_real_, NA_real_) else stats::binom.test(a, b, conf.level = CILevel)$conf.int
  sensitivity <- if (tp + fn) tp / (tp + fn) else NA_real_; specificity <- if (tn + fp) tn / (tn + fp) else NA_real_; ppv <- if (tp + fp) tp / (tp + fp) else NA_real_; npv <- if (tn + fn) tn / (tn + fn) else NA_real_
  lr_pos <- if (is.na(sensitivity) || is.na(specificity) || specificity == 1) NA_real_ else sensitivity / (1 - specificity); lr_neg <- if (is.na(sensitivity) || is.na(specificity) || specificity == 0) NA_real_ else (1 - sensitivity) / specificity
  LRCI <- function(lr, a, b, c, d, negative = FALSE) { if (!is.finite(lr) || lr <= 0 || any(c(a,b,c,d) == 0)) return(c(NA_real_,NA_real_)); se <- if (negative) sqrt(1/b - 1/(a+b) + 1/d - 1/(c+d)) else sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)); exp(log(lr) + c(-1,1) * stats::qnorm(1-(1-CILevel)/2) * se) }
  sens_ci <- RatioCI(tp,tp+fn); spec_ci <- RatioCI(tn,tn+fp); ppv_ci <- RatioCI(tp,tp+fp); npv_ci <- RatioCI(tn,tn+fn); lrp_ci <- LRCI(lr_pos,tp,fn,fp,tn); lrn_ci <- LRCI(lr_neg,tp,fn,fp,tn,TRUE)
  tibble::tibble(Sensitivity=sensitivity,Sensitivity_Lower=sens_ci[1],Sensitivity_Upper=sens_ci[2],Specificity=specificity,Specificity_Lower=spec_ci[1],Specificity_Upper=spec_ci[2],PPV=ppv,PPV_Lower=ppv_ci[1],PPV_Upper=ppv_ci[2],NPV=npv,NPV_Lower=npv_ci[1],NPV_Upper=npv_ci[2],LRPositive=lr_pos,LRPositive_Lower=lrp_ci[1],LRPositive_Upper=lrp_ci[2],LRNegative=lr_neg,LRNegative_Lower=lrn_ci[1],LRNegative_Upper=lrn_ci[2],TP=tp,FP=fp,TN=tn,FN=fn)
}

.BiomarkerThreshold <- function(y, score, method, value, CILevel) {
  if (length(unique(score)) < 2) return(NULL)
  roc <- tryCatch(pROC::roc(y, score, levels = c(0,1), direction = "<", quiet = TRUE), error = function(e) NULL)
  if (is.null(roc)) return(NULL)
  threshold <- if (method == "custom") value else { args <- list(roc, x = if (method == "youden") "best" else value, ret = "threshold", transpose = FALSE); if (method == "youden") args$best.method <- "youden" else args$input <- method; as.numeric(do.call(pROC::coords,args)$threshold[1]) }
  if (is.null(threshold) || !is.finite(threshold)) return(NULL)
  positive <- if (roc$direction == "<") score >= threshold else score <= threshold
  list(Threshold=threshold,Direction=if(roc$direction == "<")">=" else "<=",Metrics=.BiomarkerDiagnostic(y,positive,CILevel))
}

.BiomarkerValidation <- function(df_analysis, formulas, y, binary, Validation, CVFolds, BootstrapR, Seed, apparent) {
  empty <- list(Table=tibble::tibble(), Predictions=tibble::tibble())
  if (Validation == "none") return(empty)
  Metric <- function(observed,predicted) { if(binary) c(AUC=if(length(unique(observed)) == 2 && length(unique(predicted)) > 1) as.numeric(pROC::auc(pROC::roc(observed,predicted,levels=c(0,1),direction="<",quiet=TRUE))) else NA_real_,Brier=mean((observed-predicted)^2)) else { sst <- sum((observed-mean(observed))^2); c(R2=if(sst)1-sum((observed-predicted)^2)/sst else NA_real_,RMSE=sqrt(mean((observed-predicted)^2)),MAE=mean(abs(observed-predicted))) } }
  metrics <- if(binary)c("AUC","Brier") else c("R2","RMSE","MAE")
  if (Validation == "cross_validation") {
    if (CVFolds < 2 || CVFolds > nrow(df_analysis)) return(list(Table=tibble::tibble(Model=names(formulas),Validation=Validation,Metric=metrics,Apparent=NA_real_,Validated=NA_real_,Optimism=NA_real_,Status="unavailable",Note="CVFolds must be between 2 and N."),Predictions=tibble::tibble()))
    set.seed(Seed); fold <- integer(nrow(df_analysis)); if(binary) for(level in unique(y)) { index <- which(y == level); fold[index] <- sample(rep(seq_len(CVFolds),length.out=length(index))) } else fold <- sample(rep(seq_len(CVFolds),length.out=nrow(df_analysis)))
    tables <- list(); predictions <- list()
    for(name in names(formulas)) { oof <- rep(NA_real_,nrow(df_analysis)); for(k in seq_len(CVFolds)) { fit <- .BiomarkerFit(formulas[[name]],df_analysis[fold != k,,drop=FALSE],binary)$Model; if(!is.null(fit)) oof[fold == k] <- suppressWarnings(as.numeric(stats::predict(fit,df_analysis[fold == k,,drop=FALSE],type="response"))) }; valid <- is.finite(oof); values <- if(any(valid)) Metric(y[valid],oof[valid]) else stats::setNames(rep(NA_real_,length(metrics)),metrics); tables[[name]] <- tibble::tibble(Model=name,Validation=Validation,Metric=metrics,Apparent=unname(apparent[[name]][metrics]),Validated=unname(values[metrics]),Optimism=NA_real_,Status=if(all(valid))"ok" else "partial",Note=if(all(valid))NA_character_ else "Some folds did not yield predictions."); predictions[[name]] <- tibble::tibble(Row=df_analysis$.scidr_row,Fold=fold,Model=name,Observed=y,OutOfFoldPrediction=oof,Usable=valid) }
    return(list(Table=dplyr::bind_rows(tables),Predictions=dplyr::bind_rows(predictions)))
  }
  set.seed(Seed); tables <- list()
  for(name in names(formulas)) { opt <- matrix(NA_real_,BootstrapR,length(metrics),dimnames=list(NULL,metrics)); success <- 0L; for(i in seq_len(BootstrapR)) { index <- sample.int(nrow(df_analysis),replace=TRUE); fit <- .BiomarkerFit(formulas[[name]],df_analysis[index,,drop=FALSE],binary)$Model; if(is.null(fit))next; pb <- suppressWarnings(as.numeric(stats::predict(fit,df_analysis[index,,drop=FALSE],type="response"))); po <- suppressWarnings(as.numeric(stats::predict(fit,df_analysis,type="response"))); opt[i,] <- Metric(y[index],pb)[metrics]-Metric(y,po)[metrics]; success <- success+1L }; optimism <- colMeans(opt,na.rm=TRUE); optimism[!is.finite(optimism)] <- NA_real_; tables[[name]] <- tibble::tibble(Model=name,Validation="bootstrap",Metric=metrics,Apparent=unname(apparent[[name]][metrics]),Validated=unname(apparent[[name]][metrics]-optimism),Optimism=unname(optimism),Status=if(success >= 10)"ok" else "unavailable",Note=paste0(success," successful bootstrap resamples; Validated is optimism-corrected.")) }
  list(Table=dplyr::bind_rows(tables),Predictions=tibble::tibble())
}

#' Evaluate biomarker performance
#'
#' Evaluates a continuous or categorical biomarker against a binary or continuous outcome.
#' Binary analyses use ordinary logistic regression. Models with separation,
#' non-convergence, or aliased coefficients return stable unavailable metrics.
#'
#' @param data A data frame.
#' @param outcome_var One outcome variable name.
#' @param biomarker_var One biomarker variable name.
#' @param covariates Optional covariate variable names.
#' @param PositiveLevel Positive binary outcome level, or `NULL` to use the second observed level.
#' @param OutcomeType One of `"auto"`, `"binary"`, or `"continuous"`.
#' @param ThresholdMethod One of `"youden"`, `"sensitivity"`, `"specificity"`, or `"custom"`.
#' @param ThresholdValue Target sensitivity or specificity.
#' @param RawThresholdValue Custom raw-biomarker threshold.
#' @param ProbabilityThresholdValue Custom predicted-probability threshold.
#' @param Validation One of `"none"`, `"bootstrap"`, or `"cross_validation"`.
#' @param BootstrapR Number of bootstrap optimism-correction resamples.
#' @param CVFolds Number of cross-validation folds.
#' @param CIBootstrapR Number of bootstrap confidence-interval resamples.
#' @param CILevel Confidence level.
#' @param CalibrationGroups Maximum grouped-calibration bins.
#' @param Seed Random seed.
#' @param Relabel Use codebook labels, then label attributes, for presentation.
#' @param codebook Optional data frame with `Variable` and `Label`.
#' @param Verbose Print positive-level information.
#' @return A stable named list containing models, performance, thresholds, predictions,
#'   calibration, validation, plots, and metadata.
#' @examples
#' \dontrun{EvaluateBiomarkerPerformance(df, "DiseaseCohort", "NfL", c("Age", "Sex"))}
#' @export
EvaluateBiomarkerPerformance <- function(data, outcome_var, biomarker_var, covariates = NULL, PositiveLevel = NULL, OutcomeType = c("auto","binary","continuous"), ThresholdMethod = c("youden","sensitivity","specificity","custom"), ThresholdValue = NULL, RawThresholdValue = NULL, ProbabilityThresholdValue = NULL, Validation = c("none","bootstrap","cross_validation"), BootstrapR = 500, CVFolds = 10, CIBootstrapR = 500, CILevel = .95, CalibrationGroups = 10, Seed = 123, Relabel = TRUE, codebook = NULL, Verbose = TRUE) {
  OutcomeType <- match.arg(OutcomeType); ThresholdMethod <- match.arg(ThresholdMethod); Validation <- match.arg(Validation)
  if(!is.data.frame(data) || !is.character(outcome_var) || length(outcome_var)!=1 || !is.character(biomarker_var) || length(biomarker_var)!=1) stop("data must be a data frame and outcome_var and biomarker_var must each be one character string.",call.=FALSE)
  vars_model <- unique(c(outcome_var,biomarker_var,covariates)); missing_vars <- setdiff(vars_model,names(data)); if(length(missing_vars))stop("Required variables not found in data: ",paste(missing_vars,collapse=", "),call.=FALSE)
  if(!is.null(codebook) && (!is.data.frame(codebook)||!all(c("Variable","Label")%in%names(codebook))))stop("codebook must be a data frame with Variable and Label columns.",call.=FALSE)
  if(!is.numeric(CILevel)||length(CILevel)!=1||CILevel<=0||CILevel>=1)stop("CILevel must be between 0 and 1.",call.=FALSE)
  if(ThresholdMethod%in%c("sensitivity","specificity")&&(is.null(ThresholdValue)||!is.numeric(ThresholdValue)||ThresholdValue<=0||ThresholdValue>=1))stop("ThresholdValue must be between 0 and 1 for sensitivity or specificity.",call.=FALSE)
  complete <- stats::complete.cases(data[,vars_model,drop=FALSE]); df_analysis <- data[complete,vars_model,drop=FALSE]; df_analysis$.scidr_row <- which(complete); if(nrow(df_analysis)<5)stop("Too few pairwise complete observations (N < 5).",call.=FALSE)
  for(name in c(biomarker_var,covariates)) if(is.character(df_analysis[[name]])||is.logical(df_analysis[[name]]))df_analysis[[name]]<-factor(df_analysis[[name]])
  y_original <- df_analysis[[outcome_var]]; if(OutcomeType=="auto")OutcomeType<-if(is.numeric(y_original)&&length(unique(y_original))>2)"continuous" else "binary"; if(OutcomeType=="continuous"&&!is.numeric(y_original))stop("OutcomeType is continuous but outcome_var is not numeric.",call.=FALSE)
  binary <- OutcomeType=="binary"; positive_level <- NULL; negative_level <- NULL
  if(binary){ levels_outcome <- if(is.factor(y_original))levels(droplevels(y_original))else sort(unique(as.character(y_original))); if(length(levels_outcome)!=2)stop("Binary outcomes must have exactly two observed levels.",call.=FALSE); positive_level <- if(is.null(PositiveLevel))levels_outcome[2] else as.character(PositiveLevel); if(!positive_level%in%levels_outcome)stop("PositiveLevel was not found among observed outcome levels: ",paste(levels_outcome,collapse=", "),call.=FALSE); negative_level<-setdiff(levels_outcome,positive_level); df_analysis$.scidr_y<-as.integer(as.character(y_original)==positive_level);if(Verbose)message(outcome_var,': treating "',positive_level,'" as the positive outcome level.') }
  biomarker_type <- if(is.numeric(df_analysis[[biomarker_var]]))"continuous" else "categorical"; if(biomarker_type=="categorical")df_analysis[[biomarker_var]]<-droplevels(factor(df_analysis[[biomarker_var]]));if(length(unique(df_analysis[[biomarker_var]]))<2)stop("biomarker_var has fewer than two observed values.",call.=FALSE)
  outcome_label <- .BiomarkerLabel(data,outcome_var,codebook,Relabel);biomarker_label <- .BiomarkerLabel(data,biomarker_var,codebook,Relabel); response <- if(binary)".scidr_y" else outcome_var
  formulas <- list(Biomarker=.BiomarkerFormula(response,biomarker_var),Covariates=if(length(covariates)).BiomarkerFormula(response,covariates)else NULL,Adjusted=.BiomarkerFormula(response,c(biomarker_var,covariates))); fits <- lapply(formulas,function(formula)if(is.null(formula))list(Model=NULL,Status="not_applicable",Note="No covariates supplied.")else .BiomarkerFit(formula,df_analysis,binary));if(!length(covariates))fits$Adjusted<-fits$Biomarker;models<-lapply(fits,`[[`,"Model");y<-if(binary)df_analysis$.scidr_y else df_analysis[[outcome_var]]
  PerformanceTable <- tibble::tibble(Model=c("Biomarker","Covariates","Adjusted"),Status=vapply(fits,`[[`,character(1),"Status"),Note=vapply(fits,`[[`,character(1),"Note"),AUC=NA_real_,AUC_Lower=NA_real_,AUC_Upper=NA_real_,Brier=NA_real_,Brier_Lower=NA_real_,Brier_Upper=NA_real_,CalibrationIntercept=NA_real_,CalibrationSlope=NA_real_,ObservedExpectedRatio=NA_real_,R2=NA_real_,AdjustedR2=NA_real_,R2_Lower=NA_real_,R2_Upper=NA_real_,RMSE=NA_real_,RMSE_Lower=NA_real_,RMSE_Upper=NA_real_,MAE=NA_real_,MAE_Lower=NA_real_,MAE_Upper=NA_real_,DeltaAUC_vs_Covariates=NA_real_,DeltaR2_vs_Covariates=NA_real_,N=nrow(df_analysis),NPositive=if(binary)sum(y)else NA_integer_,NNegative=if(binary)sum(!y)else NA_integer_,NMissing=sum(!complete),Prevalence=if(binary)mean(y)else NA_real_)
  predictions<-list(); apparent<-list();for(name in names(models)){model<-models[[name]];if(is.null(model)||fits[[name]]$Status!="ok"){predictions[[name]]<-rep(NA_real_,nrow(df_analysis));next};predictions[[name]]<-suppressWarnings(as.numeric(stats::predict(model,type="response")));if(binary){roc<-tryCatch(pROC::roc(y,predictions[[name]],levels=c(0,1),direction="<",quiet=TRUE),error=function(e)NULL);auc<-if(is.null(roc)||length(unique(predictions[[name]]))<2)NA_real_ else as.numeric(pROC::auc(roc));auc_ci<-if(is.null(roc))c(NA_real_,NA_real_)else tryCatch(as.numeric(pROC::ci.auc(roc,conf.level=CILevel))[c(1,3)],error=function(e)c(NA_real_,NA_real_));brier<-mean((y-predictions[[name]])^2);lp<-stats::qlogis(pmin(pmax(predictions[[name]],1e-6),1-1e-6));PerformanceTable[PerformanceTable$Model==name,c("AUC","AUC_Lower","AUC_Upper","Brier","Brier_Lower","Brier_Upper","CalibrationIntercept","CalibrationSlope","ObservedExpectedRatio")]<-list(auc,auc_ci[1],auc_ci[2],brier,.BiomarkerCI(y,predictions[[name]],function(a,b)mean((a-b)^2),CIBootstrapR,CILevel,Seed,TRUE)[1],.BiomarkerCI(y,predictions[[name]],function(a,b)mean((a-b)^2),CIBootstrapR,CILevel,Seed,TRUE)[2],tryCatch(stats::coef(stats::glm(y~1,family=stats::binomial(),offset=lp))[1],error=function(e)NA_real_),tryCatch(stats::coef(stats::glm(y~lp,family=stats::binomial()))[2],error=function(e)NA_real_),if(mean(predictions[[name]]))mean(y)/mean(predictions[[name]])else NA_real_);apparent[[name]]<-c(AUC=auc,Brier=brier)}else{sm<-summary(model);r2<-sm$r.squared;rmse<-sqrt(mean((y-predictions[[name]])^2));mae<-mean(abs(y-predictions[[name]]));r2ci<-.BiomarkerCI(y,predictions[[name]],function(a,b){sst<-sum((a-mean(a))^2);if(sst)1-sum((a-b)^2)/sst else NA_real_},CIBootstrapR,CILevel,Seed);PerformanceTable[PerformanceTable$Model==name,c("R2","AdjustedR2","R2_Lower","R2_Upper","RMSE","MAE")]<-list(r2,sm$adj.r.squared,r2ci[1],r2ci[2],rmse,mae);apparent[[name]]<-c(R2=r2,RMSE=rmse,MAE=mae)}}
  if(length(covariates)){if(binary)PerformanceTable$DeltaAUC_vs_Covariates<-PerformanceTable$AUC-PerformanceTable$AUC[PerformanceTable$Model=="Covariates"]else PerformanceTable$DeltaR2_vs_Covariates<-PerformanceTable$R2-PerformanceTable$R2[PerformanceTable$Model=="Covariates"]}
  ThresholdTable<-tibble::tibble();RawBiomarkerThreshold<-NULL;ROCData<-tibble::tibble();CalibrationData<-tibble::tibble();MarginalPredictionData<-tibble::tibble();Plots<-list(ROC=NULL,Probability=NULL,Calibration=NULL,CalibrationGrouped=NULL,BiomarkerDistribution=NULL,PredictedObserved=NULL,BiomarkerRelationship=NULL,Residuals=NULL)
  if(binary&&fits$Adjusted$Status=="ok"){threshold_rows<-list();for(name in names(predictions))if(all(is.finite(predictions[[name]]))){result<-.BiomarkerThreshold(y,predictions[[name]],ThresholdMethod,ProbabilityThresholdValue,CILevel);if(!is.null(result))threshold_rows[[name]]<-dplyr::bind_cols(tibble::tibble(Source="Model probability",Model=name,Threshold=result$Threshold,ThresholdScale="Predicted probability",Direction=result$Direction),result$Metrics)};if(biomarker_type=="continuous"){result<-.BiomarkerThreshold(y,df_analysis[[biomarker_var]],ThresholdMethod,RawThresholdValue,CILevel);if(!is.null(result))RawBiomarkerThreshold<-dplyr::bind_cols(tibble::tibble(Source="Raw biomarker",Model="Biomarker",Threshold=result$Threshold,ThresholdScale=biomarker_label,Direction=result$Direction),result$Metrics)}else if(nlevels(df_analysis[[biomarker_var]])==2){level<-levels(df_analysis[[biomarker_var]])[2];RawBiomarkerThreshold<-dplyr::bind_cols(tibble::tibble(Source="Raw biomarker",Model="Biomarker",Threshold=NA_real_,ThresholdScale=paste0(biomarker_label," = ",level),Direction=paste0("Positive level: ",level)),.BiomarkerDiagnostic(y,df_analysis[[biomarker_var]]==level,CILevel))};ThresholdTable<-dplyr::bind_rows(c(threshold_rows,if(is.null(RawBiomarkerThreshold))list()else list(RawBiomarkerThreshold)))
    roc_rows<-lapply(names(predictions),function(name){if(!all(is.finite(predictions[[name]]))||length(unique(predictions[[name]]))<2)return(NULL);roc<-pROC::roc(y,predictions[[name]],levels=c(0,1),direction="<",quiet=TRUE);coords<-pROC::coords(roc,"all",ret=c("specificity","sensitivity","threshold"),transpose=FALSE);tibble::tibble(Model=name,Specificity=as.numeric(coords$specificity),Sensitivity=as.numeric(coords$sensitivity),Threshold=as.numeric(coords$threshold),FalsePositiveRate=1-as.numeric(coords$specificity))});ROCData<-dplyr::bind_rows(roc_rows);if(nrow(ROCData))Plots$ROC<-ggplot2::ggplot(ROCData,ggplot2::aes(.data$FalsePositiveRate,.data$Sensitivity,color=.data$Model))+ggplot2::geom_line()+ggplot2::geom_abline(intercept=0,slope=1,linetype=2)+ggplot2::coord_equal()+ggplot2::theme_minimal()+ggplot2::labs(x="1 - Specificity",y="Sensitivity",title=paste0("ROC curves for ",biomarker_label))
    values<-if(biomarker_type=="continuous")seq(stats::quantile(df_analysis[[biomarker_var]],.02),stats::quantile(df_analysis[[biomarker_var]],.98),length.out=100)else levels(df_analysis[[biomarker_var]]);MarginalPredictionData<-dplyr::bind_rows(lapply(values,function(value){newdata<-df_analysis;newdata[[biomarker_var]]<-if(biomarker_type=="continuous")value else factor(value,levels=levels(df_analysis[[biomarker_var]]));tibble::tibble(BiomarkerValue=value,Probability=mean(stats::predict(models$Adjusted,newdata,type="response")))}));ProbabilityPlot<-ggplot2::ggplot(MarginalPredictionData,ggplot2::aes(.data$BiomarkerValue,.data$Probability));if(biomarker_type=="continuous")ProbabilityPlot<-ProbabilityPlot+ggplot2::geom_line()else ProbabilityPlot<-ProbabilityPlot+ggplot2::geom_point(size=3);Plots$Probability<-ProbabilityPlot+ggplot2::theme_minimal()+ggplot2::labs(x=biomarker_label,y=paste0("Adjusted probability of ",positive_level))
    groups<-min(max(1,CalibrationGroups),length(unique(predictions$Adjusted)));CalibrationData<-tibble::tibble(Observed=y,Predicted=predictions$Adjusted,Group=if(groups>1)dplyr::ntile(predictions$Adjusted,groups)else 1L)%>%dplyr::group_by(.data$Group)%>%dplyr::summarise(MeanPredicted=mean(.data$Predicted),ObservedRate=mean(.data$Observed),N=dplyr::n(),.groups="drop");calibration_raw<-tibble::tibble(Observed=y,Predicted=predictions$Adjusted);Plots$Calibration<-ggplot2::ggplot(calibration_raw,ggplot2::aes(.data$Predicted,.data$Observed))+ggplot2::geom_abline(intercept=0,slope=1,linetype=2)+ggplot2::geom_smooth(method="loess",formula=y~x,se=TRUE)+ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_minimal()+ggplot2::labs(x="Predicted probability",y="Observed probability",title="Smooth calibration");Plots$CalibrationGrouped<-ggplot2::ggplot(CalibrationData,ggplot2::aes(.data$MeanPredicted,.data$ObservedRate))+ggplot2::geom_abline(intercept=0,slope=1,linetype=2)+ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::coord_equal(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_minimal()+ggplot2::labs(x="Mean predicted probability",y="Observed event rate",title="Grouped calibration")}
  if(!binary && fits$Adjusted$Status == "ok") { PlotData <- tibble::tibble(Observed = y, Fitted = predictions$Adjusted, Residual = stats::residuals(models$Adjusted), Biomarker = df_analysis[[biomarker_var]]); Plots$PredictedObserved <- ggplot2::ggplot(PlotData, ggplot2::aes(.data$Fitted, .data$Observed)) + ggplot2::geom_point(alpha = .7) + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) + ggplot2::theme_minimal() + ggplot2::labs(x = "Adjusted model prediction", y = outcome_label, title = "Predicted versus observed outcome"); Plots$Residuals <- ggplot2::ggplot(PlotData, ggplot2::aes(.data$Fitted, .data$Residual)) + ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_point(alpha = .7) + ggplot2::theme_minimal() + ggplot2::labs(x = "Fitted value", y = "Residual", title = "Adjusted model residuals"); if(biomarker_type == "continuous") Plots$BiomarkerRelationship <- ggplot2::ggplot(PlotData, ggplot2::aes(.data$Biomarker, .data$Observed)) + ggplot2::geom_point(alpha = .7) + ggplot2::geom_smooth(method = "loess", formula = y ~ x) + ggplot2::theme_minimal() + ggplot2::labs(x = biomarker_label, y = outcome_label) else Plots$BiomarkerRelationship <- ggplot2::ggplot(PlotData, ggplot2::aes(.data$Biomarker, .data$Observed)) + ggplot2::geom_boxplot() + ggplot2::theme_minimal() + ggplot2::labs(x = biomarker_label, y = outcome_label) }
  Predictions<-tibble::tibble(Row=df_analysis$.scidr_row,Observed=y,ObservedLevel=if(binary)as.character(y_original)else NA_character_,Biomarker=df_analysis[[biomarker_var]],BiomarkerPrediction=predictions$Biomarker,CovariatePrediction=predictions$Covariates,AdjustedPrediction=predictions$Adjusted);formulas_validation<-formulas[!vapply(formulas,is.null,logical(1))];if(!length(covariates))formulas_validation$Adjusted<-NULL;ValidationResult<-.BiomarkerValidation(df_analysis,formulas_validation,y,binary,Validation,CVFolds,BootstrapR,Seed,apparent)
  Metadata<-list(OutcomeVariable=outcome_var,OutcomeLabel=outcome_label,OutcomeType=OutcomeType,PositiveLevel=positive_level,NegativeLevel=negative_level,BiomarkerVariable=biomarker_var,BiomarkerLabel=biomarker_label,BiomarkerType=biomarker_type,Covariates=covariates,Formulas=formulas,ThresholdMethod=ThresholdMethod,ValidationMethod=Validation,CILevel=CILevel,Prevalence=if(binary)mean(y)else NA_real_,N=nrow(df_analysis),NPositive=if(binary)sum(y)else NA_integer_,NNegative=if(binary)sum(!y)else NA_integer_,NMissing=sum(!complete),Relabeled=Relabel,ModelStatus=tibble::tibble(Model=names(fits),Status=vapply(fits,`[[`,character(1),"Status"),Note=vapply(fits,`[[`,character(1),"Note")),MarginalPredictionData=MarginalPredictionData)
  list(Models=models,PerformanceTable=PerformanceTable,RegressionTable=tibble::tibble(),ThresholdTable=ThresholdTable,RawBiomarkerThreshold=RawBiomarkerThreshold,Predictions=Predictions,ROCData=ROCData,CalibrationData=CalibrationData,Validation=ValidationResult$Table,ValidationPredictions=ValidationResult$Predictions,Plots=Plots,Metadata=Metadata)
}
