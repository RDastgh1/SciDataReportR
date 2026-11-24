ExplanationFixture <- function(...) {
  set.seed(42)
  df_Data <- data.frame(x=rnorm(40), y=rnorm(40), Category=rep(c("A","B"),20), Date=as.Date("2020-01-01")+seq_len(40), Constant=1)
  df_Data$y[1] <- 12
  DetectDataAnomalies(df_Data, methods="isolation", seed=27,
    method_args=list(isolation=list(ntrees=20,nthreads=1)), ...)
}

test_that("explanations reconstruct scores, preserve RNG and match the backend adapter", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("isotree")
  args <- list(method="shap",rows=c(1,4), background_n=20, nsim=8, repeats=2)
  result <- ExplanationFixture(explanation_args=args)
  before <- .Random.seed
  original <- result$ReviewData
  result2 <- DetectDataAnomalies(as.data.frame(original[setdiff(names(original),".Row")]),methods="isolation",seed=27,
    method_args=list(isolation=list(ntrees=20,nthreads=1)),explanation_args=args)
  expect_identical(.Random.seed,before)
  expect_equal(result$FeatureContributions,result2$FeatureContributions)
  summary <- result$ExplanationSummary
  explained <- summary[summary$Status=="Explained",]
  expect_setequal(explained$.Row,c(1,4))
  expect_true(all(abs(explained$ReconstructionResidual)<1e-6))
  expect_true(all(summary$Status[!summary$.Row %in% c(1,4)]=="Not requested"))
  expect_false("Constant" %in% result$FeatureContributions$Variable)
  expect_s3_class(result$ReviewData$Date,"Date")
  metadata <- result$ExplanationMetadata$Strata[[1]]
  df_Model <- as.data.frame(result$ReviewData[metadata$ModelingVariables])
  df_Model$Date <- ScidrAnomalyDateDays(df_Model$Date)-result$Metadata$DateConversions$Date$ReferenceDays
  df_Model$Category <- factor(df_Model$Category)
  fit <- isotree::isolation.forest(df_Model,seed=27,ntrees=20,nthreads=1)
  reference <- match(metadata$BackgroundRows,result$ReviewData$.Row)
  baseline <- mean(predict(fit,df_Model[reference,],type="score"))
  contributions <- vapply(metadata$Seeds[["1"]],function(seed) {
    set.seed(seed)
    as.numeric(fastshap::explain(fit,X=df_Model[reference,],newdata=df_Model[1,,drop=FALSE],nsim=8,adjust=TRUE,baseline=baseline,
      pred_wrapper=function(object,newdata) as.numeric(predict(object,newdata,type="score")),parallel=FALSE)[1,metadata$ModelingVariables])
  },numeric(ncol(df_Model)))
  expect_equal(result$FeatureContributions$Contribution[result$FeatureContributions$.Row==1],rowMeans(contributions),tolerance=1e-10)
  disabled <- ExplanationFixture(explain_isolation=FALSE)
  expect_equal(result$Results,disabled$Results)
  expect_equal(ScidrReviewSnapshot(result),ScidrReviewSnapshot(disabled))
})

test_that("selection, unavailable records and explanation failures remain distinct", {
  skip_if_not_installed("fastshap")
  result <- ExplanationFixture(explanation_args=list(method="shap",top_n=1,nsim=4,repeats=1))
  expect_equal(sum(result$ExplanationSummary$Status=="Explained"),1)
  expect_true(all(is.na(result$FeatureContributions$SimulationSD)))
  expect_error(ExplanationFixture(explanation_args=list(method="shap",nsim=1)),"nsim")
  expect_error(ExplanationFixture(explanation_args=list(method="shap",rows=c(0,1))),"row references")
  expect_error(ExplanationFixture(explanation_args=list(method="shap",unknown=1)),"Unknown")
  testthat::local_mocked_bindings(explain=function(...) stop("Synthetic explanation failure"),.package="fastshap")
  failed <- ExplanationFixture(explanation_args=list(method="shap",top_n=1,nsim=4,repeats=1))
  expect_true(all(failed$Results$Status=="evaluated"))
  expect_equal(sum(failed$ExplanationSummary$Status=="Unavailable"),1)
  expect_true("explanation_failure" %in% failed$Diagnostics$Issue)
})

test_that("contribution workbooks protect a separate explanation fingerprint and legacy imports", {
  skip_if_not_installed("fastshap")
  skip_if_not_installed("openxlsx")
  result <- ExplanationFixture(explanation_args=list(method="shap",top_n=1,nsim=4,repeats=2))
  queue <- CreateAnomalyReviewQueue(result,top_n=1)
  original <- queue
  meta <- ScidrReviewExplanationMeta(queue)
  expect_null(names(meta[[1]]$values[[1]]$Observed))
  expect_type(meta[[1]]$values[[1]]$Observed,"character")
  file <- tempfile(fileext=".xlsx")
  on.exit(unlink(file),add=TRUE)
  ExportAnomalyReviewQueue(queue,file, layout = "full")
  expect_true("Contributions" %in% openxlsx::getSheetNames(file))
  expect_identical(ImportAnomalyReviewDecisions(queue,file),queue)
  EditReviewWorkbook(file,function(d) {d[[ScidrReviewHeaders("Contribution")]][1] <- "999";d},"Contributions")
  expect_error(ImportAnomalyReviewDecisions(queue,file),"Contributions.*changed")
  ExportAnomalyReviewQueue(queue,file,overwrite=TRUE, layout = "full")
  wb <- openxlsx::loadWorkbook(file);openxlsx::removeWorksheet(wb,"Contributions");openxlsx::saveWorkbook(wb,file,overwrite=TRUE)
  EditReviewWorkbook(file,function(d) d[setdiff(names(d),"ExplanationFingerprint")],"Manifest")
  expect_identical(ImportAnomalyReviewDecisions(queue,file),queue)
  expect_identical(queue,original)
  changed <- queue;changed$Sources$Analysis$FeatureContributions$Contribution[1] <- 999
  expect_silent(ScidrReviewValidate(changed))
  expect_error(ScidrReviewContributions(changed),"fingerprint changed")
  legacy <- queue;legacy$Sources$Analysis$FeatureContributions <- NULL;legacy$Sources$Analysis$ExplanationSummary <- NULL;legacy$Sources$Analysis$ExplanationMetadata <- NULL
  meta <- ScidrReviewExplanationMeta(legacy)
  expect_match(meta[[1]]$status,"rerun detection")
})

test_that("stratum backgrounds, missing inputs and original date classes are retained", {
  skip_if_not_installed("fastshap")
  set.seed(80)
  df_Data <- data.frame(ID=rep(1:20,2), Group=rep(c("A","B"),each=20),
    x=rnorm(40),y=rnorm(40),Category=rep(c("yes","no"),20),Date=as.POSIXct("2020-01-01",tz="America/New_York")+seq_len(40)*86400)
  df_Data$x[2] <- NA
  df_Data[1,c("x","y","Category","Date")] <- NA
  before <- df_Data
  result <- DetectDataAnomalies(df_Data,id_var="ID",stratify_by="Group",methods="isolation",
    method_args=list(isolation=list(ntrees=20,nthreads=1)),explanation_args=list(method="shap",rows=c(1,2,22),nsim=4,repeats=2))
  expect_identical(df_Data,before)
  expect_identical(result$ReviewData$Date,df_Data$Date)
  expect_equal(result$ExplanationSummary$Status[1],"Unavailable")
  expect_false(any(c("ID","Group") %in% result$FeatureContributions$Variable))
  expect_setequal(result$ExplanationSummary$.Row[result$ExplanationSummary$Status=="Explained"],c(2,22))
  for (metadata in result$ExplanationMetadata$Strata) {
    expect_equal(length(unique(df_Data$Group[metadata$BackgroundRows])),1)
    expect_false(1 %in% metadata$BackgroundRows)
    expect_setequal(unique(df_Data$Group[metadata$SelectedRows]), unique(df_Data$Group[metadata$BackgroundRows]))
  }
})


test_that("unrevalued tibbles retain nonsyntactic feature names in explanations", {
  skip_if_not_installed("fastshap")
  set.seed(11)
  df_Raw <- tibble::tibble(`measurement x`=rnorm(30),`measurement y`=rnorm(30),`response group`=rep(c("A","B"),15))
  before <- df_Raw
  result <- DetectDataAnomalies(df_Raw,methods="isolation",method_args=list(isolation=list(ntrees=20,nthreads=1)),
    explanation_args=list(method="shap",top_n=1,nsim=4,repeats=2))
  expect_equal(sum(result$ExplanationSummary$Status=="Explained"),1)
  expect_setequal(result$FeatureContributions$Variable,names(df_Raw))
  expect_identical(df_Raw,before)
})
