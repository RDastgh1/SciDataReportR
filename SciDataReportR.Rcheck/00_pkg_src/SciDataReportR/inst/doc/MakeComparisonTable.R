## -----------------------------------------------------------------------------
library(SciDataReportR)
library(dplyr)


## -----------------------------------------------------------------------------
data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
) 
df_Revalued <- RevaluedObj$RevaluedData


## -----------------------------------------------------------------------------
tbl_basic <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  )
)

tbl_basic


## -----------------------------------------------------------------------------
tbl_effectsize <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddEffectSize = TRUE
)

tbl_effectsize


## -----------------------------------------------------------------------------
tbl_nonparametric <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "age",
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  Parametric = FALSE,
  AddEffectSize = TRUE
)

tbl_nonparametric


## -----------------------------------------------------------------------------
tbl_covariate <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Diagnosis",
  Variables = c(
    "sex",
    "Genotype",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  Covariates = "age",
  AddEffectSize = TRUE
)

tbl_covariate


## -----------------------------------------------------------------------------
#| column: screen
tbl_pairwise <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Genotype",
  Variables = c(
    "Diagnosis",
    "age",
    "sex",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddPairwise = TRUE,
  AddEffectSize = TRUE
)

tbl_pairwise


## -----------------------------------------------------------------------------
tbl_pairwise_referent <- MakeComparisonTable(
  df_Revalued,
  CompVariable = "Genotype",
  Variables = c(
    "Diagnosis",
    "age",
    "sex",
    "AXL",
    "Calbindin",
    "Ferritin",
    "MMP7"
  ),
  AddPairwise = TRUE,
  AddEffectSize = TRUE,
  Referent = "E3E3"
)

tbl_pairwise_referent


## -----------------------------------------------------------------------------
sessionInfo()

