## -----------------------------------------------------------------------------
#| label: LoadPackages
library(SciDataReportR)
library(dplyr)
library(plotly)


## -----------------------------------------------------------------------------
#| label: LoadData
data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData


## -----------------------------------------------------------------------------
#| label: DefineMiningVariables
vars_mining <- c(
  "Diagnosis",
  "age",
  "sex",
  "AXL",
  "Calbindin",
  "Ferritin",
  "MMP7"
)

length(vars_mining)


## -----------------------------------------------------------------------------
#| label: PlotTargetedMiningMatrix
#| column: screen
#| fig-height: 5
#| fig-width: 7
MiningObj <- PlotMiningMatrix(
  data = df_Revalued,
  outcome_vars = vars_mining
)

MiningObj$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 5
#| fig-width: 7
#| 
MiningParametric <- PlotMiningMatrix(
  data = df_Revalued,
  outcome_vars = vars_mining,
  Parametric = TRUE
)

MiningParametric$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 5
#| fig-width: 7


MiningNonParametric <- PlotMiningMatrix(
  data = df_Revalued,
  outcome_vars = vars_mining,
  Parametric = FALSE
)

MiningNonParametric$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 5
#| fig-width: 7
MiningObj$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 5
#| fig-width: 7
MiningObj$FDRCorrected$plot


## -----------------------------------------------------------------------------
head(
  MiningObj$Unadjusted$PvalTable
)


## -----------------------------------------------------------------------------
#| eval: false
#| column: screen
#| fig-height: 8
#| fig-width: 12

# # Not evaluated in the shipped vignette to keep it lightweight - run interactively
# ggplotly(
#   MiningObj$Unadjusted$plot
# )


## -----------------------------------------------------------------------------
selected_vars <- c(
  "Diagnosis",
  "age",
  "sex",
  "AXL",
  "Calbindin",
  "Ferritin",
  "MMP7"
)

MiningSubset <- PlotMiningMatrix(
  Data = df_Revalued,
  OutcomeVars = selected_vars
)

MiningSubset$Unadjusted$plot


## -----------------------------------------------------------------------------
sessionInfo()

