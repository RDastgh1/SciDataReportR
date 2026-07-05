## -----------------------------------------------------------------------------
library(SciDataReportR)
library(dplyr)
library(plotly)


## -----------------------------------------------------------------------------
data("SampleData")
data("SampleVariableTypes")

RevaluedObj <- RevalueData(
  SampleData,
  SampleVariableTypes
)

df_Revalued <- RevaluedObj$RevaluedData


## -----------------------------------------------------------------------------
vars_mining <- names(df_Revalued)

length(vars_mining)


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 8
#| fig-width: 12
MiningObj <- PlotMiningMatrix(
  Data = df_Revalued,
  OutcomeVars = vars_mining
)

MiningObj$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 8
#| fig-width: 12
#| 
MiningParametric <- PlotMiningMatrix(
  Data = df_Revalued,
  OutcomeVars = vars_mining,
  Parametric = TRUE
)

MiningParametric$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 8
#| fig-width: 12


MiningNonParametric <- PlotMiningMatrix(
  Data = df_Revalued,
  OutcomeVars = vars_mining,
  Parametric = FALSE
)

MiningNonParametric$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 8
#| fig-width: 12
MiningObj$Unadjusted$plot


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 8
#| fig-width: 12
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

