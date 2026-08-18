## -----------------------------------------------------------------------------
#| label: LoadPackages
#| echo: false
library(SciDataReportR)
library(dplyr)


## -----------------------------------------------------------------------------
#| label: LoadData
#| echo: false
set.seed(20260802)

df_Raw <- data.frame(
  Group = factor(rep(c("Control", "Disease"), each = 45)),
  SymptomSeverity = c(
    sample(0:3, 45, replace = TRUE, prob = c(0.45, 0.30, 0.18, 0.07)),
    sample(0:3, 45, replace = TRUE, prob = c(0.12, 0.28, 0.35, 0.25))
  ),
  Age = round(stats::rnorm(90, mean = 48, sd = 12))
)

VariableTypes <- data.frame(
  Variable = c("Group", "SymptomSeverity", "Age"),
  Label = c("Participant group", "Symptom Severity", "Age (years)"),
  Type = c("Categorical", "Ordinal", "Continuous"),
  Recode = c("No", "Yes", "No"),
  Code = c(NA, "0=None;1=Mild;2=Moderate;3=Severe", NA)
)

df_Revalued <- RevalueData(df_Raw, VariableTypes)$RevaluedData


## -----------------------------------------------------------------------------
#| label: CategoricalComparisonTable
tbl_Categorical <- MakeComparisonTable(
  data = df_Revalued,
  group_var = "Group",
  variables = "SymptomSeverity",
  TreatOrdinalAs = "Categorical"
)

tbl_Categorical


## -----------------------------------------------------------------------------
#| label: ContinuousComparisonTable
tbl_Continuous <- MakeComparisonTable(
  data = df_Revalued,
  group_var = "Group",
  variables = "SymptomSeverity",
  TreatOrdinalAs = "Continuous"
)

tbl_Continuous


## -----------------------------------------------------------------------------
#| label: BothComparisonTable
tbl_Both <- MakeComparisonTable(
  data = df_Revalued,
  group_var = "Group",
  variables = "SymptomSeverity",
  TreatOrdinalAs = "Both"
)

tbl_Both


## -----------------------------------------------------------------------------
#| label: CategoricalOrdinalDistribution
PlotCategoricalDistributions(
  data = df_Revalued,
  variables = "SymptomSeverity",
  TreatOrdinalAs = "Categorical"
)


## -----------------------------------------------------------------------------
#| label: ContinuousOrdinalDistribution
PlotContinuousDistributions(
  data = df_Revalued,
  variables = "SymptomSeverity",
  TreatOrdinalAs = "Continuous"
)


## -----------------------------------------------------------------------------
#| label: OrdinalMiningMatrix
#| column: screen
vars_Mining <- c("Group", "SymptomSeverity", "Age")

MiningObj <- PlotMiningMatrix(
  data = df_Revalued,
  outcome_vars = vars_Mining,
  TreatOrdinalAs = "Both",
  Relabel = TRUE
)

MiningObj$Unadjusted$plot


## -----------------------------------------------------------------------------
#| label: Reproducibility
#| echo: false
#| output: false
# save.image()
print(sessionInfo())

