# Current API compatibility guide

Use current names in all new code. The generated API reference identifies every argument explicitly marked deprecated in its Rd source; this guide highlights the legacy interfaces most likely to appear in older examples.

| Current function | Use now | Do not use in new code |
| --- | --- | --- |
| `RevalueData()` | `data`, `codebook` | `DatatoRevalue`, `VarTypes` |
| `MakeComparisonTable()` | `data`, `group_var`, `variables`, `covariates`, `value_digits`, `p_digits`, `effect_size_digits` | `DataFrame`, `CompVariable`, `Variables`, `Covariates`, `ValueDigits`, `pDigits`, `EffectSizeDigits` |
| `PlotCorrelationsHeatmap()` | `data`, `predictor_vars`, `outcome_vars`, `covariates`, `TreatOrdinalAs` | `Data`, `xVars`, `yVars`, `covars`, `Ordinal` |
| `MakeUnivariateRegressionTable()` | `data`, `outcome_vars`, `predictor_vars`, `covariates` | `Data`, `OutcomeVars`, `PredictorVars`, `Covars` |
| `MultivariableRegressionTable()` | `data`, `outcome_vars`, `predictor_vars`, `covariates` | `Data`, `OutcomeVars`, `PredictorVars`, `Covars` |
| `PlotContinuousDistributions()` / `PlotMissingData()` | `data`, `variables` | `DataFrame`, `Variables` |
| `CreateSummaryTable()` | `data`, `variables`, `digits`, `TreatOrdinalAs` | `Data`, `Variables`, `numdecimals`, `Ordinal` |
| `CreatePCAObject()` | the current constructor | `CreatePCATable()` compatibility alias |
| `CreateClusterModel_SOM_MClust()` | the current constructor | `Pipeline_SOMClust()` deprecated alias |
| `PlotForestFromTable()` | the current plot function | `plotForestFromTable()` compatibility alias |

Older code may remain runnable for compatibility. Preserve it unless the task calls for modernization, but never extend it with deprecated interfaces.
