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
biomarkers <- c(
  "AXL",
  "Calbindin",
  "Ferritin",
  "MMP7",
  "MMP10",
  "NFL",
  "GFAP",
  "IL6",
  "VCAM1",
  "ICAM1",
  "CRP",
  "TNFa",
  "MCP1",
  "YKL40",
  "Tau"
)


## -----------------------------------------------------------------------------
CorrObj <- PlotCorrelationsHeatmap(
  Data = df_Revalued,
  xVars = biomarkers
)

CorrObj$Unadjusted$plot


## -----------------------------------------------------------------------------
clinical_vars <- c(
  "age",
  "Insulin"
)


## -----------------------------------------------------------------------------
CorrObj_Clinical <- PlotCorrelationsHeatmap(
  Data = df_Revalued,
  xVars = biomarkers,
  yVars = clinical_vars
)

CorrObj_Clinical$Unadjusted$plot


## -----------------------------------------------------------------------------
PearsonObj <- PlotCorrelationsHeatmap(
  Data = df_Revalued,
  xVars = biomarkers,
  method = "pearson"
)

PearsonObj$Unadjusted$plot


## -----------------------------------------------------------------------------
SpearmanObj <- PlotCorrelationsHeatmap(
  Data = df_Revalued,
  xVars = biomarkers,
  method = "spearman"
)

SpearmanObj$Unadjusted$plot


## -----------------------------------------------------------------------------
CorrObj$Unadjusted$plot


## -----------------------------------------------------------------------------
CorrObj$FDRCorrected$plot


## -----------------------------------------------------------------------------
add_r_and_stars(
  CorrObj,
  star_from = "raw"
)


## -----------------------------------------------------------------------------
add_r_and_stars(
  CorrObj,
  star_from = "fdr"
)


## -----------------------------------------------------------------------------
CorrObj_AgeAdjusted <- PlotCorrelationsHeatmap(
  Data = df_Revalued,
  xVars = biomarkers,
  covars = "age"
)

CorrObj_AgeAdjusted$FDRCorrected$plot


## -----------------------------------------------------------------------------
#| eval: false

# # Not evaluated in the shipped vignette to keep it lightweight - run interactively
# ggplotly(
#   CorrObj$FDRCorrected$plot
# )


## -----------------------------------------------------------------------------
SigPlots <- plotSigCorrelations(
  DataFrame = df_Revalued,
  CorrelationHeatmapObject = CorrObj
)


## -----------------------------------------------------------------------------
length(SigPlots)


## -----------------------------------------------------------------------------
#| column: screen
#| fig-height: 9
#| fig-width: 9
AssemblePlots(
  SigPlots,
  LegendPosition = "none",
  ncol = 3
)


## -----------------------------------------------------------------------------
head(
  CorrObj$Unadjusted$r
)


## -----------------------------------------------------------------------------
head(
  CorrObj$FDRCorrected$p
)


## -----------------------------------------------------------------------------
sessionInfo()

