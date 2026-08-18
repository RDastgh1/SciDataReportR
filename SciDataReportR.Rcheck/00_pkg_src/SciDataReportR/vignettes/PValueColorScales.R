## -----------------------------------------------------------------------------
#| label: VignetteSetup
#| include: false

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)


## -----------------------------------------------------------------------------
#| label: LoadPackages
#| message: false

library(SciDataReportR)
library(ggplot2)
library(dplyr)


## -----------------------------------------------------------------------------
#| label: BasicColorScale

df_ColorExample <- data.frame(
  Estimate = seq(-2, 2, length.out = 80),
  PValue = 10^seq(0, -8, length.out = 80)
)

ggplot(
  df_ColorExample,
  aes(x = Estimate, y = 1, color = PValue)
) +
  geom_point(size = 3) +
  scale_color_pvalue() +
  theme_minimal() +
  labs(x = "Estimate", y = NULL)


## -----------------------------------------------------------------------------
#| label: BasicFillScale

df_FillExample <- expand.grid(
  Outcome = paste("Outcome", 1:4),
  Marker = paste("Marker", 1:5)
) %>%
  dplyr::mutate(PValue = 10^seq(0, -6, length.out = dplyr::n()))

ggplot(
  df_FillExample,
  aes(x = Marker, y = Outcome, fill = PValue)
) +
  geom_tile(color = "white") +
  scale_fill_pvalue() +
  theme_minimal() +
  labs(x = NULL, y = NULL)


## -----------------------------------------------------------------------------
#| label: VolcanoPlot

set.seed(2026)

df_Volcano <- data.frame(
  Feature = paste0("Feature_", seq_len(250)),
  log2FoldChange = stats::rnorm(250),
  PValue = 10^stats::runif(250, min = -7, max = 0)
)

ggplot(
  df_Volcano,
  aes(
    x = log2FoldChange,
    y = -log10(PValue),
    color = PValue
  )
) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_pvalue() +
  theme_minimal() +
  labs(
    x = expression(log[2] * " fold change"),
    y = expression(-log[10] * "(P-value)")
  )


## -----------------------------------------------------------------------------
#| label: EnrichmentBubblePlot

df_Enrichment <- data.frame(
  Pathway = paste("Pathway", LETTERS[1:10]),
  EnrichmentRatio = seq(1.2, 3.8, length.out = 10),
  GeneCount = c(8, 12, 15, 7, 20, 11, 18, 9, 14, 22),
  PValue = c(0.32, 0.12, 0.049, 0.025, 0.011, 0.006, 0.0018, 0.0007, 1e-4, 1e-6)
)

ggplot(
  df_Enrichment,
  aes(
    x = EnrichmentRatio,
    y = reorder(Pathway, EnrichmentRatio),
    size = GeneCount,
    color = PValue
  )
) +
  geom_point() +
  scale_color_pvalue() +
  theme_minimal() +
  labs(x = "Enrichment ratio", y = NULL, size = "Gene count")


## -----------------------------------------------------------------------------
#| label: SupportedPalettes

ggplot(
  df_ColorExample,
  aes(x = Estimate, y = PValue, color = PValue)
) +
  geom_point(size = 3) +
  scale_color_pvalue(palette = "viridis") +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "P-value")


## -----------------------------------------------------------------------------
#| label: AdvancedScaleOptions

ggplot(
  df_Enrichment,
  aes(
    x = EnrichmentRatio,
    y = reorder(Pathway, EnrichmentRatio),
    fill = PValue
  )
) +
  geom_point(shape = 21, size = 6, color = "grey20") +
  scale_fill_pvalue(
    palette = "inferno",
    direction = 1,
    limits = c(1e-6, 1),
    breaks = c(1, 0.05, 0.01, 0.001, 1e-6),
    labels = c("1", "0.05", "0.01", "0.001", "1e-6")
  ) +
  theme_minimal() +
  labs(x = "Enrichment ratio", y = NULL)


## -----------------------------------------------------------------------------
#| label: IncorrectTransformedColor
#| eval: false

# # Incorrect: the p-value scale would receive transformed evidence values.
# ggplot(
#   df_Volcano,
#   aes(
#     x = log2FoldChange,
#     y = -log10(PValue),
#     color = -log10(PValue)
#   )
# ) +
#   geom_point() +
#   scale_color_pvalue()


## -----------------------------------------------------------------------------
#| label: CorrectRawPValueColor
#| eval: false

# # Correct: color receives a raw p-value between 0 and 1.
# ggplot(
#   df_Volcano,
#   aes(
#     x = log2FoldChange,
#     y = -log10(PValue),
#     color = PValue
#   )
# ) +
#   geom_point() +
#   scale_color_pvalue()

