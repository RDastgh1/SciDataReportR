# PlotMiningMatrix

Generate a matrix of statistical relationships between variables.

## Usage

``` r
PlotMiningMatrix(
  data,
  outcome_vars,
  predictor_vars = NULL,
  covariates = NULL,
  Relabel = TRUE,
  TreatOrdinalAs = "Categorical",
  Parametric = TRUE,
  Data = lifecycle::deprecated(),
  OutcomeVars = lifecycle::deprecated(),
  PredictorVars = lifecycle::deprecated(),
  fdr_scope = c("matrix", "per_outcome", "per_predictor"),
  IncludeTrendSymbol = FALSE,
  Covariates = lifecycle::deprecated(),
  Pthresh = 0.05,
  PVar = "PValue"
)
```

## Arguments

- data:

  A data frame.

- outcome_vars:

  Outcome variables.

- predictor_vars:

  Predictor variables. If NULL, uses OutcomeVars.

- covariates:

  Optional variables to exclude from the tested axes.

- Relabel:

  Use labels instead of names.

- TreatOrdinalAs:

  How ordinal variables are handled: `"Categorical"`, `"Continuous"`,
  `"Both"`, or `"Exclude"`.

- Parametric:

  Use parametric tests.

- Data:

  **Deprecated** (since 19.15.0). Use `data` instead.

- OutcomeVars:

  **Deprecated** (since 19.15.0). Use `outcome_vars` instead.

- PredictorVars:

  **Deprecated** (since 19.15.0). Use `predictor_vars` instead.

- fdr_scope:

  One of `"matrix"`, `"per_outcome"`, or `"per_predictor"`. Correction
  uses the symmetrized table before selecting the displayed grid.
  `"matrix"` corrects all p-values together (historical behavior).
  `"per_outcome"` groups by x-axis outcome (`XVar`), correcting over
  predictors; `"per_predictor"` groups by y-axis predictor (`YVar`),
  correcting over outcomes.

- IncludeTrendSymbol:

  Logical; add a plus sign for `0.05 < p <= 0.10`. Defaults to FALSE;
  the trend legend entry is omitted when disabled.

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

- Pthresh:

  Significance threshold, default `0.05`. An explicitly supplied value
  takes precedence over a legacy significance threshold argument.

- PVar:

  P-value column used for significance decisions. Defaults to raw
  `"PValue"`; select `"PValueFDR"` to use Benjamini-Hochberg adjusted
  values. Explicit legacy p-column names remain accepted. Existing,
  explicitly named FDR plots and adjusted-value columns are retained.

  Explicit `Pthresh` takes precedence over `Alpha`, `FDRAlpha`,
  `label_q`, or `adjusted_significance_threshold`. Pairwise
  `star_p = "adjusted"` retains its legacy adjustment method unless
  `PVar` is supplied; `PVar = "PValueFDR"` always selects BH.
  `star_p = "none"` suppresses legacy pairwise stars. Stars retain
  conventional p-value tiers within the requested cutoff.

## Value

List with tables and plots.

## Details

Numeric pairs use Pearson correlation when `Parametric = TRUE`, and
Spearman otherwise. Categorical/numeric pairs use ANOVA or
Kruskal-Wallis in either orientation; categorical pairs use chi-squared
tests. Logical columns are converted locally to FALSE/TRUE factors,
preserving labels. Missing effect sizes retain their p-values and plot
in neutral grey. Colors represent absolute effect sizes on a fixed 0–1
scale with a square-root transformation to emphasize smaller effects.
Values outside these limits are squished to the endpoints for color
mapping only; legend labels and returned estimates remain in original
effect-size units. Significance uses the p-values selected by `PVar`:
circles for p \> 0.05, triangles for 0.01 \< p \<= 0.05, squares for
0.001 \< p \<= 0.01, and diamonds for p \<= 0.001. With
`IncludeTrendSymbol = TRUE`, plus signs represent 0.05 \< p \<= 0.10 and
circles represent p \> 0.10. Marker sizes are fixed across plots: circle
2, plus 2.7, triangle 3.5, square 4, and diamond 6. The legend includes
all enabled bands at these sizes. `Unadjusted$Missing` lists requested
non-self pairs without a finite p-value.

## Statistical results and provenance

See
[inference-results](https://rdastgh1.github.io/SciDataReportR/reference/inference-results.md)
for the shared `Results` and `Metadata` contract.

## Examples

``` r
data(SampleData)
data(SampleVariableTypes)

# Attach labels and factor levels for readable axes
Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData

# A mining matrix over 11 mixed-type variables (categorical + continuous)
result <- PlotMiningMatrix(
  Labelled,
  outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
  predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                     "Apolipoprotein_B", "C_Reactive_Protein",
                     "Cortisol", "Insulin")
)

# The relationship plot (point shape/size encode significance from raw p)
result$Unadjusted$plot


# The p-value table carries both unadjusted (p) and FDR-adjusted (p_adj)
# p-values so results can be inspected with and without FDR correction
result$Unadjusted$PvalTable[, c("XVar", "YVar", "p", "p_adj", "Test")]
#> # A tibble: 30 × 5
#>    XVar        YVar                         p    p_adj Test       
#>    <chr>       <chr>                    <dbl>    <dbl> <chr>      
#>  1 AXL         Alpha_1_Antitrypsin   1.12e- 1 1.87e- 1 Correlation
#>  2 AXL         Alpha_2_Macroglobulin 8.35e-15 8.35e-14 Correlation
#>  3 AXL         Apolipoprotein_B      6.02e- 1 6.95e- 1 Correlation
#>  4 AXL         C_Reactive_Protein    4.82e- 1 6.48e- 1 Correlation
#>  5 AXL         Cortisol              1.85e- 2 4.07e- 2 Correlation
#>  6 AXL         Insulin               2.29e-14 1.71e-13 Correlation
#>  7 Adiponectin Alpha_1_Antitrypsin   1.21e-15 1.82e-14 Correlation
#>  8 Adiponectin Alpha_2_Macroglobulin 3.24e-12 1.95e-11 Correlation
#>  9 Adiponectin Apolipoprotein_B      2.00e-18 6.00e-17 Correlation
#> 10 Adiponectin C_Reactive_Protein    1.02e- 1 1.80e- 1 Correlation
#> # ℹ 20 more rows

# Per-outcome FDR correction instead of matrix-wide
result_perout <- PlotMiningMatrix(
  Labelled,
  outcome_vars   = c("Diagnosis", "sex", "age", "AXL", "Adiponectin"),
  predictor_vars = c("Alpha_1_Antitrypsin", "Alpha_2_Macroglobulin",
                     "Apolipoprotein_B", "C_Reactive_Protein",
                     "Cortisol", "Insulin"),
  fdr_scope = "per_outcome"
)

# An interactive version of the matrix
if (requireNamespace("plotly", quietly = TRUE)) {
  plotly::ggplotly(result$Unadjusted$plot)
}

{"x":{"data":[{"x":[4,4,4,5,1,3,3,2,3,1,3,2,3,2,3],"y":[6,4,3,3,1,6,5,5,4,3,3,3,2,2,1],"text":["XLabel: AXL receptor tyrosine kinase<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.29541085<br />stars: ns<br />stars: ns","XLabel: AXL receptor tyrosine kinase<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.16923963<br />stars: ns<br />stars: ns","XLabel: AXL receptor tyrosine kinase<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.19659438<br />stars: ns<br />stars: ns","XLabel: Adiponectin<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.29959453<br />stars: ns<br />stars: ns","XLabel: Diagnosis<br />YLabel: Insulin<br />EffectSizeAbs: 0.06324555<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.14230234<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.22679854<br />stars: ns<br />stars: ns","XLabel: Sex<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.10488088<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.18860032<br />stars: ns<br />stars: ns","XLabel: Diagnosis<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.03162278<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.22880690<br />stars: ns<br />stars: ns","XLabel: Sex<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.00717635<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: Cortisol<br />EffectSizeAbs: 0.05298856<br />stars: ns<br />stars: ns","XLabel: Sex<br />YLabel: Cortisol<br />EffectSizeAbs: 0.01220656<br />stars: ns<br />stars: ns","XLabel: Age<br />YLabel: Insulin<br />EffectSizeAbs: 0.18325466<br />stars: ns<br />stars: ns"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(197,162,165,1)","rgba(206,185,187,1)","rgba(204,180,182,1)","rgba(197,161,165,1)","rgba(213,205,206,1)","rgba(208,190,192,1)","rgba(202,175,177,1)","rgba(210,197,198,1)","rgba(205,182,184,1)","rgba(215,211,211,1)","rgba(202,174,177,1)","rgba(217,216,216,1)","rgba(214,207,208,1)","rgba(216,215,215,1)","rgba(205,183,185,1)"],"opacity":1,"size":7.559055118110237,"symbol":"circle","line":{"width":3.7795275590551185,"color":["rgba(197,162,165,1)","rgba(206,185,187,1)","rgba(204,180,182,1)","rgba(197,161,165,1)","rgba(213,205,206,1)","rgba(208,190,192,1)","rgba(202,175,177,1)","rgba(210,197,198,1)","rgba(205,182,184,1)","rgba(215,211,211,1)","rgba(202,174,177,1)","rgba(217,216,216,1)","rgba(214,207,208,1)","rgba(216,215,215,1)","rgba(205,183,185,1)"]}},"hoveron":"points","name":"ns","legendgroup":"ns","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[4,1,1,2],"y":[2,5,4,4],"text":["XLabel: AXL receptor tyrosine kinase<br />YLabel: Cortisol<br />EffectSizeAbs: 0.35927053<br />stars: *<br />stars: *","XLabel: Diagnosis<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.12247449<br />stars: *<br />stars: *","XLabel: Diagnosis<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.13784049<br />stars: *<br />stars: *","XLabel: Sex<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.12649111<br />stars: *<br />stars: *"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(192,150,154,1)","rgba(209,194,195,1)","rgba(208,191,193,1)","rgba(209,193,195,1)"],"opacity":1,"size":13.228346456692915,"symbol":"triangle-up","line":{"width":3.7795275590551185,"color":["rgba(192,150,154,1)","rgba(209,194,195,1)","rgba(208,191,193,1)","rgba(209,193,195,1)"]}},"hoveron":"points","name":"*","legendgroup":"*","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[5],"y":[2],"text":"XLabel: Adiponectin<br />YLabel: Cortisol<br />EffectSizeAbs: 0.37771007<br />stars: **<br />stars: **","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(191,147,151,1)","opacity":1,"size":15.118110236220474,"symbol":"square","line":{"width":3.7795275590551185,"color":"rgba(191,147,151,1)"}},"hoveron":"points","name":"**","legendgroup":"**","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[4,4,5,5,5,5,1,2,1,2],"y":[5,1,6,5,4,1,6,6,2,1],"text":["XLabel: AXL receptor tyrosine kinase<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.63892354<br />stars: ***<br />stars: ***","XLabel: AXL receptor tyrosine kinase<br />YLabel: Insulin<br />EffectSizeAbs: 0.63406697<br />stars: ***<br />stars: ***","XLabel: Adiponectin<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.64787064<br />stars: ***<br />stars: ***","XLabel: Adiponectin<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.60791833<br />stars: ***<br />stars: ***","XLabel: Adiponectin<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.67458199<br />stars: ***<br />stars: ***","XLabel: Adiponectin<br />YLabel: Insulin<br />EffectSizeAbs: 0.42973618<br />stars: ***<br />stars: ***","XLabel: Diagnosis<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.20736441<br />stars: ***<br />stars: ***","XLabel: Sex<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.22583180<br />stars: ***<br />stars: ***","XLabel: Diagnosis<br />YLabel: Cortisol<br />EffectSizeAbs: 0.19235384<br />stars: ***<br />stars: ***","XLabel: Sex<br />YLabel: Insulin<br />EffectSizeAbs: 0.19235384<br />stars: ***<br />stars: ***"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(168,99,109,1)","rgba(169,100,110,1)","rgba(167,98,107,1)","rgba(171,105,114,1)","rgba(165,93,103,1)","rgba(186,137,143,1)","rgba(203,178,180,1)","rgba(202,175,177,1)","rgba(204,181,183,1)","rgba(204,181,183,1)"],"opacity":1,"size":22.677165354330711,"symbol":"diamond","line":{"width":3.7795275590551185,"color":["rgba(168,99,109,1)","rgba(169,100,110,1)","rgba(167,98,107,1)","rgba(171,105,114,1)","rgba(165,93,103,1)","rgba(186,137,143,1)","rgba(203,178,180,1)","rgba(202,175,177,1)","rgba(204,181,183,1)","rgba(204,181,183,1)"]}},"hoveron":"points","name":"***","legendgroup":"***","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1],"y":[1],"name":"8185922e5ca2f57137399da5e5aab9bc","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#D9D9D9"],[0.0033444816053511705,"#D9D8D8"],[0.006688963210702341,"#D9D8D8"],[0.010033444816053512,"#D8D7D7"],[0.013377926421404682,"#D8D6D7"],[0.016722408026755852,"#D8D6D6"],[0.020066889632107024,"#D8D5D5"],[0.023411371237458192,"#D8D5D5"],[0.026755852842809364,"#D7D4D4"],[0.030100334448160536,"#D7D3D4"],[0.033444816053511704,"#D7D3D3"],[0.036789297658862873,"#D7D2D2"],[0.040133779264214048,"#D7D1D2"],[0.043478260869565216,"#D6D1D1"],[0.046822742474916385,"#D6D0D1"],[0.05016722408026756,"#D6D0D0"],[0.053511705685618728,"#D6CFCF"],[0.056856187290969896,"#D5CECF"],[0.060200668896321072,"#D5CECE"],[0.06354515050167224,"#D5CDCE"],[0.066889632107023408,"#D5CCCD"],[0.070234113712374577,"#D5CCCC"],[0.073578595317725745,"#D4CBCC"],[0.076923076923076927,"#D4CACB"],[0.080267558528428096,"#D4CACB"],[0.083612040133779264,"#D4C9CA"],[0.086956521739130432,"#D4C9C9"],[0.090301003344481601,"#D3C8C9"],[0.093645484949832769,"#D3C7C8"],[0.096989966555183951,"#D3C7C8"],[0.10033444816053512,"#D3C6C7"],[0.10367892976588629,"#D2C5C7"],[0.10702341137123746,"#D2C5C6"],[0.11036789297658862,"#D2C4C5"],[0.11371237458193979,"#D2C4C5"],[0.11705685618729096,"#D2C3C4"],[0.12040133779264214,"#D1C2C4"],[0.12374581939799331,"#D1C2C3"],[0.12709030100334448,"#D1C1C2"],[0.13043478260869565,"#D1C0C2"],[0.13377926421404682,"#D0C0C1"],[0.13712374581939799,"#D0BFC1"],[0.14046822742474915,"#D0BFC0"],[0.14381270903010032,"#D0BEBF"],[0.14715719063545149,"#CFBDBF"],[0.15050167224080269,"#CFBDBE"],[0.15384615384615385,"#CFBCBE"],[0.15719063545150502,"#CFBBBD"],[0.16053511705685619,"#CFBBBD"],[0.16387959866220736,"#CEBABC"],[0.16722408026755853,"#CEBABB"],[0.1705685618729097,"#CEB9BB"],[0.17391304347826086,"#CEB8BA"],[0.17725752508361203,"#CDB8BA"],[0.1806020066889632,"#CDB7B9"],[0.18394648829431437,"#CDB7B8"],[0.18729096989966554,"#CDB6B8"],[0.19063545150501671,"#CCB5B7"],[0.1939799331103679,"#CCB5B7"],[0.19732441471571907,"#CCB4B6"],[0.20066889632107024,"#CCB3B6"],[0.20401337792642141,"#CCB3B5"],[0.20735785953177258,"#CBB2B4"],[0.21070234113712374,"#CBB2B4"],[0.21404682274247491,"#CBB1B3"],[0.21739130434782608,"#CBB0B3"],[0.22073578595317725,"#CAB0B2"],[0.22408026755852842,"#CAAFB2"],[0.22742474916387959,"#CAAEB1"],[0.23076923076923075,"#CAAEB0"],[0.23411371237458192,"#C9ADB0"],[0.23745819397993312,"#C9ADAF"],[0.24080267558528429,"#C9ACAF"],[0.24414715719063546,"#C9ABAE"],[0.24749163879598662,"#C8ABAD"],[0.25083612040133779,"#C8AAAD"],[0.25418060200668896,"#C8AAAC"],[0.25752508361204013,"#C8A9AC"],[0.2608695652173913,"#C7A8AB"],[0.26421404682274247,"#C7A8AB"],[0.26755852842809363,"#C7A7AA"],[0.2709030100334448,"#C7A6A9"],[0.27424749163879597,"#C6A6A9"],[0.27759197324414714,"#C6A5A8"],[0.28093645484949831,"#C6A5A8"],[0.28428093645484948,"#C6A4A7"],[0.28762541806020064,"#C5A3A7"],[0.29096989966555181,"#C5A3A6"],[0.29431438127090298,"#C5A2A5"],[0.29765886287625415,"#C5A2A5"],[0.30100334448160537,"#C4A1A4"],[0.30434782608695654,"#C4A0A4"],[0.30769230769230771,"#C4A0A3"],[0.31103678929765888,"#C49FA3"],[0.31438127090301005,"#C39EA2"],[0.31772575250836121,"#C39EA1"],[0.32107023411371238,"#C39DA1"],[0.32441471571906355,"#C39DA0"],[0.32775919732441472,"#C29CA0"],[0.33110367892976589,"#C29B9F"],[0.33444816053511706,"#C29B9F"],[0.33779264214046822,"#C29A9E"],[0.34113712374581939,"#C19A9D"],[0.34448160535117056,"#C1999D"],[0.34782608695652173,"#C1989C"],[0.3511705685618729,"#C1989C"],[0.35451505016722407,"#C0979B"],[0.35785953177257523,"#C0969B"],[0.3612040133779264,"#C0969A"],[0.36454849498327757,"#C0959A"],[0.36789297658862874,"#BF9599"],[0.37123745819397991,"#BF9498"],[0.37458193979933108,"#BF9398"],[0.37792642140468224,"#BF9397"],[0.38127090301003341,"#BE9297"],[0.38461538461538458,"#BE9296"],[0.38795986622073581,"#BE9196"],[0.39130434782608697,"#BD9095"],[0.39464882943143814,"#BD9095"],[0.39799331103678931,"#BD8F94"],[0.40133779264214048,"#BD8F93"],[0.40468227424749165,"#BC8E93"],[0.40802675585284282,"#BC8D92"],[0.41137123745819398,"#BC8D92"],[0.41471571906354515,"#BC8C91"],[0.41806020066889632,"#BB8B91"],[0.42140468227424749,"#BB8B90"],[0.42474916387959866,"#BB8A8F"],[0.42809364548494983,"#BA8A8F"],[0.43143812709030099,"#BA898E"],[0.43478260869565216,"#BA888E"],[0.43812709030100333,"#BA888D"],[0.4414715719063545,"#B9878D"],[0.44481605351170567,"#B9878C"],[0.44816053511705684,"#B9868C"],[0.451505016722408,"#B9858B"],[0.45484949832775917,"#B8858B"],[0.45819397993311034,"#B8848A"],[0.46153846153846151,"#B88489"],[0.46488294314381268,"#B78389"],[0.46822742474916385,"#B78288"],[0.47157190635451507,"#B78288"],[0.47491638795986624,"#B78187"],[0.47826086956521741,"#B68187"],[0.48160535117056857,"#B68086"],[0.48494983277591974,"#B67F86"],[0.48829431438127091,"#B57F85"],[0.49163879598662208,"#B57E84"],[0.49498327759197325,"#B57E84"],[0.49832775919732442,"#B57D83"],[0.50167224080267558,"#B47C83"],[0.50501672240802675,"#B47C82"],[0.50836120401337792,"#B47B82"],[0.51170568561872909,"#B37A81"],[0.51505016722408026,"#B37A81"],[0.51839464882943143,"#B37980"],[0.52173913043478259,"#B37980"],[0.52508361204013376,"#B2787F"],[0.52842809364548493,"#B2777E"],[0.5317725752508361,"#B2777E"],[0.53511705685618727,"#B1767D"],[0.53846153846153844,"#B1767D"],[0.5418060200668896,"#B1757C"],[0.54515050167224077,"#B1747C"],[0.54849498327759194,"#B0747B"],[0.55183946488294311,"#B0737B"],[0.55518394648829428,"#B0737A"],[0.55852842809364545,"#AF727A"],[0.56187290969899661,"#AF7179"],[0.56521739130434778,"#AF7179"],[0.56856187290969895,"#AF7078"],[0.57190635451505012,"#AE7077"],[0.57525083612040129,"#AE6F77"],[0.57859531772575246,"#AE6E76"],[0.58193979933110362,"#AD6E76"],[0.58528428093645479,"#AD6D75"],[0.58862876254180596,"#AD6D75"],[0.59197324414715713,"#AC6C74"],[0.5953177257525083,"#AC6B74"],[0.59866220735785947,"#AC6B73"],[0.60200668896321075,"#AC6A73"],[0.60535117056856191,"#AB6A72"],[0.60869565217391308,"#AB6972"],[0.61204013377926425,"#AB6871"],[0.61538461538461542,"#AA6871"],[0.61872909698996659,"#AA6770"],[0.62207357859531776,"#AA666F"],[0.62541806020066892,"#A9666F"],[0.62876254180602009,"#A9656E"],[0.63210702341137126,"#A9656E"],[0.63545150501672243,"#A9646D"],[0.6387959866220736,"#A8636D"],[0.64214046822742477,"#A8636C"],[0.64548494983277593,"#A8626C"],[0.6488294314381271,"#A7626B"],[0.65217391304347827,"#A7616B"],[0.65551839464882944,"#A7606A"],[0.65886287625418061,"#A6606A"],[0.66220735785953178,"#A65F69"],[0.66555183946488294,"#A65F69"],[0.66889632107023411,"#A55E68"],[0.67224080267558528,"#A55D68"],[0.67558528428093645,"#A55D67"],[0.67892976588628762,"#A55C67"],[0.68227424749163879,"#A45C66"],[0.68561872909698995,"#A45B65"],[0.68896321070234112,"#A45A65"],[0.69230769230769229,"#A35A64"],[0.69565217391304346,"#A35964"],[0.69899665551839463,"#A35863"],[0.7023411371237458,"#A25863"],[0.70568561872909696,"#A25762"],[0.70903010033444813,"#A25762"],[0.7123745819397993,"#A15661"],[0.71571906354515047,"#A15561"],[0.71906354515050164,"#A15560"],[0.72240802675585281,"#A05460"],[0.72575250836120397,"#A0545F"],[0.72909698996655514,"#A0535F"],[0.73244147157190631,"#9F525E"],[0.73578595317725748,"#9F525E"],[0.73913043478260865,"#9F515D"],[0.74247491638795982,"#9E505D"],[0.74581939799331098,"#9E505C"],[0.74916387959866215,"#9E4F5C"],[0.75250836120401332,"#9E4F5B"],[0.75585284280936449,"#9D4E5B"],[0.75919732441471566,"#9D4D5A"],[0.76254180602006683,"#9D4D5A"],[0.76588628762541799,"#9C4C59"],[0.76923076923076916,"#9C4B59"],[0.77257525083612044,"#9C4B58"],[0.77591973244147161,"#9B4A58"],[0.77926421404682278,"#9B4A57"],[0.78260869565217395,"#9B4957"],[0.78595317725752512,"#9A4856"],[0.78929765886287628,"#9A4856"],[0.79264214046822745,"#9A4755"],[0.79598662207357862,"#994655"],[0.79933110367892979,"#994654"],[0.80267558528428096,"#994554"],[0.80602006688963213,"#984553"],[0.80936454849498329,"#984453"],[0.81270903010033446,"#984352"],[0.81605351170568563,"#974352"],[0.8193979933110368,"#974251"],[0.82274247491638797,"#974151"],[0.82608695652173914,"#964150"],[0.8294314381270903,"#964050"],[0.83277591973244147,"#963F4F"],[0.83612040133779264,"#953F4F"],[0.83946488294314381,"#953E4E"],[0.84280936454849498,"#953D4E"],[0.84615384615384615,"#943D4D"],[0.84949832775919731,"#943C4D"],[0.85284280936454848,"#943C4C"],[0.85618729096989965,"#933B4C"],[0.85953177257525082,"#933A4B"],[0.86287625418060199,"#933A4B"],[0.86622073578595316,"#92394A"],[0.86956521739130432,"#92384A"],[0.87290969899665549,"#913849"],[0.87625418060200666,"#913749"],[0.87959866220735783,"#913648"],[0.882943143812709,"#903648"],[0.88628762541806017,"#903547"],[0.88963210702341133,"#903447"],[0.8929765886287625,"#8F3346"],[0.89632107023411367,"#8F3346"],[0.89966555183946484,"#8F3245"],[0.90301003344481601,"#8E3145"],[0.90635451505016718,"#8E3144"],[0.90969899665551834,"#8E3044"],[0.91304347826086951,"#8D2F43"],[0.91638795986622068,"#8D2F43"],[0.91973244147157185,"#8D2E42"],[0.92307692307692302,"#8C2D42"],[0.92642140468227419,"#8C2C41"],[0.92976588628762535,"#8C2C41"],[0.93311036789297652,"#8B2B40"],[0.93645484949832769,"#8B2A40"],[0.93979933110367886,"#8A2A3F"],[0.94314381270903014,"#8A293F"],[0.94648829431438131,"#8A283F"],[0.94983277591973247,"#89273E"],[0.95317725752508364,"#89273E"],[0.95652173913043481,"#89263D"],[0.95986622073578598,"#88253D"],[0.96321070234113715,"#88243C"],[0.96655518394648832,"#88233C"],[0.96989966555183948,"#87233B"],[0.97324414715719065,"#87223B"],[0.97658862876254182,"#87213A"],[0.97993311036789299,"#86203A"],[0.98327759197324416,"#861F39"],[0.98662207357859533,"#851F39"],[0.98996655518394649,"#851E38"],[0.99331103678929766,"#851D38"],[0.99665551839464883,"#841C37"],[1,"#841B37"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.039999999999996,"title":"Effect size","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0.0016666666666666668,0.49999999999999994,0.70641642524925907,0.86480531910515712,0.99833333333333329],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":97.148440237213208,"l":133.69863013698634},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,5.5999999999999996],"tickmode":"array","ticktext":["Diagnosis","Sex","Age","AXL receptor tyrosine kinase","Adiponectin"],"tickvals":[1,2,3,4,5],"categoryorder":"array","categoryarray":["Diagnosis","Sex","Age","AXL receptor tyrosine kinase","Adiponectin"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,6.5999999999999996],"tickmode":"array","ticktext":["Insulin","Cortisol","C-reactive protein","Apolipoprotein B","Alpha-2-macroglobulin","Alpha-1-antitrypsin"],"tickvals":[1,2,3,3.9999999999999996,5,6],"categoryorder":"array","categoryarray":["Insulin","Cortisol","C-reactive protein","Apolipoprotein B","Alpha-2-macroglobulin","Alpha-1-antitrypsin"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"Effect size<br />Significance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"24d3709b212":{"x":{},"y":{},"colour":{},"size":{},"shape":{},"type":"scatter"}},"cur_data":"24d3709b212","visdat":{"24d3709b212":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
