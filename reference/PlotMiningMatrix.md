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
  Covariates = lifecycle::deprecated()
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

  Optional covariates (reserved for future use).

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

  Either `"matrix"` (default) or `"per_outcome"`, passed to
  [`ApplyFDRCorrection()`](https://rdastgh1.github.io/SciDataReportR/reference/ApplyFDRCorrection.md).
  `"matrix"` corrects across all pairwise p-values at once (historical
  behavior, computed on the symmetrized pair table). `"per_outcome"`
  corrects separately within each x-axis variable (`XVar`, ordered by
  `outcome_vars`).

- Covariates:

  **Deprecated** (since 19.15.0). Use `covariates` instead.

## Value

List with tables and plots.

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
#> Warning: Removed 5 rows containing missing values or values outside the scale range
#> (`geom_point()`).


# The p-value table carries both unadjusted (p) and FDR-adjusted (p_adj)
# p-values so results can be inspected with and without FDR correction
result$Unadjusted$PvalTable[, c("XVar", "YVar", "p", "p_adj", "Test")]
#> # A tibble: 18 × 5
#>    XVar        YVar                         p    p_adj Test       
#>    <chr>       <chr>                    <dbl>    <dbl> <chr>      
#>  1 AXL         Alpha_1_Antitrypsin   1.12e- 1 2.01e- 1 Correlation
#>  2 AXL         Alpha_2_Macroglobulin 8.35e-15 5.01e-14 Correlation
#>  3 AXL         Apolipoprotein_B      6.02e- 1 6.78e- 1 Correlation
#>  4 AXL         C_Reactive_Protein    4.82e- 1 6.58e- 1 Correlation
#>  5 AXL         Cortisol              1.85e- 2 4.15e- 2 Correlation
#>  6 AXL         Insulin               2.29e-14 1.03e-13 Correlation
#>  7 Adiponectin Alpha_1_Antitrypsin   1.21e-15 1.09e-14 Correlation
#>  8 Adiponectin Alpha_2_Macroglobulin 3.24e-12 1.17e-11 Correlation
#>  9 Adiponectin Apolipoprotein_B      2.00e-18 3.60e-17 Correlation
#> 10 Adiponectin C_Reactive_Protein    1.02e- 1 2.01e- 1 Correlation
#> 11 Adiponectin Cortisol              9.14e- 3 2.35e- 2 Correlation
#> 12 Adiponectin Insulin               7.08e- 4 2.12e- 3 Correlation
#> 13 age         Alpha_1_Antitrypsin   7.17e- 1 7.60e- 1 Correlation
#> 14 age         Alpha_2_Macroglobulin 3.58e- 1 5.36e- 1 Correlation
#> 15 age         Apolipoprotein_B      5.25e- 1 6.58e- 1 Correlation
#> 16 age         C_Reactive_Protein    3.49e- 1 5.36e- 1 Correlation
#> 17 age         Cortisol              9.60e- 1 9.60e- 1 Correlation
#> 18 age         Insulin               5.48e- 1 6.58e- 1 Correlation

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

{"x":{"data":[{"x":[2],"y":[2],"text":"XLabel: AXL receptor tyrosine kinase<br />YLabel: Cortisol<br />EffectSizeAbs: 0.129075315<br />size_val:  3<br />stars: *","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(253,194,41,1)","opacity":1,"size":14.105388879000483,"symbol":"triangle-up","line":{"width":1.8897637795275593,"color":"rgba(253,194,41,1)"}},"hoveron":"points","name":"*","legendgroup":"*","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[3],"y":[2],"text":"XLabel: Adiponectin<br />YLabel: Cortisol<br />EffectSizeAbs: 0.142664896<br />size_val:  4<br />stars: **","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(254,188,42,1)","opacity":1,"size":16.816969106582089,"symbol":"square","line":{"width":1.8897637795275593,"color":"rgba(254,188,42,1)"}},"hoveron":"points","name":"**","legendgroup":"**","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[3],"y":[1],"text":"XLabel: Adiponectin<br />YLabel: Insulin<br />EffectSizeAbs: 0.184673186<br />size_val:  5<br />stars: ***","type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(253,172,51,1)","opacity":1,"size":18.897637795275593,"symbol":"diamond","line":{"width":1.8897637795275593,"color":"rgba(253,172,51,1)"}},"hoveron":"points","name":"***","legendgroup":"***","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[2,2,3,3,3],"y":[5,1,6,5,4],"text":["XLabel: AXL receptor tyrosine kinase<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.408223294<br />size_val: NA<br />stars: ****","XLabel: AXL receptor tyrosine kinase<br />YLabel: Insulin<br />EffectSizeAbs: 0.402040919<br />size_val: NA<br />stars: ****","XLabel: Adiponectin<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.419736372<br />size_val: NA<br />stars: ****","XLabel: Adiponectin<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.369564701<br />size_val: NA<br />stars: ****","XLabel: Adiponectin<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.455060864<br />size_val: NA<br />stars: ****"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(223,98,99,1)","rgba(224,99,99,1)","rgba(221,94,102,1)","rgba(231,110,91,1)","rgba(213,84,110,1)"],"opacity":1,"size":[null,null,null,null,null],"symbol":[null,null,null,null,null],"line":{"width":1.8897637795275593,"color":["rgba(223,98,99,1)","rgba(224,99,99,1)","rgba(221,94,102,1)","rgba(231,110,91,1)","rgba(213,84,110,1)"]}},"hoveron":"points","name":"****","legendgroup":"****","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[2,2,2,3,1,1,1,1,1,1],"y":[6,4,3,3,6,5,4,3,2,1],"text":["XLabel: AXL receptor tyrosine kinase<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.087267572<br />size_val:  2<br />stars: ns","XLabel: AXL receptor tyrosine kinase<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.028642053<br />size_val:  2<br />stars: ns","XLabel: AXL receptor tyrosine kinase<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.038649351<br />size_val:  2<br />stars: ns","XLabel: Adiponectin<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.089756883<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: Alpha-1-antitrypsin<br />EffectSizeAbs: 0.020249956<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: Alpha-2-macroglobulin<br />EffectSizeAbs: 0.051437580<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: Apolipoprotein B<br />EffectSizeAbs: 0.035570080<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: C-reactive protein<br />EffectSizeAbs: 0.052352600<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: Cortisol<br />EffectSizeAbs: 0.002807788<br />size_val:  2<br />stars: ns","XLabel: Age<br />YLabel: Insulin<br />EffectSizeAbs: 0.033582272<br />size_val:  2<br />stars: ns"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(251,211,36,1)","rgba(244,236,39,1)","rgba(246,232,38,1)","rgba(252,210,37,1)","rgba(243,240,39,1)","rgba(247,226,37,1)","rgba(245,233,38,1)","rgba(247,226,37,1)","rgba(240,248,35,1)","rgba(245,234,38,1)"],"opacity":1,"size":7.559055118110237,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(251,211,36,1)","rgba(244,236,39,1)","rgba(246,232,38,1)","rgba(252,210,37,1)","rgba(243,240,39,1)","rgba(247,226,37,1)","rgba(245,233,38,1)","rgba(247,226,37,1)","rgba(240,248,35,1)","rgba(245,234,38,1)"]}},"hoveron":"points","name":"ns","legendgroup":"ns","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1],"y":[1],"name":"1594b806e12ab0879231f558ad12d7ae","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#F0F921"],[0.0033444816053511705,"#F0F724"],[0.006688963210702341,"#F1F625"],[0.010033444816053512,"#F1F426"],[0.013377926421404682,"#F1F326"],[0.016722408026755852,"#F2F127"],[0.020066889632107024,"#F3F027"],[0.023411371237458192,"#F3EE27"],[0.026755852842809364,"#F4ED27"],[0.030100334448160536,"#F5EC27"],[0.033444816053511704,"#F5EA26"],[0.036789297658862873,"#F5E926"],[0.040133779264214048,"#F6E826"],[0.043478260869565216,"#F6E626"],[0.046822742474916385,"#F7E425"],[0.05016722408026756,"#F7E225"],[0.053511705685618728,"#F8E125"],[0.056856187290969896,"#F8E025"],[0.060200668896321072,"#F8DE25"],[0.06354515050167224,"#F9DD25"],[0.066889632107023408,"#F9DC24"],[0.070234113712374577,"#FADA24"],[0.073578595317725745,"#FAD824"],[0.076923076923076927,"#FBD724"],[0.080267558528428096,"#FBD624"],[0.083612040133779264,"#FBD424"],[0.086956521739130432,"#FBD324"],[0.090301003344481601,"#FCD225"],[0.093645484949832769,"#FCD025"],[0.096989966555183951,"#FCCF25"],[0.10033444816053512,"#FCCD25"],[0.10367892976588629,"#FCCC25"],[0.10702341137123746,"#FDCB26"],[0.11036789297658862,"#FDCA26"],[0.11371237458193979,"#FDC827"],[0.11705685618729096,"#FDC627"],[0.12040133779264214,"#FDC527"],[0.12374581939799331,"#FDC428"],[0.12709030100334448,"#FDC328"],[0.13043478260869565,"#FDC129"],[0.13377926421404682,"#FEC029"],[0.13712374581939799,"#FEBE2A"],[0.14046822742474915,"#FEBD2A"],[0.14381270903010032,"#FEBC2B"],[0.14715719063545149,"#FEBA2C"],[0.15050167224080269,"#FEB92C"],[0.15384615384615385,"#FEB82C"],[0.15719063545150502,"#FEB72D"],[0.16053511705685619,"#FDB52E"],[0.16387959866220736,"#FDB42F"],[0.16722408026755853,"#FDB32F"],[0.1705685618729097,"#FDB22F"],[0.17391304347826086,"#FDB030"],[0.17725752508361203,"#FDAF31"],[0.1806020066889632,"#FDAE32"],[0.18394648829431437,"#FDAC33"],[0.18729096989966554,"#FDAB33"],[0.19063545150501671,"#FCAA34"],[0.1939799331103679,"#FCA934"],[0.19732441471571907,"#FCA735"],[0.20066889632107024,"#FCA636"],[0.20401337792642141,"#FCA537"],[0.20735785953177258,"#FCA338"],[0.21070234113712374,"#FBA238"],[0.21404682274247491,"#FBA139"],[0.21739130434782608,"#FBA039"],[0.22073578595317725,"#FB9F3A"],[0.22408026755852842,"#FA9E3B"],[0.22742474916387959,"#FA9C3C"],[0.23076923076923075,"#FA9B3D"],[0.23411371237458192,"#F99A3E"],[0.23745819397993312,"#F9993E"],[0.24080267558528429,"#F9983E"],[0.24414715719063546,"#F9963F"],[0.24749163879598662,"#F89540"],[0.25083612040133779,"#F89441"],[0.25418060200668896,"#F79342"],[0.25752508361204013,"#F79243"],[0.2608695652173913,"#F79044"],[0.26421404682274247,"#F79044"],[0.26755852842809363,"#F68F44"],[0.2709030100334448,"#F68D45"],[0.27424749163879597,"#F58C46"],[0.27759197324414714,"#F58B47"],[0.28093645484949831,"#F48A48"],[0.28428093645484948,"#F48948"],[0.28762541806020064,"#F48849"],[0.29096989966555181,"#F3874A"],[0.29431438127090298,"#F3854B"],[0.29765886287625415,"#F2844B"],[0.30100334448160537,"#F1834C"],[0.30434782608695654,"#F1824D"],[0.30769230769230771,"#F1814D"],[0.31103678929765888,"#F0804E"],[0.31438127090301005,"#F07F4F"],[0.31772575250836121,"#EF7E50"],[0.32107023411371238,"#EF7C51"],[0.32441471571906355,"#EE7B51"],[0.32775919732441472,"#ED7A52"],[0.33110367892976589,"#ED7A52"],[0.33444816053511706,"#ED7853"],[0.33779264214046822,"#EC7754"],[0.34113712374581939,"#EB7655"],[0.34448160535117056,"#EB7556"],[0.34782608695652173,"#EA7457"],[0.3511705685618729,"#E97357"],[0.35451505016722407,"#E97257"],[0.35785953177257523,"#E97158"],[0.3612040133779264,"#E87059"],[0.36454849498327757,"#E76F5A"],[0.36789297658862874,"#E76E5B"],[0.37123745819397991,"#E66D5C"],[0.37458193979933108,"#E56B5D"],[0.37792642140468224,"#E56B5D"],[0.38127090301003341,"#E56A5D"],[0.38461538461538458,"#E4695E"],[0.38795986622073581,"#E3685F"],[0.39130434782608697,"#E26660"],[0.39464882943143814,"#E26561"],[0.39799331103678931,"#E26561"],[0.40133779264214048,"#E16462"],[0.40468227424749165,"#E06363"],[0.40802675585284282,"#DF6263"],[0.41137123745819398,"#DE6164"],[0.41471571906354515,"#DE5F65"],[0.41806020066889632,"#DD5E66"],[0.42140468227424749,"#DD5E66"],[0.42474916387959866,"#DC5D67"],[0.42809364548494983,"#DB5C68"],[0.43143812709030099,"#DA5B69"],[0.43478260869565216,"#DA5A6A"],[0.43812709030100333,"#D9596A"],[0.4414715719063545,"#D8576B"],[0.44481605351170567,"#D8576B"],[0.44816053511705684,"#D7566C"],[0.451505016722408,"#D6556D"],[0.45484949832775917,"#D5546E"],[0.45819397993311034,"#D5536F"],[0.46153846153846151,"#D45270"],[0.46488294314381268,"#D35171"],[0.46822742474916385,"#D35071"],[0.47157190635451507,"#D24F71"],[0.47491638795986624,"#D14E72"],[0.47826086956521741,"#D04D73"],[0.48160535117056857,"#CF4C74"],[0.48494983277591974,"#CE4B75"],[0.48829431438127091,"#CD4A76"],[0.49163879598662208,"#CD4A76"],[0.49498327759197325,"#CC4977"],[0.49832775919732442,"#CC4778"],[0.50167224080267558,"#CB4679"],[0.50501672240802675,"#CA457A"],[0.50836120401337792,"#C9447A"],[0.51170568561872909,"#C9447A"],[0.51505016722408026,"#C8437B"],[0.51839464882943143,"#C7427C"],[0.52173913043478259,"#C6417D"],[0.52508361204013376,"#C5407E"],[0.52842809364548493,"#C43F7F"],[0.5317725752508361,"#C33D80"],[0.53511705685618727,"#C33D80"],[0.53846153846153844,"#C23C81"],[0.5418060200668896,"#C13B82"],[0.54515050167224077,"#C03A83"],[0.54849498327759194,"#BF3984"],[0.55183946488294311,"#BE3885"],[0.55518394648829428,"#BD3786"],[0.55852842809364545,"#BD3686"],[0.56187290969899661,"#BC3587"],[0.56521739130434778,"#BB3488"],[0.56856187290969895,"#BA3388"],[0.57190635451505012,"#B83289"],[0.57525083612040129,"#B7318A"],[0.57859531772575246,"#B6308B"],[0.58193979933110362,"#B6308B"],[0.58528428093645479,"#B52F8C"],[0.58862876254180596,"#B42E8D"],[0.59197324414715713,"#B32C8E"],[0.5953177257525083,"#B22B8F"],[0.59866220735785947,"#B12A90"],[0.60200668896321075,"#B02991"],[0.60535117056856191,"#AF2991"],[0.60869565217391308,"#AE2892"],[0.61204013377926425,"#AD2793"],[0.61538461538461542,"#AC2694"],[0.61872909698996659,"#AB2494"],[0.62207357859531776,"#AA2395"],[0.62541806020066892,"#A92395"],[0.62876254180602009,"#A82296"],[0.63210702341137126,"#A72197"],[0.63545150501672243,"#A62098"],[0.6387959866220736,"#A51F99"],[0.64214046822742477,"#A41E9A"],[0.64548494983277593,"#A21D9A"],[0.6488294314381271,"#A21C9A"],[0.65217391304347827,"#A11B9B"],[0.65551839464882944,"#A01A9C"],[0.65886287625418061,"#9E199D"],[0.66220735785953178,"#9D189D"],[0.66555183946488294,"#9C179E"],[0.66889632107023411,"#9B169F"],[0.67224080267558528,"#9A169F"],[0.67558528428093645,"#99159F"],[0.67892976588628762,"#9814A0"],[0.68227424749163879,"#9613A1"],[0.68561872909698995,"#9511A1"],[0.68896321070234112,"#9410A2"],[0.69230769230769229,"#930FA3"],[0.69565217391304346,"#920FA3"],[0.69899665551839463,"#910EA3"],[0.7023411371237458,"#8F0DA4"],[0.70568561872909696,"#8E0CA4"],[0.70903010033444813,"#8D0BA5"],[0.7123745819397993,"#8C0AA5"],[0.71571906354515047,"#8A09A5"],[0.71906354515050164,"#8909A5"],[0.72240802675585281,"#8808A6"],[0.72575250836120397,"#8707A6"],[0.72909698996655514,"#8606A6"],[0.73244147157190631,"#8405A7"],[0.73578595317725748,"#8305A7"],[0.73913043478260865,"#8205A7"],[0.74247491638795982,"#8104A7"],[0.74581939799331098,"#8004A8"],[0.74916387959866215,"#7E03A8"],[0.75250836120401332,"#7D03A8"],[0.75585284280936449,"#7C02A8"],[0.75919732441471566,"#7A02A8"],[0.76254180602006683,"#7902A8"],[0.76588628762541799,"#7801A8"],[0.76923076923076916,"#7701A8"],[0.77257525083612044,"#7501A8"],[0.77591973244147161,"#7401A8"],[0.77926421404682278,"#7301A8"],[0.78260869565217395,"#7100A8"],[0.78595317725752512,"#7000A8"],[0.78929765886287628,"#6F00A8"],[0.79264214046822745,"#6E00A8"],[0.79598662207357862,"#6C00A8"],[0.79933110367892979,"#6A00A8"],[0.80267558528428096,"#6900A8"],[0.80602006688963213,"#6800A8"],[0.80936454849498329,"#6700A8"],[0.81270903010033446,"#6600A7"],[0.81605351170568563,"#6400A7"],[0.8193979933110368,"#6300A7"],[0.82274247491638797,"#6100A7"],[0.82608695652173914,"#6001A6"],[0.8294314381270903,"#5F01A6"],[0.83277591973244147,"#5D01A6"],[0.83612040133779264,"#5C01A6"],[0.83946488294314381,"#5B01A5"],[0.84280936454849498,"#5901A5"],[0.84615384615384615,"#5801A4"],[0.84949832775919731,"#5701A4"],[0.85284280936454848,"#5601A4"],[0.85618729096989965,"#5402A4"],[0.85953177257525082,"#5302A3"],[0.86287625418060199,"#5102A3"],[0.86622073578595316,"#5002A2"],[0.86956521739130432,"#4F02A2"],[0.87290969899665549,"#4D02A1"],[0.87625418060200666,"#4C02A1"],[0.87959866220735783,"#4A03A1"],[0.882943143812709,"#4903A0"],[0.88628762541806017,"#48039F"],[0.88963210702341133,"#46039F"],[0.8929765886287625,"#45039E"],[0.89632107023411367,"#43039E"],[0.89966555183946484,"#42039E"],[0.90301003344481601,"#40049D"],[0.90635451505016718,"#3F049C"],[0.90969899665551834,"#3E049C"],[0.91304347826086951,"#3C049B"],[0.91638795986622068,"#3B049A"],[0.91973244147157185,"#39049A"],[0.92307692307692302,"#38049A"],[0.92642140468227419,"#370499"],[0.92976588628762535,"#350498"],[0.93311036789297652,"#330597"],[0.93645484949832769,"#310597"],[0.93979933110367886,"#300596"],[0.94314381270903014,"#2E0595"],[0.94648829431438131,"#2D0595"],[0.94983277591973247,"#2C0594"],[0.95317725752508364,"#2A0593"],[0.95652173913043481,"#280592"],[0.95986622073578598,"#260591"],[0.96321070234113715,"#250691"],[0.96655518394648832,"#230691"],[0.96989966555183948,"#210690"],[0.97324414715719065,"#1F068F"],[0.97658862876254182,"#1D068E"],[0.97993311036789299,"#1B068D"],[0.98327759197324416,"#1A068C"],[0.98662207357859533,"#17078B"],[0.98996655518394649,"#15078A"],[0.99331103678929766,"#120789"],[0.99665551839464883,"#100788"],[1,"#0D0887"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.039999999999996,"title":"Effect Size","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0.0016666666666666668,0.2508333333333333,0.49999999999999994,0.74916666666666665,0.99833333333333329],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":97.148440237213208,"l":133.69863013698634},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,3.6000000000000001],"tickmode":"array","ticktext":["Age","AXL receptor tyrosine kinase","Adiponectin"],"tickvals":[1,2,3],"categoryorder":"array","categoryarray":["Age","AXL receptor tyrosine kinase","Adiponectin"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-45,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.40000000000000002,6.5999999999999996],"tickmode":"array","ticktext":["Insulin","Cortisol","C-reactive protein","Apolipoprotein B","Alpha-2-macroglobulin","Alpha-1-antitrypsin"],"tickvals":[1,2,3,3.9999999999999996,5,6],"categoryorder":"array","categoryarray":["Insulin","Cortisol","C-reactive protein","Apolipoprotein B","Alpha-2-macroglobulin","Alpha-1-antitrypsin"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(235,235,235,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"","font":{"color":null,"family":null,"size":0}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"Effect Size<br />Significance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"22ec649c1017":{"x":{},"y":{},"colour":{},"size":{},"shape":{},"type":"scatter"}},"cur_data":"22ec649c1017","visdat":{"22ec649c1017":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
