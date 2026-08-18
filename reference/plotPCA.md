# Plot PCA scores

Create an interactive 2D or 3D Plotly scatter plot from a PCA object
created by
[`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).

## Usage

``` r
plotPCA(
  PCAObj,
  Var = NULL,
  t = "NULL",
  HoverVar = NULL,
  HoverVars = NULL,
  Components = NULL,
  Mode = c("auto", "3D", "2D"),
  ColorType = c("auto", "factor", "continuous"),
  Relabel = TRUE,
  Title = NULL
)
```

## Arguments

- PCAObj:

  A PCA object returned by
  [`CreatePCAObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreatePCAObject.md).
  Must contain `Scores` and `CombinedData`.

- Var:

  Optional character string naming a variable in `PCAObj$CombinedData`
  used to color points. Default is `NULL`.

- t:

  Deprecated compatibility argument for color type. Use `ColorType`
  instead. If set to `"Factor"`, `Var` is treated as categorical. Any
  other value treats `Var` as continuous. Default is `"NULL"`.

- HoverVar:

  Optional character string naming one variable in `PCAObj$CombinedData`
  to display in hover text. Retained for backward compatibility. Prefer
  `HoverVars` for new code. Default is `NULL`.

- HoverVars:

  Optional character vector naming one or more variables in
  `PCAObj$CombinedData` to display in hover text. If `NULL`, row number
  is shown. Default is `NULL`.

- Components:

  Optional character or numeric vector specifying which score columns to
  plot. Supply two components for 2D or three components for 3D. If
  `NULL`, the first three score columns are used.

- Mode:

  Character. Either `"auto"`, `"3D"`, or `"2D"`. If `"auto"`, the plot
  dimension is inferred from `Components`. If `Components = NULL`,
  `"auto"` defaults to 3D. Default is `"auto"`.

- ColorType:

  Character. Either `"auto"`, `"factor"`, or `"continuous"`. `"auto"`
  treats character, factor, logical, and labelled variables as
  categorical and numeric variables as continuous. Default is `"auto"`.

- Relabel:

  Logical. If `TRUE`, labels attached to hover variables are used in
  hover text when available. If `FALSE`, raw variable names are used.
  Default is `TRUE`.

- Title:

  Optional plot title. Default is `NULL`, which produces no title.

## Value

A Plotly htmlwidget.

## Details

By default, the function plots the first three score columns in
`PCAObj$Scores` using a 3D plot. These are not assumed to be named
`RC1`, `RC2`, and `RC3`; the function uses the current column order in
`PCAObj$Scores`. Users may manually specify which components to plot
using `Components`.

If `Mode = "auto"`, the function chooses 2D when two components are
supplied and 3D when three components are supplied. If
`Components = NULL`, it defaults to the first three score columns and
creates a 3D plot.

The function can optionally color points by a variable in
`PCAObj$CombinedData` and customize hover text using one or more
variables. When `Relabel = TRUE`, labels attached to hover variables are
used in the hover display when available.

## Examples

``` r
# \donttest{
# Build a PCA object to plot
PCAObj <- CreatePCAObject(
  data = mtcars,
  VarsToReduce = names(mtcars),
  numComponents = 3
)

# Default 3D scatter of the first three components
plotPCA(PCAObj)

{"x":{"visdat":{"224a73c5796d":["function () ","plotlyVisDat"]},"cur_data":"224a73c5796d","attrs":{"224a73c5796d":{"x":{},"y":{},"z":{},"text":{},"mode":"markers","marker":{"size":5},"hovertemplate":"%{text}<extra><\/extra>","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter3d"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"scene":{"xaxis":{"title":"RC2"},"yaxis":{"title":"RC3"},"zaxis":{"title":"RC1"}},"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0.92162135764966635,0.83654752015361356,0.6715640272032527,-0.9473444044640259,-0.86858778262704162,-1.1952350431488712,-0.6571925880363213,-0.078457169754194228,-0.23063658339644941,0.089878756509196536,0.017487502182628803,-0.89993200241266835,-0.84817398066541705,-0.92391601353367159,-1.3877655587673279,-1.3600441480792314,-1.13054556988882,0.88984275733773988,1.4272302152406233,0.99728143316896634,-0.41432829253919035,-1.0380278968493832,-0.90287464168293863,-0.4736967464125994,-0.99144757688655494,0.89079046419059815,1.5836893464014594,1.4934371121592538,1.1563875236275094,1.4066511894245342,1.2152064650335972,0.75059032886206767],"y":[0.89334758098706435,0.67268136905612741,-0.47441024082279309,-0.76370007110780858,1.0265815644063512,-1.0754144809041573,1.0396628228503126,-1.6504682290817025,-2.5159947654107064,-1.2889656924879762,-1.4829644329586469,0.54086116744142498,0.57915775252908341,0.41910233388396401,0.085250185935608971,0.07402516184738879,0.2294667755115328,-0.62370884269568305,-0.53989407851424076,-0.66763362935902026,-1.2481696161632096,1.0105243244327455,0.82207458432889036,0.98192243659931133,0.97795549663346504,-0.47982461139618698,0.73543384656869826,0.067606437227274371,1.5830022932515149,0.78470544051490332,1.0454249305484409,-0.7576378136519708],"z":[-0.27375080772972815,-0.098413532211183652,-0.74438413415473781,-0.7764831095713649,-0.97921932637951947,-0.54199690484488172,0.23052822406057222,0.22321175479300009,0.90038252567557153,1.3500228702341586,1.501626878898777,-0.16329855724796757,-0.30505876583003566,-0.12999523881353564,1.0927003010180036,1.198103961293189,1.044906508815288,-1.1434511729191681,-0.75080097168764226,-1.2458768671151788,-0.69452493093472234,-1.042595068043672,-0.86548987145415335,0.49545382458171539,-0.87403984302667237,-1.0825571064160551,-0.75959103139237072,-0.67697788166596351,0.96922270217709861,1.2011565455937656,2.8831771712741476,0.058011853023265669],"text":["Row: 1<br>RC2: 0.922<br>RC3: 0.893<br>RC1: -0.274","Row: 2<br>RC2: 0.837<br>RC3: 0.673<br>RC1: -0.098","Row: 3<br>RC2: 0.672<br>RC3: -0.474<br>RC1: -0.744","Row: 4<br>RC2: -0.947<br>RC3: -0.764<br>RC1: -0.776","Row: 5<br>RC2: -0.869<br>RC3: 1.027<br>RC1: -0.979","Row: 6<br>RC2: -1.195<br>RC3: -1.075<br>RC1: -0.542","Row: 7<br>RC2: -0.657<br>RC3: 1.04<br>RC1: 0.231","Row: 8<br>RC2: -0.078<br>RC3: -1.65<br>RC1: 0.223","Row: 9<br>RC2: -0.231<br>RC3: -2.516<br>RC1: 0.9","Row: 10<br>RC2: 0.09<br>RC3: -1.289<br>RC1: 1.35","Row: 11<br>RC2: 0.017<br>RC3: -1.483<br>RC1: 1.502","Row: 12<br>RC2: -0.9<br>RC3: 0.541<br>RC1: -0.163","Row: 13<br>RC2: -0.848<br>RC3: 0.579<br>RC1: -0.305","Row: 14<br>RC2: -0.924<br>RC3: 0.419<br>RC1: -0.13","Row: 15<br>RC2: -1.388<br>RC3: 0.085<br>RC1: 1.093","Row: 16<br>RC2: -1.36<br>RC3: 0.074<br>RC1: 1.198","Row: 17<br>RC2: -1.131<br>RC3: 0.229<br>RC1: 1.045","Row: 18<br>RC2: 0.89<br>RC3: -0.624<br>RC1: -1.143","Row: 19<br>RC2: 1.427<br>RC3: -0.54<br>RC1: -0.751","Row: 20<br>RC2: 0.997<br>RC3: -0.668<br>RC1: -1.246","Row: 21<br>RC2: -0.414<br>RC3: -1.248<br>RC1: -0.695","Row: 22<br>RC2: -1.038<br>RC3: 1.011<br>RC1: -1.043","Row: 23<br>RC2: -0.903<br>RC3: 0.822<br>RC1: -0.865","Row: 24<br>RC2: -0.474<br>RC3: 0.982<br>RC1: 0.495","Row: 25<br>RC2: -0.991<br>RC3: 0.978<br>RC1: -0.874","Row: 26<br>RC2: 0.891<br>RC3: -0.48<br>RC1: -1.083","Row: 27<br>RC2: 1.584<br>RC3: 0.735<br>RC1: -0.76","Row: 28<br>RC2: 1.493<br>RC3: 0.068<br>RC1: -0.677","Row: 29<br>RC2: 1.156<br>RC3: 1.583<br>RC1: 0.969","Row: 30<br>RC2: 1.407<br>RC3: 0.785<br>RC1: 1.201","Row: 31<br>RC2: 1.215<br>RC3: 1.045<br>RC1: 2.883","Row: 32<br>RC2: 0.751<br>RC3: -0.758<br>RC1: 0.058"],"mode":"markers","marker":{"color":"rgba(31,119,180,1)","size":5,"line":{"color":"rgba(31,119,180,1)"}},"hovertemplate":["%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>"],"type":"scatter3d","error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
# 2D scatter colored by a variable in the data
plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  Mode = "2D",
  Var = "cyl",
  ColorType = "factor"
)

{"x":{"visdat":{"224a288c3ca4":["function () ","plotlyVisDat"]},"cur_data":"224a288c3ca4","attrs":{"224a288c3ca4":{"x":{},"y":{},"text":{},"mode":"markers","hovertemplate":"%{text}<extra><\/extra>","color":["6","6","4","6","8","6","8","4","4","6","6","8","8","8","8","8","8","4","4","4","4","8","8","8","8","4","4","4","8","6","8","4"],"colors":["#0B1F5E","#E87400","#16835F"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true,"title":"RC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"RC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.74438413415473781,0.22321175479300009,0.90038252567557153,-1.1434511729191681,-0.75080097168764226,-1.2458768671151788,-0.69452493093472234,-1.0825571064160551,-0.75959103139237072,-0.67697788166596351,0.058011853023265669],"y":[0.6715640272032527,-0.078457169754194228,-0.23063658339644941,0.88984275733773988,1.4272302152406233,0.99728143316896634,-0.41432829253919035,0.89079046419059815,1.5836893464014594,1.4934371121592538,0.75059032886206767],"text":["Row: 3<br>RC1: -0.744<br>RC2: 0.672","Row: 8<br>RC1: 0.223<br>RC2: -0.078","Row: 9<br>RC1: 0.9<br>RC2: -0.231","Row: 18<br>RC1: -1.143<br>RC2: 0.89","Row: 19<br>RC1: -0.751<br>RC2: 1.427","Row: 20<br>RC1: -1.246<br>RC2: 0.997","Row: 21<br>RC1: -0.695<br>RC2: -0.414","Row: 26<br>RC1: -1.083<br>RC2: 0.891","Row: 27<br>RC1: -0.76<br>RC2: 1.584","Row: 28<br>RC1: -0.677<br>RC2: 1.493","Row: 32<br>RC1: 0.058<br>RC2: 0.751"],"mode":"markers","hovertemplate":["%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>"],"type":"scatter","name":"4","marker":{"color":"rgba(11,31,94,1)","line":{"color":"rgba(11,31,94,1)"}},"textfont":{"color":"rgba(11,31,94,1)"},"error_y":{"color":"rgba(11,31,94,1)"},"error_x":{"color":"rgba(11,31,94,1)"},"line":{"color":"rgba(11,31,94,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.27375080772972815,-0.098413532211183652,-0.7764831095713649,-0.54199690484488172,1.3500228702341586,1.501626878898777,1.2011565455937656],"y":[0.92162135764966635,0.83654752015361356,-0.9473444044640259,-1.1952350431488712,0.089878756509196536,0.017487502182628803,1.4066511894245342],"text":["Row: 1<br>RC1: -0.274<br>RC2: 0.922","Row: 2<br>RC1: -0.098<br>RC2: 0.837","Row: 4<br>RC1: -0.776<br>RC2: -0.947","Row: 6<br>RC1: -0.542<br>RC2: -1.195","Row: 10<br>RC1: 1.35<br>RC2: 0.09","Row: 11<br>RC1: 1.502<br>RC2: 0.017","Row: 30<br>RC1: 1.201<br>RC2: 1.407"],"mode":"markers","hovertemplate":["%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>"],"type":"scatter","name":"6","marker":{"color":"rgba(232,116,0,1)","line":{"color":"rgba(232,116,0,1)"}},"textfont":{"color":"rgba(232,116,0,1)"},"error_y":{"color":"rgba(232,116,0,1)"},"error_x":{"color":"rgba(232,116,0,1)"},"line":{"color":"rgba(232,116,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.97921932637951947,0.23052822406057222,-0.16329855724796757,-0.30505876583003566,-0.12999523881353564,1.0927003010180036,1.198103961293189,1.044906508815288,-1.042595068043672,-0.86548987145415335,0.49545382458171539,-0.87403984302667237,0.96922270217709861,2.8831771712741476],"y":[-0.86858778262704162,-0.6571925880363213,-0.89993200241266835,-0.84817398066541705,-0.92391601353367159,-1.3877655587673279,-1.3600441480792314,-1.13054556988882,-1.0380278968493832,-0.90287464168293863,-0.4736967464125994,-0.99144757688655494,1.1563875236275094,1.2152064650335972],"text":["Row: 5<br>RC1: -0.979<br>RC2: -0.869","Row: 7<br>RC1: 0.231<br>RC2: -0.657","Row: 12<br>RC1: -0.163<br>RC2: -0.9","Row: 13<br>RC1: -0.305<br>RC2: -0.848","Row: 14<br>RC1: -0.13<br>RC2: -0.924","Row: 15<br>RC1: 1.093<br>RC2: -1.388","Row: 16<br>RC1: 1.198<br>RC2: -1.36","Row: 17<br>RC1: 1.045<br>RC2: -1.131","Row: 22<br>RC1: -1.043<br>RC2: -1.038","Row: 23<br>RC1: -0.865<br>RC2: -0.903","Row: 24<br>RC1: 0.495<br>RC2: -0.474","Row: 25<br>RC1: -0.874<br>RC2: -0.991","Row: 29<br>RC1: 0.969<br>RC2: 1.156","Row: 31<br>RC1: 2.883<br>RC2: 1.215"],"mode":"markers","hovertemplate":["%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>"],"type":"scatter","name":"8","marker":{"color":"rgba(22,131,95,1)","line":{"color":"rgba(22,131,95,1)"}},"textfont":{"color":"rgba(22,131,95,1)"},"error_y":{"color":"rgba(22,131,95,1)"},"error_x":{"color":"rgba(22,131,95,1)"},"line":{"color":"rgba(22,131,95,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
# Add hover information
plotPCA(
  PCAObj,
  Components = c("RC1", "RC2"),
  Mode = "2D",
  HoverVars = c("mpg", "hp")
)

{"x":{"visdat":{"224ad3ad76e":["function () ","plotlyVisDat"]},"cur_data":"224ad3ad76e","attrs":{"224ad3ad76e":{"x":{},"y":{},"text":{},"mode":"markers","hovertemplate":"%{text}<extra><\/extra>","alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true,"title":"RC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"RC2"},"hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.27375080772972815,-0.098413532211183652,-0.74438413415473781,-0.7764831095713649,-0.97921932637951947,-0.54199690484488172,0.23052822406057222,0.22321175479300009,0.90038252567557153,1.3500228702341586,1.501626878898777,-0.16329855724796757,-0.30505876583003566,-0.12999523881353564,1.0927003010180036,1.198103961293189,1.044906508815288,-1.1434511729191681,-0.75080097168764226,-1.2458768671151788,-0.69452493093472234,-1.042595068043672,-0.86548987145415335,0.49545382458171539,-0.87403984302667237,-1.0825571064160551,-0.75959103139237072,-0.67697788166596351,0.96922270217709861,1.2011565455937656,2.8831771712741476,0.058011853023265669],"y":[0.92162135764966635,0.83654752015361356,0.6715640272032527,-0.9473444044640259,-0.86858778262704162,-1.1952350431488712,-0.6571925880363213,-0.078457169754194228,-0.23063658339644941,0.089878756509196536,0.017487502182628803,-0.89993200241266835,-0.84817398066541705,-0.92391601353367159,-1.3877655587673279,-1.3600441480792314,-1.13054556988882,0.88984275733773988,1.4272302152406233,0.99728143316896634,-0.41432829253919035,-1.0380278968493832,-0.90287464168293863,-0.4736967464125994,-0.99144757688655494,0.89079046419059815,1.5836893464014594,1.4934371121592538,1.1563875236275094,1.4066511894245342,1.2152064650335972,0.75059032886206767],"text":["mpg: 21<br>hp: 110<br>RC1: -0.274<br>RC2: 0.922","mpg: 21<br>hp: 110<br>RC1: -0.098<br>RC2: 0.837","mpg: 22.8<br>hp: 93<br>RC1: -0.744<br>RC2: 0.672","mpg: 21.4<br>hp: 110<br>RC1: -0.776<br>RC2: -0.947","mpg: 18.7<br>hp: 175<br>RC1: -0.979<br>RC2: -0.869","mpg: 18.1<br>hp: 105<br>RC1: -0.542<br>RC2: -1.195","mpg: 14.3<br>hp: 245<br>RC1: 0.231<br>RC2: -0.657","mpg: 24.4<br>hp: 62<br>RC1: 0.223<br>RC2: -0.078","mpg: 22.8<br>hp: 95<br>RC1: 0.9<br>RC2: -0.231","mpg: 19.2<br>hp: 123<br>RC1: 1.35<br>RC2: 0.09","mpg: 17.8<br>hp: 123<br>RC1: 1.502<br>RC2: 0.017","mpg: 16.4<br>hp: 180<br>RC1: -0.163<br>RC2: -0.9","mpg: 17.3<br>hp: 180<br>RC1: -0.305<br>RC2: -0.848","mpg: 15.2<br>hp: 180<br>RC1: -0.13<br>RC2: -0.924","mpg: 10.4<br>hp: 205<br>RC1: 1.093<br>RC2: -1.388","mpg: 10.4<br>hp: 215<br>RC1: 1.198<br>RC2: -1.36","mpg: 14.7<br>hp: 230<br>RC1: 1.045<br>RC2: -1.131","mpg: 32.4<br>hp: 66<br>RC1: -1.143<br>RC2: 0.89","mpg: 30.4<br>hp: 52<br>RC1: -0.751<br>RC2: 1.427","mpg: 33.9<br>hp: 65<br>RC1: -1.246<br>RC2: 0.997","mpg: 21.5<br>hp: 97<br>RC1: -0.695<br>RC2: -0.414","mpg: 15.5<br>hp: 150<br>RC1: -1.043<br>RC2: -1.038","mpg: 15.2<br>hp: 150<br>RC1: -0.865<br>RC2: -0.903","mpg: 13.3<br>hp: 245<br>RC1: 0.495<br>RC2: -0.474","mpg: 19.2<br>hp: 175<br>RC1: -0.874<br>RC2: -0.991","mpg: 27.3<br>hp: 66<br>RC1: -1.083<br>RC2: 0.891","mpg: 26<br>hp: 91<br>RC1: -0.76<br>RC2: 1.584","mpg: 30.4<br>hp: 113<br>RC1: -0.677<br>RC2: 1.493","mpg: 15.8<br>hp: 264<br>RC1: 0.969<br>RC2: 1.156","mpg: 19.7<br>hp: 175<br>RC1: 1.201<br>RC2: 1.407","mpg: 15<br>hp: 335<br>RC1: 2.883<br>RC2: 1.215","mpg: 21.4<br>hp: 109<br>RC1: 0.058<br>RC2: 0.751"],"mode":"markers","hovertemplate":["%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>","%{text}<extra><\/extra>"],"type":"scatter","marker":{"color":"rgba(31,119,180,1)","line":{"color":"rgba(31,119,180,1)"}},"error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}# }
```
