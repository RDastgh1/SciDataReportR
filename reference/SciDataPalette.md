# SciDataReportR qualitative color palette

Return colors from the SciDataReportR qualitative color system. The
palette is designed for categorical scientific graphics and begins with
the package anchor colors, Navy and Orange.

## Usage

``` r
SciDataPalette(n = NULL, names = TRUE)
```

## Arguments

- n:

  Number of colors to return. If `NULL`, return the full palette.

- names:

  Logical indicating whether to preserve the color names.

## Value

A character vector of hexadecimal color values.

## Examples

``` r
SciDataPalette(3)
#>      Navy    Orange   Emerald 
#> "#0B1F5E" "#E87400" "#16835F" 
SciDataPalette(8)
#>        Navy      Orange     Emerald   Cranberry        Plum    DeepTeal 
#>   "#0B1F5E"   "#E87400"   "#16835F"   "#C2185B"   "#641B68"   "#0C5F68" 
#>    Burgundy ForestGreen 
#>   "#841B37"   "#285F3D" 
SciDataPalette()
#>        Navy      Orange     Emerald   Cranberry        Plum    DeepTeal 
#>   "#0B1F5E"   "#E87400"   "#16835F"   "#C2185B"   "#641B68"   "#0C5F68" 
#>    Burgundy ForestGreen    Mulberry   Aubergine    Espresso        Blue 
#>   "#841B37"   "#285F3D"   "#9B315F"   "#49305C"   "#593B2B"   "#1769D2" 
#>        Gold      Violet       Coral   JewelBlue   LeafGreen     Magenta 
#>   "#F2B134"   "#7758A3"   "#F2645A"   "#4D73B5"   "#72B65D"   "#D62976" 
#>        Lime BurntOrange       Ochre       Mauve        Rust   SteelBlue 
#>   "#9EBB1F"   "#E65D0A"   "#B77B22"   "#B07A9E"   "#B44D2E"   "#6E8FA8" 
#>       Olive      Sienna        Mint        Rose  PowderBlue  Periwinkle 
#>   "#687548"   "#A85B32"   "#8FD8B8"   "#D9A0AE"   "#B8D8E8"   "#A6AFE8" 
#>    Lavender        Sage       Peach    Cashmere 
#>   "#C7B4D9"   "#B5C7A3"   "#F2B29A"   "#E7C6B2" 
```
