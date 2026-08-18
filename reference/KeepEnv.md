# Keep selected objects in an environment and remove everything else

Keeps only the objects specified in `Keep` within `Env` and removes all
other objects in that environment. Optionally returns a summary of what
was removed.

## Usage

``` r
KeepEnv(
  Keep,
  Env = parent.frame(),
  DryRun = FALSE,
  Invert = FALSE,
  Quiet = FALSE
)
```

## Arguments

- Keep:

  Character vector of object names to keep.

- Env:

  Environment to clean. Defaults to the calling environment.

- DryRun:

  If TRUE, does not remove anything and only reports what would be
  removed.

- Invert:

  If TRUE, removes only `Keep` and keeps everything else.

- Quiet:

  If TRUE, suppresses messages.

## Value

Invisible list with `kept` and `removed` vectors.

## Details

A long analysis session accumulates intermediates - raw imports,
temporary merges, loop variables, half-finished plots - and only a few
of those objects are worth carrying forward. Saving the whole workspace
preserves all of it, producing `.RData` files that are large, slow to
load, and full of objects whose provenance nobody remembers.

This is the inverse of writing out a long
[`rm()`](https://rdrr.io/r/base/rm.html) call: name the handful of
objects that matter and everything else goes. Because the list is
allow-list rather than deny-list, objects created after the code was
written are cleaned up too, instead of surviving because nobody
remembered to add them to the [`rm()`](https://rdrr.io/r/base/rm.html).

Run it with `DryRun = TRUE` first. It reports exactly what would be
removed without touching anything, which is worth doing whenever the
environment holds something expensive to recompute.

## Examples

``` r
# An analysis environment: a few results among many intermediates
env_Analysis <- new.env()
local({
  df_Raw <- data.frame(id = 1:5, value = rnorm(5))
  df_Clean <- df_Raw[!is.na(df_Raw$value), ]
  tmp_merge <- df_Clean
  i <- 3
  scratch_vector <- 1:100
  model_Final <- lm(value ~ id, data = df_Clean)
  df_Results <- data.frame(term = "id", estimate = coef(model_Final)[2])
}, envir = env_Analysis)

ls(env_Analysis)
#> [1] "df_Clean"       "df_Raw"         "df_Results"     "i"             
#> [5] "model_Final"    "scratch_vector" "tmp_merge"     

# Check what would go, before anything is removed
preview <- KeepEnv(
  Keep = c("df_Results", "model_Final"),
  Env = env_Analysis,
  DryRun = TRUE
)
#> Dry run: would remove 5 object(s).
preview$removed
#> [1] "df_Clean"       "df_Raw"         "i"              "scratch_vector"
#> [5] "tmp_merge"     

# Then do it for real, and save a workspace holding only what matters
KeepEnv(c("df_Results", "model_Final"), Env = env_Analysis)
#> Removing 5 object(s).
ls(env_Analysis)
#> [1] "df_Results"  "model_Final"

save(list = ls(env_Analysis), envir = env_Analysis,
     file = file.path(tempdir(), "analysis_results.RData"))

# `Invert = TRUE`: drop the objects named, keep the rest
env_Other <- new.env()
assign("df_Huge", data.frame(x = 1:10), envir = env_Other)
assign("df_Small", data.frame(x = 1:2), envir = env_Other)
KeepEnv("df_Huge", Env = env_Other, Invert = TRUE)
#> Removing 1 object(s).
ls(env_Other)
#> [1] "df_Small"
```
