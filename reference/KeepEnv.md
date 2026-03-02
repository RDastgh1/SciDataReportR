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
