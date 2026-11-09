# Choose a SOM grid size from a grid screen

Applies the map-only selection rule to a
[`ScreenSOMGrid()`](https://rdastgh1.github.io/SciDataReportR/reference/ScreenSOMGrid.md)
result. Thresholds can be changed and the choice redone without
refitting any maps.

1.  Floor: `min_nodes_per_variable` nodes per clustering variable
    (default 10), so the mixture step has enough nodes to estimate
    spread in every dimension.

2.  Ceiling: the smallest of

    - the sparsity ceiling: the largest side before the smoothed
      one-person share exceeds `max_one_person_prop` or the smoothed
      empty share exceeds `max_empty_prop`;

    - the noise ceiling: the first side whose quantization error is
      below the expected distance between a person and a retest,
      `sqrt(sum(noise_sd^2))`; larger maps only separate people by
      differences smaller than measurement noise;

    - the cap: `max_nodes_per_variable` nodes per variable (default 30),
      or `max_nodes` if smaller.

3.  Choice: the 5 \* sqrt(n) target (Vesanto et al., 2000), moved into
    the allowed range. If the floor is above the ceiling, the target is
    kept between the two and a warning is given.

The chosen map is also checked: a warning is given if its topographic
error is high compared with the other sides screened (a folded map), or
if its trustworthiness is below 0.80.

## Usage

``` r
SelectSOMGrid(
  x,
  min_nodes_per_variable = 10,
  max_nodes_per_variable = 30,
  max_nodes = NULL,
  max_one_person_prop = 0.1,
  max_empty_prop = 0.05,
  noise_sd = x$settings$jitter_sd,
  select = c("target", "largest_feasible")
)
```

## Arguments

- x:

  A `SOMGridScreen` object.

- min_nodes_per_variable:

  Floor, in nodes per clustering variable.

- max_nodes_per_variable:

  Cap, in nodes per clustering variable.

- max_nodes:

  Optional absolute cap on the number of nodes.

- max_one_person_prop:

  Largest allowed smoothed share of one-person nodes.

- max_empty_prop:

  Largest allowed smoothed share of empty nodes.

- noise_sd:

  Measurement noise SD in z-score units (one value or one per variable)
  for the noise ceiling. Defaults to the jitter used in the screen.

- select:

  `"target"` (default) moves the 5 \* sqrt(n) target into the allowed
  range; `"largest_feasible"` takes the top of the allowed range.

## Value

The screen with updated `table` pass/fail columns and `selection`.
