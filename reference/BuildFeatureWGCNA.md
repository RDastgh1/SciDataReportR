# Build a weighted feature correlation network

Build a weighted feature correlation network

## Usage

``` r
BuildFeatureWGCNA(
  data,
  variables,
  sample_id = NULL,
  network_type = c("signed", "unsigned", "signed_hybrid"),
  soft_power = NULL,
  sft_rsq = 0.85,
  min_module_size = 30,
  deep_split = 2,
  merge_cut_height = 0.25,
  keep_tom = FALSE,
  seed = NULL,
  n_modules = NULL
)
```

## Arguments

- data:

  A data frame with samples in rows and features in columns.

- variables:

  Explicit feature columns to include.

- sample_id:

  Optional unique sample-ID column retained for trait joins.

- network_type:

  One of `"signed"`, `"unsigned"`, or `"signed_hybrid"`.

- soft_power:

  Optional soft-thresholding power.

- sft_rsq:

  Target scale-free topology fit for automatic power selection.

- min_module_size:

  Minimum dynamic-tree-cut module size.

- deep_split:

  Dynamic tree-cut sensitivity from 0 through 4.

- merge_cut_height:

  Eigengene dissimilarity for module merging when `n_modules = NULL`.

- keep_tom:

  Store the topological-overlap matrix.

- seed:

  Optional random seed.

- n_modules:

  Optional positive whole number of non-grey modules to return. `NULL`
  retains height-based merging. Otherwise, the dynamic tree cut runs
  first, then its non-grey eigengenes are clustered using Pearson
  dissimilarity (`1 - cor`) and average linkage and cut into `n_modules`
  groups. Grey remains unassigned and is not counted. If too few modules
  are detected, raise `deep_split` or lower `min_module_size`. Final
  groups use the standard WGCNA colors in decreasing size order (largest
  turquoise, second blue), with ties ordered by their alphabetically
  first original module color. `merge_cut_height` is ignored for this
  merging procedure.

## Value

A `feature_wgcna_obj` with module assignments, eigengenes, diagnostics,
hub statistics, and the matrix used for fitting. `params$n_modules`
stores the requested count (or `NULL`), and `n_modules` stores the final
non-grey module count. Eigengenes and memberships reflect the final
assignments.
