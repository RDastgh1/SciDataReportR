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
  seed = NULL
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

  Eigengene dissimilarity for module merging.

- keep_tom:

  Store the topological-overlap matrix.

- seed:

  Optional random seed.

## Value

A `feature_wgcna_obj` with module assignments, eigengenes, diagnostics,
hub statistics, and the matrix used for fitting.
