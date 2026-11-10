# SOM + latent profile clustering pipeline (with AHP and distance baselines)

End-to-end pipeline to:

- Standardize variables using SciDataReportR::CreateZScoreObject() or a
  supplied Z-score object.

- Fit a Self-Organizing Map (SOM; kohonen) on complete cases.

- Generate aweSOM visualizations (Circular, Line, Cloud) with optional
  relabeling using variable labels from the original data frame.

- Cluster SOM codebook vectors using latent profile analysis (tidyLPA /
  mclust backend).

- In `method = "exploratory"`, fit a grid of models and select a
  recommended solution using an Analytic Hierarchy Process (AHP)-style
  index combining AIC, BIC, and Entropy.

- In `method = "finalize"`, fit a user-specified model and number of
  profiles.

- Map node-level clusters and posterior probabilities back to
  individuals.

- Store training variable summaries used later to quantify whether
  projected cohorts fall outside the original training range.

This supports a train once, project many clinical phenotyping workflow:
learn phenotype structure in a training cohort, then project new cohorts
into the fixed phenotype space without reclustering.

Ideal use: correlated continuous clinical or biomarker measures where a
topology-preserving map is clinically informative before model-based
profiles.

Missing data:

- SOM and clustering are fit only on rows with complete Z-scores.

- The returned `DataWithClusters` has exactly the original rows and
  columns plus canonical-cluster and hard/soft confidence columns; rows
  not used in SOM/LPA get NA for every assignment diagnostic.

- The returned `ProbFit$individual` is also full length, preserving one
  row per input row with NA posterior probabilities for rows excluded
  from SOM/LPA.

Z-score behavior:

- `ZScoreType = "Center and Scale"/"Center Only"/"Scale Only"` computes
  Z-scores from `df` via
  [`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md).

- `ZScoreType = "ZScoreObj"` projects Z-scores using an external
  `ZScoreObj` via
  [`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md).

- `ZScoreType = "PreZScored"` uses existing Z-score columns in `df`
  as-is and does not re-zscore.

Readable SOM + Mclust workflow wrapper for
`CreateClusterModel_SOM_MClust()`.

Compatibility wrapper for `Pipeline_SOM_MClust()`.

Deprecated alias for `CreateClusterModel_SOM_MClust()`.

## Usage

``` r
CreateClusterModel_SOM_MClust(
  data,
  variables = NULL,
  method = c("exploratory", "finalize", "explore"),
  k_range = 2:10,
  models = c(1, 2, 3, 6),
  final_k = NULL,
  final_model = NULL,
  ClusterVariableName = "Cluster",
  ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj",
    "PreZScored"),
  ZScoreObject = NULL,
  som_xdim = NULL,
  som_ydim = NULL,
  map_from = NULL,
  grid_screen_args = list(),
  som_topo = "hexagonal",
  som_neigh = "gaussian",
  som_epochs = 100L,
  node_weighting = c("uniform", "occupancy"),
  assignment = c("soft", "hard"),
  soft_n_nodes = 3L,
  soft_kernel = c("inverse", "gaussian"),
  soft_bandwidth = NULL,
  boundary_margin_threshold = 0.2,
  boundary_entropy_threshold = 0.6,
  seed_som = 934521L,
  seed_mixture = 93421L,
  Relabel = TRUE,
  ZScorePrefix = "Z_",
  ZScoreVars = NULL,
  id_var = NULL,
  mixture_progress = FALSE,
  mixture_em_itmax = 100L,
  mixture_em_tol = 1e-05,
  mixture_bootstrap_n = 0L,
  mixture_timeout_seconds = 120,
  mixture_drop_zero_sd = TRUE,
  mixture_zero_sd_tol = 1e-08,
  skip_model_after_n_failures = 2L,
  slow_fit_seconds = 120,
  candidate_min_nodes_per_cluster = 10,
  min_nodes_per_parameter = 3,
  min_participants_per_node = 10,
  mixture_min_node_people = 0L,
  min_cluster_prop = NULL,
  criteria_weights = NULL,
  far_distance_quantile = 0.95,
  clear_membership_threshold = 0.7,
  stability_resamples = 0L,
  stability_type = c("subsample", "seed"),
  stability_candidates = c("all", "best_per_k"),
  stability_top_n = NULL,
  seed_stability = 934522L,
  stability_fractions = 0.9,
  stability_progress = FALSE,
  stability_cores = NULL,
  df = lifecycle::deprecated(),
  id_col = lifecycle::deprecated(),
  .NodeClusterFn = NULL,
  .CandidatePairs = NULL,
  seed_lpa = lifecycle::deprecated(),
  stability_seed = lifecycle::deprecated(),
  lpa_em_itmax = lifecycle::deprecated(),
  lpa_em_tol = lifecycle::deprecated(),
  lpa_timeout_seconds = lifecycle::deprecated(),
  lpa_progress = lifecycle::deprecated(),
  lpa_bootstrap_n = lifecycle::deprecated(),
  lpa_drop_zero_sd = lifecycle::deprecated(),
  lpa_zero_sd_tol = lifecycle::deprecated(),
  lpa_min_node_n = lifecycle::deprecated(),
  min_nodes_per_cluster = lifecycle::deprecated(),
  high_dist_quantile = lifecycle::deprecated(),
  low_prob_threshold = lifecycle::deprecated()
)

Pipeline_SOM_MClust(...)

Pipeline_SOMClust(...)

CreateSOMClusterModel(...)
```

## Arguments

- data:

  Data frame containing the variables to be used in SOM and clustering.

- variables:

  Optional character vector of variable names. If NULL, numeric
  variables are auto-detected using
  `SciDataReportR::getNumVars(df, Ordinal = FALSE)`. In
  `ZScoreType = "PreZScored"`, this can also be NULL if you supply
  `ZScoreVars` or if Z-score columns can be auto-detected by prefix.

- method:

  One of `"exploratory"` (default) or `"finalize"`. In `"exploratory"`,
  a grid of models is fit and AHP chooses the recommended solution. In
  `"finalize"`, the user must specify `final_k` and `final_model`.

- k_range:

  Integer vector of numbers of clusters/profiles to consider in
  exploratory mode. Default `2:10`.

- models:

  Case-insensitive Mclust names (EEI, VVI, EEE, VVV) or IDs. Integer
  vector of model specifications for tidyLPA's mclust backend. Model 1
  uses equal variance and zero covariance; model 2 uses varying variance
  and zero covariance; model 3 uses equal variance and equal covariance;
  and model 6 uses varying variance and varying covariance.
  Zero-covariance models assume conditional independence between
  variables within each cluster. Equal parameters are shared across
  clusters; varying parameters are cluster-specific. Supported values
  and the default are `c(1, 2, 3, 6)`. Models 4 and 5 require OpenMx and
  are intentionally unsupported.

- final_k:

  Integer; number of profiles for `method = "finalize"`.

- final_model:

  Case-insensitive Mclust name or integer ID. Integer; model
  specification for `method = "finalize"` (should be one of `models`).

- ClusterVariableName:

  Name of the cluster column in the output. Defaults to `"Cluster"`. If
  this column already exists in `df`, it is overwritten (with a
  message).

- ZScoreType:

  One of:

  - `"Center and Scale"` (default)

  - `"Center Only"`

  - `"Scale Only"`

  - `"ZScoreObj"` (use an existing ZScore object)

  - `"PreZScored"` (use existing Z-score columns in df as-is)

- ZScoreObject:

  Optional ZScoreObj (from
  [`CreateZScoreObject()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateZScoreObject.md)
  or
  [`ProjectZScore()`](https://rdastgh1.github.io/SciDataReportR/reference/ProjectZScore.md))
  to use when `ZScoreType = "ZScoreObj"`.

- som_xdim, som_ydim:

  Optional integers for SOM grid dimensions. If NULL (and no
  `map_from`), the grid is chosen from the map alone by
  [`ScreenSOMGrid()`](https://rdastgh1.github.io/SciDataReportR/reference/ScreenSOMGrid.md)
  before any mixture is fitted, so no value of k and no covariance model
  influences it: a floor of 10 nodes per clustering variable, a ceiling
  from sparse nodes, measurement noise and a cap of 30 nodes per
  variable, and the 5 \* sqrt(n) target moved into that range. The
  screen is returned in `ModelInfo_SOM$grid_screen`; review it with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and override
  the choice by supplying both dimensions. Explicit dimensions are never
  resized; a warning is given when they fall below the
  10-nodes-per-variable floor.

- map_from:

  Optional exploratory `Pipeline_SOM_MClust` object. In a finalized fit,
  reuse only its resolved dimensions, topology, neighbourhood, and SOM
  seed. Clustering-variable names and order must match. Training
  preprocessing is deliberately recomputed on `data`; pass
  `ZScoreType = "ZScoreObj"` to reuse scaling explicitly.

- grid_screen_args:

  Optional named list of settings passed to the automatic grid screen
  (see
  [`ScreenSOMGrid()`](https://rdastgh1.github.io/SciDataReportR/reference/ScreenSOMGrid.md)
  and
  [`SelectSOMGrid()`](https://rdastgh1.github.io/SciDataReportR/reference/SelectSOMGrid.md)),
  e.g. `list(jitter_sd = 0.3, max_nodes = 400, refits = 3)`. Ignored
  when dimensions or `map_from` are supplied.

- som_topo:

  SOM topology for
  [`kohonen::somgrid()`](https://rdrr.io/pkg/kohonen/man/unit.distances.html),
  default `"hexagonal"`.

- som_neigh:

  SOM neighbourhood function, default `"gaussian"`.

- som_epochs:

  Positive integer batch SOM training epochs, also used for grid
  screening and stability refits. Default 100.

- node_weighting:

  Whether codebook nodes enter the mixture with equal map-space weight
  (`"uniform"`, default) or their training participant occupancy
  (`"occupancy"`). Occupancy weighting supports Models 1–3 only.

- assignment:

  Participant phenotype definition: `"hard"` maps the Mclust class of
  the best-matching SOM node; `"soft"` sums the Mclust node posterior
  probabilities across nearby codebook nodes and chooses the largest
  participant-level probability. Both outputs are always retained; the
  default `"soft"` controls the canonical `Cluster`.

- soft_n_nodes:

  Number of nearest codebook nodes used by soft projection. Their
  inverse-squared standardized-feature distances are normalized to
  proximity weights. Default `3`.

- soft_kernel:

  Soft assignment kernel: inverse squared distance (default) or
  Gaussian. Inverse is fuzzy k-nearest neighbours over prototypes with
  fuzzifier 2 (Keller, Gray and Givens, 1985).

- soft_bandwidth:

  Positive Gaussian bandwidth. NULL uses the median training distance to
  the closest node, with 1 used for a zero median. The resolved
  bandwidth is stored in Specification\$projection and reused by
  projection; refits resolve their own bandwidth.

- boundary_margin_threshold:

  Minimum difference between the largest and second-largest soft class
  probabilities before a participant is flagged as a boundary
  participant. Default `0.20`.

- boundary_entropy_threshold:

  Maximum normalized soft-membership entropy before a participant is
  flagged as a boundary participant. Default `0.60`. A flag is set if
  either boundary criterion is met.

- seed_som, seed_mixture:

  Integer seeds for SOM and LPA steps (defaults 934521 and 93421).

- Relabel:

  Logical; if TRUE (default), aweSOM plots are relabeled using variable
  labels from the *original* `df` (via Hmisc or sjlabelled when
  available) by stripping the Z-score prefix.

- ZScorePrefix:

  Character prefix used for Z-score columns when
  `ZScoreType = "PreZScored"`. Default `"Z_"`.

- ZScoreVars:

  Optional character vector of Z-score column names to use when
  `ZScoreType = "PreZScored"`. If NULL, the function attempts to infer
  them from `variables` or by detecting columns starting with
  `ZScorePrefix`.

- id_var:

  Optional character scalar. If provided and present in `df`, this
  column is carried into `ProbFit$individual` for convenience.

- mixture_progress:

  Logical; if TRUE, print short progress messages while fitting
  model/profile combinations.

- mixture_em_itmax:

  Integer; maximum number of EM iterations passed to
  [`mclust::emControl()`](https://mclust-org.github.io/mclust/reference/emControl.html).
  Use NULL to leave mclust defaults unchanged.

- mixture_em_tol:

  Numeric; EM convergence tolerance passed to
  [`mclust::emControl()`](https://mclust-org.github.io/mclust/reference/emControl.html).
  Use NULL to leave mclust defaults unchanged.

- mixture_bootstrap_n:

  Bootstrap resamples tidyLPA runs for each candidate fit: a bootstrap
  likelihood-ratio test (BLRT) against k - 1 profiles and
  MclustBootstrap. Default `0` skips them, which removes most of the
  fitting time; `BLRTStatistic` and `BLRTPValue` are then `NA`. Set to,
  e.g., `100` to report the BLRT. Stability refits never run them.

- mixture_timeout_seconds:

  Optional timeout in seconds for individual LPA fits. Use NULL to
  disable timeouts.

- mixture_drop_zero_sd:

  Logical; if TRUE, remove SOM code dimensions with near-zero standard
  deviation before LPA.

- mixture_zero_sd_tol:

  Numeric tolerance used when `mixture_drop_zero_sd = TRUE`.

- skip_model_after_n_failures:

  Optional integer; skip a model family after this many failures.

- slow_fit_seconds:

  Optional runtime threshold used to flag slow LPA fits in diagnostics.

- candidate_min_nodes_per_cluster:

  Optional minimum average SOM nodes per cluster. Default `10`. The grid
  is fixed before k is considered, so candidate profile counts with
  fewer nodes per cluster than this on the chosen map are still fitted
  but trigger one warning listing them. Use NULL to turn the warning
  off.

- min_nodes_per_parameter:

  Deprecated; only labels fit-table support. Positive number of codebook
  nodes per estimated free Gaussian-mixture parameter used to label
  candidates in the fit table (`SupportStatus`). Default `3`. It no
  longer sizes the map; candidates with fewer nodes than parameters are
  marked ineligible.

- min_participants_per_node:

  Number of complete participants per node used only for the
  `OccupancyCompliant` column of the fit table. Default `10`. It no
  longer sizes the map; sparse nodes are handled by the grid screen.

- mixture_min_node_people:

  Minimum number of training participants a SOM node must hold to enter
  the node-level mixture fit. Default `0` fits every node. Use `1` to
  leave empty nodes out, so unoccupied map regions do not count as
  observations in AIC/BIC. Nodes left out are classified afterwards by
  the fitted mixture model, so every node still carries a phenotype and
  posterior probabilities for projection.

- min_cluster_prop:

  Optional minimum share of training participants (0-1) in a candidate's
  smallest profile. Candidates below it are marked `Eligible = FALSE` in
  the fit table and are not recommended. Use it to keep only phenotypes
  large enough to describe and interpret.

- criteria_weights:

  Optional `Boundary` weights favour fewer boundary participants;
  Boundary is excluded from default weights. Weights for the composite
  rank index. `NULL` (default) weights the criteria equally. A named
  non-negative vector over `"AIC"`, `"BIC"`, `"Entropy"`, and
  `"Reproducibility"` sets weights directly. A square reciprocal
  pairwise-comparison matrix with those row and column names derives
  Analytic Hierarchy Process weights from its principal eigenvector and
  reports Saaty's consistency ratio (warning above 0.10).
  `"Reproducibility"` requires `stability_resamples > 0`.

- far_distance_quantile:

  Numeric value between 0 and 1 used to define high SOM-distance flags
  from the training distance distribution. Default is `0.95`.

- clear_membership_threshold:

  Numeric posterior probability threshold used to flag uncertain
  phenotype membership. Default is `0.70`.

- stability_resamples:

  Number of 90% participant subsample refits used to assess
  reproducibility for every successful exploratory candidate. Subsamples
  are drawn without replacement and reuse the reference model's resolved
  SOM grid dimensions and the reference SOM seed. Each refit maps the
  complete reference cohort before scoring participant assignments.
  Defaults to `0` (disabled); use `30` for an exploratory stability
  screen.

- stability_type:

  `subsample` (default), `seed`, or both. Seed refits retain every
  participant, vary SOM and mixture seeds from seed_stability, and keep
  the grid fixed. With both types, AHP reproducibility is the mean of
  their two mean ARIs (subsampling uses the primary retained fraction).
  Stability tables include StabilityType; cluster_recovery_summary holds
  per-cluster MeanJaccard by type. Only mean ARI is summarized.

- stability_candidates:

  Which candidates receive subsample stability: `"all"` (default) or
  `"best_per_k"`, the highest-ranked candidate (by the composite index
  without reproducibility, eligible first) at each number of profiles,
  so stability can inform the choice of k across the whole range.
  Candidates without stability keep `NA` reproducibility and are not
  recommended.

- stability_top_n:

  Optional integer. When supplied, stability is estimated only for the
  top `stability_top_n` candidates ranked by the composite index without
  reproducibility (eligible candidates first); the others keep `NA`
  stability and are not recommended. Use it to screen a large candidate
  grid, e.g. `stability_top_n = 5`.

- seed_stability:

  Integer seed for participant subsampling.

- stability_fractions:

  Unique participant-retention fractions for an opt-in
  deletion-stability curve. The default, `0.90`, preserves the ordinary
  90% refit screen. Each fraction refits the complete pipeline and
  projects the full reference cohort before all-participant ARI is
  scored.

- stability_progress:

  Logical; if TRUE, print subsample progress messages.

- stability_cores:

  Number of workers for stability refits. `NULL` uses a conservative
  estimate based on CPU, refit count, and memory: the smaller of 50%
  physical RAM and 75% currently available RAM. Explicit values are
  retained up to the scheduler/CPU limit and warn if their estimated
  memory exceeds that budget. Each worker holds a refit in memory.

- df:

  **Deprecated** (since 19.15.0). Use `data` instead.

- id_col:

  **Deprecated** (since 19.15.0). Use `id_var` instead.

- .NodeClusterFn:

  Internal. A function taking the SOM codebook matrix and returning a
  list with a `node_cluster` integer vector (one label per SOM node)
  and, optionally, `fit_table`, `ahp_best_row`, `recommendation`,
  `best_fit_name`, and `fit_plot`. When supplied, the SOM codebook is
  clustered by that function and the latent-profile grid is not fitted.
  Used by
  [`CreateClusterModel_SOM_HDBSCAN()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_HDBSCAN.md);
  not part of the user-facing API.

- .CandidatePairs:

  Internal. A data frame of `Model` and `Classes` pairs restricting the
  exploratory grid; used by stability refits so each resample fits only
  the shortlisted candidates.

- seed_lpa:

  Deprecated alias for `seed_mixture`.

- stability_seed:

  Deprecated alias for `seed_stability`.

- lpa_em_itmax:

  Deprecated alias for `mixture_em_itmax`.

- lpa_em_tol:

  Deprecated alias for `mixture_em_tol`.

- lpa_timeout_seconds:

  Deprecated alias for `mixture_timeout_seconds`.

- lpa_progress:

  Deprecated alias for `mixture_progress`.

- lpa_bootstrap_n:

  Deprecated alias for `mixture_bootstrap_n`.

- lpa_drop_zero_sd:

  Deprecated alias for `mixture_drop_zero_sd`.

- lpa_zero_sd_tol:

  Deprecated alias for `mixture_zero_sd_tol`.

- lpa_min_node_n:

  Deprecated alias for `mixture_min_node_people`.

- min_nodes_per_cluster:

  Deprecated alias for `candidate_min_nodes_per_cluster`.

- high_dist_quantile:

  Deprecated alias for `far_distance_quantile`.

- low_prob_threshold:

  Deprecated alias for `clear_membership_threshold`.

- ...:

  Arguments passed to `CreateClusterModel_SOM_MClust()`.

## Value

A list of class `"Pipeline_SOM_MClust"` with components:

- `method`, `vars_used`, `ZScoreType`, `ZScoreObject`, `ZScoreVars`,
  `ClusterVariableName`

- `DataWithClusters`: original `df` with canonical cluster, hard/soft
  assignment, soft membership, and boundary diagnostic columns appended

- `ProbFit$individual$SOM_Distance_LOO`: distance after removing each
  complete person from their own node mean, allowing the nearest other
  node instead. No map refit. ProjectionReference\$loo_distances stores
  these distances with SOM_Node and Cluster for arbitrary quantiles.
  ProjectionReference\$calibration stores q_star, training_rate,
  training_n, far_n, target_quantile, overall_cutoff and cluster_cutoffs
  for the combined overall OR cluster LOO rule.
  ProjectionReference\$cluster_profiles stores mean training z-score
  vectors by assigned cluster for group profile comparisons.

- `ProbFit$boundary_summary`: BoundaryReason, N and Share among complete
  participants. `ProbFit$boundary_threshold_grid` gives margin, entropy
  and BoundaryPct for 15 threshold combinations. Defaults 0.20 / 0.60
  are review conventions, not estimates.

- `ModelInfo_SOM$diagnostics`: map_quality (one-row tibble of
  quantization error, variance captured, topographic error and sparse
  node counts/shares), u_matrix (NodeID and U), and border_alignment
  (border/interior mean U, difference and one-sided 999-permutation
  p-value, conditional on permutations with both border and interior
  nodes; ValidPermutations reports that count).

- `ModelInfo_MClust$ahp_sensitivity`: weight schemes, normalized
  weights, top candidate and rank of the current recommendation.
  `BoundaryPct` in the fit table is the share of complete boundary
  participants (0–1).

- `fit_plot`: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model
  (plus reproducibility when subsample stability is enabled)

- `ModelInfo_SOM`: list with `som_model`, `som_codes`, `som_grid`,
  `som_grid_info`, `grid_screen` (the
  [`ScreenSOMGrid()`](https://rdastgh1.github.io/SciDataReportR/reference/ScreenSOMGrid.md)
  result when the grid was chosen automatically, otherwise NULL),
  `training_variable_summary`, `SOMFit` (distance diagnostics,
  baselines, and per-cluster flags), `plots` (aweSOM plots)

- `ModelInfo_MClust`: list with `lpa_models`, `fit_table`, `AHP`
  information, and `diagnostics` for LPA warnings, failures, runtimes,
  and preprocessing

- `ModelInfo_MClust$Stability`: subsample replicate, cluster recovery,
  and summary tables when `stability_resamples > 0`

- `ProbFit`: list with `node` (node-level posterior probabilities),
  `individual` (full-length per-person mapping and probabilities), and
  probability plots

## Details

`ModelInfo_MClust$fit_table` uses the fit indices returned by tidyLPA's
mclust backend. AIC, AWE, BIC, CAIC, CLC, KIC, and SABIC are
likelihood/information criteria on tidyLPA's -2 log-likelihood scale,
for which lower values are preferred when comparing candidate fits
estimated on the same data. ICL follows mclust's convention instead
(`-BIC` plus twice the summed log maximum posterior, matching
[`mclust::icl()`](https://mclust-org.github.io/mclust/reference/icl.html)),
so higher ICL values are preferred: ICL and BIC point in opposite
directions, and
[`PlotClusterFitReview()`](https://rdastgh1.github.io/SciDataReportR/reference/PlotClusterFitReview.md)
labels each panel with its preferred direction. ICL is reported for
review only and is not part of the composite index. Entropy summarizes
classification separation (higher is better). `MinProfileNodeN` and
`MaxProfileNodeN` are integer SOM-node counts in the smallest and
largest profile; the corresponding `*Proportion` fields divide those
counts by the total number of SOM nodes. `BLRTStatistic` and
`BLRTPValue` compare a k-profile model with k - 1 profiles under the
same covariance model. A small p-value supports the added profile's
statistical fit, but does not establish clinical meaning or
reproducibility; it is unavailable for k = 1. The raw `BLRTPValue`
remains in the table. Its candidate-review panel displays
`-log10(BLRTPValue)` with a dashed `0.05` p-value reference at
`-log10(0.05)`.

Lifecycle settings are intentionally strict: `k_range` and `models` are
exploratory-only, while `final_k` and `final_model` are finalize-only.
Supplying settings from both modes is an error, even when an exploratory
setting equals its default.

Stability is assessed by refitting the full SOM/LPA pipeline on
independent participant subsamples drawn without replacement, then
mapping the complete reference cohort onto each refit.
`StabilityARI_Mean` is the primary reproducibility score: the mean
adjusted Rand index across these full-cohort participant partitions. It
can be negative when agreement is worse than chance.
`stability_fractions` can request an explicit deletion curve; the
existing 90% fit-table fields remain the primary review metrics, while
`Stability$deletion_curve` records every requested retention level.
Co-membership, split/merge, and profile-recovery diagnostics describe
whether instability is localized despite imperfect exact membership;
they never replace all-participant ARI. `StabilityJaccard_Mean` is the
mean label-matched, per-profile Jaccard recovery and is retained as a
descriptive diagnostic. `ReproducibilityScore` equals
`StabilityARI_Mean`; `StabilitySuccessRate` is reported separately.
`Reproducibility_scaled` is the candidate-table min-max rescaling used
only by the AHP index and is not an independently interpretable
reproducibility measure.

The composite rank index (stored as `ahp_index`) is computed by:

1.  Z-scaling AIC, BIC, Entropy, and (when stability is run)
    reproducibility across candidate solutions (AIC/BIC are negated so
    that higher scaled scores are preferred).

2.  Taking their weighted mean. Weights are equal unless
    `criteria_weights` supplies them, either directly or as an Analytic
    Hierarchy Process pairwise-comparison matrix.

3.  Recommending the highest-scoring candidate among those marked
    `Eligible` (see `min_cluster_prop`).

The chosen weights, their method, and any consistency ratio are returned
in `ModelInfo_MClust$AHP`.

`MinProfileParticipantN` and `MinProfileParticipantProportion` give the
size of the smallest profile in training participants rather than SOM
nodes.

Phenotypes are learned by Mclust on the SOM codebook. Hard assignment
uses the best-matching node's Mclust class. Soft assignment uses
normalized inverse-squared codebook-distance weights across the nearest
nodes to combine their Mclust posterior probabilities. `SoftEntropy` is
the normalized entropy of that participant-level probability vector
(entropy divided by `log(k)`) and `SoftMargin` is the largest minus
second-largest probability. `prob_*` columns describe a SOM node's
Mclust posterior, whereas `SoftProb_*` columns describe a participant's
proximity-weighted membership. Boundary flags are confidence diagnostics
only: they do not remove, relabel, or suppress a participant during
training, stability scoring, or projection.

Changing `assignment` does not change SOM training or the fitted node
mixture. With the same seeds and selected model, node posteriors and
likelihood-based fit indices are identical. Both modes return a single
participant class; NodeN counts complete training participants on the
closest node, and Flag_OnePersonNode is NodeN == 1 (missing for
incomplete rows). Both are returned in ProbFit\$individual and
DataWithClusters. Probabilities measure ambiguity between clusters;
distance to node measures fit to the map. soft mode takes the largest
`SoftProb_*`. Hard and soft classes can therefore agree even when their
probability vectors differ, especially when nearby nodes belong to the
same phenotype. Exact soft-probability ties select the first class
deterministically. `Flag_NodeDisagreement` compares the canonical
participant class with the best-matching node's class; unmapped rows are
`NA`.

Each new fit uses one deterministic, seeded full-cohort batch SOM. Its
initialization codebook, preprocessing, grid, topology, neighbourhood,
epochs, codebooks, and node-level mixture are retained so projection is
completely frozen. Earlier fitted objects remain projectable because
projection uses their frozen map and mixture rather than refitting
either component.

`SOM_Distance` is the Euclidean distance, in Z-score units, from a
participant to their best-matching unit. Its square is the sum of the
squared per-variable residuals, returned in
`ModelInfo_SOM$SOMFit$residuals` (`Resid_<variable>`).
`ModelInfo_SOM$SOMFit$variable_contribution` summarizes each variable's
share of training misfit, and each participant's largest contributor is
recorded in `Top_Distance_Variable` and `Top_Distance_Share`.

LPA model/profile combinations are fit one at a time so that failed or
warning-producing solutions are captured in diagnostics instead of
blocking the entire pipeline. Successful fits are retained and failed
fits are listed in `ModelInfo_MClust$diagnostics`.

## Model review and refinement

Start with `ModelInfo_MClust$fit_table`, the AHP recommendation, and
`fit_plot` to compare candidate model specifications and numbers of
profiles. AHP is a useful starting point, not an automatic final
decision: retain solutions whose fit, cluster size, and clinical
interpretability are all reasonable.

The Circular and Line SOM widgets show how the analysis variables vary
across the map; the Cloud widget shows the observations within map
cells. Use them to assess whether the candidate solution represents
coherent, interpretable phenotype patterns. Inspect
`SOMFit$node_occupancy` for empty or sparse cells, and the SOM-distance
plots for clusters with systematically poor map representation. The
posterior-probability plots in `ProbFit$plots` identify uncertain node-
or individual-level assignments.

If these diagnostics suggest a weak solution, reconsider the variables,
SOM grid dimensions, model family, or candidate `k` range and refit in
exploratory mode. Once a solution is selected, refit it with
`method = "finalize"`, `final_k`, and `final_model` to create the
reusable model for projection.

Set `stability_resamples` to a positive value to assess internal
reproducibility. Each 90% subsample refits the same candidate SOM/LPA
solution and projects the original participants back into it. Mean
adjusted Rand index summarizes whole-partition agreement; cluster-wise
Jaccard recovery identifies phenotypes that dissolve despite good
overall agreement. This is internal reproducibility, not
independent-cohort validation.

## Reviewing the candidates

Exploratory mode fits every combination of `k` and model family, so the
comparison is made once, on one table, rather than by refitting by hand.
`ModelInfo_MClust$fit_table` holds every candidate solution with its
information criteria, entropy, subsample reproducibility, and the
combined AHP index side by side.

A candidate missing from that table is never a mystery:
`ModelInfo_MClust$diagnostics$lpa_fit_diagnostics` records each fit's
status, runtime, warnings, and errors, so model families that fail to
converge or exceed the fit timeout show up there.

Once a solution is chosen, refit it on its own with
`method = "finalize"` to get the reusable, projectable model.

## Stability output

Stability assesses internal reproducibility by full-pipeline 90%
participant subsampling without replacement. For each replicate, 90% of
complete participants are selected once, all preprocessing, any
reduction (PCA, MCA, or SOM), and the selected clustering method are
refit. Most pipelines compare assignments for sampled participants.
SOM + Mclust instead maps the complete reference cohort onto every refit
and compares the resulting full participant partition.
`Stability$settings$comparison_scope` records the policy used. It is an
internal sensitivity analysis, not independent-cohort validation.

`Stability$settings` records the analysis provenance:

- `resamples`: requested number of 90% subsample refits.

- `seed`: seed used to derive reproducible sampling and model seeds.

- `refit_scope`: whether agreement is calculated directly in the
  refitted sample or after refit projection onto the reference cohort.

- `comparison_scope`: the participants entering replicate agreement
  metrics; most pipelines use `"sampled_participants"`, while SOM +
  Mclust uses `"all_reference_participants"`.

- `resample_type`: `"subsample_without_replacement"` for the primary
  stability analysis.

- `resample_fraction`: the retained participant fraction, `0.90`.

- `resample_fractions`: every requested retention fraction. The default
  is `0.90`; additional values create an opt-in deletion-stability curve
  without replacing the 90% primary screen.

- `coassignment_limit`: maximum number of complete training participants
  (2,000) for which the full pairwise co-assignment matrix is
  calculated.

- `noise_policy`: whether the method has noise labels. For ordinary
  methods it is `"all clusters included"`; HDBSCAN variants retain noise
  in ARI but exclude it from phenotype-level Jaccard recovery.

- `requested_cores`, `resolved_cores`, and `worker_policy`: the
  requested and resolved refit workers and the reason for the
  resolution. Automatic resolution is capped by CPU, refit count, and
  the smaller of 50% physical RAM and 75% currently available RAM.
  Pairwise-memory methods are run serially automatically.

- `estimated_bytes_per_worker`, `estimated_total_bytes`,
  `physical_memory_bytes`, `available_memory_bytes`, and
  `safe_memory_budget_bytes`: the conservative memory audit used to
  resolve automatic workers. An explicit worker request is preserved
  (within the CPU/scheduler cap) but warns if it exceeds the estimated
  safe budget.

`Stability$replicates` has one row per requested refit. `Model` and
`Classes` identify the selected candidate (for HDBSCAN, `Classes` is the
data-derived extracted count); `Replicate` is its sequence number;
`Status` is `"success"` or a failure status; `SamplingSeed` and
`ModelSeed` identify the participant draw and refit initialization
policy; and `Error` contains the error message for an unsuccessful
refit. `ARI` is the primary successful-refit metric: it is invariant to
cluster numbering, higher is better, and it can be negative when
agreement is worse than chance.

`Stability$cluster_recovery` has one row for each reference `Cluster` in
each successful `Replicate`. `Jaccard` is the recovery of that reference
cluster after matching it to the subsample cluster with the largest
Jaccard overlap; higher is better and one is exact recovery. `Model` and
`Classes` again identify the fitted candidate.

`Stability$summary` combines successful subsample refits:
`StabilitySuccessRate` is successful refits divided by requested refits,
an operational reliability measure that does not enter the
reproducibility score. `StabilityARI_Mean` and `StabilityARI_P05` are
the mean and fifth percentile of ARI. `StabilityJaccard_Mean` and
`StabilityJaccard_Min` are respectively the mean cluster-level recovery
and the minimum of the cluster-specific mean recoveries across
replicates. `ReproducibilityScore` does not include success rate. For
most pipelines it is the mean of finite `StabilityARI_Mean` and
`StabilityJaccard_Mean` values; SOM + Mclust uses `StabilityARI_Mean`
directly because it selects on agreement of the final participant
partition.

SOM + Mclust additionally returns `StabilityARI_NonBoundary`,
`StabilityARI_NonBoundary_P05`, and non-boundary Jaccard summaries.
These are secondary diagnostics calculated only on the fixed
reference-model mask with `SoftMargin` below its threshold or
`SoftEntropy` above its threshold excluded. The primary all-participant
ARI always retains these boundary participants. `boundary_summary`
reports their count and proportion overall and within each reference
cluster. Earlier posterior-based `Core` fields are retained for
compatibility as descriptive diagnostics; they are not the
boundary-aware sensitivity measure.

`Stability$failures` repeats the replicate columns for unsuccessful
refits, making fit failures auditable without mixing them with
successful metrics.

`Stability$participant_inclusion` is one row per complete reference
participant. `.row_id` identifies the merge-safe input row, `Cluster`
its reference assignment, `SuccessfulRefits` the number of usable
refits, and `InclusionProbability` the proportion of those refits in
which the participant returned to that cluster's label-matched subsample
cluster. `Model` and `Classes` identify the candidate. Higher inclusion
is better.

`Stability$cluster_inclusion` summarizes inclusion within each reference
`Cluster`: `MeanInclusion`, `P05Inclusion`, and `MinInclusion` are the
mean, fifth percentile, and minimum participant inclusion probabilities;
`Model` and `Classes` identify the candidate. Higher values indicate
that all, not only the average, of a cluster is recovered consistently.

`Stability$coassignment` is available only when the complete training
cohort has at most `coassignment_limit` participants. Each candidate
entry has a `status` of `"available"`, `"skipped"`, or
`"not_available"`; `reason` explains a non-available result; `matrix` is
the pairwise probability that two complete reference participants are
assigned together across successful subsample refits; and `.row_id` maps
matrix rows and columns to the same identifiers returned in
`DataWithClusters` and `ProbFit$individual`. Higher matrix values mean
more consistent pairwise co-membership. Where more than one candidate is
summarized, entries are named by its `Model_Classes` key. The matrix is
diagnostic only and is never used for selection.

SOM + Mclust additionally returns memory-safe final-partition
diagnostics for every cohort size. `co_membership_by_cluster` reports
the proportion of within-reference-cluster participant pairs retained
together (`Cohesion`) and the probability of co-assignment with
outsiders (`CrossMixing`). `co_membership_between_clusters` reports
cross-reference-cluster co-assignment probabilities.
`split_merge_reference` and `split_merge_refit` identify dominant
sources, shares, and effective numbers of groups, while
`profile_recovery` reports matched-cluster correlation and RMSE for mean
clustering-variable profiles on the frozen reference scale. These
explain the nature of membership changes; they do not supersede ARI.

`deletion_curve` summarizes ARI, VI, NMI, and Fowlkes–Mallows agreement
at each requested participant-retention fraction. It is useful for
describing sensitivity to increasing cohort perturbation, not for
selecting the least demanding fraction.

`Stability$plots` contains `cluster_recovery` (per-cluster Jaccard),
`partition_metrics` (the ARI distribution), and `cluster_inclusion`; it
also contains a co-assignment heatmap when the matrix is available.
These diagnostics complement rather than replace ARI and Jaccard: none
can turn a poorly reproducible cluster into a stable phenotype.

Metric sources: Hubert and Arabie (1985) define ARI; Jaccard (1901)
defines the overlap coefficient; and Monti et al. (2003) describe
resampling-based consensus co-assignment.

## References

Saaty TL. *The Analytic Hierarchy Process*. McGraw-Hill, 1980. Kohonen
T. Self-organized formation of topologically correct feature maps.
*Biological Cybernetics*. 1982;43:59-69. Scrucca L, Fop M, Murphy TB,
Raftery AE. mclust 5. *J Stat Softw*. 2016;71(11):1-29. Celeux G,
Soromenho G. An entropy criterion for assessing the number of clusters
in a mixture model. *J Classification*. 1996;13:195-212. Cavanaugh JE. A
large-sample model selection criterion based on Kullback's symmetric
divergence. *Stat Methodol*. 1999;61:165-180. Sclove SL. Application of
model-selection criteria to some problems in multivariate analysis.
*Psychometrika*. 1987;52:333-343. Nylund KL, Asparouhov T, Muthen BO.
Deciding on the number of classes in latent class analysis and growth
mixture modeling. *Struct Equ Modeling*. 2007;14:535-569. Akogul S,
Erisoglu M. An approach for determining the number of clusters in a
model-based cluster analysis. *Entropy*. 2017;19:452.

## Examples

``` r
# \donttest{
data("SimulatedPhenotypeData")
df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
vars_Numeric <- paste0("Var", 1:12)

# Example-only display helper (not exported by SciDataReportR): round the
# numeric columns and add a frozen header.
ShowFitTable <- function(x, height = "320px") {
  x <- dplyr::mutate(
    dplyr::as_tibble(x),
    dplyr::across(dplyr::where(is.numeric), \(v) round(v, 3))
  )
  htmltools::browsable(htmltools::HTML(as.character(
    FreezeTableHeader(x, height = height, full_width = TRUE)
  )))
}

# Exploratory mode: every combination of k and model family. Model 2 is
# omitted here because it is the slowest to fit and, on this data, reaches
# the fit timeout at every k without contributing a candidate.
review <- CreateClusterModel_SOM_MClust(
  data = df_Training,
  variables = vars_Numeric,
  method = "exploratory",
  k_range = 2:5,
  models = c(1, 3),
  som_xdim = 5,
  som_ydim = 5,
  stability_resamples = 2,
  seed_stability = 20260805,
  Relabel = FALSE,
  candidate_min_nodes_per_cluster = NULL
)
#> Warning: Explicit SOM dimensions (5 x 5 = 25 nodes) are below the floor of 10 nodes per clustering variable (120 nodes for 12 variables); diagnostics are retained for review.
#> Warning: The following SOM mixture candidates have fewer mixture nodes than free parameters and are ineligible for recommendation: Model 1, k = 2; Model 3, k = 2; Model 1, k = 3; Model 3, k = 3; Model 1, k = 4; Model 3, k = 4; Model 1, k = 5; Model 3, k = 5
#> Warning: No candidate is eligible for the composite recommendation; an internal fallback is retained only to construct the review object.

# Every candidate solution
ShowFitTable(review$ModelInfo_MClust$fit_table)


 Model 
```
