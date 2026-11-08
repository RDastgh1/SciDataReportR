#' SOM + latent profile clustering pipeline (with AHP and distance baselines)
#'
#' @description
#' End-to-end pipeline to:
#' - Standardize variables using SciDataReportR::CreateZScoreObject() or a supplied
#'   Z-score object.
#' - Fit a Self-Organizing Map (SOM; kohonen) on complete cases.
#' - Generate aweSOM visualizations (Circular, Line, Cloud) with optional
#'   relabeling using variable labels from the original data frame.
#' - Cluster SOM codebook vectors using latent profile analysis (tidyLPA /
#'   mclust backend).
#' - In \code{method = "exploratory"}, fit a grid of models and select a
#'   recommended solution using an Analytic Hierarchy Process (AHP)-style
#'   index combining AIC, BIC, and Entropy.
#' - In \code{method = "finalize"}, fit a user-specified model and number
#'   of profiles.
#' - Map node-level clusters and posterior probabilities back to individuals.
#' - Store training variable summaries used later to quantify whether projected
#'   cohorts fall outside the original training range.
#'
#' This supports a train once, project many clinical phenotyping workflow: learn
#' phenotype structure in a training cohort, then project new cohorts into the
#' fixed phenotype space without reclustering.
#'
#' Ideal use: correlated continuous clinical or biomarker measures where a
#' topology-preserving map is clinically informative before model-based profiles.
#'
#' Missing data:
#' - SOM and clustering are fit only on rows with complete Z-scores.
#' - The returned \code{DataWithClusters} has exactly the original rows and
#'   columns plus canonical-cluster and hard/soft confidence columns; rows not
#'   used in SOM/LPA get NA for every assignment diagnostic.
#' - The returned \code{ProbFit$individual} is also full length, preserving
#'   one row per input row with NA posterior probabilities for rows excluded
#'   from SOM/LPA.
#'
#' Z-score behavior:
#' - \code{ZScoreType = "Center and Scale"/"Center Only"/"Scale Only"} computes
#'   Z-scores from \code{df} via \code{CreateZScoreObject()}.
#' - \code{ZScoreType = "ZScoreObj"} projects Z-scores using an external
#'   \code{ZScoreObj} via \code{ProjectZScore()}.
#' - \code{ZScoreType = "PreZScored"} uses existing Z-score columns in \code{df}
#'   as-is and does not re-zscore.
#'
#' @param data Data frame containing the variables to be used in SOM and
#'   clustering.
#' @param variables Optional character vector of variable names. If NULL,
#'   numeric variables are auto-detected using
#'   \code{SciDataReportR::getNumVars(df, Ordinal = FALSE)}. In
#'   \code{ZScoreType = "PreZScored"}, this can also be NULL if you supply
#'   \code{ZScoreVars} or if Z-score columns can be auto-detected by prefix.
#' @param method One of \code{"exploratory"} (default) or \code{"finalize"}.
#'   In \code{"exploratory"}, a grid of models is fit and AHP chooses the
#'   recommended solution. In \code{"finalize"}, the user must specify
#'   \code{final_k} and \code{final_model}.
#' @param k_range Integer vector of numbers of clusters/profiles to consider
#'   in exploratory mode. Default \code{2:10}.
#' @param models Case-insensitive Mclust names (EEI, VVI, EEE, VVV) or IDs. Integer vector of model specifications for tidyLPA's mclust
#'   backend. Model 1 uses equal variance and zero covariance; model 2 uses
#'   varying variance and zero covariance; model 3 uses equal variance and
#'   equal covariance; and model 6 uses varying variance and varying
#'   covariance. Zero-covariance models assume conditional independence
#'   between variables within each cluster. Equal parameters are shared across
#'   clusters; varying parameters are cluster-specific. Supported values and
#'   the default are \code{c(1, 2, 3, 6)}. Models 4 and 5 require OpenMx and
#'   are intentionally unsupported.
#' @param final_k Integer; number of profiles for \code{method = "finalize"}.
#' @param final_model Case-insensitive Mclust name or integer ID. Integer; model specification for \code{method = "finalize"}
#'   (should be one of \code{models}).
#' @param ClusterVariableName Name of the cluster column in the output. Defaults to
#'   \code{"Cluster"}. If this column already exists in \code{df}, it is
#'   overwritten (with a message).
#' @param ZScoreType One of:
#'   \itemize{
#'     \item \code{"Center and Scale"} (default)
#'     \item \code{"Center Only"}
#'     \item \code{"Scale Only"}
#'     \item \code{"ZScoreObj"} (use an existing ZScore object)
#'     \item \code{"PreZScored"} (use existing Z-score columns in df as-is)
#'   }
#' @param ZScoreObject Optional ZScoreObj (from \code{CreateZScoreObject()} or
#'   \code{ProjectZScore()}) to use when \code{ZScoreType = "ZScoreObj"}.
#' @param som_xdim,som_ydim Optional integers for SOM grid dimensions. If NULL
#'   (and no \code{map_from}), the grid is chosen from the map alone by
#'   [ScreenSOMGrid()] before any mixture is fitted, so no value of k and no
#'   covariance model influences it: a floor of 10 nodes per clustering
#'   variable, a ceiling from sparse nodes, measurement noise and a cap of 30
#'   nodes per variable, and the 5 * sqrt(n) target moved into that range. The
#'   screen is returned in \code{ModelInfo_SOM$grid_screen}; review it with
#'   \code{plot()} and override the choice by supplying both dimensions.
#'   Explicit dimensions are never resized; a warning is given when they fall
#'   below the 10-nodes-per-variable floor.
#' @param grid_screen_args Optional named list of settings passed to the
#'   automatic grid screen (see [ScreenSOMGrid()] and [SelectSOMGrid()]), e.g.
#'   \code{list(jitter_sd = 0.3, max_nodes = 400, refits = 3)}. Ignored when
#'   dimensions or \code{map_from} are supplied.
#' @param map_from Optional exploratory `Pipeline_SOM_MClust` object. In a
#'   finalized fit, reuse only its resolved dimensions, topology, neighbourhood,
#'   and SOM seed. Clustering-variable names and order must match. Training
#'   preprocessing is deliberately recomputed on `data`; pass
#'   `ZScoreType = "ZScoreObj"` to reuse scaling explicitly.
#' @param som_topo SOM topology for \code{kohonen::somgrid()}, default
#'   \code{"hexagonal"}.
#' @param som_epochs Positive integer batch SOM training epochs, also used
#'   for grid screening and stability refits. Default 100.
#' @param som_neigh SOM neighbourhood function, default \code{"gaussian"}.
#' @param node_weighting Whether codebook nodes enter the mixture with equal
#'   map-space weight (\code{"uniform"}, default) or their training
#'   participant occupancy (\code{"occupancy"}). Occupancy weighting supports
#'   Models 1--3 only.
#' @param assignment Participant phenotype definition: `"hard"` maps the
#'   Mclust class of the best-matching SOM node; `"soft"` sums the Mclust
#'   node posterior probabilities across nearby codebook nodes and chooses the
#'   largest participant-level probability. Both outputs are always retained;
#'   the default `"soft"` controls the canonical `Cluster`.
#' @param soft_kernel Soft assignment kernel: inverse squared distance (default)
#'   or Gaussian. Inverse is fuzzy k-nearest neighbours over prototypes with
#'   fuzzifier 2 (Keller, Gray and Givens, 1985).
#' @param soft_bandwidth Positive Gaussian bandwidth. NULL uses the median
#'   training distance to the closest node, with 1 used for a zero median.
#'   The resolved bandwidth is stored in Specification$projection and reused
#'   by projection; refits resolve their own bandwidth.
#' @param soft_n_nodes Number of nearest codebook nodes used by soft projection.
#'   Their inverse-squared standardized-feature distances are normalized to
#'   proximity weights. Default `3`.
#' @param boundary_margin_threshold Minimum difference between the largest and
#'   second-largest soft class probabilities before a participant is flagged as
#'   a boundary participant. Default `0.20`.
#' @param boundary_entropy_threshold Maximum normalized soft-membership entropy
#'   before a participant is flagged as a boundary participant. Default
#'   `0.60`. A flag is set if either boundary criterion is met.
#' @param seed_som,seed_mixture Integer seeds for SOM and LPA steps
#'   (defaults 934521 and 93421).
#' @param Relabel Logical; if TRUE (default), aweSOM plots are relabeled
#'   using variable labels from the *original* \code{df} (via Hmisc or
#'   sjlabelled when available) by stripping the Z-score prefix.
#' @param ZScorePrefix Character prefix used for Z-score columns when
#'   \code{ZScoreType = "PreZScored"}. Default \code{"Z_"}.
#' @param ZScoreVars Optional character vector of Z-score column names to use
#'   when \code{ZScoreType = "PreZScored"}. If NULL, the function attempts to
#'   infer them from \code{variables} or by detecting columns starting with
#'   \code{ZScorePrefix}.
#' @param id_var Optional character scalar. If provided and present in \code{df},
#'   this column is carried into \code{ProbFit$individual} for convenience.
#' @param mixture_progress Logical; if TRUE, print short progress messages while
#'   fitting model/profile combinations.
#' @param mixture_em_itmax Integer; maximum number of EM iterations passed to
#'   \code{mclust::emControl()}. Use NULL to leave mclust defaults unchanged.
#' @param mixture_em_tol Numeric; EM convergence tolerance passed to
#'   \code{mclust::emControl()}. Use NULL to leave mclust defaults unchanged.
#' @param mixture_drop_zero_sd Logical; if TRUE, remove SOM code dimensions with
#'   near-zero standard deviation before LPA.
#' @param mixture_zero_sd_tol Numeric tolerance used when \code{mixture_drop_zero_sd = TRUE}.
#' @param mixture_bootstrap_n Bootstrap resamples tidyLPA runs for each candidate
#'   fit: a bootstrap likelihood-ratio test (BLRT) against k - 1 profiles and
#'   MclustBootstrap. Default \code{0} skips them, which removes most of the
#'   fitting time; \code{BLRTStatistic} and \code{BLRTPValue} are then
#'   \code{NA}. Set to, e.g., \code{100} to report the BLRT. Stability refits
#'   never run them.
#' @param mixture_timeout_seconds Optional timeout in seconds for individual LPA
#'   fits. Use NULL to disable timeouts.
#' @param skip_model_after_n_failures Optional integer; skip a model family
#'   after this many failures.
#' @param slow_fit_seconds Optional runtime threshold used to flag slow LPA
#'   fits in diagnostics.
#' @param candidate_min_nodes_per_cluster Optional minimum average SOM nodes per cluster.
#'   Default \code{10}. The grid is fixed before k is considered, so candidate
#'   profile counts with fewer nodes per cluster than this on the chosen map
#'   are still fitted but trigger one warning listing them. Use NULL to turn
#'   the warning off.
#' @param min_nodes_per_parameter Deprecated; only labels fit-table support. Positive number of codebook nodes per
#'   estimated free Gaussian-mixture parameter used to label candidates in the
#'   fit table (\code{SupportStatus}). Default \code{3}. It no longer sizes the
#'   map; candidates with fewer nodes than parameters are marked ineligible.
#' @param min_participants_per_node Number of complete participants per node
#'   used only for the \code{OccupancyCompliant} column of the fit table.
#'   Default \code{10}. It no longer sizes the map; sparse nodes are handled by
#'   the grid screen.
#' @param mixture_min_node_people Minimum number of training participants a SOM node
#'   must hold to enter the node-level mixture fit. Default \code{0} fits every
#'   node. Use \code{1} to leave empty nodes out, so unoccupied map regions do
#'   not count as observations in AIC/BIC. Nodes left out are classified
#'   afterwards by the fitted mixture model, so every node still carries a
#'   phenotype and posterior probabilities for projection.
#' @param min_cluster_prop Optional minimum share of training participants
#'   (0-1) in a candidate's smallest profile. Candidates below it are marked
#'   \code{Eligible = FALSE} in the fit table and are not recommended. Use it
#'   to keep only phenotypes large enough to describe and interpret.
#' @param criteria_weights Optional `Boundary` weights favour fewer boundary
#'   participants; Boundary is excluded from default weights. Weights for the composite rank index. \code{NULL}
#'   (default) weights the criteria equally. A named non-negative vector over
#'   \code{"AIC"}, \code{"BIC"}, \code{"Entropy"}, and
#'   \code{"Reproducibility"} sets weights directly. A square reciprocal
#'   pairwise-comparison matrix with those row and column names derives
#'   Analytic Hierarchy Process weights from its principal eigenvector and
#'   reports Saaty's consistency ratio (warning above 0.10).
#'   \code{"Reproducibility"} requires \code{stability_resamples > 0}.
#' @param far_distance_quantile Numeric value between 0 and 1 used to define
#'   high SOM-distance flags from the training distance distribution. Default
#'   is \code{0.95}.
#' @param clear_membership_threshold Numeric posterior probability threshold used to
#'   flag uncertain phenotype membership. Default is \code{0.70}.
#' @param stability_type `subsample` (default), `seed`, or both. Seed refits
#'   retain every participant, vary SOM and mixture seeds from seed_stability,
#'   and keep the grid fixed. With both types, AHP reproducibility is the mean
#'   of their two mean ARIs (subsampling uses the primary retained fraction).
#'   Stability tables include StabilityType; cluster_recovery_summary holds
#'   per-cluster MeanJaccard by type. Only mean ARI is summarized.
#' @param stability_resamples Number of 90% participant subsample refits used
#'   to assess reproducibility for every successful exploratory candidate.
#'   Subsamples are drawn without replacement and reuse the reference model's
#'   resolved SOM grid dimensions and the reference SOM seed. Each refit maps
#'   the complete reference cohort before scoring participant assignments.
#'   Defaults to \code{0} (disabled); use \code{30} for an exploratory
#'   stability screen.
#' @param stability_candidates Which candidates receive subsample stability:
#'   \code{"all"} (default) or \code{"best_per_k"}, the highest-ranked
#'   candidate (by the composite index without reproducibility, eligible
#'   first) at each number of profiles, so stability can inform the choice of
#'   k across the whole range. Candidates without stability keep \code{NA}
#'   reproducibility and are not recommended.
#' @param stability_top_n Optional integer. When supplied, stability is
#'   estimated only for the top \code{stability_top_n} candidates ranked by the
#'   composite index without reproducibility (eligible candidates first); the
#'   others keep \code{NA} stability and are not recommended. Use it to screen
#'   a large candidate grid, e.g. \code{stability_top_n = 5}.
#' @param seed_stability Integer seed for participant subsampling.
#' @param stability_fractions Unique participant-retention fractions for an
#'   opt-in deletion-stability curve. The default, \code{0.90}, preserves the
#'   ordinary 90% refit screen. Each fraction refits the complete pipeline and
#'   projects the full reference cohort before all-participant ARI is scored.
#' @param stability_progress Logical; if TRUE, print subsample progress messages.
#' @param stability_cores Number of workers for stability refits. `NULL` uses
#'   a conservative estimate based on CPU, refit count, and memory: the smaller
#'   of 50% physical RAM and 75% currently available RAM. Explicit values are
#'   retained up to the scheduler/CPU limit and warn if their estimated memory
#'   exceeds that budget. Each worker holds a refit in memory.
#'
#' @details
#' `ModelInfo_MClust$fit_table` uses the fit indices returned by tidyLPA's
#' mclust backend. AIC, AWE, BIC, CAIC, CLC, KIC, and SABIC are
#' likelihood/information criteria on tidyLPA's -2 log-likelihood scale, for
#' which lower values are preferred when comparing candidate fits estimated on
#' the same data. ICL follows mclust's convention instead (\code{-BIC} plus
#' twice the summed log maximum posterior, matching \code{mclust::icl()}), so
#' higher ICL values are preferred: ICL and BIC point in opposite directions,
#' and [PlotClusterFitReview()] labels each panel with its preferred
#' direction. ICL is reported for review only and is not part of the composite
#' index. Entropy summarizes classification separation (higher is better).
#' `MinProfileNodeN` and `MaxProfileNodeN` are integer SOM-node counts in the
#' smallest and largest profile; the corresponding `*Proportion` fields divide
#' those counts by the total number of SOM nodes. `BLRTStatistic` and
#' `BLRTPValue` compare a k-profile model with k - 1 profiles under the same
#' covariance model. A small p-value supports the added profile's statistical
#' fit, but does not establish clinical meaning or reproducibility; it is
#' unavailable for k = 1.
#' The raw \code{BLRTPValue} remains in the table. Its candidate-review panel
#' displays \code{-log10(BLRTPValue)} with a dashed \code{0.05} p-value
#' reference at \code{-log10(0.05)}.
#'
#' Lifecycle settings are intentionally strict: \code{k_range} and
#' \code{models} are exploratory-only, while \code{final_k} and
#' \code{final_model} are finalize-only. Supplying settings from both modes is
#' an error, even when an exploratory setting equals its default.
#'
#' Stability is assessed by refitting the full SOM/LPA pipeline on independent
#' participant subsamples drawn without replacement, then mapping the
#' complete reference cohort onto each refit. `StabilityARI_Mean` is the
#' primary reproducibility score: the mean adjusted Rand index across these
#' full-cohort participant partitions. It can be negative when agreement is
#' worse than chance.
#' `stability_fractions` can request an explicit deletion curve; the existing
#' 90% fit-table fields remain the primary review metrics, while
#' `Stability$deletion_curve` records every requested retention level.
#' Co-membership, split/merge, and profile-recovery diagnostics describe
#' whether instability is localized despite imperfect exact membership; they
#' never replace all-participant ARI.
#' `StabilityJaccard_Mean` is the mean label-matched, per-profile Jaccard
#' recovery and is retained as a descriptive diagnostic. `ReproducibilityScore`
#' equals `StabilityARI_Mean`; `StabilitySuccessRate` is reported separately.
#' `Reproducibility_scaled` is
#' the candidate-table min-max rescaling used only by the AHP index and is not
#' an independently interpretable reproducibility measure.
#'
#' The composite rank index (stored as \code{ahp_index}) is computed by:
#' \enumerate{
#'   \item Z-scaling AIC, BIC, Entropy, and (when stability is run)
#'         reproducibility across candidate solutions (AIC/BIC are negated so
#'         that higher scaled scores are preferred).
#'   \item Taking their weighted mean. Weights are equal unless
#'         \code{criteria_weights} supplies them, either directly or as an
#'         Analytic Hierarchy Process pairwise-comparison matrix.
#'   \item Recommending the highest-scoring candidate among those marked
#'         \code{Eligible} (see \code{min_cluster_prop}).
#' }
#' The chosen weights, their method, and any consistency ratio are returned in
#' \code{ModelInfo_MClust$AHP}.
#'
#' \code{MinProfileParticipantN} and \code{MinProfileParticipantProportion}
#' give the size of the smallest profile in training participants rather than
#' SOM nodes.
#'
#' Phenotypes are learned by Mclust on the SOM codebook. Hard assignment uses
#' the best-matching node's Mclust class. Soft assignment uses normalized
#' inverse-squared codebook-distance weights across the nearest nodes to
#' combine their Mclust posterior probabilities. `SoftEntropy` is the
#' normalized entropy of that participant-level probability vector
#' (entropy divided by \code{log(k)}) and `SoftMargin` is the largest minus
#' second-largest probability. `prob_*` columns describe a SOM node's Mclust
#' posterior, whereas `SoftProb_*` columns describe a participant's
#' proximity-weighted membership. Boundary flags are confidence diagnostics
#' only: they do not remove, relabel, or suppress a participant during
#' training, stability scoring, or projection.
#'
#' Changing `assignment` does not change SOM training or the fitted node mixture.
#' With the same seeds and selected model, node posteriors and likelihood-based
#' fit indices are identical. Both modes return a single participant class;
#' NodeN counts complete training participants on the closest node, and
#' Flag_OnePersonNode is NodeN == 1 (missing for incomplete rows). Both are
#' returned in ProbFit$individual and DataWithClusters. Probabilities measure
#' ambiguity between clusters; distance to node measures fit to the map.
#' soft mode takes the largest `SoftProb_*`. Hard and soft classes can therefore
#' agree even when their probability vectors differ, especially when nearby
#' nodes belong to the same phenotype. Exact soft-probability ties select the
#' first class deterministically. `Flag_NodeDisagreement` compares the canonical
#' participant class with the best-matching node's class; unmapped rows are `NA`.
#'
#' Each new fit uses one deterministic, seeded full-cohort batch SOM. Its
#' initialization codebook, preprocessing, grid, topology, neighbourhood,
#' epochs, codebooks, and node-level mixture are retained so projection is
#' completely frozen. Earlier fitted objects remain projectable because
#' projection uses their frozen map and mixture rather than refitting either
#' component.
#'
#' \code{SOM_Distance} is the Euclidean distance, in Z-score units, from a
#' participant to their best-matching unit. Its square is the sum of the
#' squared per-variable residuals, returned in \code{ModelInfo_SOM$SOMFit$residuals}
#' (\code{Resid_<variable>}). \code{ModelInfo_SOM$SOMFit$variable_contribution}
#' summarizes each variable's share of training misfit, and each participant's
#' largest contributor is recorded in \code{Top_Distance_Variable} and
#' \code{Top_Distance_Share}.
#'
#' LPA model/profile combinations are fit one at a time so that failed or
#' warning-producing solutions are captured in diagnostics instead of blocking
#' the entire pipeline. Successful fits are retained and failed fits are listed
#' in \code{ModelInfo_MClust$diagnostics}.
#'
#' @section Model review and refinement:
#' Start with \code{ModelInfo_MClust$fit_table}, the AHP recommendation, and
#' \code{fit_plot} to compare candidate model specifications and numbers of
#' profiles. AHP is a useful starting point, not an automatic final decision:
#' retain solutions whose fit, cluster size, and clinical interpretability are
#' all reasonable.
#'
#' The Circular and Line SOM widgets show how the analysis variables vary across
#' the map; the Cloud widget shows the observations within map cells. Use them
#' to assess whether the candidate solution represents coherent, interpretable
#' phenotype patterns. Inspect \code{SOMFit$node_occupancy} for empty or sparse
#' cells, and the SOM-distance plots for clusters with systematically poor map
#' representation. The posterior-probability plots in \code{ProbFit$plots}
#' identify uncertain node- or individual-level assignments.
#'
#' If these diagnostics suggest a weak solution, reconsider the variables, SOM
#' grid dimensions, model family, or candidate \code{k} range and refit in
#' exploratory mode. Once a solution is selected, refit it with
#' \code{method = "finalize"}, \code{final_k}, and \code{final_model} to create
#' the reusable model for projection.
#'
#' Set \code{stability_resamples} to a positive value to assess internal
#' reproducibility. Each 90% subsample refits the same candidate SOM/LPA solution
#' and projects the original participants back into it. Mean adjusted Rand
#' index summarizes whole-partition agreement; cluster-wise Jaccard recovery
#' identifies phenotypes that dissolve despite good overall agreement. This is
#' internal reproducibility, not independent-cohort validation.
#'
#' @section Reviewing the candidates:
#' Exploratory mode fits every combination of `k` and model family, so the
#' comparison is made once, on one table, rather than by refitting by hand.
#' `ModelInfo_MClust$fit_table` holds every candidate solution with its
#' information criteria, entropy, subsample reproducibility, and the combined
#' AHP index side by side.
#'
#' A candidate missing from that table is never a mystery:
#' `ModelInfo_MClust$diagnostics$lpa_fit_diagnostics` records each fit's
#' status, runtime, warnings, and errors, so model families that fail to
#' converge or exceed the fit timeout show up there.
#'
#' Once a solution is chosen, refit it on its own with `method = "finalize"` to
#' get the reusable, projectable model.
#'
#' @inheritSection cluster-stability-output Stability output
#' @return A list of class \code{"Pipeline_SOM_MClust"} with components:
#'   \itemize{
#'     \item \code{method}, \code{vars_used}, \code{ZScoreType},
#'           \code{ZScoreObject}, \code{ZScoreVars}, \code{ClusterVariableName}
#'     \item \code{DataWithClusters}: original \code{df} with canonical
#'           cluster, hard/soft assignment, soft membership, and boundary
#'           diagnostic columns appended
#'     \item \code{ProbFit$individual$SOM_Distance_LOO}: distance after removing
#'       each complete person from their own node mean, allowing the nearest
#'       other node instead. No map refit. ProjectionReference$loo_distances
#'       stores these distances with SOM_Node and Cluster for arbitrary quantiles.
#'       ProjectionReference$calibration stores q_star, training_rate, training_n,
#'       far_n, target_quantile, overall_cutoff and cluster_cutoffs for the
#'       combined overall OR cluster LOO rule.
#'       ProjectionReference$cluster_profiles stores mean training z-score
#'       vectors by assigned cluster for group profile comparisons.
#'     \item \code{ProbFit$boundary_summary}: BoundaryReason, N and Share
#'       among complete participants. \code{ProbFit$boundary_threshold_grid}
#'       gives margin, entropy and BoundaryPct for 15 threshold combinations.
#'       Defaults 0.20 / 0.60 are review conventions, not estimates.
#'     \item \code{ModelInfo_SOM$diagnostics}: map_quality (one-row tibble
#'       of quantization error, variance captured, topographic error and sparse
#'       node counts/shares), u_matrix (NodeID and U), and border_alignment
#'       (border/interior mean U, difference and one-sided 999-permutation p-value,
#'       conditional on permutations with both border and interior nodes;
#'       ValidPermutations reports that count).
#'     \item \code{ModelInfo_MClust$ahp_sensitivity}: weight schemes, normalized weights,
#'       top candidate and rank of the current recommendation. `BoundaryPct` in
#'       the fit table is the share of complete boundary participants (0--1).
#'     \item \code{fit_plot}: ggplot of AIC/BIC/Entropy/BLRT p-value vs k and model
#'           (plus reproducibility when subsample stability is enabled)
#'     \item \code{ModelInfo_SOM}: list with \code{som_model},
#'           \code{som_codes}, \code{som_grid}, \code{som_grid_info},
#'           \code{grid_screen} (the [ScreenSOMGrid()] result when the grid was
#'           chosen automatically, otherwise NULL), \code{training_variable_summary},
#'           \code{SOMFit} (distance diagnostics, baselines, and per-cluster
#'           flags), \code{plots} (aweSOM plots)
#'     \item \code{ModelInfo_MClust}: list with \code{lpa_models},
#'           \code{fit_table}, \code{AHP} information, and \code{diagnostics}
#'           for LPA warnings, failures, runtimes, and preprocessing
#'     \item \code{ModelInfo_MClust$Stability}: subsample replicate, cluster
#'           recovery, and summary tables when \code{stability_resamples > 0}
#'     \item \code{ProbFit}: list with \code{node} (node-level posterior
#'           probabilities), \code{individual} (full-length per-person mapping
#'           and probabilities), and probability plots
#'   }
#'
#' @references
#' Saaty TL. \emph{The Analytic Hierarchy Process}. McGraw-Hill, 1980.
#' Kohonen T. Self-organized formation of topologically correct feature maps.
#' \emph{Biological Cybernetics}. 1982;43:59-69. Scrucca L, Fop M, Murphy TB,
#' Raftery AE. mclust 5. \emph{J Stat Softw}. 2016;71(11):1-29.
#' Celeux G, Soromenho G. An entropy criterion for assessing the number of
#' clusters in a mixture model. \emph{J Classification}. 1996;13:195-212.
#' Cavanaugh JE. A large-sample model selection criterion based on Kullback's
#' symmetric divergence. \emph{Stat Methodol}. 1999;61:165-180. Sclove SL.
#' Application of model-selection criteria to some problems in multivariate
#' analysis. \emph{Psychometrika}. 1987;52:333-343. Nylund KL, Asparouhov T,
#' Muthen BO. Deciding on the number of classes in latent class analysis and
#' growth mixture modeling. \emph{Struct Equ Modeling}. 2007;14:535-569.
#' Akogul S, Erisoglu M. An approach for determining the number of clusters in
#' a model-based cluster analysis. \emph{Entropy}. 2017;19:452.
#'
#' @param df \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param id_col \strong{Deprecated} (since 19.15.0). Use \code{id_var} instead.
#' @param .CandidatePairs Internal. A data frame of \code{Model} and
#'   \code{Classes} pairs restricting the exploratory grid; used by stability
#'   refits so each resample fits only the shortlisted candidates.
#' @param .NodeClusterFn Internal. A function taking the SOM codebook matrix and
#'   returning a list with a \code{node_cluster} integer vector (one label per
#'   SOM node) and, optionally, \code{fit_table}, \code{ahp_best_row},
#'   \code{recommendation}, \code{best_fit_name}, and \code{fit_plot}. When
#'   supplied, the SOM codebook is clustered by that function and the
#'   latent-profile grid is not fitted. Used by
#'   [CreateClusterModel_SOM_HDBSCAN()]; not part of the user-facing API.
#' @examples
#' \donttest{
#' data("SimulatedPhenotypeData")
#' df_Training <- dplyr::filter(SimulatedPhenotypeData, .data$Cohort == "Training")
#' vars_Numeric <- paste0("Var", 1:12)
#'
#' # Example-only display helper (not exported by SciDataReportR): round the
#' # numeric columns and add a frozen header.
#' ShowFitTable <- function(x, height = "320px") {
#'   x <- dplyr::mutate(
#'     dplyr::as_tibble(x),
#'     dplyr::across(dplyr::where(is.numeric), \(v) round(v, 3))
#'   )
#'   htmltools::browsable(htmltools::HTML(as.character(
#'     FreezeTableHeader(x, height = height, full_width = TRUE)
#'   )))
#' }
#'
#' # Exploratory mode: every combination of k and model family. Model 2 is
#' # omitted here because it is the slowest to fit and, on this data, reaches
#' # the fit timeout at every k without contributing a candidate.
#' review <- CreateClusterModel_SOM_MClust(
#'   data = df_Training,
#'   variables = vars_Numeric,
#'   method = "exploratory",
#'   k_range = 2:5,
#'   models = c(1, 3),
#'   som_xdim = 5,
#'   som_ydim = 5,
#'   stability_resamples = 2,
#'   seed_stability = 20260805,
#'   Relabel = FALSE,
#'   candidate_min_nodes_per_cluster = NULL
#' )
#'
#' # Every candidate solution
#' ShowFitTable(review$ModelInfo_MClust$fit_table)
#'
#' # Status, runtime, warnings, and errors for every attempted fit
#' ShowFitTable(
#'   review$ModelInfo_MClust$diagnostics$lpa_fit_diagnostics,
#'   height = "260px"
#' )
#'
#' review$ModelInfo_MClust$AHP$recommendation
#' ShowFitTable(review$Stability$summary, height = "260px")
#' review$fit_plot
#' review$fit_plot
#'
#' # Refit the chosen solution to get the projectable model
#' model <- CreateClusterModel_SOM_MClust(
#'   data = df_Training, variables = vars_Numeric, method = "finalize",
#'   final_k = 4, final_model = 1, som_xdim = 5, som_ydim = 5,
#'   stability_resamples = 2, Relabel = FALSE,
#'   candidate_min_nodes_per_cluster = NULL
#' )
#' ShowFitTable(model$ModelInfo_MClust$fit_table, height = "160px")
#' }
#' @param seed_lpa Deprecated alias for `seed_mixture`.
#' @param stability_seed Deprecated alias for `seed_stability`.
#' @param lpa_em_itmax Deprecated alias for `mixture_em_itmax`.
#' @param lpa_em_tol Deprecated alias for `mixture_em_tol`.
#' @param lpa_timeout_seconds Deprecated alias for `mixture_timeout_seconds`.
#' @param lpa_progress Deprecated alias for `mixture_progress`.
#' @param lpa_bootstrap_n Deprecated alias for `mixture_bootstrap_n`.
#' @param lpa_drop_zero_sd Deprecated alias for `mixture_drop_zero_sd`.
#' @param lpa_zero_sd_tol Deprecated alias for `mixture_zero_sd_tol`.
#' @param lpa_min_node_n Deprecated alias for `mixture_min_node_people`.
#' @param min_nodes_per_cluster Deprecated alias for `candidate_min_nodes_per_cluster`.
#' @param high_dist_quantile Deprecated alias for `far_distance_quantile`.
#' @param low_prob_threshold Deprecated alias for `clear_membership_threshold`.
#' @export
CreateClusterModel_SOM_MClust <- function(data,
    variables = NULL,
    method = c("exploratory", "finalize", "explore"),
    k_range = 2:10,
    models = c(1, 2, 3, 6),
    final_k = NULL,
    final_model = NULL,
    ClusterVariableName = "Cluster",
    ZScoreType = c("Center and Scale", "Center Only", "Scale Only", "ZScoreObj", "PreZScored"),
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
    boundary_margin_threshold = 0.20,
    boundary_entropy_threshold = 0.60,
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
    stability_fractions = 0.90,
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
    low_prob_threshold = lifecycle::deprecated()) {
  if (lifecycle::is_present(seed_lpa)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(seed_lpa)", "CreateClusterModel_SOM_MClust(seed_mixture)")
    seed_mixture <- seed_lpa
  }
  if (lifecycle::is_present(stability_seed)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(stability_seed)", "CreateClusterModel_SOM_MClust(seed_stability)")
    seed_stability <- stability_seed
  }
  if (lifecycle::is_present(lpa_em_itmax)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_em_itmax)", "CreateClusterModel_SOM_MClust(mixture_em_itmax)")
    mixture_em_itmax <- lpa_em_itmax
  }
  if (lifecycle::is_present(lpa_em_tol)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_em_tol)", "CreateClusterModel_SOM_MClust(mixture_em_tol)")
    mixture_em_tol <- lpa_em_tol
  }
  if (lifecycle::is_present(lpa_timeout_seconds)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_timeout_seconds)", "CreateClusterModel_SOM_MClust(mixture_timeout_seconds)")
    mixture_timeout_seconds <- lpa_timeout_seconds
  }
  if (lifecycle::is_present(lpa_progress)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_progress)", "CreateClusterModel_SOM_MClust(mixture_progress)")
    mixture_progress <- lpa_progress
  }
  if (lifecycle::is_present(lpa_bootstrap_n)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_bootstrap_n)", "CreateClusterModel_SOM_MClust(mixture_bootstrap_n)")
    mixture_bootstrap_n <- lpa_bootstrap_n
  }
  if (lifecycle::is_present(lpa_drop_zero_sd)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_drop_zero_sd)", "CreateClusterModel_SOM_MClust(mixture_drop_zero_sd)")
    mixture_drop_zero_sd <- lpa_drop_zero_sd
  }
  if (lifecycle::is_present(lpa_zero_sd_tol)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_zero_sd_tol)", "CreateClusterModel_SOM_MClust(mixture_zero_sd_tol)")
    mixture_zero_sd_tol <- lpa_zero_sd_tol
  }
  if (lifecycle::is_present(lpa_min_node_n)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(lpa_min_node_n)", "CreateClusterModel_SOM_MClust(mixture_min_node_people)")
    mixture_min_node_people <- lpa_min_node_n
  }
  if (lifecycle::is_present(min_nodes_per_cluster)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(min_nodes_per_cluster)", "CreateClusterModel_SOM_MClust(candidate_min_nodes_per_cluster)")
    candidate_min_nodes_per_cluster <- min_nodes_per_cluster
  }
  if (lifecycle::is_present(high_dist_quantile)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(high_dist_quantile)", "CreateClusterModel_SOM_MClust(far_distance_quantile)")
    far_distance_quantile <- high_dist_quantile
  }
  if (lifecycle::is_present(low_prob_threshold)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(low_prob_threshold)", "CreateClusterModel_SOM_MClust(clear_membership_threshold)")
    clear_membership_threshold <- low_prob_threshold
  }

  if (!missing(min_nodes_per_parameter)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(min_nodes_per_parameter)",
      details = "This argument only labels map support; eligibility always requires at least one node per mixture parameter.")
  }
  if (!missing(min_participants_per_node)) {
    lifecycle::deprecate_warn("21.18.0", "CreateClusterModel_SOM_MClust(min_participants_per_node)",
      details = "This argument only labels map support and no longer sizes the grid.")
  }
  if (!is.null(.CandidatePairs)) .CandidatePairs$Model <- .ResolveMclustModels(.CandidatePairs$Model, ".CandidatePairs$Model")
  supplied <- list(k_range = !missing(k_range), models = !missing(models),
    final_k = !missing(final_k), final_model = !missing(final_model))
  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(df)) {
    lifecycle::deprecate_warn("19.15.0", "CreateClusterModel_SOM_MClust(df)", "CreateClusterModel_SOM_MClust(data)")
    data <- df
  }
  if (!missing(data)) df <- data
  if (lifecycle::is_present(id_col)) {
    lifecycle::deprecate_warn("19.15.0", "CreateClusterModel_SOM_MClust(id_col)", "CreateClusterModel_SOM_MClust(id_var)")
    id_var <- id_col
  }
  id_col <- id_var


  method     <- .ClusterMethod(method)
  .ValidateClusterLifecycle(method, supplied[c("k_range", "models")],
    supplied[c("final_k", "final_model")])
  models <- .ResolveMclustModels(models)
  if (method == "finalize") {
    final_model <- .ResolveMclustModels(final_model, "final_model")
  }
  ZScoreType <- match.arg(ZScoreType)
  if (!is.numeric(som_epochs) || length(som_epochs) != 1L ||
      !is.finite(som_epochs) || som_epochs < 1 || som_epochs != as.integer(som_epochs)) {
    stop("som_epochs must be a positive integer.")
  }
  node_weighting <- match.arg(node_weighting)
  assignment <- match.arg(assignment)
  soft_kernel <- match.arg(soft_kernel)
  soft_bandwidth_requested <- soft_bandwidth
  if (!is.null(soft_bandwidth) && (!is.numeric(soft_bandwidth) || length(soft_bandwidth) != 1L ||
      !is.finite(soft_bandwidth) || soft_bandwidth <= 0)) stop("soft_bandwidth must be NULL or one positive finite number.")
  if (!is.numeric(soft_n_nodes) || length(soft_n_nodes) != 1L || is.na(soft_n_nodes) || soft_n_nodes < 1L || soft_n_nodes != as.integer(soft_n_nodes)) stop("soft_n_nodes must be a positive integer.")
  if (!is.numeric(boundary_margin_threshold) || length(boundary_margin_threshold) != 1L || is.na(boundary_margin_threshold) || boundary_margin_threshold < 0 || boundary_margin_threshold > 1) stop("boundary_margin_threshold must be between 0 and 1.")
  if (!is.numeric(boundary_entropy_threshold) || length(boundary_entropy_threshold) != 1L || is.na(boundary_entropy_threshold) || boundary_entropy_threshold < 0 || boundary_entropy_threshold > 1) stop("boundary_entropy_threshold must be between 0 and 1.")
  if (!requireNamespace("SciDataReportR", quietly = TRUE)) {
    stop("SciDataReportR must be installed.")
  }
  if (!requireNamespace("kohonen", quietly = TRUE)) {
    stop("Package 'kohonen' is required.")
  }
  if (!requireNamespace("aweSOM", quietly = TRUE)) {
    stop("Package 'aweSOM' is required.")
  }
  if (!suppressPackageStartupMessages(
    requireNamespace("tidyLPA", quietly = TRUE))) {
    stop("Package 'tidyLPA' is required.")
  }
  if ((!is.null(mixture_em_itmax) || !is.null(mixture_em_tol)) &&
      !requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required when mixture_em_itmax or mixture_em_tol is supplied.")
  }
  if (!is.null(mixture_timeout_seconds) &&
      !requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required when mixture_timeout_seconds is supplied. Install it or set mixture_timeout_seconds = NULL.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Packages 'dplyr', 'tidyr', and 'ggplot2' are required.")
  }

  if (!is.logical(mixture_progress) || length(mixture_progress) != 1) {
    stop("mixture_progress must be TRUE or FALSE.")
  }
  if (!is.numeric(stability_resamples) || length(stability_resamples) != 1 ||
      is.na(stability_resamples) || stability_resamples < 0 ||
      stability_resamples != as.integer(stability_resamples)) {
    stop("stability_resamples must be a single non-negative integer.")
  }
  if (!is.numeric(seed_stability) || length(seed_stability) != 1 || is.na(seed_stability)) {
    stop("seed_stability must be a single numeric value.")
  }
  stability_type <- if (missing(stability_type)) "subsample" else
    match.arg(stability_type, c("subsample", "seed"), several.ok = TRUE)
  stability_fractions <- .ValidateStabilityFractions(stability_fractions)
  primary_stability_fraction <- if (identical(stability_type, "seed")) 1 else if (0.90 %in% stability_fractions) 0.90 else stability_fractions[[1]]
  if (!is.logical(stability_progress) || length(stability_progress) != 1) {
    stop("stability_progress must be TRUE or FALSE.")
  }
  if (!is.null(stability_top_n) &&
      (!is.numeric(stability_top_n) || length(stability_top_n) != 1 ||
       is.na(stability_top_n) || stability_top_n < 1)) {
    stop("stability_top_n must be NULL or a single positive integer.")
  }
  stability_candidates <- match.arg(stability_candidates)
  .ValidateClusterStabilityCores(stability_cores)
  if (!is.numeric(mixture_bootstrap_n) || length(mixture_bootstrap_n) != 1 ||
      is.na(mixture_bootstrap_n) || mixture_bootstrap_n < 0) {
    stop("mixture_bootstrap_n must be a single non-negative integer.")
  }
  if (!is.logical(mixture_drop_zero_sd) || length(mixture_drop_zero_sd) != 1) {
    stop("mixture_drop_zero_sd must be TRUE or FALSE.")
  }
  if (!is.numeric(mixture_zero_sd_tol) || length(mixture_zero_sd_tol) != 1 || mixture_zero_sd_tol < 0) {
    stop("mixture_zero_sd_tol must be a single non-negative numeric value.")
  }
  if (!is.null(mixture_em_itmax) &&
      (!is.numeric(mixture_em_itmax) || length(mixture_em_itmax) != 1 || mixture_em_itmax < 1)) {
    stop("mixture_em_itmax must be NULL or a single positive integer.")
  }
  if (!is.null(mixture_em_tol) &&
      (!is.numeric(mixture_em_tol) || length(mixture_em_tol) != 1 || mixture_em_tol <= 0)) {
    stop("mixture_em_tol must be NULL or a single positive numeric value.")
  }
  if (!is.null(mixture_timeout_seconds) &&
      (!is.numeric(mixture_timeout_seconds) ||
       length(mixture_timeout_seconds) != 1 ||
       mixture_timeout_seconds <= 0)) {
    stop("mixture_timeout_seconds must be NULL or a single positive numeric value.")
  }
  if (!is.null(skip_model_after_n_failures) &&
      (!is.numeric(skip_model_after_n_failures) ||
       length(skip_model_after_n_failures) != 1 ||
       skip_model_after_n_failures < 1)) {
    stop("skip_model_after_n_failures must be NULL or a single positive integer.")
  }
  if (!is.null(slow_fit_seconds) &&
      (!is.numeric(slow_fit_seconds) || length(slow_fit_seconds) != 1 || slow_fit_seconds <= 0)) {
    stop("slow_fit_seconds must be NULL or a single positive numeric value.")
  }
  if (!is.null(candidate_min_nodes_per_cluster) &&
      (!is.numeric(candidate_min_nodes_per_cluster) ||
       length(candidate_min_nodes_per_cluster) != 1 ||
       candidate_min_nodes_per_cluster <= 0)) {
    stop("candidate_min_nodes_per_cluster must be NULL or a single positive numeric value.")
  }
  # The screen must train the same map as the pipeline, so seed and map design come from the main arguments
  if (is.list(grid_screen_args) && is.null(grid_screen_args$cores) && !is.null(stability_cores)) {
    grid_screen_args$cores <- stability_cores
  }
  allowed_screen_args <- setdiff(names(formals(.ScreenSOMGridMatrix)),
    c("zmat", "seed_som", "som_topo", "som_neigh", "som_epochs"))
  if (!is.list(grid_screen_args) || (length(grid_screen_args) && is.null(names(grid_screen_args))) ||
      length(setdiff(names(grid_screen_args), allowed_screen_args))) {
    stop("grid_screen_args must be a named list of grid-screen settings: ",
      paste(allowed_screen_args, collapse = ", "), ".")
  }
  if (!is.numeric(far_distance_quantile) || length(far_distance_quantile) != 1 ||
      is.na(far_distance_quantile) || far_distance_quantile <= 0 || far_distance_quantile >= 1) {
    stop("far_distance_quantile must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(clear_membership_threshold) || length(clear_membership_threshold) != 1 ||
      is.na(clear_membership_threshold) || clear_membership_threshold < 0 || clear_membership_threshold > 1) {
    stop("clear_membership_threshold must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(mixture_min_node_people) || length(mixture_min_node_people) != 1 ||
      is.na(mixture_min_node_people) || mixture_min_node_people < 0) {
    stop("mixture_min_node_people must be a single non-negative number.")
  }
  if (!is.numeric(min_nodes_per_parameter) || length(min_nodes_per_parameter) != 1L ||
      is.na(min_nodes_per_parameter) || !is.finite(min_nodes_per_parameter) ||
      min_nodes_per_parameter <= 0) {
    stop("min_nodes_per_parameter must be a single positive finite number.")
  }
  if (!is.numeric(min_participants_per_node) || length(min_participants_per_node) != 1L ||
      is.na(min_participants_per_node) || !is.finite(min_participants_per_node) ||
      min_participants_per_node <= 0) {
    stop("min_participants_per_node must be a single positive finite number.")
  }
  if (node_weighting == "occupancy" && mixture_min_node_people < 1L) {
    mixture_min_node_people <- 1L
  }
  if (!is.null(min_cluster_prop) &&
      (!is.numeric(min_cluster_prop) || length(min_cluster_prop) != 1 ||
       is.na(min_cluster_prop) || min_cluster_prop < 0 || min_cluster_prop >= 1)) {
    stop("min_cluster_prop must be NULL or a single numeric value in [0, 1).")
  }
  composite_weights <- .ResolveCompositeWeights(criteria_weights)
  if (method == "finalize" && (is.null(final_k) || is.null(final_model))) {
    stop("For method = 'finalize', final_k and final_model must be supplied.")
  }
  # Stability is estimated by refitting candidates, which only the exploratory
  # path does. Rather than silently dropping the request, refit the single
  # finalized candidate through that path and relabel the result. Recursing on
  # the evaluated arguments keeps this correct when the caller passed locals or
  # unevaluated promises.
  if (method == "finalize" && stability_resamples > 0L) {
    ## Resolve the map-design portion of map_from before recursive stability
    ## fitting. The refits intentionally receive dimensions rather than the
    ## source object, so they relearn preprocessing and the full pipeline.
    if (!is.null(map_from)) {
      if (!inherits(map_from, "Pipeline_SOM_MClust")) {
        stop("map_from must be a Pipeline_SOM_MClust object.")
      }
      if (!is.null(variables) && !identical(variables, map_from$vars_used)) {
        stop("map_from must use the same clustering variables in the same order as the finalization data.")
      }
      map_info <- map_from$ModelInfo_SOM$som_grid_info
      source_grid <- map_from$ModelInfo_SOM$som_grid
      if (is.null(map_info) || is.null(source_grid)) {
        stop("map_from does not contain a reusable SOM grid specification.")
      }
      som_xdim <- map_info$som_xdim_used
      som_ydim <- map_info$som_ydim_used
      som_topo <- source_grid$topo
      som_neigh <- as.character(source_grid$neighbourhood.fct)
      seed_som <- map_from$ModelInfo_SOM$training$seed
      map_from <- NULL
    }
    stable_model <- Recall(
      data = data, variables = variables, method = "exploratory",
      k_range = final_k, models = final_model,
      ClusterVariableName = ClusterVariableName,
      ZScoreType = ZScoreType,
      ZScoreObject = ZScoreObject, ZScorePrefix = ZScorePrefix,
      ZScoreVars = ZScoreVars, id_var = id_var,
      som_xdim = som_xdim, som_ydim = som_ydim, map_from = NULL,
      grid_screen_args = grid_screen_args, som_topo = som_topo,
      som_neigh = som_neigh, som_epochs = som_epochs,
      node_weighting = node_weighting, assignment = assignment,
      soft_n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth,
      boundary_margin_threshold = boundary_margin_threshold,
      boundary_entropy_threshold = boundary_entropy_threshold,
      seed_som = seed_som, seed_mixture = seed_mixture,
      Relabel = Relabel, mixture_progress = mixture_progress,
      mixture_em_itmax = mixture_em_itmax, mixture_em_tol = mixture_em_tol,
      mixture_bootstrap_n = mixture_bootstrap_n,
      mixture_timeout_seconds = mixture_timeout_seconds,
      mixture_drop_zero_sd = mixture_drop_zero_sd,
      mixture_zero_sd_tol = mixture_zero_sd_tol,
      skip_model_after_n_failures = skip_model_after_n_failures,
      slow_fit_seconds = slow_fit_seconds,
      candidate_min_nodes_per_cluster = candidate_min_nodes_per_cluster,
      min_nodes_per_parameter = min_nodes_per_parameter,
      min_participants_per_node = min_participants_per_node,
      mixture_min_node_people = mixture_min_node_people,
      min_cluster_prop = min_cluster_prop,
      criteria_weights = criteria_weights,
      far_distance_quantile = far_distance_quantile,
      clear_membership_threshold = clear_membership_threshold,
      stability_resamples = stability_resamples, stability_type = stability_type,
      stability_candidates = stability_candidates,
      stability_top_n = stability_top_n,
      seed_stability = seed_stability,
      stability_fractions = stability_fractions,
      stability_progress = stability_progress,
      stability_cores = stability_cores)
    # The exploratory path scores its one candidate with AHP and records no
    # selection; restate it as the user-specified finalized model.
    stable_model$method <- "finalize"
    fit_table <- stable_model$ModelInfo_MClust$fit_table
    ahp_cols <- intersect(c("AIC_scaled", "BIC_scaled", "Entropy_scaled",
      "Reproducibility_scaled", "ahp_index"), names(fit_table))
    fit_table[ahp_cols] <- NA_real_
    recommendation_txt <- paste0(
      "User-specified Model ", final_model, " with k = ", final_k, " profiles.")
    stable_model$ModelInfo_MClust$fit_table <- fit_table
    stable_model$ModelInfo_MClust$AHP <- list(
      ahp_index = fit_table$ahp_index,
      ahp_best_row = fit_table[1, ],
      recommendation = recommendation_txt)
    stable_model$ModelInfo <- stable_model$ModelInfo_MClust
    stable_model$Specification$candidate_grid <- fit_table
    stable_model$Specification$selected <- list(k = final_k, model = final_model)
    return(stable_model)
  }

  # Stable row id ----------------------------------------------------------

  df_scidr <- .AddClusterRowID(df, "df")

  if (!is.null(id_col)) {
    if (!id_col %in% names(df_scidr)) {
      stop("id_col '", id_col, "' was provided but is not in df.")
    }
  }

  # Variables --------------------------------------------------------------

  ZScoreVars_used <- NULL
  vars_used       <- NULL

  if (ZScoreType == "PreZScored") {

    if (!is.null(ZScoreVars)) {
      ZScoreVars_used <- ZScoreVars
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)

    } else if (!is.null(variables)) {

      if (all(variables %in% names(df_scidr)) &&
          all(grepl(paste0("^", ZScorePrefix), variables))) {
        ZScoreVars_used <- variables
        vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
      } else {
        ZScoreVars_used <- paste0(ZScorePrefix, variables)
        vars_used <- variables
      }

    } else {

      ZScoreVars_used <- grep(paste0("^", ZScorePrefix), names(df_scidr), value = TRUE)
      if (length(ZScoreVars_used) == 0) {
        stop("ZScoreType = 'PreZScored' requires ZScoreVars, variables, or columns starting with '",
             ZScorePrefix, "'.")
      }
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
    }

    missing_z <- setdiff(ZScoreVars_used, names(df_scidr))
    if (length(missing_z) > 0) {
      stop("The following pre-zscored variables are not in df: ",
           paste(missing_z, collapse = ", "))
    }

    z_tmp <- df_scidr[, ZScoreVars_used, drop = FALSE]
    is_num_z <- vapply(z_tmp, is.numeric, logical(1))
    if (!all(is_num_z)) {
      warning("Dropping non-numeric pre-zscored variables: ",
              paste(names(z_tmp)[!is_num_z], collapse = ", "))
      ZScoreVars_used <- ZScoreVars_used[is_num_z]
      z_tmp <- z_tmp[, is_num_z, drop = FALSE]
      vars_used <- sub(paste0("^", ZScorePrefix), "", ZScoreVars_used)
    }

    if (length(ZScoreVars_used) == 0) {
      stop("No numeric pre-zscored variables available for SOM / clustering.")
    }

  } else {

    if (is.null(variables)) {
      variables <- SciDataReportR::getNumVars(df, Ordinal = FALSE)
    }

    missing_vars <- setdiff(variables, names(df))
    if (length(missing_vars) > 0) {
      stop("The following variables are not in df: ",
           paste(missing_vars, collapse = ", "))
    }

    is_num <- vapply(df[variables], is.numeric, logical(1))
    if (!all(is_num)) {
      warning("Dropping non-numeric variables: ",
              paste(variables[!is_num], collapse = ", "))
      variables <- variables[is_num]
    }
    if (length(variables) == 0) {
      stop("No numeric variables available for SOM / clustering.")
    }

    vars_used <- variables
  }


  if (!is.null(map_from)) {
    if (method != "finalize" || !inherits(map_from, "Pipeline_SOM_MClust")) {
      stop("map_from must be an exploratory Pipeline_SOM_MClust object used with method = 'finalize'.")
    }
    if (!identical(vars_used, map_from$vars_used)) {
      stop("map_from must use the same clustering variables in the same order as the finalization data.")
    }
    map_info <- map_from$ModelInfo_SOM$som_grid_info
    source_grid <- map_from$ModelInfo_SOM$som_grid
    if (is.null(map_info) || is.null(source_grid)) {
      stop("map_from does not contain a reusable SOM grid specification.")
    }
    if ((!is.null(som_xdim) && som_xdim != map_info$som_xdim_used) ||
        (!is.null(som_ydim) && som_ydim != map_info$som_ydim_used)) {
      stop("som_xdim and som_ydim must match map_from, or be omitted.")
    }
    som_xdim <- map_info$som_xdim_used
    som_ydim <- map_info$som_ydim_used
    som_topo <- source_grid$topo
    som_neigh <- as.character(source_grid$neighbourhood.fct)
    seed_som <- map_from$ModelInfo_SOM$training$seed
  } else if (method == "finalize" && is.null(som_xdim) && is.null(som_ydim)) {
    message("No grid size or map_from supplied: the grid is chosen by the map-only screen (ScreenSOMGrid()). ",
      "It does not depend on final_k or final_model, so the same data and seed give the same grid as exploration; ",
      "supply map_from to reuse the exploratory map design exactly.")
  }

  # Training variable summary ---------------------------------------------

  training_variable_summary <- dplyr::tibble(
    Variable = vars_used
  ) %>%
    dplyr::mutate(
      present_in_training = Variable %in% names(df_scidr),
      is_numeric = vapply(Variable, function(v) {
        v %in% names(df_scidr) && is.numeric(df_scidr[[v]])
      }, logical(1)),
      n_training = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_integer_ else length(df_scidr[[v]])
      }, integer(1)),
      n_missing = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_integer_ else sum(is.na(df_scidr[[v]]))
      }, integer(1)),
      prop_missing = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr)) NA_real_ else mean(is.na(df_scidr[[v]]))
      }, numeric(1)),
      mean = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]])) NA_real_ else mean(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      sd = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]])) NA_real_ else stats::sd(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      min = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]]) || all(is.na(df_scidr[[v]]))) NA_real_ else min(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1)),
      max = vapply(Variable, function(v) {
        if (!v %in% names(df_scidr) || !is.numeric(df_scidr[[v]]) || all(is.na(df_scidr[[v]]))) NA_real_ else max(df_scidr[[v]], na.rm = TRUE)
      }, numeric(1))
    )

  # Z-scores ---------------------------------------------------------------

  if (ZScoreType == "PreZScored") {

    z_df <- df_scidr[, ZScoreVars_used, drop = FALSE]
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows in pre-zscored columns; cannot fit SOM.")
    }

    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreObject_used <- NULL

  } else if (ZScoreType == "ZScoreObj") {

    if (is.null(ZScoreObject) || !"ZScoreObj" %in% class(ZScoreObject)) {
      stop("ZScoreType = 'ZScoreObj' requires a valid ZScoreObj.")
    }

    z_res <- SciDataReportR::ProjectZScore(
      data                 = df_scidr,
      variables          = vars_used,
      parameters         = ZScoreObject,
      ParameterInputType = "ZScoreObj",
      names_prefix       = "Z_"
    )
    ZScoreObject_used <- z_res

    z_df <- ZScoreObject_used$ZScores
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows after Z-score step; cannot fit SOM.")
    }

    z_df <- z_df[, names(z_df), drop = FALSE]
    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreVars_used <- names(z_df)

  } else {

    center_flag <- ZScoreType %in% c("Center and Scale", "Center Only")
    scale_flag  <- ZScoreType %in% c("Center and Scale", "Scale Only")

    ZScoreObject_used <- SciDataReportR::CreateZScoreObject(
      data           = df_scidr,
      variables    = vars_used,
      names_prefix = "Z_",
      center       = center_flag,
      scale        = scale_flag
    )

    z_df <- ZScoreObject_used$ZScores
    complete_rows <- stats::complete.cases(z_df)
    if (!any(complete_rows)) {
      stop("No complete rows after Z-score step; cannot fit SOM.")
    }

    zmat <- as.matrix(z_df[complete_rows, , drop = FALSE])
    ZScoreVars_used <- names(z_df)
  }

  # SOM fitting ------------------------------------------------------------

  n_complete <- nrow(zmat)
  som_grid_initial_xdim <- som_xdim
  som_grid_initial_ydim <- som_ydim
  grid_design_k <- if (method == "finalize") as.integer(final_k) else max(as.integer(k_range))
  grid_nodes_per_profile <- NA_integer_   # retired k-based sizing; kept as a reported field
  grid_design_policy <- "explicit_dimensions"
  som_training_epochs <- as.integer(som_epochs)

  FitSOMForGrid <- function(xdim, ydim) {
    som_grid <- kohonen::somgrid(
      xdim = xdim, ydim = ydim, topo = som_topo,
      neighbourhood.fct = som_neigh)
    # The same seed gives each sizing iteration a deterministic initialization
    # appropriate to its own resolved grid.
    set.seed(seed_som)
    init_codes <- zmat[sample.int(nrow(zmat), nrow(som_grid$pts),
      replace = nrow(zmat) < nrow(som_grid$pts)), , drop = FALSE]
    som_model <- kohonen::supersom(list(zmat), grid = som_grid,
      init = list(init_codes), mode = "batch", rlen = som_training_epochs,
      keep.data = TRUE, cores = 1L)
    codes_raw <- som_model$codes
    som_codes <- if (is.list(codes_raw)) as.matrix(codes_raw[[1]]) else as.matrix(codes_raw)
    som_node_occupancy <- dplyr::tibble(NodeID = seq_len(nrow(som_codes))) %>%
      dplyr::left_join(
        as.data.frame(table(som_model$unit.classif)) %>%
          dplyr::transmute(NodeID = as.integer(as.character(Var1)),
            n = as.integer(Freq)),
        by = "NodeID") %>%
      dplyr::mutate(n = dplyr::if_else(is.na(n), 0L, n))
    lpa_nodes <- if (is.null(.NodeClusterFn)) {
      which(som_node_occupancy$n >= mixture_min_node_people)
    } else seq_len(nrow(som_codes))
    if (length(lpa_nodes) < 2L) {
      stop("Fewer than two SOM nodes hold at least mixture_min_node_people = ",
        mixture_min_node_people, " participants; lower mixture_min_node_people.")
    }
    mixture_codes <- som_codes[lpa_nodes, , drop = FALSE]
    dropped_lpa_vars <- character(0)
    if (is.null(.NodeClusterFn) && mixture_drop_zero_sd) {
      lpa_sd <- apply(mixture_codes, 2, stats::sd, na.rm = TRUE)
      keep_lpa_vars <- is.finite(lpa_sd) & lpa_sd > mixture_zero_sd_tol
      dropped_lpa_vars <- colnames(mixture_codes)[!keep_lpa_vars]
      mixture_codes <- mixture_codes[, keep_lpa_vars, drop = FALSE]
    }
    if (!ncol(mixture_codes)) {
      stop("No SOM code dimensions with non-zero variance are available for LPA.")
    }
    list(som_grid = som_grid, init_codes = init_codes, som_model = som_model,
      som_codes = som_codes, som_node_occupancy = som_node_occupancy,
      lpa_nodes = lpa_nodes, mixture_codes = mixture_codes,
      dropped_lpa_vars = dropped_lpa_vars)
  }

  requested_parameters <- function(p) {
    if (!is.null(.NodeClusterFn)) return(NA_integer_)
    parameter_grid <- if (method == "finalize") {
      data.frame(Model = as.integer(final_model), Classes = as.integer(final_k))
    } else {
      expand.grid(Model = as.integer(models), Classes = as.integer(k_range),
        KEEP.OUT.ATTRS = FALSE)
    }
    max(vapply(seq_len(nrow(parameter_grid)), function(i) {
      .SOMMclustParameterCount(parameter_grid$Model[[i]], p,
        parameter_grid$Classes[[i]])
    }, integer(1)))
  }

  ## Grid size is chosen from the map alone, before any mixture is fitted: no k and no covariance
  ## model enters the choice (ScreenSOMGrid()). The screen's first seed refit is the same map as
  ## FitSOMForGrid() trains below with seed_som.
  sizing_iterations <- 1L
  grid_screen <- NULL
  occupancy_ceiling_side <- NA_integer_
  automatic_dimensions <- is.null(som_xdim) || is.null(som_ydim)
  binding_constraint <- if (is.null(map_from)) "explicit_dimensions" else "map_from"
  if (automatic_dimensions) {
    screen_args <- utils::modifyList(list(seed_som = seed_som, som_topo = som_topo,
      som_neigh = som_neigh, som_epochs = som_training_epochs), grid_screen_args)
    grid_screen <- do.call(.ScreenSOMGridMatrix, c(list(zmat = zmat), screen_args))
    som_xdim <- grid_screen$selection$side
    som_ydim <- grid_screen$selection$side
    grid_design_policy <- "map_only_screen"
    binding_constraint <- if (som_xdim == grid_screen$selection$target) {
      "screen_target"
    } else if (som_xdim < grid_screen$selection$target) {
      paste0("screen_", grid_screen$selection$ceiling_binding, "_ceiling")
    } else {
      "screen_variable_floor"
    }
    message(sprintf("SOM grid chosen automatically: %d x %d. %s ", som_xdim, som_ydim,
      grid_screen$selection$reason),
      "Review with plot(<model>$ModelInfo_SOM$grid_screen); supply som_xdim and som_ydim to override.")
  }
  som_fit_for_grid <- FitSOMForGrid(som_xdim, som_ydim)

  som_grid <- som_fit_for_grid$som_grid
  init_codes <- som_fit_for_grid$init_codes
  som_model <- som_fit_for_grid$som_model
  som_codes <- som_fit_for_grid$som_codes
  som_node_occupancy <- som_fit_for_grid$som_node_occupancy
  lpa_nodes <- som_fit_for_grid$lpa_nodes
  X <- som_fit_for_grid$mixture_codes
  lpa_node_weights <- som_node_occupancy$n[lpa_nodes]
  dropped_lpa_vars <- som_fit_for_grid$dropped_lpa_vars
  if (length(dropped_lpa_vars)) {
    warning("Dropping near-zero SD SOM code dimensions before LPA: ",
      paste(dropped_lpa_vars, collapse = ", "))
  }
  mixture_dimension <- ncol(X)
  # Parameter support is reported for the requested candidates; it no longer sizes the map.
  max_mixture_parameters <- requested_parameters(mixture_dimension)
  parameter_floor_nodes <- if (is.na(max_mixture_parameters)) NA_real_ else
    min_nodes_per_parameter * max_mixture_parameters
  parameter_floor_side <- if (is.na(parameter_floor_nodes)) NA_integer_ else
    as.integer(ceiling(sqrt(parameter_floor_nodes)))
  if (!automatic_dimensions && is.null(map_from) &&
      som_xdim * som_ydim < 10 * ncol(zmat)) {
    warning(sprintf(paste0("Explicit SOM dimensions (%d x %d = %d nodes) are below the floor of 10 nodes per ",
      "clustering variable (%d nodes for %d variables); diagnostics are retained for review."),
      som_xdim, som_ydim, som_xdim * som_ydim, 10L * ncol(zmat), ncol(zmat)), call. = FALSE)
  }
  stability_core_settings <- .ResolveClusterStabilityCores(
    stability_cores, stability_resamples, n = n_complete,
    p = mixture_dimension, method = "som_mclust", som_nodes = nrow(som_codes),
    subsample_fraction = max(stability_fractions))

  # SOM distances (distance to BMU) ----------------------------------------

  SOM_Node_full <- rep(NA_integer_, nrow(df_scidr))
  SOM_Node_full[complete_rows] <- som_model$unit.classif

  # kohonen reports squared Euclidean distances; SOM_Distance is reported on
  # the Euclidean scale, in Z-score units.
  SOM_Dist_full <- rep(NA_real_, nrow(df_scidr))
  SOM_Dist_full[complete_rows] <- sqrt(som_model$distances)

  if (soft_kernel == "gaussian" && is.null(soft_bandwidth)) {
    soft_bandwidth <- stats::median(SOM_Dist_full[complete_rows])
    if (!is.finite(soft_bandwidth) || soft_bandwidth <= 0) soft_bandwidth <- 1
  }
  train_quant_err <- mean(SOM_Dist_full[complete_rows], na.rm = TRUE)

  som_residuals <- .SOMResiduals(zmat, som_codes, som_model$unit.classif,
    vars_used, df_scidr$.row_id[complete_rows])
  som_variable_contribution <- .SOMVariableContribution(som_residuals, vars_used, df)

  som_occupancy_summary <- dplyr::tibble(
    n_nodes = nrow(som_node_occupancy),
    n_empty_nodes = sum(som_node_occupancy$n == 0),
    n_singleton_nodes = sum(som_node_occupancy$n == 1),
    min_occupancy = min(som_node_occupancy$n),
    median_occupancy = stats::median(som_node_occupancy$n),
    max_occupancy = max(som_node_occupancy$n)
  )

  # AweSOM plots -----------------------------------------------------------

  CircPlot  <- aweSOM::aweSOMplot(som = som_model, type = "Circular", palvar = "viridis")
  LinePlot  <- aweSOM::aweSOMplot(som = som_model, type = "Line",     palvar = "viridis")
  CloudPlot <- aweSOM::aweSOMplot(som = som_model, type = "Cloud")

  if (Relabel) {
    relabel_fun <- function(w) {
      if (!is.null(w$x$label)) {
        zvars     <- w$x$label
        base_vars <- sub(paste0("^", ZScorePrefix), "", zvars)
        base_vars <- sub("^Z_", "", base_vars)

        get_lab <- NULL
        if (requireNamespace("Hmisc", quietly = TRUE)) {
          get_lab <- function(v) Hmisc::label(df[[v]])
        } else if (requireNamespace("sjlabelled", quietly = TRUE)) {
          get_lab <- function(v) sjlabelled::get_label(df[[v]])
        }

        if (!is.null(get_lab)) {
          new_labels <- vapply(base_vars, function(bv) {
            if (!bv %in% names(df)) return(bv)
            lab <- get_lab(bv)
            if (is.null(lab) || !nzchar(lab)) bv else as.character(lab)
          }, character(1))
        } else {
          new_labels <- base_vars
        }

        w$x$label <- unname(new_labels)
      }
      w
    }

    CircPlot  <- relabel_fun(CircPlot)
    LinePlot  <- relabel_fun(LinePlot)
  }

  som_plots <- list(
    Circular = CircPlot,
    Line     = LinePlot,
    Cloud    = CloudPlot
  )

  # Clustering the SOM codebook ---------------------------------------------

  # The default node clustering is the latent-profile grid below. A caller that
  # clusters the codebook another way supplies .NodeClusterFn and the grid is
  # never fitted: SOM + HDBSCAN, for instance, used to pay for every LPA
  # candidate and then discard all of them in favour of its own node labels.
  if (is.null(.NodeClusterFn)) {
    # LPA / mclust on SOM codes ----------------------------------------------

    if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
    X <- as.data.frame(X)
    set.seed(seed_mixture)

    lpa_control <- NULL
    if (!is.null(mixture_em_itmax) || !is.null(mixture_em_tol)) {
      lpa_control <- mclust::emControl(
        itmax = if (is.null(mixture_em_itmax)) NULL else as.integer(mixture_em_itmax),
        tol   = if (is.null(mixture_em_tol)) NULL else mixture_em_tol
      )
    }

    # This helper is intentionally kept inside the function because it isolates
    # warning capture, runtime tracking, and failed-fit handling for one fragile step.
    safe_lpa_fit <- function(X, k, model, lpa_control = NULL, mixture_timeout_seconds = NULL) {

      warnings_vec <- character(0)
      error_msg    <- NA_character_
      start_time   <- Sys.time()

      fit_expr <- function() {
        if (node_weighting == "occupancy") {
          return(.FitWeightedSOMMclust(
            X, weights = lpa_node_weights, k = k, model = model,
            itmax = if (is.null(mixture_em_itmax)) 100L else mixture_em_itmax,
            tol = if (is.null(mixture_em_tol)) 1e-05 else mixture_em_tol))
        }
        if (is.null(lpa_control)) {
          .FitTidyLPAMclust(X, k, model, nboot = mixture_bootstrap_n)
        } else {
          .FitTidyLPAMclust(X, k, model, control = lpa_control, nboot = mixture_bootstrap_n)
        }
      }

      fit <- withCallingHandlers(
        tryCatch(
          suppressMessages({
            if (is.null(mixture_timeout_seconds)) {
              fit_expr()
            } else {
              R.utils::withTimeout(
                fit_expr(),
                timeout = mixture_timeout_seconds,
                onTimeout = "error"
              )
            }
          }),
          TimeoutException = function(e) {
            error_msg <<- paste0("Timed out after ", mixture_timeout_seconds, " seconds.")
            NULL
          },
          error = function(e) {
            error_msg <<- conditionMessage(e)
            NULL
          }
        ),
        warning = function(w) {
          warnings_vec <<- c(warnings_vec, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      runtime_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      slow_fit <- !is.null(slow_fit_seconds) && runtime_seconds > slow_fit_seconds

      fit_info <- NULL
      fit_error_msg <- NA_character_
      if (!is.null(fit)) {
        fit_info <- if (inherits(fit, "WeightedSOMMclust")) {
          fit$fit_info
        } else {
          tryCatch(
            tidyLPA::get_fit(fit),
            error = function(e) {
              fit_error_msg <<- conditionMessage(e)
              NULL
            }
          )
        }
      }

      # When mclust cannot estimate a solution, tidyLPA still returns fit
      # indices but attaches no class assignments. Such a fit cannot assign
      # SOM nodes, so it is recorded as failed rather than as a candidate.
      if (!is.null(fit_info)) {
        has_assignments <- if (inherits(fit, "WeightedSOMMclust")) {
          length(fit$classification) == nrow(X)
        } else tryCatch({
          profile <- if (inherits(fit, "tidyLPA")) fit[[1]] else fit
          "Class" %in% names(tidyLPA::get_data(profile))
        }, error = function(e) FALSE)
        if (!has_assignments) {
          fit_info <- NULL
          fit_error_msg <- paste0(
            "mclust could not estimate model ", model, " with k = ", k,
            "; no class assignments were returned.")
        }
      }

      if (is.null(fit_info)) {
        if (is.na(error_msg) && !is.na(fit_error_msg)) {
          error_msg <- fit_error_msg
        }
        status <- "failed"
      } else {
        status <- "success"
      }

      diagnostics <- dplyr::tibble(
        Model           = as.integer(model),
        Classes         = as.integer(k),
        status          = status,
        runtime_seconds = runtime_seconds,
        slow_fit        = slow_fit,
        n_warnings      = length(warnings_vec),
        warnings        = paste(unique(warnings_vec), collapse = " | "),
        error           = error_msg
      )

      if (!is.null(fit_info)) {
        fit_info <- as.data.frame(fit_info)
        fit_info$Model <- as.integer(model)
        fit_info$Classes <- as.integer(k)
        # Without bootstrap resamples the BLRT is not estimated
        if (mixture_bootstrap_n < 1) {
          fit_info[intersect(c("BLRT_val", "BLRT_p"), names(fit_info))] <- NA_real_
        }
      }

      list(
        fit         = fit,
        fit_info    = fit_info,
        diagnostics = diagnostics
      )
    }

    # This helper prevents one available model from getting all of the AHP weight
    # simply because there is no variation in a fit index.
    scale_for_ahp <- function(x, higher_is_better = TRUE) {
      if (!higher_is_better) {
        x <- -x
      }

      finite_x <- x[is.finite(x)]
      if (length(finite_x) < 2 || stats::sd(finite_x) == 0) {
        return(rep(0, length(x)))
      }
      as.numeric(scale(x))
    }

    # Class and posterior for every SOM node. Nodes in the mixture fit take
    # tidyLPA's assignments; nodes left out by mixture_min_node_people are classified
    # by the fitted mclust model.
    get_node_posteriors <- function(lpa_fit) {
      if (inherits(lpa_fit, "WeightedSOMMclust")) {
        cluster <- rep(NA_integer_, nrow(som_codes))
        prob <- matrix(NA_real_, nrow(som_codes), lpa_fit$k)
        cluster[lpa_nodes] <- lpa_fit$classification
        prob[lpa_nodes, ] <- lpa_fit$z
        other_nodes <- setdiff(seq_len(nrow(som_codes)), lpa_nodes)
        if (length(other_nodes)) {
          predicted <- .WeightedSOMPosterior(lpa_fit,
            som_codes[other_nodes, colnames(X), drop = FALSE])
          cluster[other_nodes] <- predicted$classification
          prob[other_nodes, ] <- predicted$z
        }
        return(list(cluster = cluster, prob = prob))
      }
      profile <- if (inherits(lpa_fit, "tidyLPA")) {
        if (length(lpa_fit) != 1) {
          stop("An LPA candidate must contain exactly one profile.")
        }
        lpa_fit[[1]]
      } else {
        lpa_fit
      }

      profile_data <- tidyLPA::get_data(profile)
      if (!"Class" %in% names(profile_data)) {
        stop("The LPA candidate does not contain class assignments.")
      }
      prob_cols <- grep("^CPROB", names(profile_data), value = TRUE)
      k <- length(prob_cols)
      cluster <- rep(NA_integer_, nrow(som_codes))
      prob <- matrix(NA_real_, nrow(som_codes), k)
      cluster[lpa_nodes] <- as.integer(profile_data$Class)
      prob[lpa_nodes, ] <- as.matrix(profile_data[, prob_cols, drop = FALSE])

      other_nodes <- setdiff(seq_len(nrow(som_codes)), lpa_nodes)
      if (length(other_nodes) > 0) {
        predicted <- stats::predict(profile$model, newdata = as.data.frame(
          som_codes[other_nodes, colnames(X), drop = FALSE]))
        cluster[other_nodes] <- as.integer(predicted$classification)
        prob[other_nodes, ] <- predicted$z
      }
      list(cluster = cluster, prob = prob)
    }
    get_profile_assignments <- function(lpa_fit) {
      get_node_posteriors(lpa_fit)$cluster
    }
    candidate_max_prob <- function(fit_name) {
      node_max <- apply(get_node_posteriors(lpa_models_all[[fit_name]])$prob, 1, max)
      out <- rep(NA_real_, nrow(df_scidr))
      out[complete_rows] <- node_max[som_model$unit.classif]
      out
    }
    PosteriorSummary <- function(max_prob) {
      max_prob <- max_prob[complete_rows]
      dplyr::tibble(
        PosteriorMean = mean(max_prob, na.rm = TRUE),
        PosteriorMedian = stats::median(max_prob, na.rm = TRUE),
        PropLowPosterior = mean(max_prob < clear_membership_threshold, na.rm = TRUE))
    }
    # Participants per profile, from training node occupancy.
    profile_participant_sizes <- function(node_cluster, k) {
      counts <- vapply(seq_len(k), function(cluster) {
        sum(som_node_occupancy$n[node_cluster == cluster], na.rm = TRUE)
      }, numeric(1))
      list(min_n = min(counts), min_prop = min(counts) / sum(counts))
    }
    AddMixtureSupport <- function(table) {
      mixture_parameters <- vapply(seq_len(nrow(table)), function(i) {
        .SOMMclustParameterCount(table$Model[[i]], ncol(X), table$Classes[[i]])
      }, integer(1))
      table %>%
        dplyr::mutate(
          MixtureParameters = mixture_parameters,
          MixtureNodes = nrow(X),
          ParticipantsPerNode = n_complete / (som_xdim * som_ydim),
          OccupancyCompliant = .data$ParticipantsPerNode >= min_participants_per_node,
          NodesPerParameter = .data$MixtureNodes / .data$MixtureParameters,
          Overparameterized = .data$NodesPerParameter < 1,
          SupportStatus = dplyr::case_when(
            .data$Overparameterized ~ "over-parameterized",
            .data$NodesPerParameter < min_nodes_per_parameter ~ "under-supported",
            TRUE ~ "supported")
        )
    }

    cluster_jaccard <- function(reference, resampled) {
      reference <- as.integer(reference)
      resampled <- as.integer(resampled)
      valid <- !is.na(reference) & !is.na(resampled)
      reference <- reference[valid]
      resampled <- resampled[valid]

      if (length(reference) == 0) {
        return(dplyr::tibble(Cluster = integer(), Jaccard = numeric()))
      }

      dplyr::bind_rows(lapply(sort(unique(reference)), function(cluster) {
        reference_members <- reference == cluster
        candidates <- sort(unique(resampled))
        jaccard <- vapply(candidates, function(candidate) {
          candidate_members <- resampled == candidate
          sum(reference_members & candidate_members) /
            sum(reference_members | candidate_members)
        }, numeric(1))

        dplyr::tibble(
          Cluster = cluster,
          Jaccard = max(jaccard)
        )
      }))
    }

    Stability <- NULL
    if (method == "exploratory") {

      lpa_models <- list()
      fit_rows   <- list()
      diag_rows  <- list()

      model_grid <- expand.grid(
        Model   = models,
        Classes = k_range,
        KEEP.OUT.ATTRS = FALSE
      )
      # Stability refits fit only the shortlisted model/k pairs
      if (!is.null(.CandidatePairs)) {
        model_grid <- dplyr::semi_join(model_grid, .CandidatePairs, by = c("Model", "Classes"))
      }

      total_fits <- nrow(model_grid)
      .WarnSOMNodesPerCluster(nrow(X), unique(model_grid$Classes), candidate_min_nodes_per_cluster)
      model_problem_counts <- stats::setNames(rep(0L, length(unique(models))), unique(models))

      for (i in seq_len(total_fits)) {

        this_model <- model_grid$Model[i]
        this_k     <- model_grid$Classes[i]
        fit_name   <- paste0("model_", this_model, "_class_", this_k)


        if (!is.null(skip_model_after_n_failures) &&
            model_problem_counts[as.character(this_model)] >= skip_model_after_n_failures) {

          if (mixture_progress) {
            message(
              "Skipping LPA fit ", i, "/", total_fits,
              ": model ", this_model,
              ", k = ", this_k,
              " after repeated failed or slow fits."
            )
          }

          diag_rows[[fit_name]] <- dplyr::tibble(
            Model           = as.integer(this_model),
            Classes         = as.integer(this_k),
            status          = "skipped_model_family",
            runtime_seconds = NA_real_,
            slow_fit        = NA,
            n_warnings      = 0L,
            warnings        = NA_character_,
            error           = paste0(
              "Skipped after ", skip_model_after_n_failures,
              " failed or slow fits for model ", this_model, "."
            )
          )
          next
        }

        if (mixture_progress) {
          message(
            "LPA fit ", i, "/", total_fits,
            ": model ", this_model,
            ", k = ", this_k
          )
        }

        fit_res <- safe_lpa_fit(
          X           = X,
          k           = this_k,
          model       = this_model,
          lpa_control = lpa_control,
          mixture_timeout_seconds = mixture_timeout_seconds
        )

        diag_rows[[fit_name]] <- fit_res$diagnostics

        if (fit_res$diagnostics$status != "success" || isTRUE(fit_res$diagnostics$slow_fit)) {
          model_problem_counts[as.character(this_model)] <- model_problem_counts[as.character(this_model)] + 1L
        }

        if (!is.null(fit_res$fit) && !is.null(fit_res$fit_info)) {
          lpa_models[[fit_name]] <- fit_res$fit
          fit_rows[[fit_name]]   <- fit_res$fit_info
        }
      }

      lpa_diagnostics <- dplyr::bind_rows(diag_rows)
      failed_fits <- lpa_diagnostics %>%
        dplyr::filter(status != "success")

      warning_fits <- lpa_diagnostics %>%
        dplyr::filter(status == "success", n_warnings > 0)

      if (length(fit_rows) == 0) {
        stop(
          "No LPA models were successfully estimated. Inspect ModelInfo_MClust$diagnostics if available; ",
          "try fewer variables, smaller k_range, excluding unstable models, or increasing mixture_em_itmax."
        )
      }

      fit_table <- dplyr::bind_rows(fit_rows) %>%
        dplyr::mutate(
          Model = as.integer(Model),
          Classes = as.integer(Classes)
        )

      lpa_models_all <- lpa_models
      fit_table$BoundaryPct <- vapply(seq_len(nrow(fit_table)), function(i) {
        fit_name <- paste0("model_", fit_table$Model[i], "_class_", fit_table$Classes[i])
        projected <- .SOMSoftProjection(zmat, som_codes,
          .SOMNodePosteriorTable(get_node_posteriors(lpa_models_all[[fit_name]])),
          n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth, margin_threshold = boundary_margin_threshold,
          entropy_threshold = boundary_entropy_threshold)
        mean(projected$boundary)
      }, numeric(1))

      profile_sizes <- dplyr::bind_rows(lapply(names(lpa_models), function(fit_name) {
        assignment <- get_profile_assignments(lpa_models[[fit_name]])
        n_classes <- as.integer(sub("^.*_class_([0-9]+)$", "\\1", fit_name))
        counts <- tabulate(assignment, nbins = n_classes)
        participants <- profile_participant_sizes(assignment, n_classes)
        dplyr::tibble(
          Model = as.integer(sub("^model_([0-9]+)_class_.*$", "\\1", fit_name)),
          Classes = n_classes,
          MinProfileNodeN = min(counts), MaxProfileNodeN = max(counts),
          MinProfileNodeProportion = min(counts) / length(assignment),
          MaxProfileNodeProportion = max(counts) / length(assignment),
          MinProfileParticipantN = participants$min_n,
          MinProfileParticipantProportion = participants$min_prop
        ) %>%
          dplyr::bind_cols(PosteriorSummary(candidate_max_prob(fit_name)))
      }))
      fit_table <- fit_table %>%
        dplyr::rename_with(~ c(n_min = "MinProfileNodeProportion",
          n_max = "MaxProfileNodeProportion", BLRT_val = "BLRTStatistic",
          BLRT_p = "BLRTPValue")[.x], dplyr::any_of(c("n_min", "n_max", "BLRT_val", "BLRT_p"))) %>%
        dplyr::select(-dplyr::any_of(c("MinProfileNodeProportion", "MaxProfileNodeProportion"))) %>%
        dplyr::left_join(profile_sizes, by = c("Model", "Classes"))

      if (!"AIC" %in% names(fit_table)) fit_table$AIC <- NA_real_
      if (!"BIC" %in% names(fit_table)) fit_table$BIC <- NA_real_
      if (!"Entropy" %in% names(fit_table)) fit_table$Entropy <- NA_real_

      Stability <- NULL
      if (stability_resamples > 0) {
        reference_rows <- which(complete_rows)
        replicate_rows <- list()
        cluster_rows <- list()
        assignment_rows <- list()
        candidate_rows <- list()
        candidate_index <- 0L
        stability_backend <- "serial"
        resolved_stability_cores <- 1L

        # Shortlist: stability is estimated only for the top candidates by the
        # composite index without reproducibility, among eligible candidates.
        candidate_names <- names(lpa_models)
        if (stability_candidates == "best_per_k" ||
            (!is.null(stability_top_n) && stability_top_n < length(candidate_names))) {
          fit_weights <- composite_weights$weights
          fit_weights <- if (is.null(fit_weights)) {
            c(AIC = 1, BIC = 1, Entropy = 1)
          } else {
            fit_weights[intersect(names(fit_weights), c("AIC", "BIC", "Entropy", "Boundary"))]
          }
          preliminary <- as.numeric(cbind(
            AIC = scale_for_ahp(fit_table$AIC, higher_is_better = FALSE),
            BIC = scale_for_ahp(fit_table$BIC, higher_is_better = FALSE),
            Entropy = scale_for_ahp(fit_table$Entropy),
            Boundary = scale_for_ahp(fit_table$BoundaryPct, higher_is_better = FALSE))[, names(fit_weights), drop = FALSE] %*%
            (fit_weights / sum(fit_weights)))
          eligible <- if (is.null(min_cluster_prop)) rep(TRUE, nrow(fit_table)) else
            fit_table$MinProfileParticipantProportion >= min_cluster_prop
          ranked <- order(!eligible, -preliminary)
          # best_per_k keeps the highest-ranked candidate at each k, so
          # stability is compared across k rather than within the top fits
          if (stability_candidates == "best_per_k") {
            ranked <- ranked[!duplicated(fit_table$Classes[ranked])]
          }
          if (!is.null(stability_top_n)) {
            ranked <- ranked[seq_len(min(stability_top_n, length(ranked)))]
          }
          shortlist <- paste0("model_", fit_table$Model, "_class_", fit_table$Classes)[ranked]
          candidate_names <- intersect(candidate_names, shortlist)
        }
        candidate_info <- dplyr::tibble(
          fit_name = candidate_names,
          Model = as.integer(sub("^model_([0-9]+)_class_.*$", "\\1", candidate_names)),
          Classes = as.integer(sub("^.*_class_([0-9]+)$", "\\1", candidate_names)))
        # Core participants are retained as a descriptive diagnostic only.
        reference_core <- lapply(stats::setNames(candidate_names, candidate_names),
          function(fit_name) candidate_max_prob(fit_name)[reference_rows] >= clear_membership_threshold)
        reference_node_posteriors <- lapply(stats::setNames(candidate_names, candidate_names),
          function(fit_name) get_node_posteriors(lpa_models[[fit_name]]))
        reference_soft <- lapply(candidate_names, function(fit_name) {
          .SOMSoftProjection(zmat, som_codes,
            .SOMNodePosteriorTable(reference_node_posteriors[[fit_name]]),
            n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth,
            margin_threshold = boundary_margin_threshold,
            entropy_threshold = boundary_entropy_threshold)
        })
        names(reference_soft) <- candidate_names
        reference_boundary <- lapply(reference_soft, `[[`, "boundary")
        reference_assignments <- lapply(candidate_names, function(fit_name) {
          hard_assignment <- reference_node_posteriors[[fit_name]]$cluster[som_model$unit.classif]
          if (assignment == "soft") reference_soft[[fit_name]]$cluster else hard_assignment
        })
        names(reference_assignments) <- candidate_names

        # Each resample trains one SOM and fits every shortlisted candidate on
        # its codebook: the map does not depend on the model family or k, so
        # it is not retrained per candidate.
        stability_tasks <- expand.grid(
          Fraction = c(if ("subsample" %in% stability_type) stability_fractions,
            if ("seed" %in% stability_type) 1),
          Replicate = seq_len(stability_resamples),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE)
        if (stability_progress) {
          message("Running ", nrow(stability_tasks), " stability refits for ",
            length(candidate_names), " candidate", if (length(candidate_names) == 1L) "" else "s",
            " with ", stability_core_settings$resolved, " worker",
            if (stability_core_settings$resolved == 1L) "" else "s", ".")
        }
        # Capture the in-memory constructor. This matters during development:
        # multisession workers otherwise resolve an older installed namespace
        # instead of the version loaded by devtools::load_all().
        refit_constructor <- CreateClusterModel_SOM_MClust
        resample_results <- .ClusterParallelMap(seq_len(nrow(stability_tasks)), function(task_index) {
          retained_fraction <- stability_tasks$Fraction[[task_index]]
          replicate <- stability_tasks$Replicate[[task_index]]
          sampling_seed <- if (identical(retained_fraction, 0.90)) {
            as.integer(seed_stability + 100000L + replicate)
          } else {
            as.integer(seed_stability + 200000L + round(retained_fraction * 10000) * 100L + replicate)
          }
          model_seed <- if (retained_fraction == 1) as.integer(seed_stability + 300000L + replicate) else as.integer(seed_som)
          mixture_seed <- if (retained_fraction == 1) as.integer(seed_stability + 400000L + replicate) else seed_mixture
          set.seed(sampling_seed)
          subsample_rows <- if (retained_fraction == 1) reference_rows else sort(sample(reference_rows,
            max(2L, floor(length(reference_rows) * retained_fraction)), replace = FALSE))
          df_subsample <- df_scidr[subsample_rows, , drop = FALSE]

          subsample_result <- tryCatch(
            suppressWarnings(refit_constructor(
              data = df_subsample,
              variables = vars_used,
              method = "exploratory",
              k_range = sort(unique(candidate_info$Classes)),
              models = sort(unique(candidate_info$Model)),
              ClusterVariableName = ".scidr_stability_cluster",
              ZScoreType = ZScoreType,
              ZScoreObject = ZScoreObject,
              som_xdim = som_xdim,
              som_ydim = som_ydim,
              som_topo = som_topo,
              som_neigh = som_neigh, som_epochs = som_epochs,
              node_weighting = node_weighting,
              assignment = assignment,
              soft_n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth_requested,
              boundary_margin_threshold = boundary_margin_threshold,
              boundary_entropy_threshold = boundary_entropy_threshold,
              seed_som = model_seed,
              seed_mixture = mixture_seed,
              Relabel = FALSE,
              ZScorePrefix = ZScorePrefix,
              ZScoreVars = ZScoreVars_used,
              mixture_progress = FALSE,
              mixture_em_itmax = mixture_em_itmax,
              mixture_em_tol = mixture_em_tol,
              mixture_bootstrap_n = 0L,
              mixture_timeout_seconds = mixture_timeout_seconds,
              mixture_drop_zero_sd = mixture_drop_zero_sd,
              mixture_zero_sd_tol = mixture_zero_sd_tol,
              skip_model_after_n_failures = NULL,
              slow_fit_seconds = slow_fit_seconds,
              candidate_min_nodes_per_cluster = NULL,
              min_nodes_per_parameter = min_nodes_per_parameter,
              min_participants_per_node = min_participants_per_node,
              mixture_min_node_people = mixture_min_node_people,
              far_distance_quantile = far_distance_quantile,
              clear_membership_threshold = clear_membership_threshold,
              stability_resamples = 0L,
              .CandidatePairs = candidate_info[c("Model", "Classes")]
            )),
            error = function(e) e
          )

          lapply(seq_len(nrow(candidate_info)), function(i) {
            fit_name <- candidate_info$fit_name[[i]]
            candidate_model <- candidate_info$Model[[i]]
            candidate_k <- candidate_info$Classes[[i]]
            reference_assignment <- reference_assignments[[fit_name]]
            candidate_posterior <- if (inherits(subsample_result, "error")) NULL else
              subsample_result$ModelInfo_MClust$candidate_node_posteriors[[fit_name]]
            if (is.null(candidate_posterior)) {
              error_msg <- if (inherits(subsample_result, "error")) {
                conditionMessage(subsample_result)
              } else {
                "The candidate could not be estimated on this subsample."
              }
              return(list(
                replicate = dplyr::tibble(
                  Model = candidate_model, Classes = candidate_k,
                  Replicate = replicate, SamplingSeed = sampling_seed,
                  RetainedFraction = retained_fraction, RetainedParticipants = length(subsample_rows),
                  ModelSeed = model_seed, Status = "failed", ARI = NA_real_,
                  Error = error_msg),
                recovery = NULL, assignment = NULL))
            }

            refit_individual <- tryCatch({
              df_projection <- .AddClusterRowID(df_scidr, "data")
              z_vars <- subsample_result$ZScoreVars
              if (subsample_result$ZScoreType == "PreZScored") {
                z_projection <- df_projection[, z_vars, drop = FALSE]
              } else {
                z_projection <- SciDataReportR::ProjectZScore(
                  data = df_projection, variables = subsample_result$vars_used,
                  parameters = subsample_result$ZScoreObject,
                  ParameterInputType = "ZScoreObj", names_prefix = "Z_")$ZScores[, z_vars, drop = FALSE]
              }
              projected_complete <- stats::complete.cases(z_projection)
              projected_cluster <- rep(NA_integer_, nrow(df_projection))
              if (any(projected_complete)) {
                mapped_nodes <- kohonen::map(subsample_result$ModelInfo_SOM$som_model,
                  newdata = as.matrix(z_projection[projected_complete, , drop = FALSE]))$unit.classif
                hard_assignment <- candidate_posterior$cluster[mapped_nodes]
                if (assignment == "soft") {
                  soft_assignment <- .SOMSoftProjection(
                    as.matrix(z_projection[projected_complete, , drop = FALSE]),
                    subsample_result$ModelInfo_SOM$som_codes,
                    .SOMNodePosteriorTable(candidate_posterior),
                    n_nodes = soft_n_nodes, soft_kernel = soft_kernel,
                    soft_bandwidth = subsample_result$Specification$projection$soft_bandwidth,
                    margin_threshold = boundary_margin_threshold,
                    entropy_threshold = boundary_entropy_threshold)
                  projected_cluster[projected_complete] <- soft_assignment$cluster
                } else {
                  projected_cluster[projected_complete] <- hard_assignment
                }
              }
              dplyr::tibble(.row_id = df_projection$.row_id, Cluster = projected_cluster)
            }, error = function(e) e)
            if (inherits(refit_individual, "error")) {
              return(list(
                replicate = dplyr::tibble(
                  Model = candidate_model, Classes = candidate_k,
                  Replicate = replicate, SamplingSeed = sampling_seed,
                  RetainedFraction = retained_fraction, RetainedParticipants = length(subsample_rows),
                  ModelSeed = model_seed, Status = "failed", ARI = NA_real_,
                  Error = conditionMessage(refit_individual)),
                recovery = NULL, assignment = NULL))
            }
            comparison <- dplyr::tibble(
              .row_id = df_scidr$.row_id[reference_rows],
              ReferenceCluster = reference_assignment,
              Core = reference_core[[fit_name]],
              NonBoundary = !reference_boundary[[fit_name]]) %>%
              dplyr::inner_join(
                dplyr::tibble(
                  .row_id = refit_individual$.row_id,
                  RefitCluster = refit_individual$Cluster),
                by = ".row_id", relationship = "one-to-one")
            valid <- !is.na(comparison$ReferenceCluster) &
              !is.na(comparison$RefitCluster)
            ari <- if (sum(valid) > 1) {
              mclust::adjustedRandIndex(
                comparison$ReferenceCluster[valid], comparison$RefitCluster[valid])
            } else {
              NA_real_
            }
            partition_metrics <- .ClusterPartitionMetrics(
              comparison$ReferenceCluster[valid], comparison$RefitCluster[valid])
            jaccard <- cluster_jaccard(
              comparison$ReferenceCluster, comparison$RefitCluster) %>%
              dplyr::mutate(
                Model = candidate_model,
                Classes = candidate_k,
                Replicate = replicate,
                RetainedFraction = retained_fraction
              )
            core <- valid & comparison$Core %in% TRUE
            ari_core <- if (sum(core) > 1) {
              mclust::adjustedRandIndex(
                comparison$ReferenceCluster[core], comparison$RefitCluster[core])
            } else {
              NA_real_
            }
            jaccard_core <- cluster_jaccard(
              comparison$ReferenceCluster[core], comparison$RefitCluster[core]) %>%
              dplyr::mutate(Model = candidate_model, Classes = candidate_k,
                Replicate = replicate, RetainedFraction = retained_fraction)
            non_boundary <- valid & comparison$NonBoundary %in% TRUE
            ari_non_boundary <- if (sum(non_boundary) > 1) {
              mclust::adjustedRandIndex(
                comparison$ReferenceCluster[non_boundary], comparison$RefitCluster[non_boundary])
            } else {
              NA_real_
            }
            jaccard_non_boundary <- cluster_jaccard(
              comparison$ReferenceCluster[non_boundary], comparison$RefitCluster[non_boundary]) %>%
              dplyr::mutate(Model = candidate_model, Classes = candidate_k,
                Replicate = replicate, RetainedFraction = retained_fraction)

            co_membership <- .ClusterCoMembershipSummary(
              comparison$ReferenceCluster[valid], comparison$RefitCluster[valid])
            split_merge <- .ClusterSplitMergeDiagnostics(
              comparison$ReferenceCluster[valid], comparison$RefitCluster[valid])
            profile_recovery <- .ClusterProfileRecovery(
              comparison$ReferenceCluster[valid], comparison$RefitCluster[valid],
              # `zmat` is aligned to the complete reference cohort.  The
              # original data-row positions in `reference_rows` can include
              # omitted records, so they must not index this matrix.
              zmat[valid, , drop = FALSE])
            DecorateStabilityDiagnostic <- function(x) {
              if (!nrow(x)) return(x)
              dplyr::mutate(x, Model = candidate_model, Classes = candidate_k,
                Replicate = replicate, RetainedFraction = retained_fraction)
            }

            assignment_full <- rep(NA_integer_, length(reference_assignment))
            assignment_full[match(
              comparison$.row_id, df_scidr$.row_id[reference_rows])] <-
              comparison$RefitCluster
            list(
              replicate = dplyr::tibble(
                Model = candidate_model, Classes = candidate_k,
                Replicate = replicate, SamplingSeed = sampling_seed,
                RetainedFraction = retained_fraction, RetainedParticipants = length(subsample_rows),
                ModelSeed = model_seed, Status = "success", ARI = ari,
                ParticipantsCompared = sum(valid),
                ARI_Core = ari_core, CoreProportion = mean(core[valid]),
                ARI_NonBoundary = ari_non_boundary,
                NonBoundaryProportion = mean(non_boundary[valid]),
                VI = partition_metrics[["VI"]], NMI = partition_metrics[["NMI"]],
                FowlkesMallows = partition_metrics[["FowlkesMallows"]],
                Error = NA_character_),
              recovery = jaccard,
              recovery_core = jaccard_core,
              recovery_non_boundary = jaccard_non_boundary,
              co_membership_by_cluster = DecorateStabilityDiagnostic(co_membership$by_cluster),
              co_membership_between_clusters = DecorateStabilityDiagnostic(co_membership$between_clusters),
              split_merge_reference = DecorateStabilityDiagnostic(split_merge$reference),
              split_merge_refit = DecorateStabilityDiagnostic(split_merge$refit),
              profile_recovery = DecorateStabilityDiagnostic(profile_recovery),
              assignment = list(Model = candidate_model, Classes = candidate_k,
                RetainedFraction = retained_fraction,
                reference = reference_assignment, assignment = assignment_full))
          })
        }, cores = stability_core_settings$resolved)
        stability_backend <- attr(resample_results, "backend")
        resolved_stability_cores <- attr(resample_results, "resolved_cores")
        candidate_results_all <- do.call(c, resample_results)
        replicate_rows <- lapply(candidate_results_all, `[[`, "replicate")
        # Some legacy/custom node-clustering callbacks return the pre-panel
        # replicate schema. Keep those fits usable while making the new
        # descriptive measures explicitly unavailable rather than failing the
        # complete stability run.
        panel_numeric_fields <- c(
          "RetainedFraction", "RetainedParticipants", "ARI", "ARI_Core",
          "ARI_NonBoundary", "CoreProportion", "NonBoundaryProportion",
          "ParticipantsCompared", "VI", "NMI", "FowlkesMallows"
        )
        replicate_rows <- lapply(replicate_rows, function(df_replicate) {
          for (field_name in panel_numeric_fields) {
            if (!field_name %in% names(df_replicate)) {
              df_replicate[[field_name]] <- NA_real_
            }
          }
          if (all(is.na(df_replicate$RetainedFraction))) {
            df_replicate$RetainedFraction <- primary_stability_fraction
          }
          df_replicate
        })
        cluster_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "recovery"))
        cluster_core_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "recovery_core"))
        cluster_non_boundary_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "recovery_non_boundary"))
        co_membership_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "co_membership_by_cluster"))
        co_membership_between_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "co_membership_between_clusters"))
        split_merge_reference_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "split_merge_reference"))
        split_merge_refit_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "split_merge_refit"))
        profile_recovery_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "profile_recovery"))
        assignment_rows <- Filter(Negate(is.null), lapply(candidate_results_all, `[[`, "assignment"))

        for (i in seq_len(nrow(candidate_info))) {
          for (stability_kind in stability_type) {
          summary_fraction <- if (stability_kind == "seed") 1 else primary_stability_fraction
          candidate_model <- candidate_info$Model[[i]]
          candidate_k <- candidate_info$Classes[[i]]

          candidate_replicates <- dplyr::bind_rows(replicate_rows) %>%
            dplyr::filter(Model == candidate_model, Classes == candidate_k,
              .data$RetainedFraction == summary_fraction)
          candidate_clusters <- if (length(cluster_rows) > 0) {
            dplyr::bind_rows(cluster_rows) %>%
              dplyr::filter(Model == candidate_model, Classes == candidate_k,
                .data$RetainedFraction == summary_fraction)
          } else {
            dplyr::tibble(Cluster = integer(), Jaccard = numeric())
          }
          cluster_summary <- if (nrow(candidate_clusters) > 0) {
            candidate_clusters %>%
              dplyr::group_by(Cluster) %>%
              dplyr::summarise(MeanJaccard = mean(Jaccard, na.rm = TRUE), .groups = "drop")
          } else {
            dplyr::tibble(Cluster = integer(), MeanJaccard = numeric())
          }

          success_rate <- mean(candidate_replicates$Status == "success")
          mean_ari <- mean(candidate_replicates$ARI, na.rm = TRUE)
          mean_jaccard <- if (nrow(cluster_summary) > 0) mean(cluster_summary$MeanJaccard) else NA_real_
          min_jaccard <- if (nrow(cluster_summary) > 0) min(cluster_summary$MeanJaccard) else NA_real_
          reproducibility_all <- mean_ari

          core_summary <- if (length(cluster_core_rows) > 0) {
            dplyr::bind_rows(cluster_core_rows) %>%
              dplyr::filter(Model == candidate_model, Classes == candidate_k,
                .data$RetainedFraction == summary_fraction) %>%
              dplyr::group_by(Cluster) %>%
              dplyr::summarise(MeanJaccard = mean(Jaccard, na.rm = TRUE), .groups = "drop")
          } else {
            dplyr::tibble(Cluster = integer(), MeanJaccard = numeric())
          }
          mean_ari_core <- if ("ARI_Core" %in% names(candidate_replicates)) {
            mean(candidate_replicates$ARI_Core, na.rm = TRUE)
          } else {
            NA_real_
          }
          if (is.nan(mean_ari_core)) mean_ari_core <- NA_real_
          mean_jaccard_core <- if (nrow(core_summary) > 0) mean(core_summary$MeanJaccard) else NA_real_
          min_jaccard_core <- if (nrow(core_summary) > 0) min(core_summary$MeanJaccard) else NA_real_
          core_proportion <- if ("CoreProportion" %in% names(candidate_replicates)) {
            mean(candidate_replicates$CoreProportion, na.rm = TRUE)
          } else {
            NA_real_
          }
          non_boundary_summary <- if (length(cluster_non_boundary_rows) > 0) {
            dplyr::bind_rows(cluster_non_boundary_rows) %>%
              dplyr::filter(Model == candidate_model, Classes == candidate_k,
                .data$RetainedFraction == summary_fraction) %>%
              dplyr::group_by(Cluster) %>%
              dplyr::summarise(MeanJaccard = mean(Jaccard, na.rm = TRUE), .groups = "drop")
          } else {
            dplyr::tibble(Cluster = integer(), MeanJaccard = numeric())
          }
          mean_ari_non_boundary <- mean(candidate_replicates$ARI_NonBoundary, na.rm = TRUE)
          if (is.nan(mean_ari_non_boundary)) mean_ari_non_boundary <- NA_real_
          mean_jaccard_non_boundary <- if (nrow(non_boundary_summary) > 0) mean(non_boundary_summary$MeanJaccard) else NA_real_
          min_jaccard_non_boundary <- if (nrow(non_boundary_summary) > 0) min(non_boundary_summary$MeanJaccard) else NA_real_
          non_boundary_proportion <- mean(!reference_boundary[[candidate_info$fit_name[[i]]]])
          # ARI measures the agreement of the final participant partition and
          # is therefore the primary reproducibility score for selection.
          reproducibility_score <- mean_ari

          candidate_rows[[length(candidate_rows) + 1L]] <- dplyr::tibble(
            StabilityType = stability_kind,
            Model = candidate_model,
            Classes = candidate_k,
            StabilitySuccessRate = success_rate,
            StabilityARI_Mean = mean_ari,
            StabilityJaccard_Mean = mean_jaccard,
            StabilityJaccard_Min = min_jaccard,
            StabilityCoreProportion = core_proportion,
            StabilityARI_Core = mean_ari_core,
            StabilityJaccard_CoreMean = mean_jaccard_core,
            StabilityJaccard_CoreMin = min_jaccard_core,
            StabilityNonBoundaryProportion = non_boundary_proportion,
            StabilityARI_NonBoundary = mean_ari_non_boundary,
            StabilityJaccard_NonBoundaryMean = mean_jaccard_non_boundary,
            StabilityJaccard_NonBoundaryMin = min_jaccard_non_boundary,
            ReproducibilityScore_All = reproducibility_all,
            ReproducibilityScore = reproducibility_score
          )
        }

        }
        stability_summary <- dplyr::bind_rows(candidate_rows)
        stability_replicates <- dplyr::bind_rows(replicate_rows)
        SafeStabilityMean <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
        stability_deletion_curve <- stability_replicates %>%
          dplyr::group_by(.data$Model, .data$Classes, .data$RetainedFraction) %>%
          dplyr::summarise(
            StabilitySuccessRate = mean(.data$Status == "success"),
            StabilityARI_Mean = SafeStabilityMean(.data$ARI),
            StabilityVI_Mean = SafeStabilityMean(.data$VI),
            StabilityNMI_Mean = SafeStabilityMean(.data$NMI),
            StabilityFowlkesMallows_Mean = SafeStabilityMean(.data$FowlkesMallows),
            .groups = "drop")
        stability_clusters <- if (length(cluster_rows) > 0) {
          dplyr::bind_rows(cluster_rows)
        } else {
          dplyr::tibble(
            Cluster = integer(), Jaccard = numeric(), Model = integer(),
            Classes = integer(), Replicate = integer(), RetainedFraction = numeric()
          )
        }
        fit_table <- fit_table %>%
          dplyr::left_join(stability_summary %>%
            dplyr::group_by(Model, Classes) %>%
            dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop"),
            by = c("Model", "Classes"))
        primary_assignment_rows <- Filter(function(x) {
          identical(x$RetainedFraction, primary_stability_fraction)
        }, assignment_rows)
        diagnostics <- lapply(split(primary_assignment_rows, vapply(primary_assignment_rows,
          function(x) paste(x$Model, x$Classes, sep = "_"), character(1))), function(rows) {
          .ClusterStabilityDiagnostics(
            rows[[1]]$reference, lapply(rows, `[[`, "assignment"),
            seq_along(rows[[1]]$reference), df_scidr$.row_id[complete_rows])
        })
        diagnostic_keys <- names(diagnostics)
        participant_inclusion <- dplyr::bind_rows(lapply(seq_along(diagnostics), function(i) {
          key <- strsplit(diagnostic_keys[[i]], "_", fixed = TRUE)[[1]]
          dplyr::mutate(diagnostics[[i]]$participant_inclusion,
            Model = as.integer(key[[1]]), Classes = as.integer(key[[2]]))
        }))
        cluster_inclusion <- dplyr::bind_rows(lapply(seq_along(diagnostics), function(i) {
          key <- strsplit(diagnostic_keys[[i]], "_", fixed = TRUE)[[1]]
          dplyr::mutate(diagnostics[[i]]$cluster_inclusion,
            Model = as.integer(key[[1]]), Classes = as.integer(key[[2]]))
        }))
        co_membership_by_cluster <- if (length(co_membership_rows)) {
          dplyr::bind_rows(co_membership_rows)
        } else dplyr::tibble()
        co_membership_between_clusters <- if (length(co_membership_between_rows)) {
          dplyr::bind_rows(co_membership_between_rows)
        } else dplyr::tibble()
        split_merge_reference <- if (length(split_merge_reference_rows)) {
          dplyr::bind_rows(split_merge_reference_rows)
        } else dplyr::tibble()
        split_merge_refit <- if (length(split_merge_refit_rows)) {
          dplyr::bind_rows(split_merge_refit_rows)
        } else dplyr::tibble()
        profile_recovery <- if (length(profile_recovery_rows)) {
          dplyr::bind_rows(profile_recovery_rows)
        } else dplyr::tibble()
        boundary_summary <- dplyr::bind_rows(lapply(seq_len(nrow(candidate_info)), function(i) {
          fit_name <- candidate_info$fit_name[[i]]
          reference_assignment <- reference_assignments[[fit_name]]
          reference_flag <- reference_boundary[[fit_name]]
          dplyr::bind_rows(
            dplyr::tibble(
              Model = candidate_info$Model[[i]], Classes = candidate_info$Classes[[i]],
              Scope = "overall", Cluster = NA_integer_,
              Participants = length(reference_flag),
              BoundaryParticipants = sum(reference_flag, na.rm = TRUE),
              BoundaryProportion = mean(reference_flag, na.rm = TRUE)
            ),
            dplyr::tibble(Cluster = sort(unique(reference_assignment))) %>%
              dplyr::mutate(
                Model = candidate_info$Model[[i]], Classes = candidate_info$Classes[[i]],
                Scope = "reference_cluster",
                Participants = vapply(Cluster, function(cluster) sum(reference_assignment == cluster), integer(1)),
                BoundaryParticipants = vapply(Cluster, function(cluster) sum(reference_flag[reference_assignment == cluster]), integer(1)),
                BoundaryProportion = BoundaryParticipants / Participants
              ) %>%
              dplyr::select(Model, Classes, Scope, Cluster, Participants,
                BoundaryParticipants, BoundaryProportion)
          )
        }))
        Stability <- list(
          settings = list(
            resamples = as.integer(stability_resamples),
            seed = seed_stability,
            refit_scope = "full_pipeline_projected_to_reference_cohort",
            top_n = stability_top_n,
            candidates_assessed = candidate_names,
            comparison_scope = "all_reference_participants",
            resample_type = if (identical(stability_type, "subsample")) "subsample_without_replacement"
              else if (identical(stability_type, "seed")) "seed_refit_all_participants" else "subsample_and_seed",
            stability_type = stability_type,
            resample_fraction = primary_stability_fraction,
            resample_fractions = stability_fractions,
            primary_resample_fraction = primary_stability_fraction,
            coassignment_limit = 2000L,
            noise_policy = "all clusters included",
            requested_cores = stability_core_settings$requested,
            resolved_cores = resolved_stability_cores,
            backend = stability_backend,
            worker_policy = stability_core_settings$worker_policy,
            estimated_bytes_per_worker = stability_core_settings$bytes_per_worker,
            estimated_total_bytes = stability_core_settings$estimated_total_bytes,
            physical_memory_bytes = stability_core_settings$physical_memory_bytes,
            available_memory_bytes = stability_core_settings$available_memory_bytes,
            safe_memory_budget_bytes = stability_core_settings$safe_memory_budget_bytes,
            memory_method = stability_core_settings$method,
            som_grid = list(
              xdim = som_xdim,
              ydim = som_ydim,
              policy = "fixed_to_reference_fit",
              initialization_seed_policy = if ("seed" %in% stability_type) "varied_for_seed_refits" else "fixed_to_reference_fit"),
            som_training = list(
              mode = "batch",
              epochs = som_training_epochs,
              initialization_seed = seed_som,
              initialization_policy = if ("seed" %in% stability_type) "varied_for_seed_refits" else "fixed_to_reference_fit",
              cohort = if (identical(stability_type, "seed")) "all_participants" else "by_stability_type"),
            assignment = assignment,
            soft_n_nodes = as.integer(soft_n_nodes),
            boundary_margin_threshold = boundary_margin_threshold,
            boundary_entropy_threshold = boundary_entropy_threshold,
            non_boundary_mask = "reference_soft_margin_and_entropy"
          ),
          replicates = stability_replicates,
          cluster_recovery = stability_clusters,
          summary = stability_summary,
          deletion_curve = stability_deletion_curve,
          failures = stability_replicates %>% dplyr::filter(.data$Status != "success"),
          participant_inclusion = participant_inclusion,
          cluster_inclusion = cluster_inclusion,
          boundary_summary = boundary_summary,
          coassignment = lapply(diagnostics, `[[`, "coassignment"),
          co_membership_by_cluster = co_membership_by_cluster,
          co_membership_between_clusters = co_membership_between_clusters,
          split_merge_reference = split_merge_reference,
          split_merge_refit = split_merge_refit,
          profile_recovery = profile_recovery
        )
        Stability <- lapply(Stability, function(value) {
          if (is.data.frame(value) && "RetainedFraction" %in% names(value)) {
            value$StabilityType <- ifelse(value$RetainedFraction == 1, "seed", "subsample")
          }
          value
        })
        Stability$cluster_recovery_summary <- Stability$cluster_recovery %>%
          dplyr::group_by(Model, Classes, .data$StabilityType, Cluster) %>%
          dplyr::summarise(MeanJaccard = mean(Jaccard, na.rm = TRUE), .groups = "drop")
      }

      if (!is.null(Stability) && !is.null(Stability$boundary_summary)) {
        boundary_overall <- Stability$boundary_summary %>%
          dplyr::filter(.data$Scope == "overall") %>%
          dplyr::transmute(.data$Model, .data$Classes,
            BoundaryParticipantProportion = .data$BoundaryProportion)
        boundary_maximum <- Stability$boundary_summary %>%
          dplyr::filter(.data$Scope == "reference_cluster") %>%
          dplyr::group_by(.data$Model, .data$Classes) %>%
          dplyr::summarise(
            MaxClusterBoundaryProportion = max(.data$BoundaryProportion, na.rm = TRUE),
            .groups = "drop")
        fit_table <- fit_table %>%
          dplyr::left_join(boundary_overall, by = c("Model", "Classes")) %>%
          dplyr::left_join(boundary_maximum, by = c("Model", "Classes"))
      }

      fit_table <- AddMixtureSupport(fit_table) %>%
        dplyr::mutate(
          AIC_scaled = scale_for_ahp(AIC, higher_is_better = FALSE),
          BIC_scaled = scale_for_ahp(BIC, higher_is_better = FALSE),
          Entropy_scaled = scale_for_ahp(Entropy, higher_is_better = TRUE)
        )

      if (stability_resamples > 0) {
        fit_table <- fit_table %>%
          dplyr::mutate(
            Reproducibility_scaled = scale_for_ahp(ReproducibilityScore))
      }

      # Composite rank index: a weighted mean of the scaled criteria. Equal
      # weights by default; criteria_weights supplies explicit weights or an
      # AHP pairwise-comparison matrix.
      fit_table$Boundary_scaled <- scale_for_ahp(fit_table$BoundaryPct, higher_is_better = FALSE)
      fit_table$Posterior_scaled <- scale_for_ahp(fit_table$PosteriorMean)
      criteria_used <- c("AIC", "BIC", "Entropy",
        if (stability_resamples > 0) "Reproducibility",
        if (!is.null(composite_weights$weights) &&
            "Posterior" %in% names(composite_weights$weights)) "Posterior",
        if (!is.null(composite_weights$weights) &&
            "Boundary" %in% names(composite_weights$weights)) "Boundary")
      weights_used <- composite_weights$weights
      if (is.null(weights_used)) {
        weights_used <- stats::setNames(rep(1, length(criteria_used)), criteria_used)
      }
      missing_criteria <- setdiff(names(weights_used), criteria_used)
      if (length(missing_criteria) > 0) {
        stop("criteria_weights names criteria that are not computed: ",
          paste(missing_criteria, collapse = ", "),
          ". Reproducibility requires stability_resamples > 0.")
      }
      weights_used <- weights_used[names(weights_used) %in% criteria_used]
      weights_used <- weights_used / sum(weights_used)
      scaled_matrix <- as.matrix(fit_table[paste0(names(weights_used), "_scaled")])
      fit_table$ahp_index <- as.numeric(scaled_matrix %*% weights_used)

      overparameterized <- fit_table %>%
        dplyr::filter(.data$Overparameterized) %>%
        dplyr::transmute(label = paste0("Model ", .data$Model,
          ", k = ", .data$Classes)) %>%
        dplyr::pull(.data$label)
      if (length(overparameterized)) {
        warning("The following SOM mixture candidates have fewer mixture nodes than free parameters and are ineligible for recommendation: ",
          paste(overparameterized, collapse = "; "), call. = FALSE)
      }
      fit_table$Eligible <- !fit_table$Overparameterized & if (is.null(min_cluster_prop)) {
        TRUE
      } else {
        fit_table$MinProfileParticipantProportion >= min_cluster_prop
      }

      if (all(is.na(fit_table$ahp_index))) {
        stop("LPA models were estimated, but the composite rank index could not be computed because fit indices are missing.")
      }

      eligible_rows <- which(fit_table$Eligible & !is.na(fit_table$ahp_index))
      recommended <- length(eligible_rows) > 0L
      if (length(eligible_rows) == 0) {
        warning("No candidate is eligible for the composite recommendation; an internal fallback is retained only to construct the review object.",
          call. = FALSE)
        eligible_rows <- which(!is.na(fit_table$ahp_index))
      }
      best_idx     <- eligible_rows[which.max(fit_table$ahp_index[eligible_rows])]
      fit_table$Recommended <- FALSE
      if (recommended) fit_table$Recommended[[best_idx]] <- TRUE
      node_assignments <- lapply(lpa_models, get_profile_assignments)
      ahp_best_row <- if (recommended) fit_table[best_idx, ] else fit_table[0, , drop = FALSE]
      best_model   <- as.integer(fit_table$Model[[best_idx]])
      best_k       <- as.integer(fit_table$Classes[[best_idx]])
      best_name    <- paste0("model_", best_model, "_class_", best_k)
      best_fit_name <- best_name
      best_lpa     <- lpa_models[[best_name]]

      weight_txt <- if (is.null(composite_weights$weights)) {
        "equal weights"
      } else {
        paste0(composite_weights$method, " weights: ", paste0(names(weights_used), " ",
          sprintf("%.2f", weights_used), collapse = ", "))
      }
      recommendation_txt <- if (!recommended) {
        "No eligible SOM mixture candidate is available for composite recommendation."
      } else paste0(
        "Composite rank index (", paste(criteria_used, collapse = ", "), "; ",
        weight_txt, ")",
        if (!is.null(min_cluster_prop)) paste0(
          " among candidates whose smallest profile holds at least ",
          round(100 * min_cluster_prop, 1), "% of participants") else "",
        " recommends Model ", best_model, " with k = ", best_k, " profiles.")

      mdata_wide <- fit_table

    } else {

      if (is.null(final_k) || is.null(final_model)) {
        stop("For method = 'finalize', final_k and final_model must be supplied.")
      }

      if (mixture_progress) {
        message(
          "LPA finalize fit: model ", final_model,
          ", k = ", final_k
        )
      }

      fit_res <- safe_lpa_fit(
        X           = X,
        k           = final_k,
        model       = final_model,
        lpa_control = lpa_control
      )

      lpa_diagnostics <- fit_res$diagnostics
      failed_fits <- lpa_diagnostics %>%
        dplyr::filter(status != "success")
      warning_fits <- lpa_diagnostics %>%
        dplyr::filter(status == "success", n_warnings > 0)

      if (is.null(fit_res$fit) || is.null(fit_res$fit_info)) {
        stop(
          "The requested final LPA model could not be estimated: model ",
          unname(.MclustModelNames[as.character(final_model)]), ", k = ", final_k, "."
        )
      }

      lpa_models <- fit_res$fit
      best_lpa   <- lpa_models
      node_assignments <- stats::setNames(
        list(get_profile_assignments(best_lpa)),
        paste0("model_", final_model, "_class_", final_k))
      lpa_models_all <- stats::setNames(list(best_lpa), names(node_assignments))
      best_fit_name <- paste0("model_", final_model, "_class_", final_k)

      fit_table <- as.data.frame(fit_res$fit_info) %>%
        dplyr::mutate(
          Classes        = final_k,
          Model          = final_model,
          AIC_scaled     = NA_real_,
          BIC_scaled     = NA_real_,
          Entropy_scaled = NA_real_,
          ahp_index      = NA_real_
        )

      ahp_best_row       <- fit_table[1, ]
      best_model         <- final_model
      best_k             <- final_k
      recommendation_txt <- paste0(
        "User-specified Model ", final_model,
        " with k = ", final_k, " profiles."
      )

      mdata_wide <- fit_table
    }

    if (any(c("n_min", "n_max", "BLRT_val", "BLRT_p") %in% names(fit_table))) {
      fit_table <- fit_table %>%
        dplyr::rename_with(~ c(n_min = "MinProfileNodeProportion",
          n_max = "MaxProfileNodeProportion", BLRT_val = "BLRTStatistic",
          BLRT_p = "BLRTPValue")[.x], dplyr::any_of(c("n_min", "n_max", "BLRT_val", "BLRT_p")))
    }
    if (!all(c("MinProfileNodeN", "MaxProfileNodeN") %in% names(fit_table))) {
      best_node_assignment <- get_profile_assignments(best_lpa)
      counts <- tabulate(best_node_assignment, nbins = best_k)
      participants <- profile_participant_sizes(best_node_assignment, best_k)
      fit_table <- fit_table %>%
        dplyr::mutate(
          MinProfileNodeN = min(counts), MaxProfileNodeN = max(counts),
          MinProfileNodeProportion = min(counts) / length(best_node_assignment),
          MaxProfileNodeProportion = max(counts) / length(best_node_assignment),
          MinProfileParticipantN = participants$min_n,
          MinProfileParticipantProportion = participants$min_prop
        ) %>%
        dplyr::bind_cols(PosteriorSummary(candidate_max_prob(best_fit_name)))
    }
    if (!all(c("MixtureParameters", "MixtureNodes", "ParticipantsPerNode",
        "OccupancyCompliant", "NodesPerParameter", "Overparameterized",
        "SupportStatus") %in% names(fit_table))) {
      fit_table <- AddMixtureSupport(fit_table)
      overparameterized <- fit_table %>% dplyr::filter(.data$Overparameterized)
      if (nrow(overparameterized)) {
        warning("The requested SOM mixture has fewer mixture nodes than free parameters and is marked ineligible.",
          call. = FALSE)
      }
    }
    if (!"Eligible" %in% names(fit_table)) {
      fit_table$Eligible <- !fit_table$Overparameterized
    }
    if (method == "finalize") {
      fit_table$Recommended <- isTRUE(fit_table$Eligible[[1]])
      if (!isTRUE(fit_table$Eligible[[1]])) {
        ahp_best_row <- fit_table[0, , drop = FALSE]
        recommendation_txt <- paste0(
          "User-specified Model ", final_model, " with k = ", final_k,
          " is retained as a projectable review fit but is overparameterized and not recommended."
        )
      }
    }
    fit_table <- fit_table %>%
      dplyr::mutate(
        GridDesignPolicy = grid_design_policy,
        GridDesignK = grid_design_k,
        TargetNodesPerProfile = grid_nodes_per_profile,
        SOMNodeCount = som_xdim * som_ydim,
        ExpectedNodesPerProfile = .data$SOMNodeCount / .data$Classes
      )
    if (!"BoundaryPct" %in% names(fit_table)) {
      fit_table$BoundaryPct <- vapply(seq_len(nrow(fit_table)), function(i) {
        fit_name <- paste0("model_", fit_table$Model[i], "_class_", fit_table$Classes[i])
        projected <- .SOMSoftProjection(zmat, som_codes,
          .SOMNodePosteriorTable(get_node_posteriors(lpa_models_all[[fit_name]])),
          n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth, margin_threshold = boundary_margin_threshold,
          entropy_threshold = boundary_entropy_threshold)
        mean(projected$boundary)
      }, numeric(1))
    }
    fit_table$ModelName <- unname(.MclustModelNames[as.character(fit_table$Model)])
    # tidyLPA reports the information criteria on the -2 log-likelihood scale
    # (lower is better) but ICL on mclust's scale (higher is better).
    fit_table <- .SetClusterMetricDirection(fit_table,
      higher = c("ICL", "Entropy"),
      lower = c("AIC", "AWE", "BIC", "CAIC", "CLC", "KIC", "SABIC"))
    mdata_wide <- fit_table

    # Fit plot ---------------------------------------------------------------

    mdata_cols <- intersect(
      c("AIC", "BIC", "Entropy", "BoundaryPct", "PosteriorMean", "PropLowPosterior",
        "StabilityARI_Mean", "StabilityARI_NonBoundary",
        "BoundaryParticipantProportion",
        "MaxClusterBoundaryProportion", "ReproducibilityScore", "BLRTPValue"),
      names(mdata_wide)
    )
    # criteria not estimated for any candidate (e.g. BLRT without bootstraps)
    # are left out of the plot
    mdata_cols <- mdata_cols[vapply(mdata_cols,
      function(column) any(is.finite(mdata_wide[[column]])), logical(1))]

    if (length(mdata_cols) == 0) {
      mdata <- dplyr::tibble(
        Classes = integer(0),
        Model   = integer(0),
        name    = character(0),
        value   = numeric(0)
      )
    } else {
      mdata <- mdata_wide %>%
        tidyr::pivot_longer(
          cols      = dplyr::all_of(mdata_cols),
          names_to  = "name",
          values_to = "value"
        )
    }
    blrt_label <- "Bootstrap likelihood-ratio test p-value (\u2212log10 scale)"
    mdata <- mdata %>%
      dplyr::mutate(
        PlotValue = dplyr::if_else(
          .data$name == "BLRTPValue",
          -log10(pmax(.data$value, .Machine$double.xmin)), .data$value),
        name = dplyr::recode(.data$name, BLRTPValue = blrt_label),
        name = factor(.data$name, levels = c(
          "AIC", "BIC", "Entropy", "BoundaryPct", "PosteriorMean", "PropLowPosterior",
          "StabilityARI_Mean", "StabilityARI_NonBoundary",
          "BoundaryParticipantProportion",
          "MaxClusterBoundaryProportion", "ReproducibilityScore", blrt_label)))

    model_levels <- as.integer(names(.MclustModelLabels))
    model_levels <- model_levels[!is.na(model_levels)]
    model_labels <- unname(.MclustModelLabels[as.character(model_levels)])

    mdata$Model <- factor(
      mdata$Model,
      levels = model_levels,
      labels = model_labels
    )

    # Model family is a plain category, so it takes the package palette. Colours
    # are assigned over the full set of families, not the subset present, so a
    # family keeps the same colour whether or not the others converged.
    model_families <- model_labels
    pal_cols <- stats::setNames(
      .SciDataColorValues(length(model_families)), model_families)
    pal_cols <- pal_cols[names(pal_cols) %in% levels(mdata$Model)]

    fit_plot <- if ("SupportStatus" %in% names(mdata)) {
      ggplot2::ggplot(mdata, ggplot2::aes(
        x = Classes, y = PlotValue, color = Model, shape = SupportStatus)) +
        ggplot2::geom_point() + ggplot2::geom_line()
    } else {
      ggplot2::ggplot(mdata, ggplot2::aes(
        x = Classes, y = PlotValue, color = Model)) +
        ggplot2::geom_point() + ggplot2::geom_line()
    }
    fit_plot <- fit_plot +
      ggplot2::facet_wrap(~name, scales = "free_y", ncol = 1,
        labeller = ggplot2::labeller(name = function(values) {
          preferred <- .ClusterMetricDirection(values, fit_table)
          ifelse(is.na(preferred), values,
            paste0(values, " (", preferred, " is better)"))
        })) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top") +
      ggplot2::labs(x = "Number of clusters", y = "Fit index value")

    if (length(pal_cols) > 0) {
      fit_plot <- fit_plot + ggplot2::scale_color_manual(values = pal_cols)
    }
    if ("SupportStatus" %in% names(mdata)) {
      fit_plot <- fit_plot + ggplot2::scale_shape_manual(
        values = c("supported" = 16, "under-supported" = 1,
          "over-parameterized" = 4), name = "SOM mixture support",
        drop = FALSE)
    }
    if ("BLRTPValue" %in% mdata_cols) {
      fit_plot <- fit_plot + ggplot2::geom_hline(
        data = dplyr::tibble(name = factor(blrt_label, levels = levels(mdata$name)),
          PlotValue = -log10(0.05)),
        ggplot2::aes(yintercept = PlotValue), inherit.aes = FALSE,
        linetype = "dashed", color = "grey35")
    }

    # Node-level cluster and posteriors --------------------------------------

    # One row per SOM node: its codebook vector, class, and posteriors.
    # tidyLPA::get_data.tidyLPA() dispatches through an unqualified generic in
    # tidyLPA 2.0.2, so get_node_posteriors() extracts the tidyProfile first.
    node_posteriors <- get_node_posteriors(best_lpa)
    node_prob <- node_posteriors$prob
    colnames(node_prob) <- paste0("prob_", seq_len(ncol(node_prob)))

    node_df <- dplyr::bind_cols(
      dplyr::tibble(model_number = as.integer(best_model),
        classes_number = as.integer(best_k)),
      dplyr::as_tibble(som_codes),
      dplyr::as_tibble(node_prob),
      dplyr::tibble(
        Cluster = node_posteriors$cluster,
        NodeID = seq_len(nrow(som_codes)),
        NodeN = som_node_occupancy$n,
        InLPAFit = seq_len(nrow(som_codes)) %in% lpa_nodes)
    )
    prob_cols_new <- grep("^prob_", names(node_df), value = TRUE)

    if (length(prob_cols_new) > 0) {
      node_df <- node_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          max_prob      = max(dplyr::c_across(dplyr::all_of(prob_cols_new)), na.rm = TRUE),
          prob_assigned = dplyr::c_across(dplyr::all_of(prob_cols_new))[Cluster],
          uncertainty   = 1 - max_prob
        ) %>%
        dplyr::ungroup()
    } else {
      node_df <- node_df %>%
        dplyr::mutate(
          max_prob      = NA_real_,
          prob_assigned = NA_real_,
          uncertainty   = NA_real_
        )
    }

    node_df <- node_df %>%
      dplyr::mutate(
        NodeID = as.integer(NodeID),
        Cluster = as.integer(Cluster)
      )
  } else {

    node_clustering <- .NodeClusterFn(som_codes)
    if (!is.list(node_clustering) || is.null(node_clustering$node_cluster)) {
      stop(".NodeClusterFn must return a list with a `node_cluster` element.")
    }
    if (length(node_clustering$node_cluster) != nrow(som_codes)) {
      stop(".NodeClusterFn must return one cluster per SOM node.")
    }

    X <- as.data.frame(som_codes)
    dropped_lpa_vars <- character(0)
    lpa_nodes <- seq_len(nrow(som_codes))
    lpa_node_weights <- som_node_occupancy$n

    # No latent-profile model was fitted, so there are no posterior
    # probabilities to report. The columns are kept so that everything
    # downstream still sees one node table with one shape.
    node_df <- dplyr::tibble(
      NodeID        = seq_len(nrow(som_codes)),
      Cluster       = as.integer(node_clustering$node_cluster),
      max_prob      = NA_real_,
      prob_assigned = NA_real_,
      uncertainty   = NA_real_
    )

    lpa_models         <- NULL
    best_lpa           <- NULL
    node_assignments   <- NULL
    best_fit_name      <- if (is.null(node_clustering$best_fit_name)) {
      "supplied_node_clustering"
    } else node_clustering$best_fit_name
    fit_table          <- node_clustering$fit_table
    ahp_best_row       <- node_clustering$ahp_best_row
    recommendation_txt <- node_clustering$recommendation
    fit_plot           <- node_clustering$fit_plot
    lpa_diagnostics    <- dplyr::tibble()
    failed_fits        <- dplyr::tibble()
    warning_fits       <- dplyr::tibble()
    Stability          <- NULL
  }

  # Map nodes to individuals -----------------------------------------------

  # Only node assignments and posteriors are carried to individuals. The node
  # table also holds the codebook vectors under the Z-score names, which would
  # otherwise read as each participant's own Z-scores.
  individual_tbl <- dplyr::tibble(
    .row_id      = df_scidr$.row_id,
    SOM_Node     = SOM_Node_full,
    SOM_Distance = SOM_Dist_full
  ) %>%
    dplyr::left_join(.SOMNodeAssignments(node_df), by = c("SOM_Node" = "NodeID"))

  individual_tbl$SOM_Distance_LOO <- NA_real_
  individual_tbl$SOM_Distance_LOO[complete_rows] <- .SOMLeaveOneOutDistances(zmat, som_codes, som_model$unit.classif)
  individual_tbl$NodeN <- som_node_occupancy$n[match(individual_tbl$SOM_Node,
    seq_len(nrow(som_codes)))]
  individual_tbl$Flag_OnePersonNode <- individual_tbl$NodeN == 1L
  one_person_share <- mean(individual_tbl$Flag_OnePersonNode[complete_rows])
  if (one_person_share > 0.05) {
    message(sprintf("%.1f%% of complete participants are on one-person nodes. These people's node probabilities describe only themselves; certainty should come from the soft probabilities.",
      100 * one_person_share))
  }
  individual_tbl$NodeCluster <- individual_tbl$Cluster

  soft_projection <- .SOMSoftProjection(zmat, som_codes, node_df,
    n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth, margin_threshold = boundary_margin_threshold,
    entropy_threshold = boundary_entropy_threshold)
  individual_tbl$Cluster_Hard <- individual_tbl$Cluster
  individual_tbl$Cluster_Soft <- NA_integer_
  individual_tbl$SoftEntropy <- NA_real_
  individual_tbl$SoftMargin <- NA_real_
  individual_tbl$Flag_BoundaryParticipant <- NA
  individual_tbl$BoundaryReason <- NA_character_
  if (!is.null(soft_projection)) {
    individual_tbl$Cluster_Soft[complete_rows] <- soft_projection$cluster
    individual_tbl$SoftEntropy[complete_rows] <- soft_projection$entropy
    individual_tbl$SoftMargin[complete_rows] <- soft_projection$margin
    individual_tbl$Flag_BoundaryParticipant[complete_rows] <- soft_projection$boundary
    individual_tbl$BoundaryReason[complete_rows] <- soft_projection$reason
    for (j in seq_len(ncol(soft_projection$probability))) {
      individual_tbl[[paste0("SoftProb_", j)]] <- NA_real_
      individual_tbl[[paste0("SoftProb_", j)]][complete_rows] <- soft_projection$probability[, j]
    }
  }
  individual_tbl$Cluster <- if (assignment == "soft") individual_tbl$Cluster_Soft else individual_tbl$Cluster_Hard
  individual_tbl$Flag_NodeDisagreement <- individual_tbl$Cluster != individual_tbl$NodeCluster

  individual_tbl <- individual_tbl %>%
    dplyr::left_join(som_residuals$top, by = ".row_id")

  if (!is.null(id_col) && id_col %in% names(df_scidr)) {
    individual_tbl[[id_col]] <- df_scidr[[id_col]]
  }

  if (nrow(individual_tbl) != nrow(df_scidr)) {
    stop("Internal row alignment error: ProbFit$individual does not match the number of rows in df.")
  }
  attr(individual_tbl$.row_id, "label") <- "Row ID"

  # DataWithClusters: canonical label and assignment diagnostics ------------

  DataWithClusters <- df_scidr

  if (ClusterVariableName %in% names(DataWithClusters)) {
    message("Column '", ClusterVariableName, "' already exists and will be overwritten.")
  }

  DataWithClusters[[ClusterVariableName]] <- individual_tbl$Cluster
  projection_columns <- c("NodeN", "Flag_OnePersonNode", "Cluster_Hard", "Cluster_Soft", "SoftEntropy", "SoftMargin",
    "Flag_BoundaryParticipant", "BoundaryReason",
    grep("^SoftProb_", names(individual_tbl), value = TRUE))
  for (column in projection_columns) {
    if (column %in% names(individual_tbl)) {
      DataWithClusters[[column]] <- individual_tbl[[column]]
    }
  }

  # SOMFit diagnostics and distance baselines (training) -------------------

  som_fit_tbl    <- individual_tbl
  som_fit_non_na <- som_fit_tbl[!is.na(som_fit_tbl$SOM_Distance), ]

  dist_train <- som_fit_non_na$SOM_Distance
  dist_summary <- stats::quantile(
    dist_train,
    probs = c(0, 0.25, 0.5, 0.75, 0.95, 1),
    na.rm = TRUE
  )
  names(dist_summary) <- c("min", "q25", "median", "q75", "p95", "max")

  overall_mean <- mean(dist_train, na.rm = TRUE)
  overall_sd   <- stats::sd(dist_train, na.rm = TRUE)

  dist_by_cluster <- som_fit_non_na %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      mean = mean(SOM_Distance),
      sd   = stats::sd(SOM_Distance),
      p90  = stats::quantile(SOM_Distance, 0.90),
      p95  = stats::quantile(SOM_Distance, 0.95),
      p99  = stats::quantile(SOM_Distance, 0.99),
      .groups = "drop"
    )

  overall_p95 <- dist_summary["p95"]

  flag_by_cluster <- som_fit_non_na %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      n                 = dplyr::n(),
      prop_overall_high = mean(SOM_Distance > overall_p95),
      .groups = "drop"
    )

  overall_high_cutoff <- as.numeric(stats::quantile(
    dist_train,
    probs = far_distance_quantile,
    na.rm = TRUE
  ))

  dist_by_cluster_cutoffs <- som_fit_non_na %>%
    dplyr::filter(!is.na(.data$SOM_Distance), !is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      train_cluster_high_cutoff = as.numeric(stats::quantile(
        .data$SOM_Distance,
        probs = far_distance_quantile,
        na.rm = TRUE
      )),
      .groups = "drop"
    )

  som_fit_tbl <- som_fit_tbl %>%
    dplyr::mutate(
      SOMDist_z_overall = if (is.finite(overall_sd) && overall_sd > 0) {
        (.data$SOM_Distance - overall_mean) / overall_sd
      } else {
        NA_real_
      },
      SOMDist_percentile_train = vapply(.data$SOM_Distance, function(x) {
        if (is.na(x) || length(dist_train) == 0) {
          NA_real_
        } else {
          mean(dist_train <= x, na.rm = TRUE)
        }
      }, numeric(1)),
      Flag_SOMDist_overallHigh =
        !is.na(.data$SOM_Distance) &
        !is.na(overall_high_cutoff) &
        .data$SOM_Distance > overall_high_cutoff
    ) %>%
    dplyr::left_join(dist_by_cluster_cutoffs, by = "Cluster") %>%
    dplyr::mutate(
      Flag_SOMDist_clusterHigh =
        !is.na(.data$SOM_Distance) &
        !is.na(.data$train_cluster_high_cutoff) &
        .data$SOM_Distance > .data$train_cluster_high_cutoff,
      Projection_Fit_Class = dplyr::case_when(
        is.na(.data$SOM_Distance) ~ NA_character_,
        (.data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh) &
          (is.na(.data$max_prob) | .data$max_prob < clear_membership_threshold) ~
          "Potential novel phenotype",
        .data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh ~
          "Poor fit to training structure",
        !is.na(.data$max_prob) & .data$max_prob < clear_membership_threshold ~
          "Uncertain membership",
        TRUE ~ "Good fit"
      )
    )

  som_fit_non_na <- som_fit_tbl[!is.na(som_fit_tbl$SOM_Distance), ]

  training_cluster_fit_summary <- som_fit_non_na %>%
    dplyr::filter(!is.na(.data$Cluster)) %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data$SOM_Distance, na.rm = TRUE),
      median_distance = stats::median(.data$SOM_Distance, na.rm = TRUE),
      mean_probability = mean(.data$prob_assigned, na.rm = TRUE),
      median_probability = stats::median(.data$prob_assigned, na.rm = TRUE),
      mean_distance_percentile = mean(.data$SOMDist_percentile_train, na.rm = TRUE),
      median_distance_percentile = stats::median(.data$SOMDist_percentile_train, na.rm = TRUE),
      prop_high_distance = mean(.data$Flag_SOMDist_overallHigh | .data$Flag_SOMDist_clusterHigh, na.rm = TRUE),
      prop_low_probability = mean(.data$prob_assigned < clear_membership_threshold, na.rm = TRUE),
      prop_poor_fit = mean(.data$Projection_Fit_Class == "Poor fit to training structure", na.rm = TRUE),
      prop_potential_novel = mean(.data$Projection_Fit_Class == "Potential novel phenotype", na.rm = TRUE),
      .groups = "drop"
    )

  potential_novel_cases <- som_fit_tbl %>%
    dplyr::filter(.data$Projection_Fit_Class == "Potential novel phenotype")

  p_dist_hist <- ggplot2::ggplot(
    som_fit_non_na,
    ggplot2::aes(x = SOM_Distance)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "SOM distances (training data)",
      x     = "Euclidean distance to BMU (Z-score units)",
      y     = "Count"
    )

  ProjectionReference <- list(
    calibration = .CalibrateSOMFarRule(individual_tbl$SOM_Distance_LOO,
      individual_tbl$Cluster, far_distance_quantile),
    cluster_profiles = lapply(split(seq_len(nrow(zmat)), individual_tbl$Cluster[complete_rows]),
      function(idx) colMeans(zmat[idx, , drop = FALSE])),
    loo_distances = individual_tbl[, c("SOM_Node", "Cluster", "SOM_Distance_LOO")],
    training_variable_summary = training_variable_summary,
    distance_reference = list(
      mean = overall_mean,
      median = stats::median(dist_train, na.rm = TRUE),
      sd = overall_sd,
      p90 = as.numeric(stats::quantile(dist_train, 0.90, na.rm = TRUE)),
      p95 = as.numeric(stats::quantile(dist_train, 0.95, na.rm = TRUE)),
      p99 = as.numeric(stats::quantile(dist_train, 0.99, na.rm = TRUE)),
      far_distance_quantile = far_distance_quantile,
      high_dist_cutoff = overall_high_cutoff
    ),
    node_occupancy = som_node_occupancy,
    distance_scale = "euclidean",
    variable_contribution = som_variable_contribution,
    projection_thresholds = list(
      far_distance_quantile = far_distance_quantile,
      clear_membership_threshold = clear_membership_threshold
    ),
    MapHealth = list(
      n_training = nrow(df_scidr), n_complete = sum(complete_rows),
      n_excluded = sum(!complete_rows), prop_excluded = mean(!complete_rows),
      mean_distance = overall_mean,
      node_utilization = mean(som_node_occupancy$n > 0),
      empty_nodes = sum(som_node_occupancy$n == 0),
      singleton_nodes = sum(som_node_occupancy$n == 1)
    )
  )

  ProjectionReference$distance_reference <- .ClusterMetadataAliases(ProjectionReference$distance_reference)
  ProjectionReference$projection_thresholds <- .ClusterMetadataAliases(ProjectionReference$projection_thresholds)
  SOMFit <- list(
    train_quant_error = train_quant_err,
    distance_summary  = dist_summary,
    overall_mean      = overall_mean,
    overall_sd        = overall_sd,
    dist_by_cluster   = dist_by_cluster,
    flag_by_cluster   = flag_by_cluster,
    overall_p95       = overall_p95,
    node_occupancy    = som_node_occupancy,
    occupancy_summary = som_occupancy_summary,
    variable_contribution = som_variable_contribution,
    residuals         = som_residuals$residuals,
    table             = som_fit_tbl,
    plots             = list(distance_hist = p_dist_hist)
  )

  ModelInfo_SOM <- list(
    diagnostics = .SOMMapDiagnostics(som_model, zmat, node_df$Cluster, seed_som),
    som_model = som_model,
    som_codes = som_codes,
    som_grid  = som_grid,
    grid_screen = grid_screen,
    distance_scale = "euclidean",
    training_variable_summary = training_variable_summary,
    ProjectionReference = ProjectionReference,
    som_grid_info = list(
      som_xdim_initial = som_grid_initial_xdim,
      som_ydim_initial = som_grid_initial_ydim,
      som_xdim_used    = som_xdim,
      som_ydim_used    = som_ydim,
      n_nodes_used     = som_xdim * som_ydim,
      GridDesignPolicy = grid_design_policy,
      GridDesignK = grid_design_k,
      TargetNodesPerProfile = grid_nodes_per_profile,
      MinNodesPerParameter = min_nodes_per_parameter,
      MinParticipantsPerNode = min_participants_per_node,
      MixtureDimension = mixture_dimension,
      MaxMixtureParameters = max_mixture_parameters,
      ParameterFloorNodes = parameter_floor_nodes,
      ParameterFloorSide = parameter_floor_side,
      OccupancyCeilingSide = occupancy_ceiling_side,
      ParticipantsPerNode = n_complete / (som_xdim * som_ydim),
      OccupancyCompliant = n_complete / (som_xdim * som_ydim) >= min_participants_per_node,
      BindingConstraint = binding_constraint,
      MapFrom = !is.null(map_from),
      SizingIterations = sizing_iterations,
      ScreenFloorSide = if (is.null(grid_screen)) NA_integer_ else grid_screen$selection$floor,
      ScreenCeilingSide = if (is.null(grid_screen)) NA_integer_ else grid_screen$selection$ceiling,
      ScreenTargetSide = if (is.null(grid_screen)) NA_integer_ else grid_screen$selection$target,
      ExpectedNodesPerProfile = (som_xdim * som_ydim) / grid_design_k
    ),
    training = list(
      mode = "batch",
      epochs = som_training_epochs,
      seed = seed_som,
      initialization = "deterministic sampled training codebook",
      initial_codebook = init_codes,
      cohort = "complete training cohort"),
    SOMFit    = SOMFit,
    plots     = som_plots
  )

  # ProbFit diagnostics (probabilities) ------------------------------------

  node_non_na <- node_df

  # A node clustering that carries no posterior probabilities has nothing for
  # these two densities to show, and geom_density() on an all-NA column fails
  # at draw time rather than returning an empty panel.
  has_posterior <- any(is.finite(node_df$prob_assigned))

  node_ProbAssignedDensity <- if (!has_posterior) NULL else ggplot2::ggplot(
    node_non_na,
    ggplot2::aes(x = prob_assigned)
  ) +
    ggplot2::geom_density(fill = "grey40", alpha = 0.7) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Density of node prob_assigned by class",
      x     = "prob_assigned",
      y     = "Density"
    )

  indiv_non_na <- som_fit_tbl[!is.na(som_fit_tbl$max_prob), ]
  indiv_non_na$Cluster <- factor(indiv_non_na$Cluster)

  individual_ProbAssignedDensity <- if (!has_posterior) NULL else ggplot2::ggplot(
    indiv_non_na,
    ggplot2::aes(x = prob_assigned)
  ) +
    ggplot2::geom_density(fill = "darkblue", alpha = 0.7) +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Density of prob_assigned by class (training)",
      x     = "Posterior probability for assigned class",
      y     = "Density"
    )

  boundary_reason <- factor(ifelse(is.na(individual_tbl$BoundaryReason[complete_rows]),
    "not boundary", c(low_margin = "margin only", high_entropy = "entropy only",
      both = "both")[individual_tbl$BoundaryReason[complete_rows]]),
    levels = c("margin only", "entropy only", "both", "not boundary"))
  boundary_summary <- tibble::tibble(BoundaryReason = levels(boundary_reason),
    N = as.integer(table(boundary_reason))) %>%
    dplyr::mutate(Share = N / sum(N))
  boundary_threshold_grid <- expand.grid(margin = c(.10, .15, .20, .25, .30),
    entropy = c(.50, .60, .70)) %>% dplyr::as_tibble()
  boundary_threshold_grid$BoundaryPct <- vapply(seq_len(nrow(boundary_threshold_grid)), function(i) {
    mean(individual_tbl$SoftMargin[complete_rows] < boundary_threshold_grid$margin[i] |
      individual_tbl$SoftEntropy[complete_rows] > boundary_threshold_grid$entropy[i])
  }, numeric(1))
  ProbFit <- list(
    boundary_summary = boundary_summary,
    boundary_threshold_grid = boundary_threshold_grid,
    node       = node_df,
    individual = som_fit_tbl,
    plots      = list(
      node_ProbAssignedDensity       = node_ProbAssignedDensity,
      individual_ProbAssignedDensity = individual_ProbAssignedDensity
    )
  )

    lpa_preprocess <- list(
    n_rows_original          = nrow(df_scidr),
    n_complete_rows          = sum(complete_rows),
    n_excluded_rows          = sum(!complete_rows),
    prop_excluded_rows       = mean(!complete_rows),
    n_som_nodes              = nrow(som_codes),
    n_lpa_variables_original = ncol(som_codes),
    n_lpa_variables_used     = ncol(X),
    n_lpa_nodes_used         = length(lpa_nodes),
    min_nodes_per_parameter  = min_nodes_per_parameter,
    min_participants_per_node = min_participants_per_node,
    max_mixture_parameters   = max_mixture_parameters,
    parameter_floor_nodes    = parameter_floor_nodes,
      mixture_min_node_people           = mixture_min_node_people,
      node_weighting           = node_weighting,
      total_node_weight        = sum(lpa_node_weights),
      effective_node_weight    = if (node_weighting == "occupancy") sum(lpa_node_weights) else length(lpa_node_weights),
    dropped_zero_sd_vars     = dropped_lpa_vars,
    mixture_drop_zero_sd         = mixture_drop_zero_sd,
    mixture_zero_sd_tol          = mixture_zero_sd_tol,
    mixture_em_itmax             = mixture_em_itmax,
    mixture_em_tol               = mixture_em_tol,
    mixture_timeout_seconds      = mixture_timeout_seconds,
    skip_model_after_n_failures = skip_model_after_n_failures,
    slow_fit_seconds         = slow_fit_seconds,
    candidate_min_nodes_per_cluster    = candidate_min_nodes_per_cluster
  )

  lpa_preprocess <- .ClusterMetadataAliases(lpa_preprocess)
  ModelInfo_MClust <- list(
    lpa_models    = lpa_models,
    node_weighting = node_weighting,
    node_assignments = node_assignments,
    candidate_node_posteriors = if (exists("lpa_models_all")) {
      lapply(lpa_models_all, get_node_posteriors)
    } else {
      NULL
    },
    ahp_sensitivity = .SOMAHPWeightSensitivity(fit_table, composite_weights$weights,
      best_fit_name),
    fit_table     = fit_table,
    best_fit_name = best_fit_name,
    AHP         = list(
      ahp_index      = if ("ahp_index" %in% names(fit_table)) fit_table$ahp_index else NA_real_,
      ahp_best_row   = ahp_best_row,
      recommendation = recommendation_txt,
      weights        = composite_weights$weights,
      weight_method  = composite_weights$method,
      pairwise       = composite_weights$pairwise,
      consistency_ratio = composite_weights$consistency_ratio
    ),
    Stability = Stability,
    diagnostics = list(
      lpa_fit_diagnostics = lpa_diagnostics,
      failed_fits         = failed_fits,
      warning_fits        = warning_fits,
      lpa_preprocess      = lpa_preprocess
    )
  )

  if (!is.null(Stability)) Stability$plots <- .ClusterStabilityPlots(Stability)

  out <- list(
    method           = method,
    vars_used        = vars_used,
    ZScoreType       = ZScoreType,
    ZScoreObject     = ZScoreObject_used,
    ZScoreVars       = ZScoreVars_used,
    id_var           = id_var,
    ClusterVariableName      = ClusterVariableName,
    DataWithClusters = DataWithClusters,
    fit_plot         = fit_plot,
    ModelInfo        = ModelInfo_MClust,
    ModelInfo_SOM    = ModelInfo_SOM,
    ModelInfo_MClust = ModelInfo_MClust,
    ProbFit          = ProbFit,
    Stability        = Stability
  )

  out$Specification <- .ClusterSpecification(
    "SOM + Mclust", vars_used,
    list(seed_som = seed_som, seed_mixture = seed_mixture, seed_stability = seed_stability,
      stability_cores = stability_core_settings$requested,
      stability_cores_resolved = stability_core_settings$resolved),
    list(ZScoreType = ZScoreType, ZScoreObject = ZScoreObject_used,
      ZScoreVars = ZScoreVars_used, SOM = ModelInfo_SOM$som_grid_info,
      SOMTraining = ModelInfo_SOM$training, node_weighting = node_weighting),
    ModelInfo_MClust$fit_table,
    list(k = final_k, model = final_model), complete_rows,
    list(distance_metric = "Euclidean distance (Z-score units) to frozen SOM best-matching unit"))
  out$Specification$projection <- list(assignment = assignment, soft_n_nodes = soft_n_nodes, soft_kernel = soft_kernel, soft_bandwidth = soft_bandwidth,
    boundary_margin_threshold = boundary_margin_threshold,
    boundary_entropy_threshold = boundary_entropy_threshold)
  out$Specification$stability <- if (is.null(Stability)) NULL else Stability$settings
  out$Specification$selection <- list(criteria_weights = composite_weights,
    min_cluster_prop = min_cluster_prop, mixture_min_node_people = mixture_min_node_people,
    min_nodes_per_parameter = min_nodes_per_parameter,
    min_participants_per_node = min_participants_per_node,
    map_from = !is.null(map_from))

  out$Specification$selection <- .ClusterMetadataAliases(out$Specification$selection)
  class(out) <- c("Pipeline_SOM_MClust", class(out))
  out
}

#' @description Readable SOM + Mclust workflow wrapper for
#'   [CreateClusterModel_SOM_MClust()].
#' @rdname CreateClusterModel_SOM_MClust
#' @param ... Arguments passed to [CreateClusterModel_SOM_MClust()].
#' @export
Pipeline_SOM_MClust <- function(...) {
  CreateClusterModel_SOM_MClust(...)
}

#' @description Compatibility wrapper for [Pipeline_SOM_MClust()].
#' @rdname CreateClusterModel_SOM_MClust
#' @export
Pipeline_SOMClust <- function(...) {
  Pipeline_SOM_MClust(...)
}

#' @description Deprecated alias for [CreateClusterModel_SOM_MClust()].
#' @rdname CreateClusterModel_SOM_MClust
#' @export
CreateSOMClusterModel <- function(...) {
  lifecycle::deprecate_warn("20.25.0", "CreateSOMClusterModel()",
    "CreateClusterModel_SOM_MClust()")
  CreateClusterModel_SOM_MClust(...)
}

# Compare recommendations using the same z-score direction and eligibility rule.
.SOMAHPWeightSensitivity <- function(df_Fits, user_weights, recommendation) {
  raw_criteria <- c(AIC = "AIC", BIC = "BIC", Entropy = "Entropy",
    Reproducibility = "ReproducibilityScore", Posterior = "PosteriorMean", Boundary = "BoundaryPct")
  for (criterion in names(raw_criteria)) {
    column <- raw_criteria[[criterion]]
    scaled <- paste0(criterion, "_scaled")
    if (!scaled %in% names(df_Fits) && column %in% names(df_Fits)) {
      values <- df_Fits[[column]]
      if (criterion %in% c("AIC", "BIC", "Boundary")) values <- -values
      df_Fits[[scaled]] <- if (sum(is.finite(values)) < 2 || stats::sd(values, na.rm = TRUE) == 0) {
        rep(0, length(values))
      } else as.numeric(scale(values))
    }
  }
  criteria <- intersect(c("AIC", "BIC", "Entropy", "Reproducibility", "Posterior", "Boundary"),
    sub("_scaled$", "", grep("_scaled$", names(df_Fits), value = TRUE)))
  if (!length(criteria)) return(tibble::tibble())
  schemes <- list(equal = stats::setNames(rep(1, length(criteria)), criteria))
  for (criterion in criteria) schemes[[paste0("only_", criterion)]] <- stats::setNames(1, criterion)
  if (!is.null(user_weights)) schemes$user <- user_weights
  base <- if (is.null(user_weights)) stats::setNames(rep(1,
    length(setdiff(criteria, c("Boundary", "Posterior")))),
    setdiff(criteria, c("Boundary", "Posterior"))) else user_weights
  schemes$without_Boundary <- base[names(base) != "Boundary"]
  schemes$with_Boundary <- c(base[names(base) != "Boundary"], Boundary = 1)
  dplyr::bind_rows(lapply(names(schemes), function(scheme) {
    weights <- schemes[[scheme]]
    if (!length(weights) || sum(weights) <= 0) {
      return(tibble::tibble(scheme = scheme, weights = list(weights),
        top_model = NA_character_, top_k = NA_integer_, recommendation_rank = NA_integer_))
    }
    weights <- weights / sum(weights)
    score <- as.numeric(as.matrix(df_Fits[paste0(names(weights), "_scaled")]) %*% weights)
    eligible <- if ("Eligible" %in% names(df_Fits)) df_Fits$Eligible else rep(TRUE, nrow(df_Fits))
    ranked <- which(eligible & is.finite(score))
    ranked <- ranked[order(-score[ranked])]
    ids <- paste0("model_", df_Fits$Model, "_class_", df_Fits$Classes)
    top <- if (length(ranked)) ranked[1] else NA_integer_
    tibble::tibble(scheme = scheme, weights = list(weights),
      top_model = unname(.MclustModelNames[as.character(df_Fits$Model[top])]),
      top_k = df_Fits$Classes[top], recommendation_rank = match(recommendation, ids[ranked]))
  }))
}

.SOMMapDiagnostics <- function(som_model, zmat, clusters, seed_som) {
  quality <- .SOMGridMapQuality(som_model, zmat, low_density_n = 3L)
  codes <- as.matrix(if (is.list(som_model$codes)) som_model$codes[[1]] else som_model$codes)
  hits <- tabulate(quality$bmu, nrow(codes))
  adjacent <- as.matrix(stats::dist(som_model$grid$pts)) <= 1.01
  diag(adjacent) <- FALSE
  distances <- as.matrix(stats::dist(codes))
  u <- rowSums(distances * adjacent) / rowSums(adjacent)
  BorderDifference <- function(labels) {
    border <- vapply(seq_along(labels), function(i) any(labels[adjacent[i, ]] != labels[i]), logical(1))
    c(border_mean = if (any(border)) mean(u[border]) else NA_real_,
      interior_mean = if (any(!border)) mean(u[!border]) else NA_real_)
  }
  observed <- BorderDifference(clusters)
  difference <- unname(observed[1] - observed[2])
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed_som)
  permuted <- replicate(999L, {
    means <- BorderDifference(sample(clusters))
    unname(means[1] - means[2])
  })
  # No interior means no border contrast; report it rather than inventing a test.
  p_value <- if (is.finite(difference) && any(is.finite(permuted))) {
    (1 + sum(permuted >= difference, na.rm = TRUE)) / (1 + sum(is.finite(permuted)))
  } else NA_real_
  map_quality <- tibble::as_tibble(as.list(quality$metrics))
  map_quality$EmptyNodeCount <- sum(hits == 0L)
  map_quality$OnePersonNodeCount <- sum(hits == 1L)
  map_quality$FewerThan3NodeCount <- sum(hits < 3L)
  map_quality$FewerThan3NodeShare <- mean(hits < 3L)
  list(map_quality = map_quality,
    u_matrix = tibble::tibble(NodeID = seq_along(u), U = u),
    border_alignment = tibble::tibble(BorderMeanU = unname(observed[1]),
      InteriorMeanU = unname(observed[2]), Difference = difference,
      PermutationP = p_value, Permutations = 999L,
      ValidPermutations = sum(is.finite(permuted))))
}
