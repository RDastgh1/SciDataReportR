---
title: "Clustering and clinical phenotyping"
author: "Raha Dastgheyb"
date: "last-modified"
bibliography: ../inst/REFERENCES.bib
format:
  html:
    minimal: true
    embed-resources: true
    toc: true
    toc-location: left
    smooth-scroll: true
execute:
  echo: false
  message: false
  warning: false
  cache: false
vignette: >
  %\VignetteIndexEntry{Clustering and clinical phenotyping}
  %\VignetteEngine{quarto::html}
  %\VignetteEncoding{UTF-8}
---

<style>
table {
  width: 100%;
  table-layout: fixed;
  font-size: 0.68rem;
  line-height: 1.2;
}
table th, table td {
  padding: 0.35rem 0.45rem;
  overflow-wrap: anywhere;
}
</style>

# Purpose

Clinical phenotyping uses prespecified features to discover structure in a
training cohort, then places later participants in that **fixed** phenotype
space. Clustering is exploratory: fit indices, internal stability, and clean
projection do not establish that a phenotype is biologically real. Clinical
interpretation and independent-cohort validation remain essential.

Every pipeline in this package uses the same train-once/project-many contract:

- preprocessing and model parameters are learned from training data only;
- `DataWithClusters` preserves the original row order and columns;
- incomplete cases remain in the output with `NA` assignments;
- bootstrap stability refits the **entire pipeline**, including scaling and
  PCA/MCA when present;
- fit, membership uncertainty or distance, and projection support are
  reviewed before phenotypes are named.

`method = "exploratory"` (or the shorter alias `"explore"`) compares candidate
settings and gives an advisory AHP-style recommendation. The researcher then
creates an explicit `method = "finalize"` model for projection.

## Breaking output-name changes

The unreleased clustering interface now uses `DataWithClusters` rather than
`df_with_clusters`, and `ClusterVariableName` rather than `ClusterName`.
`ModelInfo_MClust` is the only correctly cased Mclust-specific field; the
erroneous `ModelInfo_Mclust` field has been removed. `complete_rows`,
`CandidateAudit`, and internally generated `.scidr_rowid` columns are no
longer public outputs. Use `Specification$n_complete` for aggregate training
support and the `fit_table` status/warning/error columns to audit candidates.

Every method stores its figures in the same places, so the same review sequence
works whichever approach you choose. A figure lives beside the table it
summarises rather than in one flat list:

| Location | What it shows |
|---|---|
| `fit_plot` | candidate solutions across the search grid |
| `ModelInfo$plots` | structure of the selected solution |
| `ModelInfo$FitDiagnostics$plots` | how individual training cases sit inside that structure |
| `ProbFit$plots` | membership confidence |
| `Stability$plots` | bootstrap reproducibility |
| `ProjectionFit$plots` | projected cases against the frozen training reference |

The *contents* are method-specific: K-means and Gower/PAM give you a silhouette
profile, Mclust a BIC curve and an uncertainty map, HDBSCAN a persistence
review, LCA item response profiles, and the reduction pipelines a scree plot.
The SOM pipeline keeps its own richer set, including the interactive `aweSOM`
maps in `ModelInfo_SOM$plots`.

# Choosing a method

| Method | Data and geometry | Primary diagnostics | Projection rule | Ideal situation | Main caveat |
|---|---|---|---|---|---|
| SOM + Mclust | Continuous; topology and gradients matter | map fit, posterior, SOM distance, stability | scaling → SOM node → mixture class | phenotype maps and nonlinear gradients | most complex workflow |
| Mclust | Continuous; Gaussian mixtures | BIC/ICL, entropy, posterior, stability | exact mixture posterior | probabilistic groups with different covariance | sensitive to distributional mismatch |
| PCA + Mclust | Correlated continuous measures | scree, mixture fit, posterior, stability | frozen PCA → mixture posterior | clustering a lower-dimensional signal space | PCs can obscure individual variables |
| K-means | Numeric; compact centroid geometry | WSS, silhouette, CH, margin, stability | nearest frozen centroid | simple, interpretable baseline | no probability model |
| PCA + K-means | Correlated numeric measures | scree, silhouette, distance, stability | frozen PCA → nearest centroid | fast reduced-space baseline | inherits PCA and K-means assumptions |
| HDBSCAN | Numeric; irregular or unequal-density groups | persistence, membership, noise, support | nearest eligible training core | clusters plus meaningful noise | conservative projection may label many cases noise |
| Gower + PAM | Mixed numeric/categorical/ordinal | silhouette, medoids, distance, stability | nearest frozen medoid | mixed clinical measures and exemplar participants | coding and variable selection define distance |
| LCA | Categorical response patterns | BIC, entropy, response probabilities, stability | latent-class posterior | symptom or questionnaire response phenotypes | local-independence assumption |
| MCA + Mclust | Nominal categorical measures | inertia, mixture fit, posterior, stability | frozen MCA → mixture posterior | reduced categorical response space | MCA dimensions can be abstract |

# Load the benchmark


::: {.cell}

:::



::: {.cell}

:::


::: callout-note
The data file is `SimulatedPhenotypeData`; the codebook is
`SimulatedPhenotypeVariableTypes`. The benchmark contains 320 training and 160
projection participants. `TruthCluster`, `TruthDensityGroup`, and `Cohort` are
never supplied to a clustering function. Truth is revealed only as a teaching
check after fitting.
:::


::: {.cell}

:::


The 12 numeric measures are three indicators of each of four independent
simulated axes. The default 85% rule therefore retains at least four PCs. Four
balanced truth groups use different combinations of these axes. Separate
`DensityX`/`DensityY` coordinates contain two unequal-density groups and noise;
three categorical variables contain deliberately strong class-dependent
response patterns. That separation is useful for checking implementation; it
is intentionally easier than most real categorical phenotyping problems.

# Constructor interface

All clustering constructors use `method = "exploratory"` to review candidates
and `method = "finalize"` to create a reusable model. Fixed-count methods use
`k_range` and `final_k`; density methods use `minPts_range`, epsilon settings,
and their matching finalized inputs. Numeric methods use `ZScoreType`; `Scaling`
is retained as a compatibility alias and must agree if both are supplied.

| Workflow | Candidate settings | Finalized settings | Method-specific control |
|---|---|---|---|
| Mclust, PCA + Mclust, MCA + Mclust, SOM + Mclust | `k_range`, `models` | `final_k`, `final_model` | `models`: 1 = EEI, 2 = VVI, 3 = EEE |
| K-means, PCA + K-means, Gower + PAM, latent class | `k_range` | `final_k` | `nstart` for K-means; `nrep` for LCA |
| HDBSCAN, SOM + HDBSCAN | `minPts_range`, epsilon range | `final_minPts`, final epsilon | extracted class count is data-derived |

# Continuous workflows

## Mclust

Gaussian mixture models provide posterior membership and allow clusters to
differ in covariance [@scrucca2016].

::: callout-note
**When to use it:** continuous measures with approximately mixture-shaped
groups. **Inspect before naming:** BIC/ICL, entropy, posterior confidence,
small clusters, bootstrap ARI/Jaccard, and projected uncertainty.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 11
  Model ModelName Classes     BIC     ICL Entropy MinClusterN SizeBalance
  <int> <chr>       <int>   <dbl>   <dbl>   <dbl>       <int>       <dbl>
1     1 EEI            10  -5756.  -5771.   0.975           8       0.138
2     1 EEI             9  -5779.  -5793.   0.976           8       0.138
3     1 EEI             8  -5790.  -5803.   0.972          21       0.356
4     1 EEI             7  -5908.  -5923.   0.967          22       0.275
5     1 EEI             6  -5980.  -5995.   0.968          22       0.275
6     1 EEI             4  -6009.  -6009.   1.000          80       1    
7     1 EEI             5  -6019.  -6034.   0.967          36       0.45 
8     1 EEI             3  -8630.  -8630.   1.000          80       0.5  
9     1 EEI             2 -10446. -10446.   0.995          80       0.333
# ℹ 3 more variables: ReproducibilityScore <dbl>, ahp_index <dbl>,
#   Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends Gaussian-mixture model (Model = 1, Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustReview-7.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MclustFinalize-8.png){width=672}
:::
:::


## K-means

K-means is a useful scaled baseline for compact groups [@macqueen1967].

::: callout-note
**When to use it:** continuous features with roughly compact, similarly shaped
groups. **Inspect before naming:** WSS, silhouette, Calinski--Harabasz index,
centroid distance, assignment margin, and bootstrap stability.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 8
  Classes   WSS Silhouette CalinskiHarabasz MinClusterN ReproducibilityScore
    <int> <dbl>      <dbl>            <dbl>       <int>                <dbl>
1       2 2914.      0.272             99.8          80                0.174
2       3 2073.      0.373            134.           80                0.783
3       4 1264.      0.490            214.           80                1    
4       5 1118.      0.453            191.           25                0.884
5       6  977.      0.416            183.           25                0.772
6       7  843.      0.382            185.           25                0.748
7       8  716.      0.345            194.           23                0.794
8       9  668.      0.332            184.           14                0.753
9      10  625.      0.313            177.           14                0.859
# ℹ 2 more variables: ahp_index <dbl>, Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends K-means k (Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansReview-7.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-8.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/KMeansFinalize-9.png){width=672}
:::
:::


## PCA followed by clustering

PCA concentrates correlated variation into orthogonal dimensions
[@jolliffe2016]. These pipelines refit scaling and PCA inside every stability
bootstrap, then use the frozen training transformation for projection.

::: callout-note
**When to use it:** correlated continuous variables where a reduced signal
space is clearer than the raw feature space. **Inspect before naming:** scree
plot, retained component count/loadings, downstream fit and stability, score
space separation, and projection distance or uncertainty.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 10
  Model ModelName Classes    BIC    ICL Entropy MinClusterN ReproducibilityScore
  <int> <chr>       <int>  <dbl>  <dbl>   <dbl>       <int>                <dbl>
1     1 EEI             4 -2100. -2100.   1.000          80                1    
2     1 EEI             5 -2118. -2150.   0.942          22                0.886
3     1 EEI             6 -2127. -2162.   0.939           5                0.735
4     1 EEI             7 -2141. -2199.   0.906           5                0.758
5     1 EEI             8 -2146. -2219.   0.887           9                0.890
6     1 EEI             9 -2173. -2277.   0.857           9                0.645
7     1 EEI            10 -2200. -2326.   0.836           7                0.652
8     1 EEI             3 -3105. -3106.   0.991          80                0.565
9     1 EEI             2 -3676. -3720.   0.740          80                0.587
# ℹ 2 more variables: ahp_index <dbl>, Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends PCA plus Gaussian-mixture model (Model = 1, Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] 4
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustReview-7.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAMclustFinalize-8.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 8
  Classes   WSS Silhouette CalinskiHarabasz MinClusterN ReproducibilityScore
    <int> <dbl>      <dbl>            <dbl>       <int>                <dbl>
1       2  977.      0.282             97.5          80                0.416
2       3  679.      0.407            139.           80                0.482
3       4  383.      0.543            246.           80                1    
4       5  328.      0.513            227.           24                0.892
5       6  277.      0.479            227.           24                0.785
6       7  226.      0.466            242.           24                0.964
7       8  178.      0.444            274.           23                0.948
8       9  162.      0.432            268.           13                0.801
9      10  144.      0.423            270.           13                0.760
# ℹ 2 more variables: ahp_index <dbl>, Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends PCA plus K-means k (Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] 4
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansReview-7.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-8.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PCAKMeansFinalize-9.png){width=672}
:::
:::


## HDBSCAN

HDBSCAN is designed for irregular clusters, unequal density, and noise
[@mcinnes2017]. A projected `0` is an out-of-support/noise assignment, not a
new phenotype.

::: callout-note
**When to use it:** density structure is more plausible than centroid geometry
and noise is scientifically meaningful. **Inspect before naming:** persistence,
membership, noise burden, outlier scores, stability, and nearest-core support.
All-noise candidates are valid negative results.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 6 × 10
  Classes MinPts Epsilon Persistence MeanMembershipProbability NoiseProportion
    <int>  <dbl>   <dbl>       <dbl>                     <dbl>           <dbl>
1       3      6    0           774.                     0.710          0.0688
2       2     10    0          1112.                     0.654          0.134 
3       2     14    0           911.                     0.637          0.141 
4       3      6    0.05        774.                     0.710          0.0688
5       2     10    0.05       1112.                     0.654          0.134 
6       2     14    0.05        911.                     0.637          0.141 
# ℹ 4 more variables: MinClusterN <int>, ReproducibilityScore <dbl>,
#   ahp_index <dbl>, Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends HDBSCAN density setting (Classes = 2, MinPts = 14, Epsilon = 0.05). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANReview-6.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-8.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-9.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/HDBSCANFinalize-10.png){width=672}
:::
:::


## SOM followed by Mclust

Self-organizing maps preserve local topology before mixture clustering
[@kohonen1982].

::: callout-note
**When to use it:** a map of phenotype gradients is as useful as the final
labels. **Inspect before naming:** map topology, quantization/SOM distance,
node occupancy, posterior confidence, stability, and projected distance flags.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
  Model Classes     LogLik parameters  n       AIC       AWE      BIC     CAIC
1     1       2 -329.58579         37 25 733.17158 1006.5018 778.2700 815.2700
2     1       3 -287.80626         50 25 675.61252 1045.5632 736.5563 786.5563
3     1       4 -182.88068         63 25 491.76135  958.3427 568.5505 631.5505
4     1       5 -140.28938         76 25 432.57876  995.8480 525.2133 601.2133
5     1       6 -104.86298         89 25 387.72596 1047.6860 496.2059 585.2059
6     1       7  -68.97640        102 25 341.95280 1098.6035 466.2781 568.2781
7     1       8   92.86636        115 25  44.26728  897.6087 184.4380 299.4380
8     1       9   86.29563        128 25  83.40875 1033.4415 239.4249 367.4249
9     1      10  105.45796        141 25  71.08408 1117.8079 242.9456 383.9456
        CLC      KIC     SABIC       ICL   Entropy  prob_min prob_max
1  661.0382 773.1716  663.5296 -779.1319 0.9332915 0.9301147 0.998899
2  577.5494 728.6125  581.5017 -737.0754 0.9684421 0.9582059 0.999478
3  367.7584 557.7614  373.1817 -568.5659 0.9985197 0.9989606 1.000000
4  282.5787 511.5788  289.5303 -525.2138 0.9999468 0.9999677 1.000000
5  211.7258 479.7260  220.2087 -496.2064 0.9999445 0.9999460 1.000000
6  139.9528 446.9528  149.9667 -466.2782 0.9999983 0.9999985 1.000000
7 -183.7327 162.2673 -172.1877 -184.4380 1.0000000 1.0000000 1.000000
8 -170.5918 214.4087 -157.5150 -239.4288 0.9997308 0.9980997 1.000000
9 -208.9167 215.0841 -194.3085 -242.9518 0.9995963 0.9978376 1.000000
  BLRTStatistic BLRTPValue MinProfileNodeN MaxProfileNodeN
1      54.51982 0.00990099               5              20
2      83.55906 0.00990099               5              15
3     209.85117 0.00990099               5               8
4      85.18259 0.00990099               4               7
5      70.85280 0.00990099               2               6
6      71.77316 0.00990099               2               6
7     323.68552 0.00990099               1               5
8     -13.14147 1.00000000               1               4
9      38.32467 0.00990099               1               5
  MinProfileNodeProportion MaxProfileNodeProportion StabilitySuccessRate
1                     0.20                     0.80                    1
2                     0.20                     0.60                    1
3                     0.20                     0.32                    1
4                     0.16                     0.28                    1
5                     0.08                     0.24                    1
6                     0.08                     0.24                    1
7                     0.04                     0.20                    1
8                     0.04                     0.16                    1
9                     0.04                     0.20                    1
  StabilityARI_Mean StabilityARI_P05 StabilityJaccard_Mean StabilityJaccard_Min
1         0.4082278      -0.02926955             0.6722532           0.60126050
2         0.4576730       0.45465979             0.6625065           0.49689441
3         0.9916404       0.99164035             0.9937886           0.98750000
4         0.9855618       0.98012283             0.8090879           0.08950617
5         0.9025672       0.88134138             0.6409295           0.01234568
6         0.8579067       0.85590269             0.5675309           0.01282051
7         0.6903812       0.64346240             0.4655000           0.02282051
8         0.7664925       0.62412601             0.5173137           0.01808996
9         0.6573845       0.64362723             0.5162323           0.01474359
  ReproducibilityScore  AIC_scaled  BIC_scaled Entropy_scaled
1            0.5402405 -1.45123232 -1.43753336     -2.3884933
2            0.5600897 -1.22594363 -1.24245873     -0.8768048
3            0.9927145 -0.50634208 -0.45677714      0.4167145
4            0.8973249 -0.27469881 -0.25411004      0.4780891
5            0.7717484 -0.09914298 -0.11845642      0.4779903
6            0.7127188  0.08001517  0.02150125      0.4803067
7            0.5779406  1.24516939  1.33953094      0.4803777
8            0.6419031  1.09196797  1.08238409      0.4688018
9            0.5868084  1.14020727  1.06591940      0.4630181
  Reproducibility_scaled  ahp_index
1            -0.98258933 -1.5649621
2            -0.85891595 -1.0510308
3             1.83661466  0.3225525
4             1.24227619  0.2978891
5             0.45985369  0.1800611
6             0.09206148  0.1684711
7            -0.74769334  0.5793462
8            -0.34916598  0.5734970
9            -0.69244142  0.4941758
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP (AIC, BIC, Entropy, reproducibility) recommends Model 1 with k = 8 profiles."
```


:::

::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 8
  Model Classes StabilitySuccessRate StabilityARI_Mean StabilityARI_P05
  <int>   <int>                <dbl>             <dbl>            <dbl>
1     1       2                    1             0.408          -0.0293
2     1       3                    1             0.458           0.455 
3     1       4                    1             0.992           0.992 
4     1       5                    1             0.986           0.980 
5     1       6                    1             0.903           0.881 
6     1       7                    1             0.858           0.856 
7     1       8                    1             0.690           0.643 
8     1       9                    1             0.766           0.624 
9     1      10                    1             0.657           0.644 
# ℹ 3 more variables: StabilityJaccard_Mean <dbl>, StabilityJaccard_Min <dbl>,
#   ReproducibilityScore <dbl>
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMReview-1.png){width=672}
:::

::: {.cell-output-display}

```{=html}
<h4 id="aweSOMwidget-bf70b9e3242e59d606c9-info"></h4>
<h4 id="aweSOMwidget-bf70b9e3242e59d606c9-message"></h4>
<div id="aweSOMwidget-bf70b9e3242e59d606c9" class="aweSOMwidget html-widget" style="width:400px;height:auto;display:block; margin:auto; margin-top:5px; margin-bottom:5px;"></div>
<p id="aweSOMwidget-bf70b9e3242e59d606c9-names"></p>
<svg id="aweSOMwidget-bf70b9e3242e59d606c9-placeHolder"></svg>
<script type="application/json" data-for="aweSOMwidget-bf70b9e3242e59d606c9">{"x":{"plotType":"Circular","sizeInfo":400,"gridInfo":{"nbLines":5,"nbColumns":5,"topology":"hexagonal"},"superclass":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"superclassColor":["#8DD3C7"],"cellNames":["245, 249, 256, 266, 272, 276, 278, 283, 287, 306, 311, 317, 318, 320","241, 246, 248, 251, 253, 255, 257, 258, 263, 273, 277, 279, 282, 285, 289, 290, 294, 300, 301, 305, 308, 309, 310, 312, 315","86","83, 85, 88, 91, 96, 98, 101, 103, 107, 111, 112, 115, 119, 120, 126, 139, 142, 143, 144, 148, 156, 158, 159","87, 90, 93, 94, 102, 104, 105, 108, 110, 114, 116, 121, 124, 125, 127, 135, 137, 145, 146, 151, 160","242, 244, 250, 254, 260, 262, 265, 268, 271, 280, 281, 286, 288, 292, 295, 303, 313, 316, 319","247, 252, 259, 261, 264, 267, 269, 270, 284, 293, 298, 302, 307, 314","243, 274, 275, 291, 296, 297, 299, 304","81, 89, 97, 99, 106, 118, 128, 131, 132, 133, 134, 136, 138, 140, 141, 149, 150, 153, 154, 157","82, 84, 92, 95, 100, 109, 113, 117, 122, 123, 129, 130, 147, 152, 155","","","","55","","164, 170, 172, 177, 179, 180, 189, 192, 195, 199, 207, 208, 209, 212, 213, 214, 222, 226, 227, 237, 238","165, 171, 174, 181, 187, 193, 202, 204, 215, 217, 220, 230, 234, 236","168, 176, 188","10, 11, 15, 17, 18, 20, 22, 23, 24, 25, 31, 35, 36, 38, 42, 46, 48, 49, 54, 56, 59, 60, 63, 66, 69, 75, 77, 79","2, 5, 7, 12, 14, 27, 29, 30, 34, 37, 40, 43, 64, 65, 67, 71, 74","162, 163, 166, 169, 184, 186, 194, 196, 198, 200, 203, 221, 223, 224, 225, 240","161, 167, 173, 175, 178, 182, 183, 185, 190, 191, 197, 201, 205, 206, 210, 211, 216, 218, 219, 228, 229, 231, 232, 233, 235, 239","51, 57","3, 16, 28, 44, 45, 47, 50, 52, 53, 70, 72, 76","1, 4, 6, 8, 9, 13, 19, 21, 26, 32, 33, 39, 41, 58, 61, 62, 68, 73, 78, 80"],"clustering":[25,20,24,25,20,25,20,25,25,19,19,20,25,20,19,24,19,19,25,19,25,19,19,19,19,25,20,24,20,20,19,25,25,20,19,19,20,19,25,20,25,19,20,24,24,19,24,19,19,24,23,24,24,19,14,19,23,25,19,19,25,25,19,20,20,19,20,25,19,24,20,24,25,20,19,24,19,25,19,25,9,10,4,10,4,3,5,4,9,5,4,10,5,5,10,4,9,4,9,10,4,5,4,5,5,9,4,5,10,5,4,4,10,5,4,5,10,9,4,4,5,10,10,5,5,4,5,9,10,10,9,9,9,9,5,9,5,9,4,9,9,4,4,4,5,5,10,4,9,9,5,10,9,9,10,4,9,4,4,5,22,21,21,16,17,21,22,18,21,16,17,16,22,17,22,18,16,22,16,16,17,22,22,21,22,21,17,18,16,22,22,16,17,21,16,21,22,21,16,21,22,17,21,17,22,22,16,16,16,22,22,16,16,16,17,22,17,22,22,17,21,16,21,21,21,16,16,22,22,17,22,22,22,17,22,17,16,16,22,21,2,6,8,6,1,2,7,2,1,6,2,7,2,6,2,1,2,2,7,6,7,6,2,7,6,1,7,6,7,7,6,1,2,8,8,1,2,1,2,6,6,2,1,7,2,6,1,6,2,2,8,6,7,2,6,8,8,7,8,2,2,7,6,8,2,1,7,2,2,2,1,2,6,7,2,6,1,1,6,1],"cellPop":[14,25,1,23,21,19,14,8,20,15,0,0,0,1,0,21,14,3,28,17,16,26,2,12,20],"showAxes":true,"transparency":true,"showNames":true,"legendPos":"beside","legendFontsize":14,"legendReverse":false,"label":["Z_Var1","Z_Var2","Z_Var3","Z_Var4","Z_Var5","Z_Var6","Z_Var7","Z_Var8","Z_Var9","Z_Var10","Z_Var11","Z_Var12"],"normalizedValues":[[0.05,0.09845135479254073,0.06312558905577184,0.06804589402451992,0.05,0.09562431447933992,0.9500000000000001,0.8719523815407212,0.9500000000000001,0.6125837584075988,0.6574812705985936,0.6189469366040105],[0.08857759008163235,0.105907117310156,0.05707779430786286,0.1164670975411831,0.1365144712177756,0.158128292275026,0.8693835944789301,0.9167990424720932,0.8142688420065105,0.4165458030197262,0.4067459306473301,0.4540418841331614],[0.5127482286954579,0.6944997552058161,0.7116845081424661,0.1776252714038998,0.2271494354616294,0.05,0.2904353973865388,0.09658671886773591,0.4386978566656336,0.1996928641168658,0.258927868858518,0.2240552363030068],[0.8463174876607997,0.7913454155652869,0.8137645060098385,0.08947558292678637,0.06312166436245348,0.06754147867233545,0.1785922754420079,0.05,0.08698476912207531,0.5049428448424013,0.5143411172019526,0.4741359320068767],[0.899666242582813,0.7494631208952575,0.800118658133368,0.09561683296522269,0.08329863411359048,0.117481214764859,0.1484538116315015,0.08173531005936152,0.1458747283361189,0.7422282108380349,0.7324091182742033,0.7795045842836833],[0.1723642461609303,0.15042448426648,0.09931246026264401,0.08693921281492301,0.1543322456186746,0.1148290262470307,0.8478759077654557,0.9129889488103128,0.8222246510744543,0.9500000000000001,0.9500000000000001,0.9500000000000001],[0.2200473668688029,0.1394343422715492,0.1984524160897938,0.09445511244572544,0.1557004848104969,0.1556964000138447,0.7455441604210447,0.8183432579098927,0.7815195911018031,0.6399143069541775,0.5995028462650621,0.572473244497671],[0.1566290131357191,0.1512547317439397,0.1083330581265388,0.1241016704448015,0.1504020253645396,0.146011049938514,0.7458783734097891,0.8242889528052476,0.8574662665711995,0.07062022119489191,0.1712498728175972,0.05032324957055508],[0.8928488048577458,0.7808596638964286,0.7670626007377043,0.1705459069395758,0.07163674373181438,0.09877259363757007,0.1209272091010103,0.1160784369405797,0.1344285104547597,0.05,0.05,0.05],[0.9390495977925585,0.814456653110988,0.930119813883826,0.05,0.06879258897071303,0.1142208005674119,0.05,0.1546527658267469,0.1637454554858418,0.3495911243759209,0.3530346031869097,0.4379888502891326],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[0.9500000000000001,0.9500000000000001,0.9500000000000001,0.7334793096612111,0.7590438378361971,0.6382314534017896,0.5989851653526258,0.502513047210761,0.9326748438819138,0.2773860674543602,0.3676331934633182,0.4261300789288175],[null,null,null,null,null,null,null,null,null,null,null,null],[0.107076470223123,0.1301933997994004,0.09069392929784147,0.8184476193155803,0.8085573508692021,0.7610409720035139,0.1529619772292712,0.06129777498029563,0.0556583104542358,0.8577754136127497,0.8206888711124233,0.8288170311681037],[0.1868254014170879,0.1374715803409088,0.07585036277110642,0.788165110629112,0.8976602468275806,0.8242063136666126,0.1840499175612459,0.2222890566769031,0.1960544180205192,0.6229131171744096,0.6061827507297842,0.5290955601486533],[0.2942152035592111,0.3071597250249885,0.4006346796841981,0.8764287861711313,0.9353679015032749,0.8082692254396113,0.3628598220001634,0.2351249015335853,0.3589330547006762,0.3080830218281665,0.3139747191932172,0.2693725486356411],[0.8202323809911798,0.729317630020988,0.811246449889903,0.7989183476648734,0.865801378745682,0.8185196446774852,0.8564631857274251,0.8828949747052077,0.9063217062659349,0.2341240322158078,0.1860005441047004,0.1623954148160907],[0.9068247409097622,0.7914930450640806,0.8332766236584539,0.7664472706110095,0.8161305648171556,0.7707450512339071,0.8388396191683064,0.8109263905019205,0.9016572166628178,0.5881482175986895,0.5513059513580776,0.5416143170181472],[0.07296095597137137,0.05,0.05,0.7324965704546614,0.8166365242714549,0.7485998917095468,0.1958947583535889,0.05105847956324232,0.05,0.4990938049126612,0.4927145127169753,0.5227316327582177],[0.1274149126551875,0.1227519273249815,0.08821069320029765,0.8146231109471999,0.8108189499797752,0.8281801856145649,0.1889102534326543,0.08020487648742039,0.1570414576840522,0.2183949517643701,0.2243559715027022,0.1594628718621929],[0.6085180307050093,0.6058139507770435,0.660233094764862,0.9500000000000001,0.9109032877171012,0.9500000000000001,0.7329504375955722,0.7423822176330352,0.6291617524753448,0.3209524533113836,0.2768158175039977,0.2963113639306945],[0.8599176731012064,0.739331789850999,0.8686203252012453,0.807711451749885,0.9500000000000001,0.8049947441415318,0.8776792008696508,0.9500000000000001,0.8285603010944435,0.4534310214800714,0.4547727977579033,0.4310136681163551],[0.8782461014080221,0.7305742406758964,0.7969968769385186,0.7997929674667038,0.9033880892655833,0.837145979950884,0.8279572375090215,0.9252489469209939,0.861458630972154,0.8407296949080898,0.7918050008929475,0.7867269957170437]],"realValues":[[-1.134,-0.999,-1.016,-1.034,-1.059,-0.993,1.235,0.916,1.169,0.376,0.591,0.456],[-1.036,-0.977,-1.031,-0.903,-0.845,-0.822,1.013,1.023,0.822,-0.303,-0.316,-0.108],[0.032,0.731,0.658,-0.737,-0.621,-1.118,-0.582,-0.9389999999999999,-0.138,-1.054,-0.851,-0.894],[0.873,1.012,0.921,-0.976,-1.027,-1.07,-0.89,-1.051,-1.037,0.003,0.073,-0.039],[1.007,0.891,0.886,-0.96,-0.977,-0.9330000000000001,-0.973,-0.975,-0.887,0.826,0.862,1.005],[-0.825,-0.848,-0.922,-0.983,-0.801,-0.9409999999999999,0.953,1.014,0.842,1.546,1.649,1.587],[-0.705,-0.88,-0.666,-0.963,-0.798,-0.829,0.672,0.787,0.738,0.471,0.381,0.297],[-0.865,-0.845,-0.899,-0.882,-0.8110000000000001,-0.855,0.672,0.802,0.9330000000000001,-1.501,-1.168,-1.487],[0.99,0.982,0.801,-0.756,-1.006,-0.984,-1.049,-0.893,-0.916,-1.573,-1.606,-1.488],[1.106,1.079,1.222,-1.083,-1.013,-0.9419999999999999,-1.244,-0.8,-0.841,-0.535,-0.51,-0.162],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[1.134,1.473,1.273,0.772,0.695,0.488,0.268,0.032,1.125,-0.785,-0.458,-0.203],[null,null,null,null,null,null,null,null,null,null,null,null],[-0.99,-0.906,-0.945,1.003,0.8179999999999999,0.824,-0.961,-1.024,-1.117,1.226,1.181,1.173],[-0.789,-0.885,-0.983,0.92,1.038,0.996,-0.875,-0.638,-0.758,0.412,0.405,0.149],[-0.518,-0.393,-0.145,1.16,1.131,0.952,-0.382,-0.608,-0.342,-0.679,-0.652,-0.739],[0.8070000000000001,0.832,0.915,0.95,0.959,0.98,0.977,0.9419999999999999,1.057,-0.9350000000000001,-1.114,-1.104],[1.025,1.013,0.972,0.862,0.836,0.85,0.929,0.77,1.046,0.292,0.207,0.192],[-1.076,-1.139,-1.05,0.769,0.838,0.79,-0.842,-1.048,-1.132,-0.017,-0.005,0.127],[-0.9389999999999999,-0.928,-0.951,0.992,0.823,1.007,-0.862,-0.978,-0.858,-0.989,-0.976,-1.114],[0.273,0.474,0.525,1.36,1.071,1.339,0.637,0.606,0.349,-0.634,-0.786,-0.647],[0.907,0.861,1.063,0.974,1.168,0.944,1.035,1.102,0.859,-0.175,-0.142,-0.186],[0.953,0.836,0.878,0.952,1.052,1.031,0.899,1.043,0.9429999999999999,1.167,1.076,1.029]],"isCatBarplot":false,"showSC":true,"nVars":12,"labelColor":["#440154","#482173","#433E85","#38598C","#2D708E","#25858E","#1E9B8A","#2BB07F","#51C56A","#85D54A","#C2DF23","#FDE725"]},"evals":[],"jsHooks":[]}</script>
```

:::

::: {.cell-output-display}

```{=html}
<h4 id="aweSOMwidget-26c0b86ced9b810441f5-info"></h4>
<h4 id="aweSOMwidget-26c0b86ced9b810441f5-message"></h4>
<div id="aweSOMwidget-26c0b86ced9b810441f5" class="aweSOMwidget html-widget" style="width:400px;height:auto;display:block; margin:auto; margin-top:5px; margin-bottom:5px;"></div>
<p id="aweSOMwidget-26c0b86ced9b810441f5-names"></p>
<svg id="aweSOMwidget-26c0b86ced9b810441f5-placeHolder"></svg>
<script type="application/json" data-for="aweSOMwidget-26c0b86ced9b810441f5">{"x":{"plotType":"Line","sizeInfo":400,"gridInfo":{"nbLines":5,"nbColumns":5,"topology":"hexagonal"},"superclass":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"superclassColor":["#8DD3C7"],"cellNames":["245, 249, 256, 266, 272, 276, 278, 283, 287, 306, 311, 317, 318, 320","241, 246, 248, 251, 253, 255, 257, 258, 263, 273, 277, 279, 282, 285, 289, 290, 294, 300, 301, 305, 308, 309, 310, 312, 315","86","83, 85, 88, 91, 96, 98, 101, 103, 107, 111, 112, 115, 119, 120, 126, 139, 142, 143, 144, 148, 156, 158, 159","87, 90, 93, 94, 102, 104, 105, 108, 110, 114, 116, 121, 124, 125, 127, 135, 137, 145, 146, 151, 160","242, 244, 250, 254, 260, 262, 265, 268, 271, 280, 281, 286, 288, 292, 295, 303, 313, 316, 319","247, 252, 259, 261, 264, 267, 269, 270, 284, 293, 298, 302, 307, 314","243, 274, 275, 291, 296, 297, 299, 304","81, 89, 97, 99, 106, 118, 128, 131, 132, 133, 134, 136, 138, 140, 141, 149, 150, 153, 154, 157","82, 84, 92, 95, 100, 109, 113, 117, 122, 123, 129, 130, 147, 152, 155","","","","55","","164, 170, 172, 177, 179, 180, 189, 192, 195, 199, 207, 208, 209, 212, 213, 214, 222, 226, 227, 237, 238","165, 171, 174, 181, 187, 193, 202, 204, 215, 217, 220, 230, 234, 236","168, 176, 188","10, 11, 15, 17, 18, 20, 22, 23, 24, 25, 31, 35, 36, 38, 42, 46, 48, 49, 54, 56, 59, 60, 63, 66, 69, 75, 77, 79","2, 5, 7, 12, 14, 27, 29, 30, 34, 37, 40, 43, 64, 65, 67, 71, 74","162, 163, 166, 169, 184, 186, 194, 196, 198, 200, 203, 221, 223, 224, 225, 240","161, 167, 173, 175, 178, 182, 183, 185, 190, 191, 197, 201, 205, 206, 210, 211, 216, 218, 219, 228, 229, 231, 232, 233, 235, 239","51, 57","3, 16, 28, 44, 45, 47, 50, 52, 53, 70, 72, 76","1, 4, 6, 8, 9, 13, 19, 21, 26, 32, 33, 39, 41, 58, 61, 62, 68, 73, 78, 80"],"clustering":[25,20,24,25,20,25,20,25,25,19,19,20,25,20,19,24,19,19,25,19,25,19,19,19,19,25,20,24,20,20,19,25,25,20,19,19,20,19,25,20,25,19,20,24,24,19,24,19,19,24,23,24,24,19,14,19,23,25,19,19,25,25,19,20,20,19,20,25,19,24,20,24,25,20,19,24,19,25,19,25,9,10,4,10,4,3,5,4,9,5,4,10,5,5,10,4,9,4,9,10,4,5,4,5,5,9,4,5,10,5,4,4,10,5,4,5,10,9,4,4,5,10,10,5,5,4,5,9,10,10,9,9,9,9,5,9,5,9,4,9,9,4,4,4,5,5,10,4,9,9,5,10,9,9,10,4,9,4,4,5,22,21,21,16,17,21,22,18,21,16,17,16,22,17,22,18,16,22,16,16,17,22,22,21,22,21,17,18,16,22,22,16,17,21,16,21,22,21,16,21,22,17,21,17,22,22,16,16,16,22,22,16,16,16,17,22,17,22,22,17,21,16,21,21,21,16,16,22,22,17,22,22,22,17,22,17,16,16,22,21,2,6,8,6,1,2,7,2,1,6,2,7,2,6,2,1,2,2,7,6,7,6,2,7,6,1,7,6,7,7,6,1,2,8,8,1,2,1,2,6,6,2,1,7,2,6,1,6,2,2,8,6,7,2,6,8,8,7,8,2,2,7,6,8,2,1,7,2,2,2,1,2,6,7,2,6,1,1,6,1],"cellPop":[14,25,1,23,21,19,14,8,20,15,0,0,0,1,0,21,14,3,28,17,16,26,2,12,20],"showAxes":true,"transparency":true,"showNames":true,"legendPos":"beside","legendFontsize":14,"legendReverse":false,"label":["Z_Var1","Z_Var2","Z_Var3","Z_Var4","Z_Var5","Z_Var6","Z_Var7","Z_Var8","Z_Var9","Z_Var10","Z_Var11","Z_Var12"],"normalizedValues":[[0.05,0.09845135479254073,0.06312558905577184,0.06804589402451992,0.05,0.09562431447933992,0.9500000000000001,0.8719523815407212,0.9500000000000001,0.6125837584075988,0.6574812705985936,0.6189469366040105],[0.08857759008163235,0.105907117310156,0.05707779430786286,0.1164670975411831,0.1365144712177756,0.158128292275026,0.8693835944789301,0.9167990424720932,0.8142688420065105,0.4165458030197262,0.4067459306473301,0.4540418841331614],[0.5127482286954579,0.6944997552058161,0.7116845081424661,0.1776252714038998,0.2271494354616294,0.05,0.2904353973865388,0.09658671886773591,0.4386978566656336,0.1996928641168658,0.258927868858518,0.2240552363030068],[0.8463174876607997,0.7913454155652869,0.8137645060098385,0.08947558292678637,0.06312166436245348,0.06754147867233545,0.1785922754420079,0.05,0.08698476912207531,0.5049428448424013,0.5143411172019526,0.4741359320068767],[0.899666242582813,0.7494631208952575,0.800118658133368,0.09561683296522269,0.08329863411359048,0.117481214764859,0.1484538116315015,0.08173531005936152,0.1458747283361189,0.7422282108380349,0.7324091182742033,0.7795045842836833],[0.1723642461609303,0.15042448426648,0.09931246026264401,0.08693921281492301,0.1543322456186746,0.1148290262470307,0.8478759077654557,0.9129889488103128,0.8222246510744543,0.9500000000000001,0.9500000000000001,0.9500000000000001],[0.2200473668688029,0.1394343422715492,0.1984524160897938,0.09445511244572544,0.1557004848104969,0.1556964000138447,0.7455441604210447,0.8183432579098927,0.7815195911018031,0.6399143069541775,0.5995028462650621,0.572473244497671],[0.1566290131357191,0.1512547317439397,0.1083330581265388,0.1241016704448015,0.1504020253645396,0.146011049938514,0.7458783734097891,0.8242889528052476,0.8574662665711995,0.07062022119489191,0.1712498728175972,0.05032324957055508],[0.8928488048577458,0.7808596638964286,0.7670626007377043,0.1705459069395758,0.07163674373181438,0.09877259363757007,0.1209272091010103,0.1160784369405797,0.1344285104547597,0.05,0.05,0.05],[0.9390495977925585,0.814456653110988,0.930119813883826,0.05,0.06879258897071303,0.1142208005674119,0.05,0.1546527658267469,0.1637454554858418,0.3495911243759209,0.3530346031869097,0.4379888502891326],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[0.9500000000000001,0.9500000000000001,0.9500000000000001,0.7334793096612111,0.7590438378361971,0.6382314534017896,0.5989851653526258,0.502513047210761,0.9326748438819138,0.2773860674543602,0.3676331934633182,0.4261300789288175],[null,null,null,null,null,null,null,null,null,null,null,null],[0.107076470223123,0.1301933997994004,0.09069392929784147,0.8184476193155803,0.8085573508692021,0.7610409720035139,0.1529619772292712,0.06129777498029563,0.0556583104542358,0.8577754136127497,0.8206888711124233,0.8288170311681037],[0.1868254014170879,0.1374715803409088,0.07585036277110642,0.788165110629112,0.8976602468275806,0.8242063136666126,0.1840499175612459,0.2222890566769031,0.1960544180205192,0.6229131171744096,0.6061827507297842,0.5290955601486533],[0.2942152035592111,0.3071597250249885,0.4006346796841981,0.8764287861711313,0.9353679015032749,0.8082692254396113,0.3628598220001634,0.2351249015335853,0.3589330547006762,0.3080830218281665,0.3139747191932172,0.2693725486356411],[0.8202323809911798,0.729317630020988,0.811246449889903,0.7989183476648734,0.865801378745682,0.8185196446774852,0.8564631857274251,0.8828949747052077,0.9063217062659349,0.2341240322158078,0.1860005441047004,0.1623954148160907],[0.9068247409097622,0.7914930450640806,0.8332766236584539,0.7664472706110095,0.8161305648171556,0.7707450512339071,0.8388396191683064,0.8109263905019205,0.9016572166628178,0.5881482175986895,0.5513059513580776,0.5416143170181472],[0.07296095597137137,0.05,0.05,0.7324965704546614,0.8166365242714549,0.7485998917095468,0.1958947583535889,0.05105847956324232,0.05,0.4990938049126612,0.4927145127169753,0.5227316327582177],[0.1274149126551875,0.1227519273249815,0.08821069320029765,0.8146231109471999,0.8108189499797752,0.8281801856145649,0.1889102534326543,0.08020487648742039,0.1570414576840522,0.2183949517643701,0.2243559715027022,0.1594628718621929],[0.6085180307050093,0.6058139507770435,0.660233094764862,0.9500000000000001,0.9109032877171012,0.9500000000000001,0.7329504375955722,0.7423822176330352,0.6291617524753448,0.3209524533113836,0.2768158175039977,0.2963113639306945],[0.8599176731012064,0.739331789850999,0.8686203252012453,0.807711451749885,0.9500000000000001,0.8049947441415318,0.8776792008696508,0.9500000000000001,0.8285603010944435,0.4534310214800714,0.4547727977579033,0.4310136681163551],[0.8782461014080221,0.7305742406758964,0.7969968769385186,0.7997929674667038,0.9033880892655833,0.837145979950884,0.8279572375090215,0.9252489469209939,0.861458630972154,0.8407296949080898,0.7918050008929475,0.7867269957170437]],"realValues":[[-1.134,-0.999,-1.016,-1.034,-1.059,-0.993,1.235,0.916,1.169,0.376,0.591,0.456],[-1.036,-0.977,-1.031,-0.903,-0.845,-0.822,1.013,1.023,0.822,-0.303,-0.316,-0.108],[0.032,0.731,0.658,-0.737,-0.621,-1.118,-0.582,-0.9389999999999999,-0.138,-1.054,-0.851,-0.894],[0.873,1.012,0.921,-0.976,-1.027,-1.07,-0.89,-1.051,-1.037,0.003,0.073,-0.039],[1.007,0.891,0.886,-0.96,-0.977,-0.9330000000000001,-0.973,-0.975,-0.887,0.826,0.862,1.005],[-0.825,-0.848,-0.922,-0.983,-0.801,-0.9409999999999999,0.953,1.014,0.842,1.546,1.649,1.587],[-0.705,-0.88,-0.666,-0.963,-0.798,-0.829,0.672,0.787,0.738,0.471,0.381,0.297],[-0.865,-0.845,-0.899,-0.882,-0.8110000000000001,-0.855,0.672,0.802,0.9330000000000001,-1.501,-1.168,-1.487],[0.99,0.982,0.801,-0.756,-1.006,-0.984,-1.049,-0.893,-0.916,-1.573,-1.606,-1.488],[1.106,1.079,1.222,-1.083,-1.013,-0.9419999999999999,-1.244,-0.8,-0.841,-0.535,-0.51,-0.162],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null],[1.134,1.473,1.273,0.772,0.695,0.488,0.268,0.032,1.125,-0.785,-0.458,-0.203],[null,null,null,null,null,null,null,null,null,null,null,null],[-0.99,-0.906,-0.945,1.003,0.8179999999999999,0.824,-0.961,-1.024,-1.117,1.226,1.181,1.173],[-0.789,-0.885,-0.983,0.92,1.038,0.996,-0.875,-0.638,-0.758,0.412,0.405,0.149],[-0.518,-0.393,-0.145,1.16,1.131,0.952,-0.382,-0.608,-0.342,-0.679,-0.652,-0.739],[0.8070000000000001,0.832,0.915,0.95,0.959,0.98,0.977,0.9419999999999999,1.057,-0.9350000000000001,-1.114,-1.104],[1.025,1.013,0.972,0.862,0.836,0.85,0.929,0.77,1.046,0.292,0.207,0.192],[-1.076,-1.139,-1.05,0.769,0.838,0.79,-0.842,-1.048,-1.132,-0.017,-0.005,0.127],[-0.9389999999999999,-0.928,-0.951,0.992,0.823,1.007,-0.862,-0.978,-0.858,-0.989,-0.976,-1.114],[0.273,0.474,0.525,1.36,1.071,1.339,0.637,0.606,0.349,-0.634,-0.786,-0.647],[0.907,0.861,1.063,0.974,1.168,0.944,1.035,1.102,0.859,-0.175,-0.142,-0.186],[0.953,0.836,0.878,0.952,1.052,1.031,0.899,1.043,0.9429999999999999,1.167,1.076,1.029]],"isCatBarplot":false,"showSC":true,"nVars":12,"labelColor":["#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080","#808080"]},"evals":[],"jsHooks":[]}</script>
```

:::

::: {.cell-output-display}

```{=html}
<h4 id="aweSOMwidget-5226b60e9222b46835e1-info"></h4>
<h4 id="aweSOMwidget-5226b60e9222b46835e1-message"></h4>
<div id="aweSOMwidget-5226b60e9222b46835e1" class="aweSOMwidget html-widget" style="width:400px;height:auto;display:block; margin:auto; margin-top:5px; margin-bottom:5px;"></div>
<p id="aweSOMwidget-5226b60e9222b46835e1-names"></p>
<svg id="aweSOMwidget-5226b60e9222b46835e1-placeHolder"></svg>
<script type="application/json" data-for="aweSOMwidget-5226b60e9222b46835e1">{"x":{"plotType":"Cloud","sizeInfo":400,"gridInfo":{"nbLines":5,"nbColumns":5,"topology":"hexagonal"},"superclass":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"superclassColor":["#8DD3C7"],"cellNames":["245, 249, 256, 266, 272, 276, 278, 283, 287, 306, 311, 317, 318, 320","241, 246, 248, 251, 253, 255, 257, 258, 263, 273, 277, 279, 282, 285, 289, 290, 294, 300, 301, 305, 308, 309, 310, 312, 315","86","83, 85, 88, 91, 96, 98, 101, 103, 107, 111, 112, 115, 119, 120, 126, 139, 142, 143, 144, 148, 156, 158, 159","87, 90, 93, 94, 102, 104, 105, 108, 110, 114, 116, 121, 124, 125, 127, 135, 137, 145, 146, 151, 160","242, 244, 250, 254, 260, 262, 265, 268, 271, 280, 281, 286, 288, 292, 295, 303, 313, 316, 319","247, 252, 259, 261, 264, 267, 269, 270, 284, 293, 298, 302, 307, 314","243, 274, 275, 291, 296, 297, 299, 304","81, 89, 97, 99, 106, 118, 128, 131, 132, 133, 134, 136, 138, 140, 141, 149, 150, 153, 154, 157","82, 84, 92, 95, 100, 109, 113, 117, 122, 123, 129, 130, 147, 152, 155","","","","55","","164, 170, 172, 177, 179, 180, 189, 192, 195, 199, 207, 208, 209, 212, 213, 214, 222, 226, 227, 237, 238","165, 171, 174, 181, 187, 193, 202, 204, 215, 217, 220, 230, 234, 236","168, 176, 188","10, 11, 15, 17, 18, 20, 22, 23, 24, 25, 31, 35, 36, 38, 42, 46, 48, 49, 54, 56, 59, 60, 63, 66, 69, 75, 77, 79","2, 5, 7, 12, 14, 27, 29, 30, 34, 37, 40, 43, 64, 65, 67, 71, 74","162, 163, 166, 169, 184, 186, 194, 196, 198, 200, 203, 221, 223, 224, 225, 240","161, 167, 173, 175, 178, 182, 183, 185, 190, 191, 197, 201, 205, 206, 210, 211, 216, 218, 219, 228, 229, 231, 232, 233, 235, 239","51, 57","3, 16, 28, 44, 45, 47, 50, 52, 53, 70, 72, 76","1, 4, 6, 8, 9, 13, 19, 21, 26, 32, 33, 39, 41, 58, 61, 62, 68, 73, 78, 80"],"clustering":[25,20,24,25,20,25,20,25,25,19,19,20,25,20,19,24,19,19,25,19,25,19,19,19,19,25,20,24,20,20,19,25,25,20,19,19,20,19,25,20,25,19,20,24,24,19,24,19,19,24,23,24,24,19,14,19,23,25,19,19,25,25,19,20,20,19,20,25,19,24,20,24,25,20,19,24,19,25,19,25,9,10,4,10,4,3,5,4,9,5,4,10,5,5,10,4,9,4,9,10,4,5,4,5,5,9,4,5,10,5,4,4,10,5,4,5,10,9,4,4,5,10,10,5,5,4,5,9,10,10,9,9,9,9,5,9,5,9,4,9,9,4,4,4,5,5,10,4,9,9,5,10,9,9,10,4,9,4,4,5,22,21,21,16,17,21,22,18,21,16,17,16,22,17,22,18,16,22,16,16,17,22,22,21,22,21,17,18,16,22,22,16,17,21,16,21,22,21,16,21,22,17,21,17,22,22,16,16,16,22,22,16,16,16,17,22,17,22,22,17,21,16,21,21,21,16,16,22,22,17,22,22,22,17,22,17,16,16,22,21,2,6,8,6,1,2,7,2,1,6,2,7,2,6,2,1,2,2,7,6,7,6,2,7,6,1,7,6,7,7,6,1,2,8,8,1,2,1,2,6,6,2,1,7,2,6,1,6,2,2,8,6,7,2,6,8,8,7,8,2,2,7,6,8,2,1,7,2,2,2,1,2,6,7,2,6,1,1,6,1],"cellPop":[14,25,1,23,21,19,14,8,20,15,0,0,0,1,0,21,14,3,28,17,16,26,2,12,20],"showAxes":true,"transparency":true,"showNames":true,"legendPos":"none","legendFontsize":14,"legendReverse":false,"normalizedValues":[[-0.08938189578986433,0.02225745718722899],[0.1779259962946015,0.392361114958159],[-0.3328025249433563,0.2721464677434906],[-0.2628817845356601,-0.0662219839256217],[-0.2170081724350225,-0.1482747456051056],[0.439521748901184,0.04905978886800465],[0.155888056164681,-0.0933329321298446],[-0.001881963554883455,-0.1915537777429545],[-0.1649392389431272,-0.1301998720168935],[-0.5,0.1840207908833286],[0.2656708239009459,-0.1206521040068907],[0.08843117495568029,-0.3944094774097216],[0.4635831705535347,-0.04973406947672107],[-0.2687369018712231,-0.1470689637027187],[0.2827125625805907,-0.2827929409803042],[0.1686033021972875,-0.08913111684871106],[-0.07801435819361993,0.08706668369849223],[0.01108971585778436,0.1117158263926534],[-0.2338658563360615,-0.2575598368204493],[0.3108843705574283,-0.176839617488303],[0.08747689296451898,0.1301418166630946],[-0.1847782886840978,0.2103742467115529],[0.3675698430103161,0.007765746294110243],[-0.02718571138637405,-0.3225787337307827],[0.3235660361937899,0.360353924198164],[-0.2373727345527932,-0.0698705729960712],[-0.02041268262869129,0.01815658894458303],[-0.2988342458375819,0.1323485518120052],[0.1699219856154359,0.1357631423711284],[-0.03057605846650279,-0.1851114915042828],[0.07199110054001252,-0.2499045994154777],[-0.1286903145947414,0.1351566193192346],[-0.04914326921767456,0.2987054952066347],[-0.5,0.2981657265457095],[-0.03662831011431286,-0.2705075055136492],[-0.4611163778164124,0.1032643195865759],[0.28697293643007,0.045072752582855],[0.3310559737510965,0.2077975132010809],[-0.1341235077556343,0.1952661905406892],[0.1702695568181226,0.03600629267439973],[-0.2618206829416927,0.192254658517879],[0.2474381821485624,0.1331034061534553],[0.1082175286727302,-0.08211901567867166],[0.5,0.3199271385366139],[-0.09750085854445112,0.2375169162992073],[0.2264632273187764,-0.05742654872040292],[-0.1629970797478023,-0.003327094051776843],[0.07800561296251875,-0.2927746520589867],[-0.06539763413159524,0.4299604252299882],[-0.08123756824378706,-0.1828304144434907],[-0.5,-0.07715383899221386],[0.486490323235776,-0.1381606505913332],[-0.3475357683688711,-0.09900032565908647],[-0.1634795981554602,-0.02359323361384202],[0,0],[0.03051309240137814,-0.03576190278799275],[0.3587414258300325,-0.2274917788057088],[0.1074203100451177,-0.2473062175996077],[0.3440139632814646,0.2223861461184218],[0.08187180145075637,0.01017132473052257],[-0.1069523223183983,0.09956210733902074],[0.1664099569285262,-0.1224408853061939],[-0.3141290632573182,0.1328535653352398],[-0.1918245217343945,-0.01953978634439462],[-0.1083385113732159,-0.02521497224163119],[-0.3632118039040403,-0.1002984289281991],[0.2174330249101365,0.05119795518043583],[0.5,0.09655169580173795],[-0.440707337042587,-0.1226522044404035],[0.2155987735884045,-0.1565871879802502],[0.1542298706097905,0.06814182131159466],[0.1150271098021611,0.04017551187787936],[0.2643405689980403,-0.08625925618473229],[-0.1923932819621976,0.05020599004750494],[0.04999692866307373,0.1208902998952655],[-0.1648114631377791,-0.3330777966945483],[0.08717745739652927,-0.1212172978656451],[-0.3231930378021977,-0.0582773252648426],[-0.4753722093292083,-0.144724448877974],[-0.03450604004819302,0.06046796789056315],[-0.040108446004873,-0.07793831907427437],[-0.06064438837344607,-0.02302380805553271],[0.001238202355336529,-0.2435268095180662],[-0.01363749874624288,0.4716442117665638],[-0.06497805120250506,-0.1147459853150022],[0,0],[0.02586617430301886,-0.1021587100077977],[-0.3031570167062717,0.2256900945029411],[0.07571553370117405,-0.01899964852652224],[0.03635913458171355,0.03477816357777783],[0.149001389823149,-0.3204716883650386],[0.4344847679428247,0.1686676356528224],[-0.04708773282006083,0.002026837752921653],[-0.02056167593431095,0.06988763179323736],[0.07020064300802681,-0.09410636542520201],[-0.1522312956742732,0.08915082955194387],[0.1572319553471493,0.1804398075457374],[0.2852340685050946,0.06284252586788032],[0.1874446837409529,0.008485051384593174],[0.1787671715171925,0.1557926198631452],[-0.2224877947008295,-0.1323279523587233],[0.1148674793322899,-0.05269432518477182],[-0.05088429790873922,0.0113270348340577],[0.1287468738956757,0.005479336318578643],[0.1781247303307083,-0.06065374768749909],[-0.1328097966713778,-0.1441893954654176],[-0.1939472998964609,0.007474713139084296],[-0.04117461550017203,-0.1488259007642918],[-0.05650484523243408,0.1603197073102481],[0.1406472735158469,0.0477190977532241],[-0.1427636555645252,-0.1564140483868865],[-0.08151279517263296,0.01857033369662298],[0.1590175233994827,-0.1449949327036845],[0.07480767636076563,-0.06274747241625828],[-0.1488241410030403,0.0959728939462131],[0.04210345962529045,-0.004058571915385089],[-0.1493441313212801,-0.09651714364050294],[-0.5,0.1102904413411138],[0.3025039194305784,0.02907578856842383],[0.4344878713214298,-0.03108202919876034],[0.1572252905227816,0.010970142635812],[0.2925535559826373,-0.1463814898224874],[-0.3809599552898058,-0.2221546377615003],[-0.1980611359041776,0.2003540416887248],[-5.965745535401835e-05,0.1484687878009851],[0.01755457293947833,0.207912174829761],[-0.01675151713284034,0.03329682033073033],[0.1113798034032872,0.2280425116658407],[-0.5,0.2723856115860148],[0.2429621098939951,-0.04467970588487583],[0.1149548371084284,-0.2002122945348831],[-0.02950039822559994,-0.003096369812564357],[0.0697532839085543,0.05972738362577612],[0.2121096468326938,0.004748637506602593],[-0.1927964150687709,-0.03278259390700223],[-0.2941928794186796,0.0703994694180503],[-0.5,-0.05600182626894866],[-0.01678709162833203,0.01084030429124742],[-0.0938844334048282,0.05095462022531713],[-0.08952429252639447,-0.1576840612379901],[0.02653736044693853,-0.066497713156721],[0.3686361172045781,-0.0404794741635425],[-0.1583468928439754,-0.287122681956853],[0.0512453701215966,-0.1522124667027014],[-0.01273646119680446,-0.1646667744061346],[-0.02025095573517651,-0.02819334442117554],[-0.03660573804304038,-0.2293385866708015],[-0.01738989752479119,0.01585599699690747],[0.05415439153445571,0.2175662140388391],[0.06616662013932556,-0.03093891501970307],[0.03392651110643087,0.07292434490092223],[-0.1230376786434609,-0.1827690221843024],[-0.06512105112797337,-0.03915275022555137],[0.05734024729523971,-0.07664956389702658],[-0.05725153609444823,-0.04484409402990575],[-0.2666479201638837,0.04167973248980457],[0.03525559214503109,-0.0751807898671456],[0.07988587082423619,0.5],[0.2072681092412778,0.121876397316615],[0.1168055631731453,0.0868780624263505],[-0.1532359338184639,-0.009915666181058932],[-0.04963904742949746,0.1198707141287242],[-0.3070676114026669,0.2383080833524118],[0.1591908973647137,-0.02270025950669045],[-0.2562541813839005,0.3149479678474931],[0.007126175767629328,0.2120632233319239],[-0.1552777011368947,0.2378566024362341],[0.1304368167857688,0.3030713069233165],[0.1089780067627075,-0.2528532101009047],[0.1184935695911645,0.08148040266518401],[0.1710178730662782,-0.1779535673657376],[-0.01101319913836078,-0.01369157024237466],[-0.5,-0.2013811179449495],[-0.3134045181973488,-0.2281661837893616],[0.3647001403380249,0.0737960094780907],[-0.5,-0.08334338641777943],[-0.4816395688008772,-0.01888169644004377],[-0.001621033216052492,0.2067784271685295],[-0.08371329805559284,0.04303583518033199],[-0.1072536325012596,-0.04478099913611722],[0.3974136612747673,-0.01946131074760801],[-0.2687953388096785,0.1160604814108759],[0.2193628982533596,-0.3584667094545663],[-0.09787076656260575,-0.1597593370989621],[-0.02106149651310212,0.1525974595746983],[-0.09905969293270256,0.2039136257311889],[-0.1989266071267451,0.2495472468201839],[0.3695631832142311,-0.219727920505537],[-0.09009663543063608,0.01098956139244789],[0.2518361221661125,-0.06183490734393567],[-0.01493184251308764,0.2504926766568088],[0.07493626191900307,-0.209063418923689],[0.1423573533834718,0.1015104848720622],[-0.008822786083754626,-0.1282953178079248],[0.1824067164079552,-0.02178469886095115],[0.5,0.00829498286844528],[0.3162282613684992,-0.004420672758541508],[0.1048766132799501,0.1302052558615436],[0.1857682492451339,0.1132328352662291],[-0.09428541293183292,-0.1682239454476465],[-0.1874786498360104,0.01366543040262586],[0.3139490171221729,0.1961603192953724],[-0.1060931052407268,-0.2249211582190416],[-0.2484989243054426,-0.003170003251132018],[-0.09847226976650789,-0.2221556356659909],[-0.2973366941829101,0.1638965096600576],[-0.5,0.1284995863585953],[0.08556168091898263,-0.08757220630631277],[0.09228890325643295,-0.03432695326091327],[0.1747765710951646,0.07786932154619018],[-0.07904260646340965,-0.1265193529981467],[-0.1686835065607275,-0.1250002328943316],[0.1010611296348635,0.08151793631347948],[0.1213699663466653,0.09464554114951911],[-0.2923758779439,0.1386659306775103],[0.1594257337228851,-0.07487652700229813],[0.4373123473906121,-0.5],[-0.11597243997137,-0.08211509010844031],[-0.1841364770153937,-0.2192996629073541],[0.4015218338660848,0.3838825048637826],[0.07613642359408813,0.0130878279693154],[0.0101187634097774,0.0320169973720556],[-0.1054044541103744,-0.06999006284629199],[0.1806942489409952,0.1150068706228903],[-0.1114009357376747,-0.1941916836544691],[0.009707251355108144,0.04606731073862554],[-0.02020428219368069,-0.1365266037817038],[0.07331899549297026,-0.0001007097660163893],[0.04048898669569853,-0.008606512607263318],[-0.1239140990246843,-0.02252679293007432],[0.1604172131461947,0.2396145186570048],[0.2585239092201858,0.09771392456207324],[0.1474579440305236,0.09415026767005255],[-0.1279823926038131,-0.04899474299833684],[0.2245821335672241,-0.3321431217360905],[-0.3022154855175522,-0.3844418532941534],[0.1719585498676824,0.05733699587736274],[0.1497421833636533,0.02550563703929727],[-0.3137564258539622,-0.02265594274859144],[0.001832344086465382,0.1574841313087975],[-0.2168683972345294,-0.01014773809882736],[-0.08300337911262157,-0.2489613918350222],[0.3838383299187048,-0.02712470402097772],[0.1171142333146382,0.06639926304871456],[-0.04306086798379497,0.005582677856877118],[0.02181036565126977,0.00247944730905728],[-0.1194732665664654,-0.04742101184646859],[-0.03923024066217956,0.07776833907657861],[-0.078112635372336,0.007869801996848105],[-0.05147877734233067,0.04290249670707622],[-0.2851459344804723,0.1235238756880858],[-0.180854347921886,0.09125397588963592],[0.3149959099415864,-0.001797343441180604],[-0.2660781733537448,-0.01944784680003437],[-0.05466227187829559,0.1053792747500747],[0.008178388109746944,0.06241929280146331],[0.1348203390364852,-0.03181402380840659],[-0.08715338472253291,-0.01019684120957954],[0.23817321237244,0.1318600734383538],[0.1728929634221039,-0.05180991535015687],[-0.231335842633987,-0.02917939186125514],[-0.2303918821007541,0.1979476480114949],[-0.2811936643975716,-0.05522778751536555],[-0.1205433417703569,0.08142531418092448],[0.2580914852895205,0.1041886303713239],[-0.06618188107087529,0.1758457194209726],[0.05489957277615262,0.08102720467722865],[0.07987223837205767,-0.04628361234029613],[-0.3302223542416196,0.3944648167718734],[0.06350440413195849,0.06788303324754791],[-0.01161411905136438,-0.03183489834095397],[0.1896133555695576,0.1065060398856759],[-0.1235730661170932,0.07106155719887684],[-0.08550433882229859,-0.1293811895753408],[-0.174612666915525,-0.1478470573409218],[0.2548800650737459,-0.2695544300700939],[0.2398667648080732,-0.2045332628828764],[0.1466993262037284,0.10182949496232],[0.2852130434497482,0.1655938797350217],[0.02352588134412396,0.1607713723228746],[0.5,0.02990915989482134],[-0.1491332165394827,0.002812022747462035],[0.3330122003218917,0.03143447432117661],[-0.1982272645901017,-0.1092204488091731],[-0.04720072267433473,-0.2987960271532223],[0.05503339936251514,-0.104310189326056],[-0.5,0.04100685382720781],[-0.3076634513466345,-0.1445428780661447],[0.02679257663532474,0.4026406448480872],[-0.5,0.1719266442138961],[-0.1452084490438434,-0.2944677606578869],[-0.07594877623652956,0.1197285458501745],[0.1148876453170544,-0.1224022368923045],[0.1350867770660264,-0.06783310265874243],[-0.1270159907432891,0.05424301633000384],[0.5,0.009860784730502656],[-0.2474504823079845,0.2146869279431961],[-0.04293712863314839,-0.3503088732837162],[0.06077095616363746,0.2501189155844323],[0.09901855540148885,0.1758055375764442],[0.03056200899477448,-0.397507538933194],[0.251449988193384,-0.1398749834931487],[-0.2779914391227751,-0.09226863447754934],[-0.2918333489926909,0.1241540833369958],[0.05689576010711991,-0.05631320561847207],[-0.1011618593440059,0.07210938485323252],[0.5,0.2334678196864561],[0.07882612922134409,-0.2799883440086148],[0.142874012440696,0.1985376528896103],[-0.1465401708856597,-0.3424703729444387],[-0.005112361591315916,-0.2707643053580606],[0.05838157694126285,0.1060882742083453],[-0.1907036451727603,0.1458841536434723],[0.0006787235665755946,-0.282975291705954],[0.305557249896951,0.1530084380313785],[0.1814415215093356,-0.05665918565082251],[-0.2487072109964292,-0.1886737471522243],[0.01758435809574096,-0.04124956926351005],[0.2339179109685092,-0.1258557339929205],[0.09236912298434587,0.1656383119181153]],"realValues":["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320"],"label":"","labelColor":"#440154","cloudColor":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"fullData":[]},"evals":[],"jsHooks":[]}</script>
```

:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/SOMFinalize-8.png){width=672}
:::
:::


# Mixed and categorical workflows

## Gower distance followed by PAM

Gower distance combines mixed feature types, while PAM represents clusters by
observed medoids [@gower1971; @kaufman1990]. Projection uses training ranges,
levels, and medoids rather than recalculating them in the new cohort.

::: callout-note
**When to use it:** numeric, binary, ordinal, and nominal measures must be
clustered together. **Inspect before naming:** candidate silhouette, size
balance, medoids, numeric and categorical profiles, distance/margin, and
bootstrap stability.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 8
  Classes Silhouette MeanMedoidDistance MinClusterN SizeBalance
    <int>      <dbl>              <dbl>       <dbl>       <dbl>
1       2      0.260             0.252          119       0.592
2       3      0.485             0.168           84       0.592
3       4      0.752             0.0832          80       1    
4       5      0.617             0.0797          29       0.362
5       6      0.485             0.0772          23       0.288
6       7      0.349             0.0748          23       0.288
7       8      0.217             0.0729          16       0.25 
8       9      0.192             0.0716          16       0.25 
9      10      0.168             0.0703          16       0.25 
# ℹ 3 more variables: ReproducibilityScore <dbl>, ahp_index <dbl>,
#   Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends Gower/PAM k (Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMReview-6.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-8.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/GowerPAMFinalize-9.png){width=672}
:::
:::


## Latent class analysis

LCA models class-specific categorical response probabilities [@linzer2011].

::: callout-note
**When to use it:** symptom, diagnosis, questionnaire, or assay-call items are
the phenotype definition. **Inspect before naming:** BIC, entropy, class size,
response-probability heatmap, posterior confidence, and bootstrap stability.
Low posterior confidence is evidence of weak class separation.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 9
  Classes   BIC   AIC Entropy MinClassN SizeBalance ReproducibilityScore
    <int> <dbl> <dbl>   <dbl>     <int>       <dbl>                <dbl>
1       4 1328. 1181.   0.999        80      1                     0.971
2       5 1385. 1201.   0.949         3      0.0375                0.938
3       6 1443. 1221.   0.931         2      0.025                 0.992
4       7 1501. 1241.   0.925         2      0.025                 0.881
5       8 1558. 1261.   0.824         2      0.025                 0.850
6       9 1616. 1280.   0.764         2      0.0256                0.889
7       3 1643. 1534.   0.999        80      0.5                   0.968
8      10 1674. 1300.   0.878         1      0.0127                0.816
9       2 1983. 1911.   0.999       160      1                     0.986
# ℹ 2 more variables: ahp_index <dbl>, Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends latent class count (Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassReview-6.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/LatentClassFinalize-7.png){width=672}
:::
:::


## MCA followed by Mclust

MCA is a dimension-reduction method for nominal categorical measures
[@greenacre2017]. The bootstrap refits MCA before refitting the mixture.

::: callout-note
**When to use it:** several nominal variables form a useful lower-dimensional
category space. **Inspect before naming:** inertia/scree plot, MCA dimensions,
mixture fit, score-space separation, posterior confidence, and stability.
:::


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 9 × 11
  Model ModelName Classes    BIC    ICL Entropy MinClusterN SizeBalance
  <int> <chr>       <int>  <dbl>  <dbl>   <dbl>       <int>       <dbl>
1     1 EEI            10  1516.  1516.   1.000           2      0.0263
2     1 EEI             9  1351.  1351.   1.000           3      0.0395
3     1 EEI             8  1057.  1057.   1.000           4      0.0526
4     1 EEI             7   861.   861.   1.000           5      0.0641
5     1 EEI             6   583.   583.   1.000           4      0.05  
6     1 EEI             5   400.   400.   1.000           8      0.1   
7     1 EEI             4   173.   173.   1.000          80      1     
8     1 EEI             3 -1113. -1113.   1.000          80      0.5   
9     1 EEI             2 -2486. -2487.   0.991          80      0.333 
# ℹ 3 more variables: ReproducibilityScore <dbl>, ahp_index <dbl>,
#   Recommended <lgl>
```


:::

::: {.cell-output .cell-output-stdout}

```
[1] "AHP-style review recommends MCA plus Gaussian-mixture model (Model = 1, Classes = 4). Review this advisory choice alongside the candidate plots."
```


:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustReview-7.png){width=672}
:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-1.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-2.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-3.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-4.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-5.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-6.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-7.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-8.png){width=672}
:::

::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/MCAMclustFinalize-9.png){width=672}
:::
:::


# Reveal truth for teaching only

The benchmark truth was not used during fitting. ARI is label-invariant, so
cluster numbers do not need to match the truth labels. Strong simulation
recovery confirms that the examples behave as designed; it is not external or
biological validation.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 8 × 3
  Method        TrainingARI ProjectionARI
  <chr>               <dbl>         <dbl>
1 Mclust              1             1    
2 K-means             1             1    
3 PCA + Mclust        1             1    
4 PCA + K-means       1             1    
5 Gower + PAM         1             1    
6 LCA                 1             0.872
7 MCA + Mclust        1             1    
8 SOM + Mclust        0.992         1    
```


:::
:::



::: {.cell}
::: {.cell-output .cell-output-stdout}

```
# A tibble: 2 × 5
  Cohort       ARI NoiseProportion NoiseSensitivity NoiseSpecificity
  <chr>      <dbl>           <dbl>            <dbl>            <dbl>
1 Training   0.841           0.134            0.672                1
2 Projection 0.845           0.125            0.625                1
```


:::
:::


# Compare and interpret profiles

The final label is not the phenotype interpretation. Compare occupancy and
feature distributions across plausible models, then use codebook labels when
describing the profiles.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
  MclustCluster KMeansCluster  n
1             1             1 80
2             2             4 80
3             3             2 80
4             4             3 80
```


:::
:::



::: {.cell}
::: {.cell-output-display}
![](ClusteringPhenotyping_files/figure-html/PlotClusterProfiles-1.png){width=672}
:::
:::


# Reproducibility

# Candidate metrics and SOM + Mclust workflow

All clustering functions have two modes. `method = "exploratory"` compares the
requested 2:10 candidate settings and returns the fit table, candidate audit,
and AHP advisory only. `method = "finalize"` requires the selected setting and
returns assignments, profiles, fit diagnostics, and a reusable projection
model. AHP is advisory; its selected row is retained in the table but is never
highlighted on a figure.

`SizeBalance` is the smallest extracted cluster/class size divided by the
largest (1 is perfectly balanced). `MinClusterN`/`MinClassN` is the smallest
group size. Bootstrap `ARI` measures agreement for the entire partition,
whereas per-cluster Jaccard recovery measures whether each reference phenotype
reappears after refitting. `ReproducibilityScore` is the mean of bootstrap mean
ARI and mean Jaccard recovery; bootstrap success rate is reported separately.

### Complementary stability diagnostics

Every successful full-pipeline bootstrap refit also records variation of
information (`VI`; lower is better), normalized mutual information (`NMI`; higher
is better), and Fowlkes--Mallows agreement (higher is better). These are
secondary partition diagnostics: they do not enter `ReproducibilityScore` or
the AHP advisory index. `participant_inclusion` reports the probability that a
participant returns to its reference cluster's label-matched bootstrap cluster;
`cluster_inclusion` summarizes that probability by reference cluster.

For complete-case cohorts of 2,000 or fewer participants, `coassignment`
contains the bootstrap probability that each pair is assigned together and
`Stability$plots$coassignment` displays it. Above that size, the matrix is
skipped with a recorded reason to avoid unbounded memory use. For density-based
methods, noise remains in global partition metrics but is excluded from
inclusion and co-assignment summaries. These diagnostics explain instability;
they cannot make a weak phenotype stable.

## SOM + Mclust candidate-table fields

`CreateClusterModel_SOM_MClust()` uses tidyLPA with its mclust backend to fit
latent profiles to frozen SOM node codebooks. `ShowFitTable()` in its help
example is only a local display helper: it rounds a table and displays it with
a frozen header; it is not a SciDataReportR function.

| Field | Definition and interpretation |
|---|---|
| `AIC`, `AWE`, `BIC`, `CAIC`, `CLC`, `KIC`, `SABIC`, `ICL` | Likelihood/information criteria from the LPA fit. Lower values are preferred only when comparing candidates fitted to the same SOM nodes and data. AWE, CLC, and ICL incorporate classification uncertainty [@celeux1996]; KIC is based on symmetric Kullback divergence [@cavanaugh1999], and SABIC is sample-size adjusted [@sclove1987]. `CAIC` and `CLC` are separate criteria, not one combined `CAICCLC` field. |
| `Entropy` | Classification separation summary; higher values indicate cleaner profile separation. |
| `MinProfileNodeN`, `MaxProfileNodeN` | Integer SOM-node counts in the smallest and largest profile. |
| `MinProfileNodeProportion`, `MaxProfileNodeProportion` | The corresponding node counts divided by the total SOM-node count. These replace tidyLPA's ambiguous decimal `n_min` and `n_max` fields. |
| `BLRTStatistic`, `BLRTPValue` | Bootstrap likelihood-ratio statistic and p-value for comparing the k-profile fit with k - 1 profiles [@nylund2007]. A small p-value supports the additional profile's statistical fit, but does not by itself establish stability, clinical meaning, or a final choice; it is unavailable for one profile. |
| `StabilityARI_Mean` | Mean adjusted Rand index across successful full-pipeline bootstrap refits. It may be negative: this means agreement is worse than expected by chance, not an error. |
| `StabilityJaccard_Mean` | Mean label-matched, per-profile Jaccard recovery across successful full-pipeline bootstrap refits. |
| `ReproducibilityScore` | Mean of the finite `StabilityARI_Mean` and `StabilityJaccard_Mean`; `StabilitySuccessRate` remains separate. |
| `Reproducibility_scaled` | Min-max rescaling of `ReproducibilityScore` across the current candidate table. It is used only as an AHP input [@akogul2017], not as an independently interpretable reproducibility metric. |

| Method | Candidate metrics |
|---|---|
| Mclust, PCA + Mclust, MCA + Mclust | BIC/ICL/AIC (information criteria), entropy, maximum posterior uncertainty, minimum cluster size, size balance |
| K-means, PCA + K-means | WSS, between-cluster sum of squares, silhouette width, Calinski–Harabasz index, minimum size, size balance |
| Gower + PAM | silhouette width, mean distance to assigned medoid, minimum size, size balance |
| Latent class | BIC, AIC, log likelihood, entropy, minimum class size, size balance |
| HDBSCAN and SOM + HDBSCAN | persistence, noise proportion, membership probability, minimum size, size balance; extracted class count is data-derived |
| SOM + Mclust | LPA AIC/BIC/entropy plus SOM node occupancy and bootstrap stability |

## Candidate-metric definitions for every clustering model

All exploratory `fit_table`s retain the raw metrics produced by their fitted
method. AHP is advisory only: any `*_scaled` field is a within-search
rescaling used to construct its index, not a standalone scientific result.

| Metric | Definition |
|---|---|
| `Classes` | Requested class/cluster count for fixed-count methods. HDBSCAN variants instead report the data-derived number extracted from each density setting. |
| `MinClusterN` / `MinClassN` | Smallest extracted participant cluster/class size. `SizeBalance` is smallest divided by largest extracted participant group size. SOM + Mclust reports its separate integer node-count and node-proportion fields above. |
| `StabilitySuccessRate` | Proportion of requested bootstrap refits that successfully complete. It is reported separately and is never multiplied into reproducibility. |
| `StabilityARI_Mean`, `StabilityARI_P05` | Mean and fifth percentile of the adjusted Rand index across successful full-pipeline bootstrap refits [@hubert1985]. ARI can be negative when agreement is worse than chance. |
| `StabilityJaccard_Mean`, `StabilityJaccard_Min` | Mean and minimum label-matched per-cluster Jaccard recovery across successful bootstrap refits [@jaccard1901]. |
| `ReproducibilityScore` | Mean of finite `StabilityARI_Mean` and `StabilityJaccard_Mean`; it is unavailable when neither component is estimable. |
| `AHPIndex` / `ahp_index` and `*_scaled` | Advisory within-candidate ranking, calculated from the documented method-specific raw metrics. The recommendation does not finalize a model or alter any fit plot. |

| Pipeline | Method-specific candidate metrics |
|---|---|
| Mclust, PCA + Mclust, MCA + Mclust | `BIC`, `ICL`, `AIC`, entropy, and maximum posterior uncertainty. BIC/ICL/AIC compare fitted Gaussian mixtures [@scrucca2016]; lower information criteria are preferred, entropy is higher-is-better, and uncertainty is lower-is-better. |
| K-means, PCA + K-means | Within-cluster sum of squares (`WSS`; lower), between-cluster sum of squares (`BSS`; higher), mean silhouette width (higher), and Calinski--Harabasz index (higher) [@calinski1974; @rousseeuw1987]. |
| Gower + PAM | Mean silhouette width (higher) and mean assigned-medoid Gower distance (lower). The latter is a dissimilarity, not a probability [@gower1971; @kaufman1990]. |
| Latent class analysis | Log likelihood (higher), `AIC`/`BIC` (lower), and entropy (higher) from the fitted item-response model [@linzer2011]. |
| HDBSCAN | `Persistence` (mean cluster persistence; higher), `NoiseProportion` (lower), and the data-derived extracted class count. Membership probability and outlier score are assignment diagnostics, not candidate-selection probabilities [@mcinnes2017]. |
| SOM + Mclust | The LPA criteria, node allocation, BLRT, and bootstrap metrics defined in the dedicated SOM + Mclust table below. SOM distances describe map representation, not phenotype probability [@kohonen1982]. |
| SOM + HDBSCAN | Node-codebook HDBSCAN persistence, noise proportion, extracted node-cluster count, and bootstrap stability. `minPts` and epsilon are requested settings; class count is data-derived. |

`CreateClusterModel_SOM_MClust()` is the explicit SOM + Mclust API.
`CreateSOMClusterModel()` remains available for compatibility but is deprecated.
The SOM first freezes Z-score preprocessing and learns map nodes; Mclust then
clusters those node prototypes. A SOM distance is the distance, in this frozen
scaled feature space, between a participant and their best-matching map node.
Lower values mean the participant is better represented by the map; they do
not mean stronger membership in a phenotype. Interpret it against the frozen
training distribution only: distances above the training p95 flag a case that
is poorly represented by the training phenotype space.

# How projection works

Every fitted clustering object freezes its preprocessing and records it in
`Specification`; every projected object reports participant-level diagnostics
in `ProjectionFit$individual`, cohort summaries in `ProjectionFit$summary`,
cluster summaries in `ProjectionFit$by_cluster`, support violations in
`ProjectionFit$out_of_support`, and the applied thresholds in
`ProjectionFit$policy`.

| Workflow | Frozen projection rule | Assignment evidence |
|---|---|---|
| Mclust | Gaussian-mixture `predict()` after training scaling | Posterior probabilities and assigned-component Mahalanobis distance |
| K-means | Nearest frozen centroid | Euclidean centroid distance and second-centroid margin; not a probability |
| PCA/MCA + clustering | Training scaling and frozen reduction, then the fitted downstream model | The downstream model's metric/probability |
| Gower + PAM | Nearest frozen observed medoid using training ranges and levels | Gower medoid distance and runner-up margin |
| Latent class | Frozen class prevalence and item-response probabilities | Posterior class probabilities |
| HDBSCAN | Conservative nearest in-training-cluster support rule | Nearest-core distance and support flag; **not native HDBSCAN prediction** |
| SOM + Mclust | Frozen scaling, SOM BMU, and node mixture probabilities | BMU distance and node-derived posterior probability |
| SOM + HDBSCAN | Frozen scaling, SOM BMU, and node HDBSCAN label | BMU distance and inherited node cluster/noise label; **not native HDBSCAN prediction** |

`ReproducibilityScore` is the mean of bootstrap mean ARI and mean per-cluster
Jaccard recovery. ARI summarizes agreement for the whole partition, whereas
Jaccard asks whether each individual phenotype was recovered. Bootstrap
success rate is reported separately and is not folded into the score.

<details>
<summary>Session information</summary>


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.5.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_4.0.3          dplyr_1.2.1            SciDataReportR_20.24.0

loaded via a namespace (and not attached):
  [1] ggstatsplot_1.0.0      splines_4.5.3          later_1.4.8           
  [4] tibble_3.3.1           R.oo_1.27.1            datawizard_1.3.1      
  [7] rpart_4.1.27           fastDummies_1.7.6      lifecycle_1.0.5       
 [10] gtsummary_2.5.1        rstatix_1.1.0          globals_0.19.1        
 [13] lattice_0.22-9         MASS_7.3-65            insight_1.5.2         
 [16] flashClust_1.1-4       backports_1.5.1        magrittr_2.0.5        
 [19] Hmisc_5.2-6            rmarkdown_2.31         yaml_2.3.12           
 [22] httpuv_1.6.17          otel_0.2.0             poLCA_1.6.0.2         
 [25] RColorBrewer_1.1-3     multcomp_1.4-30        abind_1.4-8           
 [28] quadprog_1.5-8         purrr_1.2.2            R.utils_2.13.0        
 [31] nnet_7.3-20            TH.data_1.1-5          sandwich_3.1-1        
 [34] labelled_2.16.0        tidyLPA_2.0.2          ggrepel_0.9.8         
 [37] irlba_2.3.7            listenv_1.0.0          kohonen_3.0.13        
 [40] correlation_0.8.8      proto_1.0.0            aweSOM_1.3            
 [43] texreg_1.39.5          parallelly_1.48.0      codetools_0.2-20      
 [46] DT_0.34.0              ggtext_0.1.2           xml2_1.6.0            
 [49] tidyselect_1.2.1       farver_2.1.2           viridis_0.6.5         
 [52] effectsize_1.0.3       stats4_4.5.3           base64enc_0.1-6       
 [55] showtext_0.9-8         jsonlite_2.0.0         progressr_1.0.0       
 [58] Formula_1.2-5          survival_3.8-6         emmeans_2.0.3         
 [61] dbscan_1.2.5           tools_4.5.3            progress_1.2.3        
 [64] Rcpp_1.1.2             glue_1.8.1             mnormt_2.1.2          
 [67] gridExtra_2.3.1        xfun_0.60              MplusAutomation_1.3   
 [70] withr_3.0.3            fastmap_1.2.0          boot_1.3-32           
 [73] digest_0.6.39          R6_2.6.1               mime_0.13             
 [76] estimability_2.0.0     colorspace_2.1-3       dichromat_2.0-0.1     
 [79] R.methodsS3_1.8.2      utf8_1.2.6             tidyr_1.3.2           
 [82] generics_0.1.4         data.table_1.18.4      prettyunits_1.2.0     
 [85] httr_1.4.8             htmlwidgets_1.6.4      parameters_0.29.2     
 [88] scatterplot3d_0.3-45   pkgconfig_2.0.3        gtable_0.3.6          
 [91] statsExpressions_2.0.0 S7_0.2.2               htmltools_0.5.9       
 [94] lavaan_0.6-21          carData_3.0-6          sysfonts_0.8.9        
 [97] multcompView_0.1-12    scales_1.4.0           leaps_3.2             
[100] snakecase_0.11.1       knitr_1.51             rstudioapi_0.19.0     
[103] coda_0.19-4.1          checkmate_2.3.4        nlme_3.1-169          
[106] showtextdb_3.0         zoo_1.9-0              stringr_1.6.0         
[109] sjlabelled_1.2.0       parallel_4.5.3         nonnest2_0.5-9        
[112] foreign_0.8-91         pillar_1.11.1          grid_4.5.3            
[115] vctrs_0.7.3            RANN_2.6.2             promises_1.5.0        
[118] car_3.1-5              xtable_1.8-8           cluster_2.1.8.2       
[121] htmlTable_2.5.0        paletteer_1.7.0        evaluate_1.0.5        
[124] pbivnorm_0.6.0         gsubfn_0.7             mvtnorm_1.4-2         
[127] cli_3.6.6              compiler_4.5.3         rlang_1.3.0           
[130] crayon_1.5.3           rstantools_2.7.0       future.apply_1.20.2   
[133] labeling_0.4.3         mclust_6.1.3           rematch2_2.1.2        
[136] plyr_1.8.9             forcats_1.0.1          stringi_1.8.7         
[139] psych_2.6.5            pander_0.6.6           viridisLite_0.4.3     
[142] CompQuadForm_1.4.4     bayestestR_0.18.1      Matrix_1.7-5          
[145] tidySEM_0.2.10         hms_1.1.4              patchwork_1.3.2       
[148] future_1.70.0          FactoMineR_2.16        shiny_1.14.0          
[151] haven_2.5.5            gridtext_0.1.6         igraph_2.3.3          
[154] broom_1.0.13           RcppParallel_6.2.0    
```


:::
:::


</details>

# References

::: {#refs}
:::

