# Shared clustering figures

Figure builders shared by every clustering pipeline in this package.
Each `Create*ClusterModel()` and `Project*Cluster()` function stores its
figures beside the object they describe, following the layout
established by
[`CreateClusterModel_SOM_MClust()`](https://rdastgh1.github.io/SciDataReportR/reference/CreateClusterModel_SOM_MClust.md):

- `fit_plot` reviews candidate solutions across the search grid.

- `ModelInfo$plots` describes the structure of the selected solution.

- `ModelInfo$FitDiagnostics$plots` describes how well individual
  training cases sit inside that structure.

- `ProbFit$plots` describes membership confidence.

- `Stability$plots` describes 80% subsample reproducibility.

- `ProjectionFit$plots` describes how projected cases compare with the
  frozen training reference.

There is deliberately no flat `plots` list: a figure lives next to the
table it summarises so that the object stays self-describing.
