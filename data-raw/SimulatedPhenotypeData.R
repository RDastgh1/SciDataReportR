# Fixed benchmark data for the clustering and phenotyping vignette.
set.seed(20260804)

rows_TrainingPerCluster <- 80L
rows_ProjectionPerCluster <- 40L
TruthCluster <- factor(rep(paste("Cluster", 1:4), each = rows_TrainingPerCluster + rows_ProjectionPerCluster))
Cohort <- factor(rep(rep(c("Training", "Projection"), c(rows_TrainingPerCluster, rows_ProjectionPerCluster)), 4))
cluster_means <- 2.5 * rbind(c(1, 1, 1), c(1, -1, -1), c(-1, 1, -1), c(-1, -1, 1))
latent_scores <- do.call(rbind, lapply(seq_len(4), function(i) {
  cbind(
    sapply(seq_len(3), function(axis) stats::rnorm(120L, mean = cluster_means[i, axis], sd = 0.5)),
    stats::rnorm(120L, mean = 0, sd = 2.5)
  )
}))
df_Continuous <- as.data.frame(do.call(cbind, lapply(seq_len(4), function(axis) {
  sapply(seq_len(3), function(indicator) latent_scores[, axis] + stats::rnorm(nrow(latent_scores), sd = 0.8))
})))
names(df_Continuous) <- paste0("Var", seq_len(12))
df_Continuous <- df_Continuous[, paste0("Var", seq_len(12)), drop = FALSE]

set.seed(20260805)
df_Categorical <- as.data.frame(lapply(seq_len(3), function(variable_i) {
  response_levels <- paste0("Item ", variable_i, ": ", LETTERS[seq_len(4)])
  factor(
    unlist(lapply(seq_len(4), function(cluster_i) {
      probabilities <- rep(0.01, 4)
      probabilities[[cluster_i]] <- 0.97
      sample(response_levels, 120L, replace = TRUE, prob = probabilities)
    })),
    levels = response_levels)
}))
names(df_Categorical) <- paste0("CatVar", seq_len(3))

set.seed(20260806)
TruthDensityGroup <- factor(c(rep("DenseA", 128), rep("DenseB", 128), rep("Noise", 64), rep("DenseA", 64), rep("DenseB", 64), rep("Noise", 32)))
DensityX <- c(rnorm(128, -2, .25), rnorm(128, 2, .7), runif(64, -5, 5), rnorm(64, -2, .25), rnorm(64, 2, .7), runif(32, -5, 5))
DensityY <- c(rnorm(128, -2, .25), rnorm(128, 2, .7), runif(64, -5, 5), rnorm(64, -2, .25), rnorm(64, 2, .7), runif(32, -5, 5))

SimulatedPhenotypeData <- cbind(data.frame(ParticipantID = sprintf("SIM%03d", seq_len(480)), Cohort, TruthCluster, TruthDensityGroup), df_Continuous, df_Categorical, DensityX, DensityY)
SimulatedPhenotypeData$Var12[SimulatedPhenotypeData$Cohort == "Projection"][1:5] <- NA_real_
for (variable_i in seq_len(12)) {
  attr(SimulatedPhenotypeData[[paste0("Var", variable_i)]], "label") <-
    paste("Var", variable_i)
}
for (variable_i in seq_len(3)) {
  attr(SimulatedPhenotypeData[[paste0("CatVar", variable_i)]], "label") <-
    paste("CatVar", variable_i)
}
attr(SimulatedPhenotypeData$DensityX, "label") <- "Density X"
attr(SimulatedPhenotypeData$DensityY, "label") <- "Density Y"
SimulatedPhenotypeVariableTypes <- data.frame(
  Variable = names(SimulatedPhenotypeData),
  Label = c(
    "Participant identifier", "Cohort", "Simulated cluster (truth only)",
    "Simulated density group (truth only)", paste("Var", 1:12),
    paste("CatVar", 1:3), "Density X", "Density Y"),
  Type = c("Categorical", "Categorical", "Categorical", "Categorical", rep("Double", 12), rep("Categorical", 3), "Double", "Double"),
  Recode = "No",
  Code = NA_character_,
  MissingCode = NA_character_
)
usethis::use_data(SimulatedPhenotypeData, SimulatedPhenotypeVariableTypes, overwrite = TRUE)
