#' Plot 3D PCA Scores
#'
#' This function creates a 3D scatter plot of PCA scores using Plotly. It allows for optional coloring by a variable.
#'
#' @param PCAObj A list containing PCA results, including a `Scores` matrix with PCA scores and a `CombinedData` dataframe with the original data.
#' @param Var An optional string specifying the variable name in `CombinedData` to color the points by. Default is `NULL`.
#' @param t An optional string specifying the type of variable for coloring: `"Factor"` for categorical variables or anything else for continuous variables. Default is `"NULL"`.
#'
#' @return A Plotly object representing the 3D PCA plot.
#'
#' @importFrom plotly plot_ly layout
#' @export
#'


plotPCA <- function(PCAObj, Var = NULL, t = "NULL") {
  if (is.null(Var)) {
    p <- plotly::plot_ly(x = PCAObj$Scores[,1],
                         y = PCAObj$Scores[,2],
                         z = PCAObj$Scores[,3]) %>%
      plotly::layout(
        scene = list(
          xaxis = list(title = colnames(PCAObj$Scores)[1]),
          yaxis = list(title = colnames(PCAObj$Scores)[2]),
          zaxis = list(title = colnames(PCAObj$Scores)[3])
        )
      )
  } else {
    if (t == "Factor") {
      p <- plotly::plot_ly(x = PCAObj$Scores[,1],
                           y = PCAObj$Scores[,2],
                           z = PCAObj$Scores[,3],
                           color = as.factor(PCAObj$CombinedData[[Var]])) %>%
        plotly::layout(
          title = Var,
          scene = list(
            xaxis = list(title = colnames(PCAObj$Scores)[1]),
            yaxis = list(title = colnames(PCAObj$Scores)[2]),
            zaxis = list(title = colnames(PCAObj$Scores)[3])
          )
        )
    } else {
      p <- plotly::plot_ly(x = PCAObj$Scores[,1],
                           y = PCAObj$Scores[,2],
                           z = PCAObj$Scores[,3],
                           color = PCAObj$CombinedData[[Var]]) %>%
        plotly::layout(
          title = Var,
          scene = list(
            xaxis = list(title = colnames(PCAObj$Scores)[1]),
            yaxis = list(title = colnames(PCAObj$Scores)[2]),
            zaxis = list(title = colnames(PCAObj$Scores)[3])
          )
        )
    }
  }
  return(p)
}
