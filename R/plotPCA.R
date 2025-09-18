#' Plot 3D PCA Scores
#'
#' Creates an interactive 3-D scatter plot of PCA scores with Plotly.
#' You can colour points by a variable **and** customise the hover label.
#'
#' @param PCAObj A list containing:
#'   * `Scores` – a matrix/data-frame with the PCA scores.
#'   * `CombinedData` – the original data used to calculate the PCA.
#' @param Var   (string) Optional variable in `CombinedData` to colour by.
#'              *Default:* `NULL` (no colouring).
#' @param t     (string)  If colouring, set to `"Factor"` for categorical;
#'              anything else is treated as continuous. *Default:* `"NULL"`.
#' @param HoverVar (string) Optional variable in `CombinedData` whose value
#'              will be displayed on hover.
#'              *Default:* `NULL`, which shows the row number.
#'
#' @return A Plotly htmlwidget.
#' @export
#'
plotPCA <- function(PCAObj, Var = NULL, t = "NULL", HoverVar = NULL) {

  ## ── 1. Build the hover text ────────────────────────────────────────────────
  hover_text <- if (is.null(HoverVar)) {
    # default: row numbers
    paste0("Row: ", seq_len(nrow(PCAObj$CombinedData)))
  } else {
    # user-supplied column
    paste0(HoverVar, ": ", PCAObj$CombinedData[[HoverVar]])
  }

  ## ── 2. Helper to create the plotly object ─────────────────────────────────
  make_plot <- function(colour_vec = NULL, colour_title = NULL) {
    plotly::plot_ly(
      x    = PCAObj$Scores[, 1],
      y    = PCAObj$Scores[, 2],
      z    = PCAObj$Scores[, 3],
      text = hover_text,
      color = colour_vec,
      hovertemplate = paste(
        "%{text}<br>",
        "PC1: %{x:.3f}<br>",
        "PC2: %{y:.3f}<br>",
        "PC3: %{z:.3f}<extra></extra>"
      )
    ) |>
      plotly::layout(
        title = colour_title,
        scene = list(
          xaxis = list(title = colnames(PCAObj$Scores)[1]),
          yaxis = list(title = colnames(PCAObj$Scores)[2]),
          zaxis = list(title = colnames(PCAObj$Scores)[3])
        )
      )
  }

  ## ── 3. Decide colouring logic & build plot ────────────────────────────────
  if (is.null(Var)) {
    p <- make_plot()
  } else if (t == "Factor") {
    p <- make_plot(
      colour_vec   = as.factor(PCAObj$CombinedData[[Var]]),
      colour_title = Var
    )
  } else {
    p <- make_plot(
      colour_vec   = PCAObj$CombinedData[[Var]],
      colour_title = Var
    )
  }

  return(p)
}
