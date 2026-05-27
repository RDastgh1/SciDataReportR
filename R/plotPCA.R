#' Plot PCA scores
#'
#' Create an interactive 2D or 3D Plotly scatter plot from a PCA object created
#' by `CreatePCAObject()`.
#'
#' By default, the function plots the first three score columns in
#' `PCAObj$Scores` using a 3D plot. These are not assumed to be named `RC1`,
#' `RC2`, and `RC3`; the function uses the current column order in
#' `PCAObj$Scores`. Users may manually specify which components to plot using
#' `Components`.
#'
#' If `Mode = "auto"`, the function chooses 2D when two components are supplied
#' and 3D when three components are supplied. If `Components = NULL`, it defaults
#' to the first three score columns and creates a 3D plot.
#'
#' The function can optionally color points by a variable in
#' `PCAObj$CombinedData` and customize hover text using one or more variables.
#' When `Relabel = TRUE`, labels attached to hover variables are used in the
#' hover display when available.
#'
#' @param PCAObj A PCA object returned by `CreatePCAObject()`. Must contain
#'   `Scores` and `CombinedData`.
#' @param Var Optional character string naming a variable in
#'   `PCAObj$CombinedData` used to color points. Default is `NULL`.
#' @param t Deprecated compatibility argument for color type. Use `ColorType`
#'   instead. If set to `"Factor"`, `Var` is treated as categorical. Any other
#'   value treats `Var` as continuous. Default is `"NULL"`.
#' @param HoverVar Optional character string naming one variable in
#'   `PCAObj$CombinedData` to display in hover text. Retained for backward
#'   compatibility. Prefer `HoverVars` for new code. Default is `NULL`.
#' @param HoverVars Optional character vector naming one or more variables in
#'   `PCAObj$CombinedData` to display in hover text. If `NULL`, row number is
#'   shown. Default is `NULL`.
#' @param Components Optional character or numeric vector specifying which score
#'   columns to plot. Supply two components for 2D or three components for 3D.
#'   If `NULL`, the first three score columns are used.
#' @param Mode Character. Either `"auto"`, `"3D"`, or `"2D"`. If `"auto"`,
#'   the plot dimension is inferred from `Components`. If `Components = NULL`,
#'   `"auto"` defaults to 3D. Default is `"auto"`.
#' @param ColorType Character. Either `"auto"`, `"factor"`, or `"continuous"`.
#'   `"auto"` treats character, factor, logical, and labelled variables as
#'   categorical and numeric variables as continuous. Default is `"auto"`.
#' @param Relabel Logical. If `TRUE`, labels attached to hover variables are
#'   used in hover text when available. If `FALSE`, raw variable names are used.
#'   Default is `TRUE`.
#' @param Title Optional plot title. Default is `NULL`, which produces no title.
#'
#' @return A Plotly htmlwidget.
#'
#' @examples
#' \dontrun{
#' plotPCA(PCAObj)
#'
#' plotPCA(
#'   PCAObj,
#'   Components = c("RC1", "RC2"),
#'   Var = "SITE",
#'   ColorType = "factor"
#' )
#'
#' plotPCA(
#'   PCAObj,
#'   Components = c("RC2", "RC1", "RC4"),
#'   Mode = "3D"
#' )
#'
#' plotPCA(
#'   PCAObj,
#'   Components = c("RC1", "RC2"),
#'   HoverVars = c("SubjectID", "SITE", "Visit")
#' )
#'
#' plotPCA(
#'   PCAObj,
#'   Components = c("RC1", "RC2"),
#'   HoverVars = c("SubjectID", "SITE", "Visit"),
#'   Relabel = FALSE
#' )
#' }
#'
#' @export
plotPCA <- function(PCAObj,
                    Var = NULL,
                    t = "NULL",
                    HoverVar = NULL,
                    HoverVars = NULL,
                    Components = NULL,
                    Mode = c("auto", "3D", "2D"),
                    ColorType = c("auto", "factor", "continuous"),
                    Relabel = TRUE,
                    Title = NULL) {

  Mode <- match.arg(Mode)
  ColorType <- match.arg(ColorType)

  # Validate inputs

  if (!is.list(PCAObj)) {
    stop("PCAObj must be a list returned by CreatePCAObject().")
  }

  if (!"Scores" %in% names(PCAObj)) {
    stop("PCAObj must contain a Scores element.")
  }

  if (!"CombinedData" %in% names(PCAObj)) {
    stop("PCAObj must contain a CombinedData element.")
  }

  if (!is.logical(Relabel) || length(Relabel) != 1) {
    stop("Relabel must be a single logical value: TRUE or FALSE.")
  }

  Scores <- as.data.frame(PCAObj$Scores, check.names = FALSE)
  CombinedData <- as.data.frame(PCAObj$CombinedData, check.names = FALSE)

  if (is.null(colnames(Scores)) || any(colnames(Scores) == "")) {
    stop("PCAObj$Scores must have valid column names.")
  }

  if (anyDuplicated(colnames(Scores)) > 0) {
    stop(
      "PCAObj$Scores has duplicated column names: ",
      paste(unique(colnames(Scores)[duplicated(colnames(Scores))]), collapse = ", "),
      ". Score column names must be unique."
    )
  }

  if (nrow(Scores) != nrow(CombinedData)) {
    stop(
      "PCAObj$Scores and PCAObj$CombinedData must have the same number of rows. ",
      "Received ", nrow(Scores), " score rows and ",
      nrow(CombinedData), " combined data rows."
    )
  }

  if (!is.null(Var) && !Var %in% colnames(CombinedData)) {
    stop(
      "Var must be NULL or a column in PCAObj$CombinedData. ",
      "Missing column: ", Var
    )
  }

  if (!is.null(HoverVar)) {
    if (!is.character(HoverVar) || length(HoverVar) != 1) {
      stop("HoverVar must be NULL or a single column name in PCAObj$CombinedData.")
    }

    if (!HoverVar %in% colnames(CombinedData)) {
      stop(
        "HoverVar must be NULL or a column in PCAObj$CombinedData. ",
        "Missing column: ", HoverVar
      )
    }
  }

  if (!is.null(HoverVars)) {
    if (!is.character(HoverVars)) {
      stop("HoverVars must be NULL or a character vector of column names in PCAObj$CombinedData.")
    }

    missing_hover_vars <- setdiff(HoverVars, colnames(CombinedData))

    if (length(missing_hover_vars) > 0) {
      stop(
        "Some HoverVars are not columns in PCAObj$CombinedData: ",
        paste(missing_hover_vars, collapse = ", ")
      )
    }
  }

  if (!is.null(HoverVar) && is.null(HoverVars)) {
    HoverVars <- HoverVar
  }

  # Select components

  if (is.null(Components)) {
    if (ncol(Scores) < 3) {
      stop(
        "Components = NULL defaults to a 3D plot and requires at least 3 score columns. ",
        "PCAObj$Scores has ", ncol(Scores), ". Supply two components or use Mode = '2D'."
      )
    }

    Components <- colnames(Scores)[seq_len(3)]
  }

  if (is.numeric(Components)) {
    if (any(Components < 1) || any(Components > ncol(Scores))) {
      stop(
        "Numeric Components must be valid column positions in PCAObj$Scores. ",
        "PCAObj$Scores has ", ncol(Scores), " columns."
      )
    }

    Components <- colnames(Scores)[Components]
  }

  if (!is.character(Components)) {
    stop("Components must be NULL, a numeric vector, or a character vector.")
  }

  if (!length(Components) %in% c(2, 3)) {
    stop(
      "Components must contain exactly 2 components for a 2D plot or 3 components for a 3D plot. ",
      "Received ", length(Components), "."
    )
  }

  missing_components <- setdiff(Components, colnames(Scores))

  if (length(missing_components) > 0) {
    stop(
      "Some Components are not columns in PCAObj$Scores: ",
      paste(missing_components, collapse = ", ")
    )
  }

  if (Mode == "auto") {
    Mode <- if (length(Components) == 2) {
      "2D"
    } else {
      "3D"
    }
  }

  required_components <- if (Mode == "3D") {
    3
  } else {
    2
  }

  if (length(Components) != required_components) {
    stop(
      "Mode = '", Mode, "' requires exactly ",
      required_components, " components. Received ",
      length(Components), "."
    )
  }

  selected_scores <- Scores[, Components, drop = FALSE]

  score_is_numeric <- vapply(selected_scores, is.numeric, logical(1))

  if (!all(score_is_numeric)) {
    stop(
      "Selected PCA score columns must be numeric. Non-numeric selected columns: ",
      paste(Components[!score_is_numeric], collapse = ", ")
    )
  }

  # Prepare hover variable labels

  hover_var_labels <- HoverVars

  if (!is.null(HoverVars) && Relabel) {
    hover_var_labels <- vapply(
      HoverVars,
      function(x) {
        this_label <- sjlabelled::get_label(CombinedData[[x]])

        if (is.null(this_label) ||
            length(this_label) == 0 ||
            is.na(this_label) ||
            this_label == "") {
          x
        } else {
          as.character(this_label)
        }
      },
      character(1)
    )
  }

  if (!is.null(HoverVars) && !Relabel) {
    hover_var_labels <- HoverVars
  }

  # Build plotting data with internal names to avoid duplicate column issues

  PlotData <- data.frame(
    .pca_x = selected_scores[[1]],
    .pca_y = selected_scores[[2]],
    stringsAsFactors = FALSE
  )

  if (Mode == "3D") {
    PlotData$.pca_z <- selected_scores[[3]]
  }

  # Build hover text

  if (is.null(HoverVars)) {
    hover_text <- paste0("Row: ", seq_len(nrow(Scores)))
  } else {
    hover_text <- purrr::map_chr(
      seq_len(nrow(CombinedData)),
      function(i) {
        paste(
          paste0(
            hover_var_labels,
            ": ",
            vapply(
              HoverVars,
              function(x) {
                as.character(CombinedData[[x]][i])
              },
              character(1)
            )
          ),
          collapse = "<br>"
        )
      }
    )
  }

  hover_text <- paste0(
    hover_text,
    "<br>",
    Components[1], ": ", round(selected_scores[[1]], 3),
    "<br>",
    Components[2], ": ", round(selected_scores[[2]], 3)
  )

  if (Mode == "3D") {
    hover_text <- paste0(
      hover_text,
      "<br>",
      Components[3], ": ", round(selected_scores[[3]], 3)
    )
  }

  PlotData$.hover_text <- hover_text

  # Decide color type

  colour_vec <- NULL

  if (!is.null(Var)) {
    colour_vec <- CombinedData[[Var]]

    if (!identical(t, "NULL")) {
      if (identical(t, "Factor")) {
        ColorType <- "factor"
      } else {
        ColorType <- "continuous"
      }
    }

    if (ColorType == "auto") {
      if (is.factor(colour_vec) ||
          is.character(colour_vec) ||
          is.logical(colour_vec) ||
          inherits(colour_vec, "labelled")) {
        ColorType <- "factor"
      } else {
        ColorType <- "continuous"
      }
    }

    if (ColorType == "factor") {
      colour_vec <- as.factor(colour_vec)
    }
  }

  # Build plot

  if (Mode == "3D") {
    p <- plotly::plot_ly(
      data = PlotData,
      x = ~.pca_x,
      y = ~.pca_y,
      z = ~.pca_z,
      text = ~.hover_text,
      color = colour_vec,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 5),
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
      plotly::layout(
        title = Title,
        scene = list(
          xaxis = list(title = Components[1]),
          yaxis = list(title = Components[2]),
          zaxis = list(title = Components[3])
        )
      )
  }

  if (Mode == "2D") {
    p <- plotly::plot_ly(
      data = PlotData,
      x = ~.pca_x,
      y = ~.pca_y,
      text = ~.hover_text,
      color = colour_vec,
      type = "scatter",
      mode = "markers",
      hovertemplate = "%{text}<extra></extra>"
    ) %>%
      plotly::layout(
        title = Title,
        xaxis = list(title = Components[1]),
        yaxis = list(title = Components[2])
      )
  }

  return(p)
}
