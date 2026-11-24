#' Create directional heatmaps across continuous & binary variables
#'
#' Combines:
#'  - Continuous~Continuous (Pearson/Spearman)
#'  - Binary~Binary (Phi; 1 == PositiveLevel)
#'  - Binary~Continuous (r_pb; 1 == PositiveLevel)
#' into a square or rectangular heatmap with raw and FDR-star overlays.
#'
#' Constant variables (no variation in the current data) are
#' automatically excluded before computing any tiles.
#'
#' @param data A dataframe.
#' @param covariates Optional character vector naming numeric or categorical
#'   adjustment variables. Covariates are excluded from the heatmap axes;
#'   explicitly requesting them on an axis produces a warning. Binary variables
#'   retain their positive-level 0/1 coding. Adjusted estimates are partial
#'   Pearson correlations, using pair-specific complete cases and residual test
#'   degrees of freedom `n - rank(design) - 1` (including the intercept).
#'   Failed adjustments warn and are marked failed in `Results` and `Metadata`;
#'   fallback estimates do not carry adjusted significance.
#'   `Unadjusted` refers to raw p-values, including when estimates adjust for
#'   covariates. Without covariates, existing estimates and tests are unchanged.
#' @param variables Character vector of variables to include (subset of `data`
#'   columns). The analysis is symmetric: every variable is related to every
#'   other, so a single variable set defines both axes. If NULL, uses all
#'   detected continuous + binary vars.
#' @param predictor_vars,outcome_vars Optional column (x-axis) and row (y-axis)
#'   variable sets, respectively, matching historical `xVars` and `yVars`.
#'   Supply both for rectangular mode; `variables` (with a warning) and
#'   `triangle` are ignored. Self-pairs are blank. A single set uses both axes.
#'   Legacy positional character vectors in arguments two and three mean
#'   columns and rows; a logical third argument still means `Relabel`.
#' @param Relabel Logical; use sjlabelled variable labels if present.
#' @param Ordinal Logical; passed to [PlotPointCorrelationsHeatmap()] for the
#'   binary~continuous block, where it controls whether ordinal variables are
#'   treated as continuous. Defaults to `TRUE`.
#' @param fdr_scope One of `"matrix"` (default), `"per_outcome"`, or
#'   `"per_predictor"`, threaded
#'   through to the three sub-analyses ([PlotCorrelationsHeatmap()],
#'   [PlotPhiHeatmap()], [PlotPointCorrelationsHeatmap()]). Correction is
#'   applied within each sub-analysis block (continuous~continuous,
#'   binary~binary, binary~continuous), matching historical behavior; each
#'   sub-function's documented outcome orientation applies within its block.
#'   In two-set mode only requested pairs enter these families; point-biserial
#'   outcomes remain the continuous variables regardless of display axis.
#' @param triangle Ignored in two-set mode. Display `"full"` (default), `"upper"`, or `"lower"` half
#'   of a symmetric directional matrix. This affects plots only and is applied
#'   only when the same variables occur on both axes.
#' @param Data \strong{Deprecated} (since 19.15.0). Use \code{data} instead.
#' @param xVars \strong{Deprecated} (since 19.15.0). Use \code{predictor_vars}
#'   instead.
#' @param yVars \strong{Deprecated} (since 19.15.0). Use \code{outcome_vars}
#'   instead.
#' @return list(Unadjusted, FDRCorrected, Relabel, BinaryMapping, Excluded)
#' @note
#' The analysis covers continuous and binary variables. Multi-level categorical
#' variables (more than two levels) are not placed on the heatmap. `Ordinal`
#' affects only the binary~continuous block, which is the one sub-analysis that
#' accepts it.
#'
#' @examples
#' data(SampleData)
#' data(SampleVariableTypes)
#'
#' # Attach labels and factor levels so binary variables are detected and
#' # axis labels are readable
#' Labelled <- RevalueData(SampleData, SampleVariableTypes)$RevaluedData
#'
#' # A mix of binary categorical (Diagnosis, sex) and continuous variables
#' result <- PlotDirectionalHeatmaps(
#'   Labelled,
#'   variables = c("Diagnosis", "sex", "age", "AXL", "Adiponectin",
#'                 "Alpha_1_Antitrypsin", "C_Reactive_Protein", "Cortisol",
#'                 "Insulin", "Leptin")
#' )
#'
#' # Raw p-value directional heatmap
#' result$Unadjusted$plot
#'
#' # FDR-adjusted directional heatmap
#' result$FDRCorrected$plot
#'
#' # Show each symmetric variable pair once.
#' upper <- PlotDirectionalHeatmaps(
#'   Labelled,
#'   variables = c("Diagnosis", "sex", "age", "AXL", "Adiponectin",
#'                 "Alpha_1_Antitrypsin", "C_Reactive_Protein", "Cortisol",
#'                 "Insulin", "Leptin"),
#'   triangle = "upper"
#' )
#' upper$FDRCorrected$plot
#'
#' # Different column and row sets (historical xVars/yVars orientation).
#' rectangular <- PlotDirectionalHeatmaps(
#'   Labelled, predictor_vars = c("Diagnosis", "age"),
#'   outcome_vars = c("sex", "AXL", "Leptin")
#' )
#' rectangular$Unadjusted$plot
#'
#' # Adjust every correlation block for age and sex.
#' adjusted <- PlotDirectionalHeatmaps(
#'   Labelled, variables = c("Diagnosis", "AXL", "Leptin"),
#'   covariates = c("age", "sex")
#' )
#' adjusted$FDRCorrected$plot
#'
#' # How binary variables were coded (which level counts as the positive one)
#' result$BinaryMapping
#' @param IncludeTrendSymbol Logical; add a plus sign for `0.05 < p <= 0.10`.
#' @param cluster_rows,cluster_columns Logical; cluster heatmap rows and columns.
#' @export
#' @inheritParams inference-results
#' @section Statistical results and provenance:
#' See [inference-results] for the shared `Results` and `Metadata` contract.
PlotDirectionalHeatmaps <- function(data,
    variables = NULL,
    Relabel = TRUE,
    Ordinal = TRUE,
    fdr_scope = c("matrix", "per_outcome", "per_predictor"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    triangle = c("full", "upper", "lower"),
    IncludeTrendSymbol = FALSE,
    Pthresh = 0.05,
    PVar = "PValue",
    predictor_vars = NULL,
    outcome_vars = NULL,
    Data = lifecycle::deprecated(),
    xVars = lifecycle::deprecated(),
    yVars = lifecycle::deprecated(),
    covariates = NULL) {
  .scidr_inference <- ScidrStartInference("PlotDirectionalHeatmaps", match.call(), environment())
  on.exit({
    if (!is.null(returnValue())) {
      return(ScidrFinishInference(returnValue(), .scidr_inference, environment()))
    }
  }, add = TRUE)

  # Deprecated argument shims (SciDataReportR 19.15.0)
  if (lifecycle::is_present(Data)) {
    lifecycle::deprecate_warn("19.15.0", "PlotDirectionalHeatmaps(Data)", "PlotDirectionalHeatmaps(data)")
    data <- Data
  }
  if (!missing(data)) Data <- data
  if (lifecycle::is_present(xVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotDirectionalHeatmaps(xVars)", "PlotDirectionalHeatmaps(predictor_vars)")
  } else {
    xVars <- NULL
  }
  if (lifecycle::is_present(yVars)) {
    lifecycle::deprecate_warn("19.15.0", "PlotDirectionalHeatmaps(yVars)", "PlotDirectionalHeatmaps(outcome_vars)")
  } else {
    yVars <- NULL
  }
  # Historical positional calls used (data, xVars, yVars).
  if (is.character(variables) && is.character(Relabel)) {
    if (is.null(predictor_vars)) predictor_vars <- variables
    if (is.null(outcome_vars)) outcome_vars <- Relabel
    variables <- NULL
    Relabel <- TRUE
  }
  if (is.null(predictor_vars)) predictor_vars <- xVars
  if (is.null(outcome_vars)) outcome_vars <- yVars
  two_set <- !is.null(predictor_vars) && !is.null(outcome_vars)
  if (two_set && !is.null(variables)) {
    warning("`variables` is ignored when both `predictor_vars` and `outcome_vars` are supplied.", call. = FALSE)
  }
  if (!is.null(predictor_vars) || !is.null(outcome_vars)) {
    xVars <- predictor_vars
    yVars <- outcome_vars
    if (is.null(xVars)) xVars <- yVars
    if (is.null(yVars)) yVars <- xVars
  } else {
    xVars <- variables
    yVars <- variables
  }
  xVars <- unique(xVars)
  yVars <- unique(yVars)
  fdr_scope <- match.arg(fdr_scope)
  triangle <- if (two_set) "full" else match.arg(triangle)
  if (!is.logical(cluster_rows) || length(cluster_rows) != 1 || is.na(cluster_rows)) {
    stop("cluster_rows must be TRUE or FALSE.")
  }
  if (!is.logical(cluster_columns) || length(cluster_columns) != 1 || is.na(cluster_columns)) {
    stop("cluster_columns must be TRUE or FALSE.")
  }


  # ---- Validate inputs ----------------------------------------------------
  if (!is.data.frame(Data)) {
    stop("Data must be a data.frame.", call. = FALSE)
  }

  covariates <- ScidrValidateVariables(Data, covariates, "covariates")
  vars_Overlap <- intersect(unique(c(xVars, yVars)), covariates)
  if (length(vars_Overlap)) {
    warning("Excluding covariates from heatmap axes: ",
      paste(vars_Overlap, collapse = ", "), call. = FALSE)
  }

  # ---- Prepare data: protect against make.names() behavior downstream -----
  orig_names <- colnames(Data)
  safe_names <- make.names(orig_names, unique = TRUE)

  # Map original -> safe and safe -> original
  map_orig_to_safe <- stats::setNames(safe_names, orig_names)
  map_safe_to_orig <- stats::setNames(orig_names, safe_names)

  Data_safe <- Data
  colnames(Data_safe) <- safe_names

  # Helper to map vectors while preserving unknowns
  to_safe <- function(v) {
    v <- as.character(v)
    out <- unname(map_orig_to_safe[v])
    out[is.na(out)] <- v[is.na(out)]
    out
  }
  to_orig <- function(v) {
    v <- as.character(v)
    out <- unname(map_safe_to_orig[v])
    out[is.na(out)] <- v[is.na(out)]
    out
  }

  # ---- helpers ------------------------------------------------------------
  uniq2 <- function(x) length(unique(stats::na.omit(x)))

  has_two_levels <- function(v_safe) {
    uniq2(Data_safe[[v_safe]]) == 2
  }

  has_variance <- function(v_safe) {
    x <- Data_safe[[v_safe]]
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    length(x) >= 2 && stats::var(x) > 0
  }

  # ---- 1) Determine candidate vars ---------------------------------------
  if (is.null(xVars)) {
    Cont_all_safe <- getNumVars(Data_safe)
    Cat_all_safe  <- getBinaryVars(Data_safe)
    xVars_safe <- colnames(Data_safe)[colnames(Data_safe) %in% c(Cat_all_safe, Cont_all_safe)]
    xVars_orig <- to_orig(xVars_safe)
  } else {
    xVars_orig <- as.character(xVars)
    xVars_safe <- to_safe(xVars_orig)
  }

  if (is.null(yVars)) {
    yVars_safe <- xVars_safe
    yVars_orig <- xVars_orig
  } else {
    yVars_orig <- as.character(yVars)
    yVars_safe <- to_safe(yVars_orig)
  }

  vars_CovariatesSafe <- to_safe(covariates)
  xVars_safe <- setdiff(xVars_safe, vars_CovariatesSafe)
  yVars_safe <- setdiff(yVars_safe, vars_CovariatesSafe)
  cand_safe <- unique(c(xVars_safe, yVars_safe))

  # Guard against truly missing vars (in either naming universe)
  missing_safe <- setdiff(cand_safe, colnames(Data_safe))
  if (length(missing_safe) > 0) {
    missing_orig <- to_orig(missing_safe)
    stop(
      "Variables not found in `data`: ",
      paste(missing_orig, collapse = ", "),
      " (supplied to `variables`).",
      call. = FALSE
    )
  }

  # initial typing from your helpers, but we'll enforce non-constant below
  Cat0_safe  <- intersect(getBinaryVars(dplyr::select(Data_safe, dplyr::all_of(cand_safe))), cand_safe)
  Cont0_safe <- intersect(getNumVars(   dplyr::select(Data_safe, dplyr::all_of(cand_safe))), cand_safe)

  # ---- 2) Drop constants BEFORE any subcalls ------------------------------
  CatVars_safe  <- Cat0_safe[vapply(Cat0_safe, has_two_levels, logical(1))]
  ContVars_safe <- Cont0_safe[vapply(Cont0_safe, has_variance, logical(1))]

  dropped_safe <- setdiff(cand_safe, c(CatVars_safe, ContVars_safe))
  if (length(dropped_safe)) {
    reasons <- vapply(dropped_safe, function(vs) {
      if (vs %in% Cont0_safe && !has_variance(vs))    return("no variance")
      if (vs %in% Cat0_safe  && !has_two_levels(vs))  return("not two levels")
      if (!(vs %in% c(Cat0_safe, Cont0_safe)))        return("not binary/continuous")
      "excluded"
    }, character(1))

    Excluded <- tibble::tibble(
      Variable = to_orig(dropped_safe),
      Reason   = reasons
    )
  } else {
    Excluded <- tibble::tibble(Variable = character(), Reason = character())
  }

  # refresh axes after exclusions
  xVars_safe <- xVars_safe[xVars_safe %in% c(CatVars_safe, ContVars_safe)]
  yVars_safe <- yVars_safe[yVars_safe %in% c(CatVars_safe, ContVars_safe)]

  xVars_orig <- to_orig(xVars_safe)
  yVars_orig <- to_orig(yVars_safe)

  # Build a single binary mapping (stable 1 == PositiveLevel)
  BinaryMapping <- NULL
  if (length(CatVars_safe) > 0) {
    BinaryMapping <- createBinaryMapping(Data_safe, CatVars_safe)

    # Best-effort: if mapping is a data frame with a Variable column, map back to original names
    if (is.data.frame(BinaryMapping) && "Variable" %in% names(BinaryMapping)) {
      BinaryMapping$Variable <- to_orig(BinaryMapping$Variable)
    }

    # Best-effort: if mapping is a named list keyed by variables, rename to original
    if (is.list(BinaryMapping) && !is.data.frame(BinaryMapping) && !is.null(names(BinaryMapping))) {
      nm <- names(BinaryMapping)
      names(BinaryMapping) <- to_orig(nm)
    }
  }

  binary_map_safe <- BinaryMapping
  if (!is.null(binary_map_safe)) binary_map_safe$Variable <- to_safe(binary_map_safe$Variable)

  # ---- 3) Collect tiles from each sub-plotter -----------------------------
  df_Combined <- tibble::tibble(
    XVar        = character(),
    YVar        = character(),
    correlation = numeric(),
    test        = character()
  )

  if (two_set) {
    # Compute each requested cell once; no union-by-union matrix is built.
    df_Pairs <- tidyr::expand_grid(XVar = xVars_safe, YVar = yVars_safe)
    if (nrow(df_Pairs)) df_Combined <- dplyr::bind_rows(lapply(seq_len(nrow(df_Pairs)), function(i) {
      column <- df_Pairs$XVar[i]
      row <- df_Pairs$YVar[i]
      column_binary <- column %in% CatVars_safe
      row_binary <- row %in% CatVars_safe
      test <- if (column_binary && row_binary) "Phi" else if (column_binary || row_binary) "Point Correlation" else "pearson"
      estimate <- p_value <- NA_real_
      n_pairs <- sum(stats::complete.cases(Data_safe[, unique(c(column, row, vars_CovariatesSafe)), drop = FALSE]))
      if (column != row) {
        if (test == "pearson") {
          pair <- PlotCorrelationsHeatmap(Data_safe, predictor_vars = column,
            outcome_vars = row, Relabel = Relabel, IncludeTrendSymbol = IncludeTrendSymbol, covariates = vars_CovariatesSafe)
          estimate <- pair$Unadjusted$r[1, 1]
          p_value <- pair$Unadjusted$p[1, 1]
          n_pairs <- pair$Unadjusted$npairs[1, 1]
        } else if (test == "Phi") {
          pair <- PlotPhiHeatmap(Data_safe, c(column, row), Relabel = Relabel,
            binary_map = binary_map_safe, IncludeTrendSymbol = IncludeTrendSymbol, covariates = vars_CovariatesSafe)
          df_Cell <- pair$Unadjusted$PvalTable
          df_Cell <- df_Cell[df_Cell$XVar == column & df_Cell$YVar == row, ]
          estimate <- df_Cell$Phi
          p_value <- df_Cell$p_value
          n_pairs <- df_Cell$nPairs
        } else {
          binary <- if (column_binary) column else row
          continuous <- if (column_binary) row else column
          pair <- PlotPointCorrelationsHeatmap(Data_safe, binary, continuous,
            Relabel = Relabel, Ordinal = Ordinal, binary_map = binary_map_safe,
            IncludeTrendSymbol = IncludeTrendSymbol, covariates = vars_CovariatesSafe)
          df_Cell <- pair$Unadjusted$PvalTable
          if (nrow(df_Cell)) {
            estimate <- df_Cell$correlation
            p_value <- df_Cell$p_value
            n_pairs <- df_Cell$nPairs
          }
        }
      }
      Label <- function(variable) {
        if (!Relabel) return(variable)
        label <- sjlabelled::get_label(Data_safe[[variable]])
        if (is.null(label) || is.na(label) || label == "") variable else label
      }
      tibble::tibble(XVar = column, YVar = row, correlation = estimate,
        test = test, P = p_value, nPairs = n_pairs,
        XLabel = Label(column), YLabel = Label(row))
    }))
    if (nrow(df_Combined)) {
      df_Combined <- df_Combined %>%
        dplyr::group_by(test) %>%
        dplyr::mutate(
          P_adj = ApplyFDRCorrection(P, fdr_scope = fdr_scope,
            outcome_ids = ifelse(test == "Point Correlation" & XVar %in% ContVars_safe, XVar, YVar),
            predictor_ids = ifelse(test == "Point Correlation" & XVar %in% ContVars_safe, YVar, XVar)),
          stars = .FormatPValueStars(ScidrInferenceP(P, P_adj), IncludeTrendSymbol = IncludeTrendSymbol),
          stars_FDR = .FormatPValueStars(P_adj, IncludeTrendSymbol = IncludeTrendSymbol)
        ) %>%
        dplyr::ungroup()
    }
  } else {
  # Continuous ~ Continuous
  if (length(ContVars_safe) > 0) {
    O_ContCont <- PlotCorrelationsHeatmap(Data_safe, ContVars_safe, Relabel = Relabel,
                                          fdr_scope = fdr_scope,
                                          IncludeTrendSymbol = IncludeTrendSymbol, covariates = vars_CovariatesSafe)
    df_ContCont <- O_ContCont$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = R, test = O_ContCont$method)

    df_Combined <- dplyr::bind_rows(df_Combined, df_ContCont)
  }

  # Binary ~ Binary (Phi)
  if (length(CatVars_safe) > 1) {
    O_CatCat <- PlotPhiHeatmap(Data_safe, CatVars_safe, Relabel = Relabel,
                               binary_map = binary_map_safe, fdr_scope = fdr_scope,
                               IncludeTrendSymbol = IncludeTrendSymbol, covariates = vars_CovariatesSafe)
    df_CatCat <- O_CatCat$Unadjusted$plot$data %>%
      dplyr::mutate(correlation = Phi, test = "Phi")

    df_Combined <- dplyr::bind_rows(df_Combined, df_CatCat)
  }

  # Binary ~ Continuous (point-biserial)
  if (length(CatVars_safe) > 0 && length(ContVars_safe) > 0) {
    O_CatCont <- PlotPointCorrelationsHeatmap(
      Data_safe, CatVars_safe, ContVars_safe,
      Relabel = Relabel, Ordinal = Ordinal, binary_map = binary_map_safe,
      fdr_scope = fdr_scope, IncludeTrendSymbol = IncludeTrendSymbol,
      covariates = vars_CovariatesSafe
    )

    df_CatCont <- O_CatCont$Unadjusted$plot$data %>%
      dplyr::mutate(stars = `p<.05`, stars_FDR = p.adj.signif) %>%
      dplyr::select(-`p<.05`, -p.adj.signif) %>%
      dplyr::mutate(xVar = CategoricalVariable, yVar = ContinuousVariable)

    # make square by duplicating with swapped axes/labels
    df_CatCont1 <- df_CatCont %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable) %>%
      dplyr::mutate(XVar = xVar, YVar = yVar) %>%
      dplyr::select(-xVar, -yVar)

    df_CatCont2 <- df_CatCont %>%
      dplyr::mutate(XVar = yVar, YVar = xVar) %>%
      dplyr::select(-xVar, -yVar) %>%
      dplyr::mutate(XL = XLabel, YL = YLabel) %>%
      dplyr::mutate(YLabel = XL, XLabel = YL) %>%
      dplyr::select(-CategoricalVariable, -ContinuousVariable, -XL, -YL)

    df_CatContSquare <- dplyr::bind_rows(df_CatCont1, df_CatCont2)
    df_CatContSquare$test <- "Point Correlation"

    df_Combined <- dplyr::bind_rows(df_Combined, df_CatContSquare)
  }

  }

  # Child adjustment failures apply to both display orientations. Suppress
  # fallback inference before drawing stars or forming correction families.
  if (length(covariates) && nrow(df_Combined)) {
    df_Failed <- dplyr::bind_rows(lapply(.scidr_inference$Children,
      function(child) child$Metadata$Failed))
    if (nrow(df_Failed)) {
      failed <- rep(FALSE, nrow(df_Combined))
      for (i in seq_len(nrow(df_Failed))) {
        failed <- failed |
          (df_Combined$XVar == df_Failed$Predictor[i] & df_Combined$YVar == df_Failed$Outcome[i]) |
          (df_Combined$YVar == df_Failed$Predictor[i] & df_Combined$XVar == df_Failed$Outcome[i])
      }
      for (column in intersect(c("P", "p_value", "P_adj", "p.adj"), names(df_Combined))) {
        df_Combined[[column]][failed] <- NA_real_
      }
      for (column in c("P_adj", "stars", "stars_FDR")) {
        if (!column %in% names(df_Combined)) {
          df_Combined[[column]] <- if (column == "P_adj") {
            rep(NA_real_, nrow(df_Combined))
          } else {
            rep("", nrow(df_Combined))
          }
        }
      }
      for (block in unique(df_Combined$test)) {
        indices <- which(df_Combined$test == block)
        df_Block <- df_Combined[indices, , drop = FALSE]
        raw_p <- if (two_set || block == "pearson") df_Block$P else df_Block$p_value
        if (!two_set && block == "pearson") {
          # Continuous square blocks correct unique symmetric pairs.
          vars_Block <- unique(c(df_Block$XVar, df_Block$YVar))
          p_matrix <- matrix(NA_real_, length(vars_Block), length(vars_Block),
            dimnames = list(vars_Block, vars_Block))
          coordinates <- cbind(match(df_Block$XVar, vars_Block), match(df_Block$YVar, vars_Block))
          p_matrix[coordinates] <- raw_p
          adjusted <- ApplyFDRCorrection(p_matrix, fdr_scope = fdr_scope)[coordinates]
        } else {
          outcome <- ifelse(block == "Point Correlation" & df_Block$XVar %in% ContVars_safe,
            df_Block$XVar, df_Block$YVar)
          predictor <- ifelse(block == "Point Correlation" & df_Block$XVar %in% ContVars_safe,
            df_Block$YVar, df_Block$XVar)
          keys <- paste(predictor, outcome, sep = "\r")
          keep <- if (!two_set && block == "Point Correlation") !duplicated(keys) else rep(TRUE, length(keys))
          adjusted <- ApplyFDRCorrection(raw_p[keep], fdr_scope = fdr_scope,
            outcome_ids = outcome[keep], predictor_ids = predictor[keep])
          if (!all(keep)) adjusted <- adjusted[match(keys, keys[keep])]
        }
        df_Combined$P_adj[indices] <- adjusted
        df_Combined$stars[indices] <- .FormatPValueStars(
          ScidrInferenceP(raw_p, adjusted), IncludeTrendSymbol = IncludeTrendSymbol)
        df_Combined$stars_FDR[indices] <- .FormatPValueStars(adjusted,
          IncludeTrendSymbol = IncludeTrendSymbol)
        if ("p.adj" %in% names(df_Combined)) df_Combined$p.adj[indices] <- adjusted
      }
    }
  }

  # ---- 4) Map back to original var names and filter to requested axes -----
  # Ensure we have XVar/YVar columns and convert them back to original names
  if ("XVar" %in% names(df_Combined)) df_Combined$XVar <- to_orig(df_Combined$XVar)
  if ("YVar" %in% names(df_Combined)) df_Combined$YVar <- to_orig(df_Combined$YVar)

  df_Combined_plot <- df_Combined %>%
    dplyr::filter(XVar %in% xVars_orig, YVar %in% yVars_orig)

  # handle fully-empty safely
  if (nrow(df_Combined_plot) == 0L) {
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::ggtitle("No valid variable pairs after excluding constants")

    Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)
    FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = empty_plot)

    return(list(
      Unadjusted   = Unadjusted,
      FDRCorrected = FDRCorrected,
      Relabel      = Relabel,
      AxisOrder    = list(
        rows = character(0), columns = character(0),
        cluster_rows = cluster_rows, cluster_columns = cluster_columns
      ),
      BinaryMapping= BinaryMapping,
      Excluded     = Excluded
    ))
  }

  # Lock label order using the original variable order
  ordered_xlabels <- vapply(xVars_orig, function(var) {
    as.character(df_Combined_plot$XLabel[df_Combined_plot$XVar == var][1])
  }, character(1))
  ordered_xlabels <- unique(ordered_xlabels[!is.na(ordered_xlabels)])

  ordered_ylabels <- vapply(yVars_orig, function(var) {
    as.character(df_Combined_plot$YLabel[df_Combined_plot$YVar == var][1])
  }, character(1))
  ordered_ylabels <- unique(ordered_ylabels[!is.na(ordered_ylabels)])

  triangle_display <- .ResolveHeatmapTriangle(
    data = df_Combined_plot, row_id = "YVar", column_id = "XVar",
    value = "correlation", row_order = yVars_orig, column_order = xVars_orig,
    cluster_rows = cluster_rows, cluster_columns = cluster_columns,
    triangle = triangle
  )
  AxisOrder <- triangle_display$AxisOrder
  df_Combined_plot_display <- triangle_display$PlotData
  if (two_set) {
    df_Combined_plot_display <- df_Combined_plot_display %>%
      dplyr::filter(XVar != YVar)
  }
  x_axis_labels <- stats::setNames(ordered_xlabels, xVars_orig)
  y_axis_labels <- stats::setNames(ordered_ylabels, yVars_orig)
  df_Combined_plot$XOrder <- factor(df_Combined_plot$XVar, levels = AxisOrder$columns)
  df_Combined_plot$YOrder <- factor(df_Combined_plot$YVar, levels = rev(AxisOrder$rows))
  df_Combined_plot_display$XOrder <- factor(df_Combined_plot_display$XVar, levels = AxisOrder$columns)
  df_Combined_plot_display$YOrder <- factor(df_Combined_plot_display$YVar, levels = rev(AxisOrder$rows))

  # ---- 5) Build layered heatmap ------------------------------------------
  p <- ggplot2::ggplot()

  if (nrow(dplyr::filter(df_Combined_plot_display, test == "pearson")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot_display, test == "pearson"),
        ggplot2::aes(x = XOrder, y = YOrder, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = if (length(covariates)) "Partial correlation" else "r") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot_display, test == "spearman")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot_display, test == "spearman"),
        ggplot2::aes(x = XOrder, y = YOrder, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1), name = "\u03C1") +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot_display, test == "Phi")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot_display, test == "Phi"),
        ggplot2::aes(x = XOrder, y = YOrder, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(
        limits = c(-1, 1),
        name   = if (length(covariates)) "Partial correlation" else expression(Phi),  # plotmath: portable under C / non-UTF8 locale
        low    = scales::muted("purple"),
        high   = scales::muted("green")
      ) +
      ggnewscale::new_scale_fill()
  }

  if (nrow(dplyr::filter(df_Combined_plot_display, test == "Point Correlation")) > 0) {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(df_Combined_plot_display, test == "Point Correlation"),
        ggplot2::aes(x = XOrder, y = YOrder, fill = correlation)
      ) +
      ggplot2::scale_fill_gradient2(
        limits = c(-1, 1),
        name   = if (length(covariates)) "Partial correlation" else expression(r[pb]),
        low    = scales::muted("#FFA500"),
        high   = scales::muted("#008080")
      ) +
      ggnewscale::new_scale_fill()
  }

  # reset limits to preserve ordering
  p <- p +
    ggplot2::scale_x_discrete(limits = AxisOrder$columns, labels = x_axis_labels) +
    ggplot2::scale_y_discrete(limits = rev(AxisOrder$rows), labels = y_axis_labels)

  # annotate with raw and FDR stars
  p_raw <- p +
    ggplot2::geom_text(
      data = df_Combined_plot_display,
      ggplot2::aes(x = XOrder, y = YOrder, label = stars),
      color = "black"
    ) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 7, color = "grey40"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text  = ggplot2::element_text(size = 8),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(5.5, 5.5, 26, 5.5)
    )

  p_FDR <- p +
    ggplot2::geom_text(
      data = df_Combined_plot_display,
      ggplot2::aes(x = XOrder, y = YOrder, label = stars_FDR),
      color = "black"
    ) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 7, color = "grey40"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text  = ggplot2::element_text(size = 8),
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = ggplot2::margin(5.5, 5.5, 26, 5.5)
    )

  # ---- 6) Return ----------------------------------------------------------
  p_raw <- p_raw + ggplot2::labs(subtitle = "Significance: raw p-values (no multiple-comparison correction)")
  p_FDR <- p_FDR + ggplot2::labs(subtitle = ScidrAdjustedSignificanceSubtitle(fdr_scope, family = "separately within each test type"))

  Unadjusted   <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_raw)
  FDRCorrected <- list(Relabel = Relabel, data = df_Combined_plot, plot = p_FDR)

  out <- list(
    Unadjusted    = Unadjusted,
    FDRCorrected  = FDRCorrected,
    Relabel       = Relabel,
    AxisOrder     = AxisOrder,
    Triangle      = triangle,
    BinaryMapping = BinaryMapping,
    Excluded      = Excluded
  )
  # Standardized p-value element aliases (old names kept)
  out$p <- out$Unadjusted
  out$p_fdr <- out$FDRCorrected
  out
}
