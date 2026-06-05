#' Plot merge validation diagnostics
#'
#' Create diagnostic plots from a `ValidateMerge()` result object. This function
#' visualizes key merge-audit outputs, including validation check status, key
#' coverage, join-variable auditing, duplicate-variable agreement, and
#' duplicate-variable conflict counts. Use this after running `ValidateMerge()`
#' to quickly inspect whether a merged dataset appears trustworthy.
#'
#' The function expects the object returned by `ValidateMerge()`. It does not
#' re-run any merge validation checks.
#'
#' @param MergeObj A list returned by `ValidateMerge()`.
#' @param Plot Character value specifying which plot to return. Options are
#'   `"All"`, `"Checks"`, `"Coverage"`, `"JoinAudit"`, `"Agreement"`, and
#'   `"Conflicts"`. Default is `"All"`.
#' @param Interactive Logical; if `TRUE`, plots are converted to interactive
#'   `plotly` objects using `plotly::ggplotly()`. Default is `TRUE`.
#'
#' @return If `Plot = "All"`, a named list of plots. Otherwise, a single plot
#' object. Plot objects are either `ggplot` objects or `plotly` htmlwidgets,
#' depending on `Interactive`.
#'
#' @export
PlotMergeValidation <- function(
    MergeObj,
    Plot = c(
      "All",
      "Checks",
      "Coverage",
      "JoinAudit",
      "Agreement",
      "Conflicts"
    ),
    Interactive = TRUE
) {

  # Validate inputs

  Plot <- match.arg(Plot)

  if (!is.list(MergeObj)) {
    stop("MergeObj must be a list returned by ValidateMerge().")
  }

  required_elements <- c(
    "Checks",
    "IDCoverage",
    "JoinAudit",
    "DuplicateVariables",
    "Summary"
  )

  missing_elements <- setdiff(
    required_elements,
    names(MergeObj)
  )

  if (length(missing_elements) > 0) {
    stop(
      "MergeObj is missing required element(s): ",
      paste(missing_elements, collapse = ", "),
      ". Please recreate the object using the latest version of ValidateMerge()."
    )
  }

  if (is.null(MergeObj$Checks)) {
    stop("MergeObj$Checks is NULL. Please recreate the object using the latest version of ValidateMerge().")
  }

  if (is.null(MergeObj$IDCoverage)) {
    stop("MergeObj$IDCoverage is NULL. Please recreate the object using the latest version of ValidateMerge().")
  }

  if (is.null(MergeObj$JoinAudit)) {
    stop("MergeObj$JoinAudit is NULL. Please recreate the object using the latest version of ValidateMerge().")
  }

  if (is.null(MergeObj$DuplicateVariables)) {
    stop("MergeObj$DuplicateVariables is NULL. Please recreate the object using the latest version of ValidateMerge().")
  }

  if (!is.logical(Interactive) || length(Interactive) != 1) {
    stop("Interactive must be TRUE or FALSE.")
  }

  if (Interactive && !requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "The plotly package is required when Interactive = TRUE. ",
      "Install it with install.packages('plotly') or set Interactive = FALSE."
    )
  }

  # Set palette

  status_colors <- c(
  "PASS" = "#2E7D32",
  "WARNING" = "#F9A825",
  "FAIL" = "#C62828"
)

coverage_colors <- c(
  "Matching" = "#2E7D32",
  "Left Only" = "#F9A825",
  "Right Only" = "#EF6C00"
)

join_role_colors <- c(
  "Specified Key" = "#1565C0",
  "Overlap Not In Keys" = "#F9A825"
)

  # Prepare checks data

  checks_df <- MergeObj$Checks %>%
    dplyr::mutate(
      Status = factor(
        Status,
        levels = c(
          "PASS",
          "WARNING",
          "FAIL"
        )
      ),
      Check = factor(
        Check,
        levels = rev(Check)
      ),
      HoverText = paste0(
        "Check: ",
        as.character(Check),
        "<br>Status: ",
        Status,
        "<br>Count: ",
        Count,
        "<br>Details: ",
        Details
      )
    )

  # Prepare coverage data

  coverage_df <- tibble::tibble(
    Category = c(
      "Matching",
      "Left Only",
      "Right Only"
    ),
    Count = c(
      nrow(MergeObj$IDCoverage$Matching),
      nrow(MergeObj$IDCoverage$LeftOnly),
      nrow(MergeObj$IDCoverage$RightOnly)
    )
  ) %>%
    dplyr::mutate(
      Category = factor(
        Category,
        levels = c(
          "Matching",
          "Left Only",
          "Right Only"
        )
      ),
      Percent = if (sum(Count) > 0) {
        round(
          100 * Count / sum(Count),
          1
        )
      } else {
        rep(
          NA_real_,
          dplyr::n()
        )
      },
      Label = paste0(
        Count,
        "\n",
        Percent,
        "%"
      ),
      HoverText = paste0(
        Category,
        "<br>Count: ",
        Count,
        "<br>Percent of displayed key groups: ",
        Percent,
        "%"
      )
    )

  # Prepare join audit data

  join_audit_df <- MergeObj$JoinAudit %>%
    dplyr::mutate(
      JoinRole = factor(
        JoinRole,
        levels = c(
          "Specified Key",
          "Overlap Not In Keys"
        )
      ),
      Variable = factor(
        Variable,
        levels = rev(Variable)
      ),
      HoverText = paste0(
        "Variable: ",
        Variable,
        "<br>Role: ",
        JoinRole,
        "<br>Is key: ",
        IsKey
      )
    )

  # Prepare duplicate variable data

  duplicate_df <- MergeObj$DuplicateVariables

  if (nrow(duplicate_df) > 0) {
    duplicate_df <- duplicate_df %>%
      dplyr::mutate(
        HoverText = paste0(
          "Variable: ",
          Variable,
          "<br>Agreement: ",
          Agreement,
          "%",
          "<br>Conflicts: ",
          Conflicts,
          "<br>Missingness conflicts: ",
          MissingnessConflicts,
          "<br>Both missing: ",
          BothMissing,
          "<br>Left class: ",
          LeftClass,
          "<br>Right class: ",
          RightClass
        )
      )
  }

  # Build checks plot

  checks_plot <- ggplot2::ggplot(
    checks_df,
    ggplot2::aes(
      x = Check,
      y = 1,
      fill = Status,
      text = HoverText
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Merge validation checks",
      x = NULL,
      y = NULL,
      fill = "Status"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    ) + ggplot2::scale_fill_manual(
  values = status_colors,
  drop = FALSE
)

  # Build coverage plot

 coverage_plot <- ggplot2::ggplot(
  coverage_df,
  ggplot2::aes(
    x = Category,
    y = Count,
    fill = Category,
    text = HoverText
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = Label
    ),
    vjust = -0.25,
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(
    values = coverage_colors,
    drop = FALSE
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Merge coverage",
    subtitle = "Key combinations by source overlap",
    x = NULL,
    y = "Key combinations",
    fill = "Coverage"
  )

  # Build join audit plot

  join_audit_plot <- ggplot2::ggplot(
    join_audit_df,
    ggplot2::aes(
      x = Variable,
      y = 1,
      fill = JoinRole,
      text = HoverText
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Join audit",
      x = NULL,
      y = NULL,
      fill = "Join role"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )  +ggplot2::scale_fill_manual(
  values = join_role_colors,
  drop = FALSE
)

  # Build agreement plot

  if (nrow(duplicate_df) > 0) {

    agreement_plot <- ggplot2::ggplot(
      duplicate_df,
      ggplot2::aes(
        x = stats::reorder(
          Variable,
          Agreement
        ),
        y = Agreement,
        text = HoverText
      )
    ) +
      ggplot2::geom_col(fill = "#1565C0") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Duplicate-variable agreement",
        x = NULL,
        y = "Agreement (%)"
      ) +
      ggplot2::ylim(
        0,
        100
      )

  } else {

    agreement_plot <- ggplot2::ggplot(
      tibble::tibble(
        x = 1,
        y = 1,
        Label = "No unresolved duplicate variables detected"
      ),
      ggplot2::aes(
        x = x,
        y = y,
        label = Label
      )
    ) +
      ggplot2::geom_text(
        size = 4
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = "Duplicate-variable agreement"
      )

  }

  # Build conflict plot

  if (nrow(duplicate_df) > 0) {

    conflict_plot <- ggplot2::ggplot(
      duplicate_df,
      ggplot2::aes(
        x = stats::reorder(
          Variable,
          Conflicts
        ),
        y = Conflicts,
        text = HoverText
      )
    ) +
      ggplot2::geom_col(fill = "#C62828") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Duplicate-variable conflicts",
        x = NULL,
        y = "Conflicts"
      )

  } else {

    conflict_plot <- ggplot2::ggplot(
      tibble::tibble(
        x = 1,
        y = 1,
        Label = "No variable conflicts detected"
      ),
      ggplot2::aes(
        x = x,
        y = y,
        label = Label
      )
    ) +
      ggplot2::geom_text(
        size = 4
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = "Duplicate-variable conflicts"
      )

  }

  # Convert to interactive plots when requested

  if (Interactive) {
    checks_plot <- plotly::ggplotly(
      checks_plot,
      tooltip = "text"
    )

    coverage_plot <- plotly::ggplotly(
      coverage_plot,
      tooltip = "text"
    )

    join_audit_plot <- plotly::ggplotly(
      join_audit_plot,
      tooltip = "text"
    )

    agreement_plot <- plotly::ggplotly(
      agreement_plot,
      tooltip = "text"
    )

    conflict_plot <- plotly::ggplotly(
      conflict_plot,
      tooltip = "text"
    )
  }

  # Build output list

  plots <- list(
    Checks = checks_plot,
    Coverage = coverage_plot,
    JoinAudit = join_audit_plot,
    Agreement = agreement_plot,
    Conflicts = conflict_plot
  )

  # Return result

  if (Plot == "All") {
    return(plots)
  }

  return(plots[[Plot]])
}
