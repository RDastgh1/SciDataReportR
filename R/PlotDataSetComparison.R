#' Plot dataset comparison diagnostics
#'
#' Create diagnostic plots from a `CompareDatasets()` result object. This
#' function visualizes dataset-version changes, including check status, summary
#' metrics, structure changes, and variable-level value changes.
#'
#' Use this function after running `CompareDatasets()` to quickly inspect what
#' changed between two versions of a dataset.
#'
#' @param CompareObj A list returned by `CompareDatasets()`.
#' @param Plot Character value specifying which plot to return. Options are
#'   `"All"`, `"Checks"`, `"SummaryMetrics"`, `"StructureChanges"`,
#'   `"VariableChanges"`, and `"TopChangedVariables"`. Default is `"All"`.
#' @param Interactive Logical; if `TRUE`, plots are converted to interactive
#'   `plotly` objects using `plotly::ggplotly()`. Default is `TRUE`.
#' @param TopN Integer number of variables or records to preview in plots and
#'   hover text. Default is `10`.
#'
#' @return If `Plot = "All"`, a named list of plots. Otherwise, a single plot
#' object. Plot objects are either `ggplot` objects or `plotly` htmlwidgets,
#' depending on `Interactive`.
#'
#' @export
PlotDatasetComparison <- function(
    CompareObj,
    Plot = c(
      "All",
      "Checks",
      "SummaryMetrics",
      "StructureChanges",
      "VariableChanges",
      "TopChangedVariables"
    ),
    Interactive = TRUE,
    TopN = 10
) {

  # Validate inputs

  Plot <- match.arg(Plot)

  if (!is.list(CompareObj)) {
    stop("CompareObj must be a list returned by CompareDatasets().")
  }

  required_elements <- c(
    "Checks",
    "Summary",
    "AddedRecords",
    "RemovedRecords",
    "StructureChanges",
    "AddedVariables",
    "RemovedVariables",
    "DuplicateKeys",
    "ComparisonKeys",
    "NameRepairAudit",
    "ComparisonVariableMap",
    "ClassAudit",
    "ModifiedValues",
    "VariableChangeSummary",
    "TopChangedVariables",
    "SuspiciousChanges"
  )

  missing_elements <- setdiff(
    required_elements,
    names(CompareObj)
  )

  if (length(missing_elements) > 0) {
    stop(
      "CompareObj is missing required element(s): ",
      paste(missing_elements, collapse = ", "),
      ". Please provide an object returned by CompareDatasets()."
    )
  }

  if (!is.logical(Interactive) || length(Interactive) != 1) {
    stop("Interactive must be TRUE or FALSE.")
  }

  if (!is.numeric(TopN) || length(TopN) != 1 || TopN < 1) {
    stop("TopN must be a single positive number.")
  }

  TopN <- as.integer(TopN)

  if (Interactive && !requireNamespace("plotly", quietly = TRUE)) {
    stop(
      "The plotly package is required when Interactive = TRUE. ",
      "Install it with install.packages('plotly') or set Interactive = FALSE."
    )
  }

  # Prepare plotting palettes

  status_colors <- c(
    "PASS" = "#2E7D32",
    "WARNING" = "#F9A825",
    "FAIL" = "#C62828"
  )

  summary_colors <- c(
    "Records" = "#1565C0",
    "Variables" = "#6A1B9A",
    "Values" = "#00897B",
    "Warnings" = "#F9A825"
  )

  structure_colors <- c(
    "Added Records" = "#2E7D32",
    "Removed Records" = "#C62828",
    "Added Variables" = "#1565C0",
    "Removed Variables" = "#EF6C00"
  )

  # Prepare hover previews for checks

  key_type_preview <- if ("KeyTypes" %in% names(CompareObj) && nrow(CompareObj$KeyTypes) > 0) {
    tmp <- CompareObj$KeyTypes %>%
      dplyr::filter(
        OldType != NewType
      )

    if (nrow(tmp) > 0) {
      paste(
        tmp %>%
          dplyr::slice_head(n = TopN) %>%
          dplyr::mutate(
            Text = paste0(Key, " (Old: ", OldType, "; New: ", NewType, ")")
          ) %>%
          dplyr::pull(Text),
        collapse = "<br>"
      )
    } else {
      "None"
    }
  } else {
    "None"
  }

  duplicate_key_preview <- if (
    "DuplicateKeys" %in% names(CompareObj) &&
      (nrow(CompareObj$DuplicateKeys$Old) > 0 || nrow(CompareObj$DuplicateKeys$New) > 0)
  ) {
    old_preview <- if (nrow(CompareObj$DuplicateKeys$Old) > 0) {
      CompareObj$DuplicateKeys$Old %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::select(where(~ !is.list(.x))) %>%
        dplyr::select(1:min(ncol(.), 5)) %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    new_preview <- if (nrow(CompareObj$DuplicateKeys$New) > 0) {
      CompareObj$DuplicateKeys$New %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::select(where(~ !is.list(.x))) %>%
        dplyr::select(1:min(ncol(.), 5)) %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    paste0(
      "Old duplicate examples:<br>",
      old_preview,
      "<br><br>New duplicate examples:<br>",
      new_preview
    )
  } else {
    "None"
  }

  added_record_preview <- if (nrow(CompareObj$AddedRecords) > 0) {
    CompareObj$AddedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  removed_record_preview <- if (nrow(CompareObj$RemovedRecords) > 0) {
    CompareObj$RemovedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  added_variable_preview <- if (nrow(CompareObj$AddedVariables) > 0) {
    CompareObj$AddedVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::pull(Variable) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  removed_variable_preview <- if (nrow(CompareObj$RemovedVariables) > 0) {
    CompareObj$RemovedVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::pull(Variable) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  name_repair_preview <- if (nrow(CompareObj$NameRepairAudit) > 0) {
    CompareObj$NameRepairAudit %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (Old: ",
          OldVariableCount,
          "; New: ",
          NewVariableCount,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  skipped_variable_preview <- if (
    "ComparisonVariableMap" %in% names(CompareObj) &&
      nrow(CompareObj$ComparisonVariableMap) > 0
  ) {
    tmp <- CompareObj$ComparisonVariableMap %>%
      dplyr::filter(!Compared)

    if (nrow(tmp) > 0) {
      tmp %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::mutate(
          Text = paste0(
            Variable,
            " [",
            SkipReason,
            "]"
          )
        ) %>%
        dplyr::pull(Text) %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }
  } else {
    "None"
  }

  class_change_preview <- if (nrow(CompareObj$ClassAudit) > 0) {
    tmp <- CompareObj$ClassAudit %>%
      dplyr::filter(
        Compared,
        ClassChanged
      )

    if (nrow(tmp) > 0) {
      tmp %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::mutate(
          Text = paste0(
            Variable,
            " (Old: ",
            OldClass,
            "; New: ",
            NewClass,
            ")"
          )
        ) %>%
        dplyr::pull(Text) %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }
  } else {
    "None"
  }

  modified_value_preview <- if (nrow(CompareObj$VariableChangeSummary) > 0) {
    CompareObj$VariableChangeSummary %>%
      dplyr::arrange(
        dplyr::desc(Changes),
        Variable
      ) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          Changes,
          "; ",
          PercentChanged,
          "%)"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  suspicious_preview <- if (nrow(CompareObj$SuspiciousChanges) > 0) {
    CompareObj$SuspiciousChanges %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " [",
          Reason,
          "]"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  excluded_key_preview <- if (
    "ComparisonKeys" %in% names(CompareObj) &&
      nrow(CompareObj$ComparisonKeys$ExcludedDuplicateKeys) > 0
  ) {
    CompareObj$ComparisonKeys$ExcludedDuplicateKeys %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  # Prepare checks data

  checks_df <- CompareObj$Checks %>%
    dplyr::mutate(
      DisplayCheck = dplyr::case_when(
        Check == "Name Repair Differences" ~ "Name Repair",
        Check == "Variables Skipped From Cell Comparison" ~ "Skipped Variables",
        Check == "Duplicate Keys Excluded From Cell Comparison" ~ "Excluded Keys",
        TRUE ~ Check
      ),
      IssuePreview = dplyr::case_when(
        Check == "Key Types" ~ key_type_preview,
        Check == "Duplicate Keys" ~ duplicate_key_preview,
        Check == "Records Added" ~ added_record_preview,
        Check == "Records Removed" ~ removed_record_preview,
        Check == "Variables Added" ~ added_variable_preview,
        Check == "Variables Removed" ~ removed_variable_preview,
        Check == "Name Repair Differences" ~ name_repair_preview,
        Check == "Variables Skipped From Cell Comparison" ~ skipped_variable_preview,
        Check == "Class Changes" ~ class_change_preview,
        Check == "Values Modified" ~ modified_value_preview,
        Check == "Suspicious Changes" ~ suspicious_preview,
        Check == "Duplicate Keys Excluded From Cell Comparison" ~ excluded_key_preview,
        TRUE ~ "None"
      ),
      Status = factor(
        Status,
        levels = c(
          "PASS",
          "WARNING",
          "FAIL"
        )
      ),
      DisplayCheck = factor(
        DisplayCheck,
        levels = rev(DisplayCheck)
      ),
      HoverText = paste0(
        "Check: ",
        Check,
        "<br>Status: ",
        Status,
        "<br>Count: ",
        Count,
        "<br>Details: ",
        Details,
        "<br><br>Examples:<br>",
        IssuePreview
      )
    )

  # Prepare summary metrics data

  summary_row <- CompareObj$Summary[1, ]

  added_records_preview_summary <- if (nrow(CompareObj$AddedRecords) > 0) {
    CompareObj$AddedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  removed_records_preview_summary <- if (nrow(CompareObj$RemovedRecords) > 0) {
    CompareObj$RemovedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  added_variables_preview_summary <- if (nrow(CompareObj$AddedVariables) > 0) {
    CompareObj$AddedVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = dplyr::if_else(
          "NewVariables" %in% names(.),
          paste0(Variable, " [", NewVariables, "]"),
          Variable
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  removed_variables_preview_summary <- if (nrow(CompareObj$RemovedVariables) > 0) {
    CompareObj$RemovedVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = dplyr::if_else(
          "OldVariables" %in% names(.),
          paste0(Variable, " [", OldVariables, "]"),
          Variable
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  modified_values_preview_summary <- if (nrow(CompareObj$VariableChangeSummary) > 0) {
    CompareObj$VariableChangeSummary %>%
      dplyr::arrange(
        dplyr::desc(Changes),
        Variable
      ) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          Changes,
          " changes; ",
          PercentChanged,
          "%)"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  suspicious_changes_preview_summary <- if (nrow(CompareObj$SuspiciousChanges) > 0) {
    CompareObj$SuspiciousChanges %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " [",
          Reason,
          "]"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  name_repair_preview_summary <- if (nrow(CompareObj$NameRepairAudit) > 0) {
    CompareObj$NameRepairAudit %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (Old raw columns: ",
          OldVariableCount,
          "; New raw columns: ",
          NewVariableCount,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  skipped_variables_preview_summary <- if (
    "ComparisonVariableMap" %in% names(CompareObj) &&
      nrow(CompareObj$ComparisonVariableMap) > 0
  ) {
    skipped_tmp <- CompareObj$ComparisonVariableMap %>%
      dplyr::filter(!Compared)

    if (nrow(skipped_tmp) > 0) {
      skipped_tmp %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::mutate(
          Text = paste0(
            Variable,
            " [",
            SkipReason,
            "]"
          )
        ) %>%
        dplyr::pull(Text) %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }
  } else {
    "None"
  }

  summary_metrics_df <- tibble::tibble(
    Metric = c(
      "Added Records",
      "Removed Records",
      "Added Variables",
      "Removed Variables",
      "Modified Values",
      "Suspicious Changes",
      "Name Repair",
      "Skipped Variables"
    ),
    Value = c(
      summary_row$AddedRecords,
      summary_row$RemovedRecords,
      summary_row$AddedVariables,
      summary_row$RemovedVariables,
      summary_row$ModifiedValues,
      summary_row$SuspiciousChanges,
      summary_row$NameRepairDifferences,
      summary_row$VariablesSkippedFromCellComparison
    ),
    Group = c(
      "Records",
      "Records",
      "Variables",
      "Variables",
      "Values",
      "Warnings",
      "Warnings",
      "Warnings"
    ),
    Examples = c(
      added_records_preview_summary,
      removed_records_preview_summary,
      added_variables_preview_summary,
      removed_variables_preview_summary,
      modified_values_preview_summary,
      suspicious_changes_preview_summary,
      name_repair_preview_summary,
      skipped_variables_preview_summary
    )
  ) %>%
    dplyr::mutate(
      Metric = factor(
        Metric,
        levels = rev(Metric)
      ),
      HoverText = paste0(
        "Metric: ",
        Metric,
        "<br>Value: ",
        Value,
        "<br>Group: ",
        Group,
        "<br><br>Examples:<br>",
        Examples
      )
    )
  # Prepare structure changes data

  structure_df <- tibble::tibble(
    Category = c(
      "Added Records",
      "Removed Records",
      "Added Variables",
      "Removed Variables"
    ),
    Count = c(
      nrow(CompareObj$AddedRecords),
      nrow(CompareObj$RemovedRecords),
      nrow(CompareObj$AddedVariables),
      nrow(CompareObj$RemovedVariables)
    ),
    Preview = c(
      added_record_preview,
      removed_record_preview,
      added_variable_preview,
      removed_variable_preview
    )
  ) %>%
    dplyr::mutate(
      Category = factor(
        Category,
        levels = c(
          "Added Records",
          "Removed Records",
          "Added Variables",
          "Removed Variables"
        )
      ),
      HoverText = paste0(
        Category,
        "<br>Count: ",
        Count,
        "<br><br>Examples:<br>",
        Preview
      )
    )

  # Prepare variable-change data

  variable_change_df <- CompareObj$VariableChangeSummary

  if (nrow(variable_change_df) > 0) {
    variable_change_df <- variable_change_df %>%
      dplyr::arrange(
        dplyr::desc(Changes),
        Variable
      ) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Variable = factor(
          Variable,
          levels = rev(Variable)
        ),
        HoverText = paste0(
          "Variable: ",
          Variable,
          "<br>Changes: ",
          Changes,
          "<br>Percent changed: ",
          PercentChanged,
          "%",
          "<br>Compared records: ",
          ComparedRecords,
          "<br>Old variable: ",
          OldVariable,
          "<br>New variable: ",
          NewVariable,
          "<br>Old class: ",
          OldClass,
          "<br>New class: ",
          NewClass
        )
      )
  }

  top_changed_df <- CompareObj$TopChangedVariables

  if (nrow(top_changed_df) > 0) {
    top_changed_df <- top_changed_df %>%
      dplyr::arrange(
        dplyr::desc(PercentChanged),
        dplyr::desc(Changes),
        Variable
      ) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Variable = factor(
          Variable,
          levels = rev(Variable)
        ),
        HoverText = paste0(
          "Variable: ",
          Variable,
          "<br>Changes: ",
          Changes,
          "<br>Percent changed: ",
          PercentChanged,
          "%",
          "<br>Compared records: ",
          ComparedRecords,
          "<br>Old variable: ",
          OldVariable,
          "<br>New variable: ",
          NewVariable,
          "<br>Old class: ",
          OldClass,
          "<br>New class: ",
          NewClass
        )
      )
  }

  # Build checks plot

  checks_plot <- ggplot2::ggplot(
    checks_df,
    ggplot2::aes(
      x = DisplayCheck,
      y = 1,
      fill = Status,
      text = HoverText
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = status_colors
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Dataset comparison checks",
      subtitle = "PASS, WARNING, and FAIL status by comparison check",
      x = NULL,
      y = NULL,
      fill = "Status"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )

  # Build summary metrics plot

  summary_metrics_plot <- ggplot2::ggplot(
    summary_metrics_df,
    ggplot2::aes(
      x = Metric,
      y = Value,
      fill = Group,
      text = HoverText
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Value
      ),
      hjust = -0.15,
      size = 3.4
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = summary_colors
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Dataset comparison summary",
      subtitle = "High-level record, variable, value, and warning counts",
      x = NULL,
      y = "Count",
      fill = "Group"
    ) +
    ggplot2::expand_limits(
      y = max(summary_metrics_df$Value, na.rm = TRUE) * 1.12
    )

  # Build structure changes plot

  structure_plot <- ggplot2::ggplot(
    structure_df,
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
        label = Count
      ),
      vjust = -0.25,
      size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      values = structure_colors
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Record and structure changes",
      subtitle = "Records and variables added or removed",
      x = NULL,
      y = "Count",
      fill = NULL
    ) +
    ggplot2::expand_limits(
      y = max(structure_df$Count, na.rm = TRUE) * 1.12
    )

  # Build variable changes plot

  if (nrow(variable_change_df) > 0) {

    variable_changes_plot <- ggplot2::ggplot(
      variable_change_df,
      ggplot2::aes(
        x = Variable,
        y = Changes,
        text = HoverText
      )
    ) +
      ggplot2::geom_col(
        fill = "#1565C0"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = Changes
        ),
        hjust = -0.15,
        size = 3.4
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Variable-level changes",
        subtitle = paste0(
          "Top ",
          min(TopN, nrow(variable_change_df)),
          " variables by modified values"
        ),
        x = NULL,
        y = "Modified values"
      ) +
      ggplot2::expand_limits(
        y = max(variable_change_df$Changes, na.rm = TRUE) * 1.12
      )

  } else {

    variable_changes_plot <- ggplot2::ggplot(
      tibble::tibble(
        x = 1,
        y = 1,
        Label = "No modified values detected"
      ),
      ggplot2::aes(
        x = x,
        y = y,
        label = Label,
        text = Label
      )
    ) +
      ggplot2::geom_text(
        size = 4
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = "Variable-level changes"
      )

  }

  # Build top changed variables plot

  if (nrow(top_changed_df) > 0) {

    top_changed_plot <- ggplot2::ggplot(
      top_changed_df,
      ggplot2::aes(
        x = Variable,
        y = PercentChanged,
        text = HoverText
      )
    ) +
      ggplot2::geom_col(
        fill = "#6A1B9A"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = paste0(
            PercentChanged,
            "%"
          )
        ),
        hjust = -0.15,
        size = 3.4
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Top changed variables",
        subtitle = paste0(
          "Top ",
          min(TopN, nrow(top_changed_df)),
          " variables by percent changed"
        ),
        x = NULL,
        y = "Percent changed"
      ) +
      ggplot2::expand_limits(
        y = max(top_changed_df$PercentChanged, na.rm = TRUE) * 1.12
      )

  } else {

    top_changed_plot <- ggplot2::ggplot(
      tibble::tibble(
        x = 1,
        y = 1,
        Label = "No changed variables detected"
      ),
      ggplot2::aes(
        x = x,
        y = y,
        label = Label,
        text = Label
      )
    ) +
      ggplot2::geom_text(
        size = 4
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = "Top changed variables"
      )

  }

  # Convert to interactive plots when requested

  if (Interactive) {
    checks_plot <- plotly::ggplotly(
      checks_plot,
      tooltip = "text"
    ) %>%
      plotly::config(
        displayModeBar = FALSE
      )

    summary_metrics_plot <- plotly::ggplotly(
      summary_metrics_plot,
      tooltip = "text"
    ) %>%
      plotly::config(
        displayModeBar = FALSE
      )

    structure_plot <- plotly::ggplotly(
      structure_plot,
      tooltip = "text"
    ) %>%
      plotly::config(
        displayModeBar = FALSE
      )

    variable_changes_plot <- plotly::ggplotly(
      variable_changes_plot,
      tooltip = "text"
    ) %>%
      plotly::config(
        displayModeBar = FALSE
      )

    top_changed_plot <- plotly::ggplotly(
      top_changed_plot,
      tooltip = "text"
    ) %>%
      plotly::config(
        displayModeBar = FALSE
      )
  }

  # Build output list

  plots <- list(
    Checks = checks_plot,
    SummaryMetrics = summary_metrics_plot,
    StructureChanges = structure_plot,
    VariableChanges = variable_changes_plot,
    TopChangedVariables = top_changed_plot
  )

  # Return result

  if (Plot == "All") {
    return(plots)
  }

  return(plots[[Plot]])
}
