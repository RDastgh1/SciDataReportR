#' Explore dataset comparison results interactively
#'
#' Create an interactive HTML dashboard from a `CompareDatasets()` result object.
#' This function is designed for data review and quality-control workflows. It
#' displays high-level summary cards, a traffic-light checks table, and an
#' expandable variable-change explorer that shows side-by-side old and new
#' values for modified cells.
#'
#' This function is intended for interactive review rather than publication
#' tables. It returns an HTML object that can be rendered in the RStudio Viewer,
#' Quarto, R Markdown, or Shiny.
#'
#' @param CompareObj A list returned by `CompareDatasets()`.
#' @param Title Character title shown at the top of the dashboard. Default is
#'   `"Dataset comparison explorer"`.
#' @param TopN Integer number of example variables or records to show in
#'   previews and expanded sections. Default is `10`.
#'
#' @return An `htmltools::tagList()` object containing an interactive dashboard.
#'
#' @export
ExploreDatasetComparison <- function(
    CompareObj,
    Title = "Dataset comparison explorer",
    TopN = 10
) {

  # Validate inputs

  if (!is.list(CompareObj)) {
    stop("CompareObj must be a list returned by CompareDatasets().")
  }

  required_elements <- c(
    "Summary",
    "Checks",
    "AddedRecords",
    "RemovedRecords",
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

  if (!is.character(Title) || length(Title) != 1) {
    stop("Title must be a single character value.")
  }

  if (!is.numeric(TopN) || length(TopN) != 1 || TopN < 1) {
    stop("TopN must be a single positive number.")
  }

  TopN <- as.integer(TopN)

  if (!requireNamespace("reactable", quietly = TRUE)) {
    stop(
      "The reactable package is required. ",
      "Install it with install.packages('reactable')."
    )
  }

  if (!requireNamespace("htmltools", quietly = TRUE)) {
    stop(
      "The htmltools package is required. ",
      "Install it with install.packages('htmltools')."
    )
  }

  # Prepare shared values

  summary_row <- CompareObj$Summary[1, ]

  status_colors <- c(
    "PASS" = "#2E7D32",
    "WARNING" = "#F9A825",
    "FAIL" = "#C62828"
  )

  status_backgrounds <- c(
    "PASS" = "#E8F5E9",
    "WARNING" = "#FFF8E1",
    "FAIL" = "#FFEBEE"
  )

  status_icons <- c(
    "PASS" = "\u25cf",
    "WARNING" = "\u25cf",
    "FAIL" = "\u25cf"
  )

  known_modified_cols <- c(
    "Variable",
    "OldVariable",
    "NewVariable",
    "OldValue",
    "NewValue",
    "OldClass",
    "NewClass",
    "ChangeType"
  )

  key_columns <- setdiff(
    names(CompareObj$ModifiedValues),
    known_modified_cols
  )

  # Build issue previews

  added_records_preview <- if (nrow(CompareObj$AddedRecords) > 0) {
    CompareObj$AddedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  removed_records_preview <- if (nrow(CompareObj$RemovedRecords) > 0) {
    CompareObj$RemovedRecords %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  added_variables_preview <- if (nrow(CompareObj$AddedVariables) > 0) {
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

  removed_variables_preview <- if (nrow(CompareObj$RemovedVariables) > 0) {
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

  duplicate_key_preview <- if (
    nrow(CompareObj$DuplicateKeys$Old) > 0 ||
      nrow(CompareObj$DuplicateKeys$New) > 0
  ) {

    old_duplicate_preview <- if (nrow(CompareObj$DuplicateKeys$Old) > 0) {
      old_tmp <- CompareObj$DuplicateKeys$Old %>%
        dplyr::slice_head(n = TopN)

      old_tmp <- old_tmp[
        ,
        seq_len(min(ncol(old_tmp), 5)),
        drop = FALSE
      ]

      old_tmp %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    new_duplicate_preview <- if (nrow(CompareObj$DuplicateKeys$New) > 0) {
      new_tmp <- CompareObj$DuplicateKeys$New %>%
        dplyr::slice_head(n = TopN)

      new_tmp <- new_tmp[
        ,
        seq_len(min(ncol(new_tmp), 5)),
        drop = FALSE
      ]

      new_tmp %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    paste0(
      "<strong>Old duplicate examples:</strong><br>",
      old_duplicate_preview,
      "<br><br><strong>New duplicate examples:</strong><br>",
      new_duplicate_preview
    )
  } else {
    "None"
  }

  name_repair_preview <- if (nrow(CompareObj$NameRepairAudit) > 0) {
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

  skipped_variable_preview <- if (nrow(CompareObj$ComparisonVariableMap) > 0) {
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

  class_change_preview <- if (nrow(CompareObj$ClassAudit) > 0) {
    class_tmp <- CompareObj$ClassAudit %>%
      dplyr::filter(
        Compared,
        ClassChanged
      )

    if (nrow(class_tmp) > 0) {
      class_tmp %>%
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

  modified_values_preview <- if (nrow(CompareObj$VariableChangeSummary) > 0) {
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
    nrow(CompareObj$ComparisonKeys$ExcludedDuplicateKeys) > 0
  ) {
    CompareObj$ComparisonKeys$ExcludedDuplicateKeys %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  key_type_preview <- if ("KeyTypes" %in% names(CompareObj) && nrow(CompareObj$KeyTypes) > 0) {
    key_tmp <- CompareObj$KeyTypes %>%
      dplyr::filter(
        OldType != NewType
      )

    if (nrow(key_tmp) > 0) {
      key_tmp %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::mutate(
          Text = paste0(
            Key,
            " (Old: ",
            OldType,
            "; New: ",
            NewType,
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

  # Prepare summary cards

  summary_cards <- tibble::tibble(
    Metric = c(
      "Added records",
      "Removed records",
      "Added variables",
      "Removed variables",
      "Modified values",
      "Suspicious changes",
      "Duplicate keys",
      "Skipped variables"
    ),
    Value = c(
      summary_row$AddedRecords,
      summary_row$RemovedRecords,
      summary_row$AddedVariables,
      summary_row$RemovedVariables,
      summary_row$ModifiedValues,
      summary_row$SuspiciousChanges,
      summary_row$DuplicateKeyGroups_Old + summary_row$DuplicateKeyGroups_New,
      summary_row$VariablesSkippedFromCellComparison
    ),
    Preview = c(
      added_records_preview,
      removed_records_preview,
      added_variables_preview,
      removed_variables_preview,
      modified_values_preview,
      suspicious_preview,
      duplicate_key_preview,
      skipped_variable_preview
    )
  ) %>%
    dplyr::mutate(
      Status = dplyr::case_when(
        Metric == "Duplicate keys" & Value > 0 ~ "FAIL",
        Value > 0 ~ "WARNING",
        TRUE ~ "PASS"
      ),
      Color = dplyr::case_when(
        Status == "PASS" ~ status_colors[["PASS"]],
        Status == "WARNING" ~ status_colors[["WARNING"]],
        Status == "FAIL" ~ status_colors[["FAIL"]],
        TRUE ~ "#546E7A"
      ),
      Background = dplyr::case_when(
        Status == "PASS" ~ status_backgrounds[["PASS"]],
        Status == "WARNING" ~ status_backgrounds[["WARNING"]],
        Status == "FAIL" ~ status_backgrounds[["FAIL"]],
        TRUE ~ "#ECEFF1"
      )
    )

  summary_card_tags <- purrr::pmap(
    summary_cards,
    function(Metric, Value, Preview, Status, Color, Background) {
      htmltools::tags$div(
        class = "sdr-card",
        style = paste0(
          "border-left: 6px solid ",
          Color,
          "; background: ",
          Background,
          ";"
        ),
        htmltools::tags$div(
          class = "sdr-card-label",
          Metric
        ),
        htmltools::tags$div(
          class = "sdr-card-value",
          Value
        ),
        htmltools::tags$div(
          class = "sdr-card-status",
          paste0(Status)
        )
      )
    }
  )

  # Prepare checks table

  checks_table_df <- CompareObj$Checks %>%
    dplyr::mutate(
      DisplayCheck = dplyr::case_when(
        Check == "Name Repair Differences" ~ "Name Repair",
        Check == "Variables Skipped From Cell Comparison" ~ "Skipped Variables",
        Check == "Duplicate Keys Excluded From Cell Comparison" ~ "Excluded Keys",
        TRUE ~ Check
      ),
      Examples = dplyr::case_when(
        Check == "Key Types" ~ key_type_preview,
        Check == "Duplicate Keys" ~ duplicate_key_preview,
        Check == "Records Added" ~ added_records_preview,
        Check == "Records Removed" ~ removed_records_preview,
        Check == "Variables Added" ~ added_variables_preview,
        Check == "Variables Removed" ~ removed_variables_preview,
        Check == "Name Repair Differences" ~ name_repair_preview,
        Check == "Variables Skipped From Cell Comparison" ~ skipped_variable_preview,
        Check == "Class Changes" ~ class_change_preview,
        Check == "Values Modified" ~ modified_values_preview,
        Check == "Suspicious Changes" ~ suspicious_preview,
        Check == "Duplicate Keys Excluded From Cell Comparison" ~ excluded_key_preview,
        TRUE ~ "None"
      )
    ) %>%
    dplyr::select(
      Status,
      Check = DisplayCheck,
      Count,
      Details,
      Examples
    )

  checks_table <- reactable::reactable(
    checks_table_df,
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    striped = TRUE,
    compact = TRUE,
    defaultPageSize = 12,
    columns = list(
      Status = reactable::colDef(
        width = 130,
        cell = function(value) {
          color <- status_colors[[value]]
          background <- status_backgrounds[[value]]
          icon <- status_icons[[value]]

          htmltools::tags$span(
            style = paste0(
              "display:inline-block;",
              "padding:4px 9px;",
              "border-radius:999px;",
              "font-weight:700;",
              "color:",
              color,
              ";background:",
              background,
              ";"
            ),
            paste(icon, value)
          )
        }
      ),
      Check = reactable::colDef(
        minWidth = 190
      ),
      Count = reactable::colDef(
        align = "center",
        width = 90
      ),
      Details = reactable::colDef(
        minWidth = 350
      ),
      Examples = reactable::colDef(
        show = FALSE,
        html = TRUE
      )
    ),
    details = function(index) {
      row <- checks_table_df[index, ]

      htmltools::tags$div(
        class = "sdr-detail-box",
        htmltools::tags$div(
          class = "sdr-detail-title",
          paste0(row$Check, " details")
        ),
        htmltools::tags$p(
          htmltools::HTML(htmltools::htmlEscape(row$Details))
        ),
        htmltools::tags$div(
          class = "sdr-detail-subtitle",
          "Examples"
        ),
        htmltools::tags$div(
          class = "sdr-detail-examples",
          htmltools::HTML(row$Examples)
        )
      )
    }
  )

  # Prepare variable changes table

  variable_changes_df <- CompareObj$VariableChangeSummary

  if (nrow(variable_changes_df) > 0) {
    variable_changes_df <- variable_changes_df %>%
      dplyr::mutate(
        PercentChangedLabel = paste0(PercentChanged, "%"),
        ClassChangeLabel = dplyr::case_when(
          ClassChanged ~ "Class changed",
          TRUE ~ "No class change"
        )
      )
  } else {
    variable_changes_df <- tibble::tibble(
      Variable = character(),
      Changes = integer(),
      PercentChanged = numeric(),
      PercentChangedLabel = character(),
      ComparedRecords = integer(),
      OldVariable = character(),
      NewVariable = character(),
      OldClass = character(),
      NewClass = character(),
      ClassChanged = logical(),
      ClassChangeLabel = character()
    )
  }

  variable_changes_table <- reactable::reactable(
    variable_changes_df,
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    striped = TRUE,
    compact = TRUE,
    defaultPageSize = 10,
    defaultSorted = "Changes",
    defaultSortOrder = "desc",
    columns = list(
      Variable = reactable::colDef(
        minWidth = 220
      ),
      Changes = reactable::colDef(
        align = "center",
        width = 100,
        style = function(value) {
          if (is.na(value) || value == 0) {
            list()
          } else {
            list(
              fontWeight = "700",
              color = "#1565C0"
            )
          }
        }
      ),
      PercentChanged = reactable::colDef(
        name = "% Changed",
        align = "center",
        width = 120,
        cell = function(value) {
          paste0(value, "%")
        },
        style = function(value) {
          if (is.na(value)) {
            list()
          } else if (value >= 10) {
            list(
              fontWeight = "700",
              color = "#C62828",
              background = "#FFEBEE"
            )
          } else if (value > 0) {
            list(
              fontWeight = "700",
              color = "#1565C0"
            )
          } else {
            list()
          }
        }
      ),
      PercentChangedLabel = reactable::colDef(
        show = FALSE
      ),
      ComparedRecords = reactable::colDef(
        name = "Compared records",
        align = "center",
        width = 140
      ),
      OldVariable = reactable::colDef(
        name = "Old variable",
        minWidth = 180
      ),
      NewVariable = reactable::colDef(
        name = "New variable",
        minWidth = 180
      ),
      OldClass = reactable::colDef(
        name = "Old class",
        width = 130
      ),
      NewClass = reactable::colDef(
        name = "New class",
        width = 130
      ),
      ClassChanged = reactable::colDef(
        show = FALSE
      ),
      ClassChangeLabel = reactable::colDef(
        name = "Class status",
        width = 140,
        cell = function(value) {
          if (is.na(value)) {
            return("")
          }

          if (value == "Class changed") {
            htmltools::tags$span(
              style = paste0(
                "display:inline-block;",
                "padding:4px 9px;",
                "border-radius:999px;",
                "font-weight:700;",
                "color:#C62828;",
                "background:#FFEBEE;"
              ),
              value
            )
          } else {
            htmltools::tags$span(
              style = paste0(
                "display:inline-block;",
                "padding:4px 9px;",
                "border-radius:999px;",
                "color:#2E7D32;",
                "background:#E8F5E9;"
              ),
              value
            )
          }
        }
      )
    ),
    details = function(index) {
      if (nrow(variable_changes_df) == 0) {
        return(NULL)
      }

      selected_variable <- as.character(variable_changes_df$Variable[index])

      detail_df <- CompareObj$ModifiedValues %>%
        dplyr::filter(
          Variable == selected_variable
        ) %>%
        dplyr::slice_head(n = TopN)

      if (nrow(detail_df) == 0) {
        return(
          htmltools::tags$div(
            class = "sdr-detail-box",
            "No modified cell examples available for this variable."
          )
        )
      }

      detail_df <- detail_df %>%
        dplyr::mutate(
          OldValueEscaped = htmltools::htmlEscape(OldValue),
          NewValueEscaped = htmltools::htmlEscape(NewValue),
          Change = paste0(
            "<span class='sdr-old-value'>",
            OldValueEscaped,
            "</span>",
            " &rarr; ",
            "<span class='sdr-new-value'>",
            NewValueEscaped,
            "</span>"
          )
        )

      detail_display <- detail_df %>%
        dplyr::select(
          dplyr::all_of(key_columns),
          ChangeType,
          Change,
          OldClass,
          NewClass
        )

      reactable::reactable(
        detail_display,
        searchable = TRUE,
        filterable = TRUE,
        highlight = TRUE,
        bordered = TRUE,
        striped = TRUE,
        compact = TRUE,
        defaultPageSize = min(TopN, 10),
        columns = list(
          Change = reactable::colDef(
            html = TRUE,
            minWidth = 260
          ),
          ChangeType = reactable::colDef(
            name = "Change type",
            minWidth = 190
          ),
          OldClass = reactable::colDef(
            name = "Old class",
            width = 120
          ),
          NewClass = reactable::colDef(
            name = "New class",
            width = 120
          )
        )
      )
    }
  )

  # Assemble dashboard

  dashboard <- htmltools::tagList(
    htmltools::tags$style(
      htmltools::HTML(
        "
        .sdr-dashboard {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          color: #263238;
        }

        .sdr-title {
          font-size: 28px;
          font-weight: 750;
          margin: 0 0 6px 0;
        }

        .sdr-subtitle {
          color: #607D8B;
          margin-bottom: 22px;
          font-size: 14px;
        }

        .sdr-card-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
          gap: 12px;
          margin-bottom: 24px;
        }

        .sdr-card {
          border-radius: 14px;
          padding: 14px 16px;
          box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        }

        .sdr-card-label {
          font-size: 13px;
          color: #546E7A;
          margin-bottom: 6px;
        }

        .sdr-card-value {
          font-size: 30px;
          font-weight: 800;
          line-height: 1;
        }

        .sdr-card-status {
          margin-top: 8px;
          font-size: 12px;
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: 0.04em;
          color: #455A64;
        }

        .sdr-section {
          margin-top: 28px;
          margin-bottom: 10px;
        }

        .sdr-section-title {
          font-size: 21px;
          font-weight: 750;
          margin-bottom: 4px;
        }

        .sdr-section-subtitle {
          font-size: 13px;
          color: #607D8B;
          margin-bottom: 12px;
        }

        .sdr-detail-box {
          background: #FAFAFA;
          border-left: 4px solid #1565C0;
          padding: 12px 14px;
          margin: 8px;
          border-radius: 8px;
        }

        .sdr-detail-title {
          font-weight: 750;
          margin-bottom: 8px;
        }

        .sdr-detail-subtitle {
          font-weight: 700;
          margin-top: 10px;
          margin-bottom: 4px;
          color: #455A64;
        }

        .sdr-detail-examples {
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
          font-size: 12px;
          line-height: 1.45;
          color: #37474F;
        }

        .sdr-old-value {
          color: #C62828;
          background: #FFEBEE;
          padding: 2px 5px;
          border-radius: 5px;
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        }

        .sdr-new-value {
          color: #2E7D32;
          background: #E8F5E9;
          padding: 2px 5px;
          border-radius: 5px;
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        }
        "
      )
    ),
    htmltools::tags$div(
      class = "sdr-dashboard",
      htmltools::tags$div(
        class = "sdr-title",
        Title
      ),
      htmltools::tags$div(
        class = "sdr-subtitle",
        "Interactive review of dataset version changes from CompareDatasets()."
      ),
      htmltools::tags$div(
        class = "sdr-card-grid",
        summary_card_tags
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Validation checks"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Search, filter, sort, and expand rows to inspect examples."
        ),
        checks_table
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Variable changes"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Expand a variable to review old and new values side by side."
        ),
        variable_changes_table
      )
    )
  )

  return(dashboard)
}
