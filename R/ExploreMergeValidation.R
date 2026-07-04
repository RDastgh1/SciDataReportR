#' Explore merge validation results interactively
#'
#' Create an interactive HTML dashboard from a `ValidateMerge()` result object.
#' This function is designed for merge quality-control workflows. It displays a
#' traffic-light checks table (with rows/columns/unique-key context folded in
#' as informational rows), a coverage explorer, an expandable
#' duplicate-variable conflict explorer, and suggested actions.
#'
#' This function is intended for interactive review rather than publication
#' tables. It returns an HTML object that can be rendered in the RStudio Viewer,
#' Quarto, R Markdown, Shiny, or saved as HTML. If needed in an interactive
#' console, wrap the result with `htmltools::browsable()`.
#'
#' @param MergeObj A list returned by `ValidateMerge()`.
#' @param Title Character title shown at the top of the dashboard. Default is
#'   `"Merge validation explorer"`.
#' @param TopN Integer number of example variables or records to show in
#'   previews and expanded sections. Default is `10`.
#' @param TableHeight Height in pixels for scrollable reactable tables.
#'   Default is `350`.
#' @param Detail Either `"Compact"` (default) or `"Full"`. In `"Compact"`
#'   mode, the coverage explorer and conflicts explorer render as collapsed
#'   click-to-expand accordion sections labeled with their item counts. In
#'   `"Full"` mode, the same sections are expanded by default. In both modes,
#'   sections with nothing to show (no unmatched keys, no duplicated
#'   variables, no suggested actions) are omitted entirely.
#'
#' @return An `htmltools::tagList()` object containing an interactive dashboard.
#'
#' @export
ExploreMergeValidation <- function(
    MergeObj,
    Title = "Merge validation explorer",
    TopN = 10,
    TableHeight = 350,
    Detail = c("Compact", "Full")
) {

  # Validate inputs

  Detail <- match.arg(Detail)

  if (!is.list(MergeObj)) {
    stop("MergeObj must be a list returned by ValidateMerge().")
  }

  required_elements <- c(
    "Summary",
    "Checks",
    "ReadyForAnalysis",
    "IDCoverage",
    "DuplicateKeys",
    "OverlappingVariables",
    "PotentialMergeRisk",
    "JoinAudit",
    "OverlapAudit",
    "DuplicateVariables",
    "SuspiciousConflicts",
    "VariableConflicts"
  )

  missing_elements <- setdiff(
    required_elements,
    names(MergeObj)
  )

  if (length(missing_elements) > 0) {
    stop(
      "MergeObj is missing required element(s): ",
      paste(missing_elements, collapse = ", "),
      ". Please provide an object returned by ValidateMerge()."
    )
  }

  if (!is.character(Title) || length(Title) != 1) {
    stop("Title must be a single character value.")
  }

  if (!is.numeric(TopN) || length(TopN) != 1 || TopN < 1) {
    stop("TopN must be a single positive number.")
  }

  if (!is.numeric(TableHeight) || length(TableHeight) != 1 || TableHeight < 100) {
    stop("TableHeight must be a single numeric value of at least 100.")
  }

  TopN <- as.integer(TopN)
  TableHeight <- as.integer(TableHeight)

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

  expanded <- Detail == "Full"

  # Prepare shared values

  summary_row <- MergeObj$Summary[1, ]

  status_colors <- c(
    "PASS" = "#2E7D32",
    "WARNING" = "#F9A825",
    "FAIL" = "#C62828",
    "INFO" = "#1565C0"
  )

  status_backgrounds <- c(
    "PASS" = "#E8F5E9",
    "WARNING" = "#FFF8E1",
    "FAIL" = "#FFEBEE",
    "INFO" = "#E3F2FD"
  )

  status_icons <- c(
    "PASS" = "●",
    "WARNING" = "●",
    "FAIL" = "●",
    "INFO" = "●"
  )

  key_columns <- names(MergeObj$IDCoverage$Matching)

  known_conflict_cols <- c(
    "Variable",
    "LeftValue",
    "RightValue",
    "ConflictType"
  )

  conflict_key_columns <- setdiff(
    names(MergeObj$VariableConflicts),
    known_conflict_cols
  )

  # Prepare dataset context values (folded into the checks table)

  left_rows <- if ("LeftRows" %in% names(summary_row)) summary_row$LeftRows else NA_integer_
  right_rows <- if ("RightRows" %in% names(summary_row)) summary_row$RightRows else NA_integer_
  merged_rows <- if ("MergedRows" %in% names(summary_row)) summary_row$MergedRows else NA_integer_

  left_columns <- if ("LeftColumns" %in% names(summary_row)) summary_row$LeftColumns else NA_integer_
  right_columns <- if ("RightColumns" %in% names(summary_row)) summary_row$RightColumns else NA_integer_
  merged_columns <- if ("MergedColumns" %in% names(summary_row)) summary_row$MergedColumns else NA_integer_

  left_unique_keys <- if ("LeftUniqueKeys" %in% names(summary_row)) summary_row$LeftUniqueKeys else NA_integer_
  right_unique_keys <- if ("RightUniqueKeys" %in% names(summary_row)) summary_row$RightUniqueKeys else NA_integer_
  merged_unique_keys <- if ("MergedUniqueKeys" %in% names(summary_row)) summary_row$MergedUniqueKeys else NA_integer_

  # Build issue previews

  matching_preview <- if (nrow(MergeObj$IDCoverage$Matching) > 0) {
    MergeObj$IDCoverage$Matching %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  left_only_preview <- if (nrow(MergeObj$IDCoverage$LeftOnly) > 0) {
    MergeObj$IDCoverage$LeftOnly %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  right_only_preview <- if (nrow(MergeObj$IDCoverage$RightOnly) > 0) {
    MergeObj$IDCoverage$RightOnly %>%
      dplyr::slice_head(n = TopN) %>%
      apply(1, paste, collapse = " | ") %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  coverage_preview <- paste0(
    "<strong>Matching examples:</strong><br>",
    matching_preview,
    "<br><br><strong>Left only examples:</strong><br>",
    left_only_preview,
    "<br><br><strong>Right only examples:</strong><br>",
    right_only_preview
  )

  duplicate_key_preview <- if (
    nrow(MergeObj$DuplicateKeys$Left) > 0 ||
      nrow(MergeObj$DuplicateKeys$Right) > 0 ||
      nrow(MergeObj$DuplicateKeys$Merged) > 0
  ) {

    left_duplicate_preview <- if (nrow(MergeObj$DuplicateKeys$Left) > 0) {
      left_tmp <- MergeObj$DuplicateKeys$Left %>%
        dplyr::slice_head(n = TopN)

      left_tmp <- left_tmp[
        ,
        seq_len(min(ncol(left_tmp), 5)),
        drop = FALSE
      ]

      left_tmp %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    right_duplicate_preview <- if (nrow(MergeObj$DuplicateKeys$Right) > 0) {
      right_tmp <- MergeObj$DuplicateKeys$Right %>%
        dplyr::slice_head(n = TopN)

      right_tmp <- right_tmp[
        ,
        seq_len(min(ncol(right_tmp), 5)),
        drop = FALSE
      ]

      right_tmp %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    merged_duplicate_preview <- if (nrow(MergeObj$DuplicateKeys$Merged) > 0) {
      merged_tmp <- MergeObj$DuplicateKeys$Merged %>%
        dplyr::slice_head(n = TopN)

      merged_tmp <- merged_tmp[
        ,
        seq_len(min(ncol(merged_tmp), 5)),
        drop = FALSE
      ]

      merged_tmp %>%
        apply(1, paste, collapse = " | ") %>%
        paste(collapse = "<br>")
    } else {
      "None"
    }

    paste0(
      "<strong>Left duplicate examples:</strong><br>",
      left_duplicate_preview,
      "<br><br><strong>Right duplicate examples:</strong><br>",
      right_duplicate_preview,
      "<br><br><strong>Merged duplicate examples:</strong><br>",
      merged_duplicate_preview
    )
  } else {
    "None"
  }

  overlap_preview <- if (nrow(MergeObj$OverlappingVariables) > 0) {
    MergeObj$OverlappingVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::pull(Variable) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  unresolved_duplicate_preview <- if (nrow(MergeObj$DuplicateVariables) > 0) {
    MergeObj$DuplicateVariables %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          XVariable,
          " / ",
          YVariable,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  conflict_preview <- if (nrow(MergeObj$DuplicateVariables) > 0) {
    MergeObj$DuplicateVariables %>%
      dplyr::arrange(
        dplyr::desc(Conflicts),
        Variable
      ) %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          Conflicts,
          " conflicts; ",
          Agreement,
          "% agreement)"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  suspicious_preview <- if (nrow(MergeObj$SuspiciousConflicts) > 0) {
    MergeObj$SuspiciousConflicts %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Variable,
          " (",
          Agreement,
          "% agreement; ",
          LeftClass,
          " vs ",
          RightClass,
          ")"
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  key_type_preview <- if ("KeyTypes" %in% names(MergeObj) && nrow(MergeObj$KeyTypes) > 0) {
    key_tmp <- MergeObj$KeyTypes %>%
      dplyr::filter(
        LeftType != RightType |
          LeftType != MergedType |
          RightType != MergedType
      )

    if (nrow(key_tmp) > 0) {
      key_tmp %>%
        dplyr::slice_head(n = TopN) %>%
        dplyr::mutate(
          Text = paste0(
            Key,
            " (Left: ",
            LeftType,
            "; Right: ",
            RightType,
            "; Merged: ",
            MergedType,
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

  row_inflation_preview <- if ("RowInflationFactor" %in% names(summary_row)) {
    paste0(
      "Row inflation factor: ",
      summary_row$RowInflationFactor
    )
  } else {
    "Not available"
  }

  ready_preview <- if (isTRUE(MergeObj$ReadyForAnalysis)) {
    "No major merge-integrity blockers detected."
  } else {
    "Major merge-integrity blockers detected. Review duplicate keys, unresolved duplicate variables, and other warnings."
  }

  # Prepare checks table
  #
  # Dataset context (rows, columns, unique keys) is folded in as informational
  # rows rather than shown as separate fingerprint/summary cards.

  context_rows_df <- tibble::tibble(
    Status = "INFO",
    Check = c(
      "Rows",
      "Columns",
      "Unique keys"
    ),
    Count = c(
      as.numeric(merged_rows),
      as.numeric(merged_columns),
      as.numeric(merged_unique_keys)
    ),
    Details = c(
      paste0(
        "Left: ", left_rows,
        "; Right: ", right_rows,
        "; Merged: ", merged_rows, "."
      ),
      paste0(
        "Left: ", left_columns,
        "; Right: ", right_columns,
        "; Merged: ", merged_columns, "."
      ),
      paste0(
        "Unique key combinations - Left: ", left_unique_keys,
        "; Right: ", right_unique_keys,
        "; Merged: ", merged_unique_keys, "."
      )
    ),
    Examples = "None"
  )

  checks_table_df <- MergeObj$Checks %>%
    dplyr::mutate(
      DisplayCheck = dplyr::case_when(
        Check == "Unresolved Duplicate Variables" ~ "Unresolved duplicates",
        Check == "Overlapping Variables" ~ "Overlapping variables",
        Check == "Variable Conflicts" ~ "Variable conflicts",
        Check == "Suspicious Conflicts" ~ "Suspicious conflicts",
        Check == "Merge Readiness" ~ "Merge readiness",
        TRUE ~ Check
      ),
      Examples = dplyr::case_when(
        Check == "Key Types" ~ key_type_preview,
        Check == "Duplicate Keys" ~ duplicate_key_preview,
        Check == "Coverage" ~ coverage_preview,
        Check == "Row Inflation" ~ row_inflation_preview,
        Check == "Overlapping Variables" ~ overlap_preview,
        Check == "Unresolved Duplicate Variables" ~ unresolved_duplicate_preview,
        Check == "Variable Conflicts" ~ conflict_preview,
        Check == "Suspicious Conflicts" ~ suspicious_preview,
        Check == "Merge Readiness" ~ ready_preview,
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

  checks_table_df <- dplyr::bind_rows(
    context_rows_df,
    checks_table_df
  )

  checks_table <- reactable::reactable(
    checks_table_df,
    width = "100%",
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    striped = TRUE,
    compact = TRUE,
    pagination = FALSE,
    height = TableHeight,
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

      if (row$Status == "INFO") {
        return(NULL)
      }

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

  # Prepare coverage table (omitted entirely when there are no unmatched keys)

  unmatched_key_count <- nrow(MergeObj$IDCoverage$LeftOnly) +
    nrow(MergeObj$IDCoverage$RightOnly)

  coverage_section <- NULL

  if (unmatched_key_count > 0) {

    coverage_table_df <- dplyr::bind_rows(
      MergeObj$IDCoverage$Matching %>%
        dplyr::mutate(CoverageStatus = "Matching"),
      MergeObj$IDCoverage$LeftOnly %>%
        dplyr::mutate(CoverageStatus = "Left only"),
      MergeObj$IDCoverage$RightOnly %>%
        dplyr::mutate(CoverageStatus = "Right only")
    ) %>%
      dplyr::select(
        CoverageStatus,
        dplyr::all_of(key_columns)
      )

    coverage_table <- reactable::reactable(
      coverage_table_df,
      width = "100%",
      searchable = TRUE,
      filterable = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      striped = TRUE,
      compact = TRUE,
      pagination = FALSE,
      height = TableHeight,
      columns = list(
        CoverageStatus = reactable::colDef(
          name = "Coverage",
          width = 130,
          cell = function(value) {
            color <- dplyr::case_when(
              value == "Matching" ~ "#2E7D32",
              value == "Left only" ~ "#F9A825",
              value == "Right only" ~ "#EF6C00",
              TRUE ~ "#546E7A"
            )

            background <- dplyr::case_when(
              value == "Matching" ~ "#E8F5E9",
              value == "Left only" ~ "#FFF8E1",
              value == "Right only" ~ "#FFF3E0",
              TRUE ~ "#ECEFF1"
            )

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
              value
            )
          }
        )
      )
    )

    coverage_section <- htmltools::tags$details(
      class = "sdr-accordion",
      open = if (expanded) NA else NULL,
      htmltools::tags$summary(
        class = "sdr-accordion-summary",
        paste0(
          "Coverage explorer (",
          unmatched_key_count,
          " unmatched)"
        )
      ),
      htmltools::tags$div(
        class = "sdr-section-subtitle",
        "Review matching, left-only, and right-only key combinations."
      ),
      coverage_table
    )
  }

  # Prepare duplicate-variable conflict table (omitted when no .x/.y pairs)

  duplicate_variable_count <- nrow(MergeObj$DuplicateVariables)

  conflicts_section <- NULL

  if (duplicate_variable_count > 0) {

    duplicate_variables_df <- MergeObj$DuplicateVariables %>%
      dplyr::mutate(
        AgreementLabel = paste0(Agreement, "%"),
        ClassComparison = paste0(LeftClass, " vs ", RightClass)
      )

    conflicts_table <- reactable::reactable(
      duplicate_variables_df,
      width = "100%",
      searchable = TRUE,
      filterable = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      striped = TRUE,
      compact = TRUE,
      pagination = FALSE,
      height = min(TableHeight, 300),
      defaultSorted = "Conflicts",
      defaultSortOrder = "desc",
      columns = list(
        Variable = reactable::colDef(
          minWidth = 220
        ),
        XVariable = reactable::colDef(
          name = ".x variable",
          minWidth = 180
        ),
        YVariable = reactable::colDef(
          name = ".y variable",
          minWidth = 180
        ),
        LeftClass = reactable::colDef(
          show = FALSE
        ),
        RightClass = reactable::colDef(
          show = FALSE
        ),
        ClassComparison = reactable::colDef(
          name = "Classes",
          minWidth = 160
        ),
        Agreement = reactable::colDef(
          name = "Agreement",
          align = "center",
          width = 110,
          cell = function(value) {
            paste0(value, "%")
          },
          style = function(value) {
            if (is.na(value)) {
              list()
            } else if (value < 75) {
              list(
                fontWeight = "700",
                color = "#C62828",
                background = "#FFEBEE"
              )
            } else if (value < 100) {
              list(
                fontWeight = "700",
                color = "#F9A825",
                background = "#FFF8E1"
              )
            } else {
              list(
                color = "#2E7D32",
                background = "#E8F5E9"
              )
            }
          }
        ),
        AgreementLabel = reactable::colDef(
          show = FALSE
        ),
        Conflicts = reactable::colDef(
          align = "center",
          width = 100,
          style = function(value) {
            if (is.na(value) || value == 0) {
              list()
            } else {
              list(
                fontWeight = "700",
                color = "#C62828"
              )
            }
          }
        ),
        MissingnessConflicts = reactable::colDef(
          name = "Missingness conflicts",
          align = "center",
          minWidth = 170
        ),
        BothMissing = reactable::colDef(
          name = "Both missing",
          align = "center",
          width = 120
        ),
        TotalRows = reactable::colDef(
          name = "Rows",
          align = "center",
          width = 90
        )
      ),
      details = function(index) {
        selected_variable <- as.character(duplicate_variables_df$Variable[index])

        detail_df <- MergeObj$VariableConflicts %>%
          dplyr::filter(
            Variable == selected_variable
          ) %>%
          dplyr::slice_head(n = TopN)

        if (nrow(detail_df) == 0) {
          return(
            htmltools::tags$div(
              class = "sdr-detail-box",
              "No conflict examples available for this variable."
            )
          )
        }

        detail_df <- detail_df %>%
          dplyr::mutate(
            LeftValueEscaped = htmltools::htmlEscape(LeftValue),
            RightValueEscaped = htmltools::htmlEscape(RightValue),
            Conflict = paste0(
              "<span class='sdr-left-value'>",
              LeftValueEscaped,
              "</span>",
              " &rarr; ",
              "<span class='sdr-right-value'>",
              RightValueEscaped,
              "</span>"
            )
          )

        detail_display <- detail_df %>%
          dplyr::select(
            dplyr::all_of(conflict_key_columns),
            ConflictType,
            Conflict
          )

        reactable::reactable(
          detail_display,
          width = "100%",
          searchable = TRUE,
          filterable = TRUE,
          highlight = TRUE,
          bordered = TRUE,
          striped = TRUE,
          compact = TRUE,
          pagination = FALSE,
          height = min(TableHeight, 260),
          columns = list(
            Conflict = reactable::colDef(
              html = TRUE,
              minWidth = 260
            ),
            ConflictType = reactable::colDef(
              name = "Conflict type",
              minWidth = 210
            )
          )
        )
      }
    )

    conflicts_section <- htmltools::tags$details(
      class = "sdr-accordion",
      open = if (expanded) NA else NULL,
      htmltools::tags$summary(
        class = "sdr-accordion-summary",
        paste0(
          "Duplicate-variable conflicts (",
          duplicate_variable_count,
          " variable",
          ifelse(duplicate_variable_count == 1, "", "s"),
          ")"
        )
      ),
      htmltools::tags$div(
        class = "sdr-section-subtitle",
        "Expand a variable to review conflicting .x and .y values side by side."
      ),
      conflicts_table
    )
  }

  # Prepare suggested actions table (omitted when empty)

  actions_section <- NULL

  if ("SuggestedActions" %in% names(MergeObj) && nrow(MergeObj$SuggestedActions) > 0) {

    actions_table <- reactable::reactable(
      MergeObj$SuggestedActions,
      width = "100%",
      searchable = TRUE,
      filterable = TRUE,
      highlight = TRUE,
      bordered = TRUE,
      striped = TRUE,
      compact = TRUE,
      pagination = FALSE,
      height = min(TableHeight, 250),
      columns = list(
        Priority = reactable::colDef(
          width = 120,
          cell = function(value) {
            color <- dplyr::case_when(
              value == "High" ~ "#C62828",
              value == "Medium" ~ "#F9A825",
              TRUE ~ "#546E7A"
            )

            background <- dplyr::case_when(
              value == "High" ~ "#FFEBEE",
              value == "Medium" ~ "#FFF8E1",
              TRUE ~ "#ECEFF1"
            )

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
              value
            )
          }
        ),
        Action = reactable::colDef(
          minWidth = 500
        )
      )
    )

    actions_section <- htmltools::tags$div(
      class = "sdr-section",
      htmltools::tags$div(
        class = "sdr-section-title",
        "Suggested actions"
      ),
      htmltools::tags$div(
        class = "sdr-section-subtitle",
        "Recommended next steps generated from the merge audit."
      ),
      actions_table
    )
  }

  # Assemble dashboard

  dashboard <- htmltools::tagList(
    htmltools::tags$style(
      htmltools::HTML(
        "
        .sdr-dashboard {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          color: #263238;
          width: 100%;
  max-width: none;
  overflow-x: auto;
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

        .sdr-accordion {
          margin-top: 28px;
          margin-bottom: 10px;
          border: 1px solid #E0E6ED;
          border-radius: 14px;
          padding: 10px 16px;
          background: #F8FAFC;
        }

        .sdr-accordion[open] {
          background: #FFFFFF;
        }

        .sdr-accordion-summary {
          font-size: 21px;
          font-weight: 750;
          cursor: pointer;
          padding: 4px 0;
        }

        .sdr-accordion-summary:hover {
          color: #1565C0;
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

        .sdr-left-value {
          color: #C62828;
          background: #FFEBEE;
          padding: 2px 5px;
          border-radius: 5px;
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        }

        .sdr-right-value {
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
        "Interactive review of merge integrity from ValidateMerge()."
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Validation checks"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Search, filter, sort, and expand rows to inspect merge-integrity examples. INFO rows summarize rows, columns, and unique keys across the source and merged datasets."
        ),
        checks_table
      ),
      coverage_section,
      conflicts_section,
      actions_section
    )
  )

  return(dashboard)
}
