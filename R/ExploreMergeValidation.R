#' Explore merge validation results interactively
#'
#' Create an interactive HTML dashboard from a `ValidateMerge()` result object.
#' This function is designed for merge quality-control workflows. It displays
#' dataset fingerprints, high-level summary cards, a traffic-light checks table,
#' coverage diagnostics, join-variable auditing, and an expandable
#' duplicate-variable conflict explorer.
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
#'
#' @return An `htmltools::tagList()` object containing an interactive dashboard.
#'
#' @export
ExploreMergeValidation <- function(
    MergeObj,
    Title = "Merge validation explorer",
    TopN = 10,
    TableHeight = 350
) {

  # Validate inputs

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

  # Prepare shared values

  summary_row <- MergeObj$Summary[1, ]

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

  # Prepare fingerprint values

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

  suggested_actions_preview <- if ("SuggestedActions" %in% names(MergeObj) && nrow(MergeObj$SuggestedActions) > 0) {
    MergeObj$SuggestedActions %>%
      dplyr::slice_head(n = TopN) %>%
      dplyr::mutate(
        Text = paste0(
          Priority,
          ": ",
          Action
        )
      ) %>%
      dplyr::pull(Text) %>%
      paste(collapse = "<br>")
  } else {
    "None"
  }

  # Prepare summary card values

  duplicate_key_total <- if ("DuplicateKeyGroups_Left" %in% names(summary_row)) {
    summary_row$DuplicateKeyGroups_Left +
      summary_row$DuplicateKeyGroups_Right +
      summary_row$DuplicateKeyGroups_Merged
  } else {
    nrow(MergeObj$DuplicateKeys$Left) +
      nrow(MergeObj$DuplicateKeys$Right) +
      nrow(MergeObj$DuplicateKeys$Merged)
  }

  coverage_issue_count <- if ("LeftOnlyKeys" %in% names(summary_row)) {
    summary_row$LeftOnlyKeys + summary_row$RightOnlyKeys
  } else {
    nrow(MergeObj$IDCoverage$LeftOnly) +
      nrow(MergeObj$IDCoverage$RightOnly)
  }

  overlapping_variable_count <- if ("OverlappingVariables" %in% names(summary_row)) {
    summary_row$OverlappingVariables
  } else {
    nrow(MergeObj$OverlappingVariables)
  }

  unresolved_duplicate_count <- if ("UnresolvedDuplicateVariables" %in% names(summary_row)) {
    summary_row$UnresolvedDuplicateVariables
  } else {
    nrow(MergeObj$DuplicateVariables)
  }

  conflict_count <- if ("VariableConflictCount" %in% names(summary_row)) {
    summary_row$VariableConflictCount
  } else {
    nrow(MergeObj$VariableConflicts)
  }

  suspicious_count <- if ("SuspiciousConflictCount" %in% names(summary_row)) {
    summary_row$SuspiciousConflictCount
  } else {
    nrow(MergeObj$SuspiciousConflicts)
  }

  row_inflation_value <- if ("RowInflationFactor" %in% names(summary_row)) {
    summary_row$RowInflationFactor
  } else {
    NA_real_
  }

  row_inflation_status <- dplyr::case_when(
    is.na(row_inflation_value) ~ "PASS",
    row_inflation_value > 2 ~ "FAIL",
    row_inflation_value > 1.05 ~ "WARNING",
    TRUE ~ "PASS"
  )

  # Build merge fingerprint cards

  fingerprint_cards <- tibble::tibble(
    Dataset = c(
      "Left dataset",
      "Right dataset",
      "Merged dataset"
    ),
    Rows = c(
      left_rows,
      right_rows,
      merged_rows
    ),
    Variables = c(
      left_columns,
      right_columns,
      merged_columns
    ),
    UniqueKeys = c(
      left_unique_keys,
      right_unique_keys,
      merged_unique_keys
    )
  ) %>%
    dplyr::mutate(
      Preview = paste0(
        "Rows: ",
        Rows,
        "\nVariables: ",
        Variables,
        "\nUnique key combinations: ",
        UniqueKeys
      )
    )

  fingerprint_card_tags <- purrr::pmap(
    fingerprint_cards,
    function(Dataset, Rows, Variables, UniqueKeys, Preview) {
      htmltools::tags$div(
        class = "sdr-fingerprint-card",
        title = Preview,
        htmltools::tags$div(
          class = "sdr-fingerprint-title",
          Dataset
        ),
        htmltools::tags$div(
          class = "sdr-fingerprint-grid-inner",
          htmltools::tags$div(
            class = "sdr-fingerprint-metric",
            htmltools::tags$div(
              class = "sdr-fingerprint-value",
              Rows
            ),
            htmltools::tags$div(
              class = "sdr-fingerprint-label",
              "Rows"
            )
          ),
          htmltools::tags$div(
            class = "sdr-fingerprint-metric",
            htmltools::tags$div(
              class = "sdr-fingerprint-value",
              Variables
            ),
            htmltools::tags$div(
              class = "sdr-fingerprint-label",
              "Variables"
            )
          ),
          htmltools::tags$div(
            class = "sdr-fingerprint-metric",
            htmltools::tags$div(
              class = "sdr-fingerprint-value",
              UniqueKeys
            ),
            htmltools::tags$div(
              class = "sdr-fingerprint-label",
              "Unique keys"
            )
          )
        )
      )
    }
  )

  # Build merge quality cards

  summary_cards <- tibble::tibble(
    Metric = c(
      "Ready for analysis",
      "Duplicate keys",
      "Coverage issues",
      "Overlapping variables",
      "Unresolved duplicates",
      "Variable conflicts",
      "Suspicious conflicts",
      "Row inflation"
    ),
    Value = c(
      ifelse(isTRUE(MergeObj$ReadyForAnalysis), "Yes", "No"),
      as.character(duplicate_key_total),
      as.character(coverage_issue_count),
      as.character(overlapping_variable_count),
      as.character(unresolved_duplicate_count),
      as.character(conflict_count),
      as.character(suspicious_count),
      as.character(row_inflation_value)
    ),
    Status = c(
      ifelse(isTRUE(MergeObj$ReadyForAnalysis), "PASS", "FAIL"),
      ifelse(duplicate_key_total > 0, "FAIL", "PASS"),
      ifelse(coverage_issue_count > 0, "WARNING", "PASS"),
      ifelse(overlapping_variable_count > 0, "WARNING", "PASS"),
      ifelse(unresolved_duplicate_count > 0, "FAIL", "PASS"),
      ifelse(conflict_count > 0, "WARNING", "PASS"),
      ifelse(suspicious_count > 0, "WARNING", "PASS"),
      row_inflation_status
    ),
    Preview = c(
      ready_preview,
      duplicate_key_preview,
      coverage_preview,
      overlap_preview,
      unresolved_duplicate_preview,
      conflict_preview,
      suspicious_preview,
      row_inflation_preview
    )
  ) %>%
    dplyr::mutate(
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
    function(Metric, Value, Status, Preview, Color, Background) {
      htmltools::tags$div(
        class = "sdr-card",
        style = paste0(
          "border-left: 6px solid ",
          Color,
          "; background: ",
          Background,
          ";"
        ),
        title = gsub("<br>", "\n", Preview, fixed = TRUE),
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
          Status
        )
      )
    }
  )

  # Prepare checks table

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

  checks_table <- reactable::reactable(
    checks_table_df,
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

  # Prepare coverage table

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

  # Prepare join audit table

  join_audit_table_df <- MergeObj$JoinAudit %>%
    dplyr::mutate(
      IsKeyLabel = ifelse(IsKey, "Yes", "No")
    )

  join_audit_table <- reactable::reactable(
    join_audit_table_df,
    searchable = TRUE,
    filterable = TRUE,
    highlight = TRUE,
    bordered = TRUE,
    striped = TRUE,
    compact = TRUE,
    pagination = FALSE,
    height = min(TableHeight, 250),
    columns = list(
      Variable = reactable::colDef(
        minWidth = 220
      ),
      InBoth = reactable::colDef(
        show = FALSE
      ),
      IsKey = reactable::colDef(
        show = FALSE
      ),
      IsKeyLabel = reactable::colDef(
        name = "Specified key",
        width = 140,
        cell = function(value) {
          if (value == "Yes") {
            htmltools::tags$span(
              style = paste0(
                "display:inline-block;",
                "padding:4px 9px;",
                "border-radius:999px;",
                "font-weight:700;",
                "color:#2E7D32;",
                "background:#E8F5E9;"
              ),
              value
            )
          } else {
            htmltools::tags$span(
              style = paste0(
                "display:inline-block;",
                "padding:4px 9px;",
                "border-radius:999px;",
                "font-weight:700;",
                "color:#F9A825;",
                "background:#FFF8E1;"
              ),
              value
            )
          }
        }
      ),
      JoinRole = reactable::colDef(
        name = "Join role",
        minWidth = 180
      )
    )
  )

  # Prepare duplicate-variable conflict table

  duplicate_variables_df <- MergeObj$DuplicateVariables

  if (nrow(duplicate_variables_df) > 0) {
    duplicate_variables_df <- duplicate_variables_df %>%
      dplyr::mutate(
        AgreementLabel = paste0(Agreement, "%"),
        ClassComparison = paste0(LeftClass, " vs ", RightClass)
      )
  } else {
    duplicate_variables_df <- tibble::tibble(
      Variable = character(),
      XVariable = character(),
      YVariable = character(),
      LeftClass = character(),
      RightClass = character(),
      Agreement = numeric(),
      AgreementLabel = character(),
      Conflicts = integer(),
      MissingnessConflicts = integer(),
      BothMissing = integer(),
      TotalRows = integer(),
      ClassComparison = character()
    )
  }

  conflicts_table <- reactable::reactable(
    duplicate_variables_df,
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
      if (nrow(duplicate_variables_df) == 0) {
        return(NULL)
      }

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

  # Prepare suggested actions table

  actions_table <- if ("SuggestedActions" %in% names(MergeObj) && nrow(MergeObj$SuggestedActions) > 0) {
    reactable::reactable(
      MergeObj$SuggestedActions,
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
  } else {
    reactable::reactable(
      tibble::tibble(
        Message = "No suggested actions available."
      ),
      bordered = TRUE,
      compact = TRUE,
      pagination = FALSE,
      height = min(TableHeight, 200)
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

        .sdr-fingerprint-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
          gap: 12px;
          margin-bottom: 22px;
        }

        .sdr-fingerprint-card {
          background: #F8FAFC;
          border: 1px solid #E0E6ED;
          border-left: 6px solid #1565C0;
          border-radius: 14px;
          padding: 14px 16px;
          box-shadow: 0 1px 4px rgba(0,0,0,0.06);
        }

        .sdr-fingerprint-title {
          font-size: 14px;
          font-weight: 750;
          color: #37474F;
          margin-bottom: 10px;
        }

        .sdr-fingerprint-grid-inner {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 8px;
        }

        .sdr-fingerprint-metric {
          background: #FFFFFF;
          border-radius: 10px;
          padding: 10px 8px;
          text-align: center;
          border: 1px solid #EEF2F5;
        }

        .sdr-fingerprint-value {
          font-size: 24px;
          font-weight: 800;
          line-height: 1;
          color: #1565C0;
        }

        .sdr-fingerprint-label {
          margin-top: 6px;
          font-size: 11px;
          color: #607D8B;
          text-transform: uppercase;
          letter-spacing: 0.04em;
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
        class = "sdr-section-title",
        "Merge fingerprint"
      ),
      htmltools::tags$div(
        class = "sdr-section-subtitle",
        "Rows, variables, and unique key combinations before and after the merge."
      ),
      htmltools::tags$div(
        class = "sdr-fingerprint-grid",
        fingerprint_card_tags
      ),
      htmltools::tags$div(
        class = "sdr-section-title",
        "Merge quality summary"
      ),
      htmltools::tags$div(
        class = "sdr-section-subtitle",
        "High-level merge integrity indicators."
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
          "Search, filter, sort, and expand rows to inspect merge-integrity examples."
        ),
        checks_table
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Coverage explorer"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Review matching, left-only, and right-only key combinations."
        ),
        coverage_table
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Join audit"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Review variables present in both source datasets and whether they were specified as keys."
        ),
        join_audit_table
      ),
      htmltools::tags$div(
        class = "sdr-section",
        htmltools::tags$div(
          class = "sdr-section-title",
          "Duplicate-variable conflicts"
        ),
        htmltools::tags$div(
          class = "sdr-section-subtitle",
          "Expand a variable to review conflicting .x and .y values side by side."
        ),
        conflicts_table
      ),
      htmltools::tags$div(
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
    )
  )

  return(dashboard)
}
