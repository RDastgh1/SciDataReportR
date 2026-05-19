# =============================================================================
# CodebookMergeApp
# =============================================================================

#' Interactive codebook harmonization dashboard
#'
#' Launch a Shiny dashboard for reviewing and harmonizing multiple
#' codebooks before deterministic merging with MergeCodebooks().
#'
#' @param codebooks Named list of codebook data frames.
#' @param VariableCol Name of variable identifier column.
#' @param auto_type_mapping Logical; normalize common type synonyms.
#' @param ignore_columns Optional metadata columns to ignore.
#'
#' @return Launches a Shiny app.
#' @export

CodebookMergeApp <- function(
    codebooks,
    VariableCol = "Variable",
    auto_type_mapping = TRUE,
    ignore_columns = NULL
) {

  # Validate inputs --------------------------------------------------------

  if (!is.list(codebooks)) {
    stop("codebooks must be a named list.")
  }

  if (is.null(names(codebooks))) {
    stop("codebooks must be a named list.")
  }

  if (any(names(codebooks) == "")) {
    stop("All codebooks must have names.")
  }

  # Preserve original order -----------------------------------------------

  original_order <- purrr::imap_dfr(

    codebooks,

    function(cb, source_name) {

      tibble::tibble(
        SourceCodebook = source_name,
        !!VariableCol := cb[[VariableCol]],
        OriginalOrder = seq_len(nrow(cb))
      )
    }
  )

  # Helper functions -------------------------------------------------------

  normalize_type <- function(x) {

    x <- stringr::str_trim(as.character(x))

    dplyr::case_when(

      x %in% c(
        "Double",
        "double",
        "numeric",
        "Numeric"
      ) ~ "Continuous",

      x %in% c(
        "Integer",
        "integer"
      ) ~ "Continuous",

      TRUE ~ x
    )
  }

  # Structure overview -----------------------------------------------------

  structure_table <- purrr::imap_dfr(

    codebooks,

    function(cb, nm) {

      tibble::tibble(
        Source = nm,
        Column = colnames(cb)
      )
    }
  ) %>%

    dplyr::mutate(
      Present = "✓"
    ) %>%

    tidyr::pivot_wider(
      names_from = Source,
      values_from = Present
    )

  # Standardize codebooks --------------------------------------------------

  all_columns <- unique(
    unlist(
      purrr::map(
        codebooks,
        colnames
      )
    )
  )

  codebooks_std <- purrr::imap(

    codebooks,

    function(cb, source_name) {

      missing_cols <- setdiff(
        all_columns,
        colnames(cb)
      )

      if (length(missing_cols) > 0) {
        cb[missing_cols] <- NA
      }

      cb <- cb[, all_columns]

      cb <- cb %>%

        dplyr::mutate(
          dplyr::across(
            everything(),
            as.character
          )
        )

      if (
        auto_type_mapping &&
        "Type" %in% colnames(cb)
      ) {

        cb <- cb %>%

          dplyr::mutate(
            Type = normalize_type(Type)
          )
      }

      cb$SourceCodebook <- source_name

      cb
    }
  )

  merged_long <- dplyr::bind_rows(codebooks_std)

  metadata_cols <- setdiff(
    colnames(merged_long),
    c(
      VariableCol,
      "SourceCodebook",
      ignore_columns
    )
  )

  # Long comparison --------------------------------------------------------

  long_compare <- merged_long %>%

    tidyr::pivot_longer(
      cols = dplyr::all_of(metadata_cols),
      names_to = "Column",
      values_to = "Value"
    ) %>%

    dplyr::filter(
      !is.na(Value),
      Value != ""
    )

  # Conflicts --------------------------------------------------------------

  conflicts <- long_compare %>%

    dplyr::group_by(
      .data[[VariableCol]],
      Column
    ) %>%

    dplyr::summarise(
      n_unique = dplyr::n_distinct(Value),
      .groups = "drop"
    ) %>%

    dplyr::filter(
      n_unique > 1
    )

  # Variable presence ------------------------------------------------------

  variable_presence <- merged_long %>%

    dplyr::distinct(
      .data[[VariableCol]],
      SourceCodebook
    ) %>%

    dplyr::mutate(
      Present = "✓"
    ) %>%

    tidyr::pivot_wider(
      names_from = SourceCodebook,
      values_from = Present
    )

  # UI ---------------------------------------------------------------------

  ui <- shiny::fluidPage(

    theme = bslib::bs_theme(
      version = 5,
      bootswatch = "minty",
      primary = "#1E6F5C"
    ),

    shinyjs::useShinyjs(),

    shiny::titlePanel(
      "Codebook Harmonization Dashboard"
    ),

    shiny::tabsetPanel(

      shiny::tabPanel(

        "Overview",

        br(),

        shiny::fluidRow(

          shiny::column(
            6,

            div(
              class = "cardbox",

              h3("Variable Presence"),

              DT::DTOutput(
                "presence_table"
              )
            )
          ),

          shiny::column(
            6,

            div(
              class = "cardbox",

              h3("Structure Comparison"),

              DT::DTOutput(
                "structure_table"
              )
            )
          )
        )
      ),

      shiny::tabPanel(

        "Harmonization",

        br(),

        shiny::fluidRow(

          shiny::column(

            width = 7,

            div(
              class = "cardbox",

              h3("Conflict Browser"),

              DT::DTOutput(
                "harmonization_table"
              )
            )
          ),

          shiny::column(

            width = 5,

            div(
              class = "cardbox",

              h3("Resolution"),

              shiny::uiOutput(
                "selected_conflict_ui"
              ),

              br(),

              shinyWidgets::pickerInput(
                "resolution_choice",
                "Choose Resolution",
                choices = NULL
              ),

              shiny::conditionalPanel(
                condition = "input.resolution_choice == 'Custom'",

                shiny::textInput(
                  "custom_resolution",
                  "Custom Resolution"
                )
              ),

              br(),

              shiny::actionButton(
                "save_resolution",
                "Save Resolution",
                class = "btn-success"
              ),

              br(),
              br(),

              DT::DTOutput(
                "saved_resolutions"
              )
            )
          )
        )
      ),

      shiny::tabPanel(

        "Export",

        br(),

        shiny::fluidRow(

          shiny::column(
            6,

            div(
              class = "cardbox",

              h3("Merge Rules"),

              shiny::fluidRow(

                shiny::column(
                  6,

                  shiny::actionButton(
                    "generate_rules",
                    "Generate Rules"
                  )
                ),

                shiny::column(
                  6,

                  shiny::actionButton(
                    "copy_rules",
                    "Copy To Clipboard",
                    class = "btn-primary"
                  )
                )
              ),

              br(),

              shiny::textAreaInput(
                "rules_text",
                NULL,
                width = "100%",
                height = "550px"
              )
            )
          ),

          shiny::column(
            6,

            div(
              class = "cardbox",

              h3("Merged Codebook Preview"),

              DT::DTOutput(
                "merged_preview"
              ),

              br(),

              shiny::downloadButton(
                "download_merged",
                "Download CSV"
              )
            )
          )
        )
      )
    )
  )

  # Server ----------------------------------------------------------------

  server <- function(
    input,
    output,
    session
  ) {

    resolutions <- shiny::reactiveVal(

      tibble::tibble(
        Variable = character(),
        Column = character(),
        Resolution = character()
      )
    )

    # Overview tables ------------------------------------------------------

    output$presence_table <- DT::renderDT({

      DT::datatable(

        variable_presence,

        options = list(
          scrollY = "500px",
          paging = FALSE,
          scrollX = TRUE
        ),

        rownames = FALSE
      )
    })

    output$structure_table <- DT::renderDT({

      DT::datatable(

        structure_table,

        options = list(
          scrollY = "500px",
          paging = FALSE,
          scrollX = TRUE
        ),

        rownames = FALSE
      )
    })

    # Harmonization display ------------------------------------------------

    harmonization_display <- shiny::reactive({

      conflict_rows <- long_compare %>%

        dplyr::filter(
          paste(
            .data[[VariableCol]],
            Column
          ) %in%
            paste(
              conflicts[[VariableCol]],
              conflicts$Column
            )
        )

      wide_conflicts <- conflict_rows %>%

        dplyr::select(
          dplyr::all_of(VariableCol),
          Column,
          SourceCodebook,
          Value
        ) %>%

        tidyr::pivot_wider(
          names_from = SourceCodebook,
          values_from = Value
        )

      source_order <- names(codebooks)

      wide_conflicts$Resolution <- purrr::pmap_chr(

        wide_conflicts[, source_order],

        function(...) {

          vals <- list(...)

          first_nonmissing <- vals[
            !is.na(vals) &
              vals != ""
          ]

          if (length(first_nonmissing) == 0) {
            return(NA_character_)
          }

          first_nonmissing[[1]]
        }
      )

      saved_rules <- resolutions()

      if (nrow(saved_rules) > 0) {

        for (i in seq_len(nrow(saved_rules))) {

          matching_row <- which(
            wide_conflicts[[VariableCol]] ==
              saved_rules$Variable[i] &
              wide_conflicts$Column ==
              saved_rules$Column[i]
          )

          if (length(matching_row) > 0) {

            wide_conflicts$Resolution[
              matching_row
            ] <- saved_rules$Resolution[i]
          }
        }
      }

      wide_conflicts
    })

    # Harmonization table --------------------------------------------------

    output$harmonization_table <- DT::renderDT({

      df <- harmonization_display()

      dt <- DT::datatable(

        df,

        selection = "single",

        options = list(
          scrollY = "700px",
          paging = FALSE,
          scrollX = TRUE
        ),

        rownames = FALSE
      )

      source_cols <- names(codebooks)

      mismatch_rows <- which(

        purrr::map_lgl(
          1:nrow(df),

          function(i) {

            row_vals <- unlist(
              df[i, source_cols]
            )

            row_vals <- row_vals[
              !is.na(row_vals)
            ]

            length(unique(row_vals)) > 1
          }
        )
      )

      for (col in source_cols) {

        dt <- DT::formatStyle(

          dt,

          columns = col,

          backgroundColor = DT::styleEqual(
            df[[col]][mismatch_rows],
            rep(
              "#FFF3CD",
              length(mismatch_rows)
            )
          ),

          fontWeight = "bold"
        )
      }

      auto_resolution <- purrr::pmap_chr(

        df[, source_cols],

        function(...) {

          vals <- list(...)

          first_nonmissing <- vals[
            !is.na(vals) &
              vals != ""
          ]

          if (length(first_nonmissing) == 0) {
            return(NA_character_)
          }

          first_nonmissing[[1]]
        }
      )

      resolution_colors <- ifelse(
        df$Resolution == auto_resolution,
        "#F1F3F5",
        "#D1F2E1"
      )

      dt <- DT::formatStyle(

        dt,

        columns = "Resolution",

        backgroundColor = DT::styleEqual(
          df$Resolution,
          resolution_colors
        ),

        fontWeight = "bold"
      )

      dt
    })

    # Selected conflict ----------------------------------------------------

    selected_conflict <- shiny::reactive({

      shiny::req(
        input$harmonization_table_rows_selected
      )

      harmonization_display()[
        input$harmonization_table_rows_selected,
      ]
    })

    # Resolution UI --------------------------------------------------------

    output$selected_conflict_ui <- shiny::renderUI({

      shiny::req(selected_conflict())

      row <- selected_conflict()

      selected_value <- ifelse(
        input$resolution_choice %in% c("", NA),
        row$Resolution,
        input$resolution_choice
      )

      values <- row %>%

        dplyr::select(
          -dplyr::all_of(c(
            VariableCol,
            "Column",
            "Resolution"
          ))
        ) %>%

        tidyr::pivot_longer(
          everything(),
          names_to = "Source",
          values_to = "Value"
        )

      shiny::tagList(

        h4(
          paste0(
            row[[VariableCol]],
            " : ",
            row$Column
          )
        ),

        br(),

        purrr::map(
          1:nrow(values),

          function(i) {

            is_selected <- identical(
              values$Value[i],
              selected_value
            )

            bg_color <- ifelse(
              is_selected,
              "#D1F2E1",
              "#F8F9FA"
            )

            border_color <- ifelse(
              is_selected,
              "#1E6F5C",
              "#D9D9D9"
            )

            div(

              style = paste0(
                "
                padding:12px;
                margin-bottom:12px;
                background-color:", bg_color, ";
                border-radius:12px;
                border-left:6px solid ", border_color, ";
                "
              ),

              tags$b(values$Source[i]),
              br(),
              values$Value[i]
            )
          }
        )
      )
    })

    observe({

      shiny::req(selected_conflict())

      row <- selected_conflict()

      choices <- unique(
        unlist(
          row %>%

            dplyr::select(
              -dplyr::all_of(c(
                VariableCol,
                "Column",
                "Resolution"
              ))
            )
        )
      )

      choices <- choices[
        !is.na(choices)
      ]

      shinyWidgets::updatePickerInput(

        session,

        "resolution_choice",

        choices = c(
          choices,
          "Custom"
        ),

        selected = row$Resolution
      )
    })

    observeEvent(

      input$save_resolution,

      {

        shiny::req(selected_conflict())

        row <- selected_conflict()

        resolution_value <- ifelse(
          input$resolution_choice == "Custom",
          input$custom_resolution,
          input$resolution_choice
        )

        new_rule <- tibble::tibble(
          Variable = row[[VariableCol]],
          Column = row$Column,
          Resolution = resolution_value
        )

        resolutions(

          dplyr::bind_rows(
            resolutions(),
            new_rule
          ) %>%

            dplyr::group_by(
              Variable,
              Column
            ) %>%

            dplyr::slice_tail(
              n = 1
            ) %>%

            dplyr::ungroup()
        )
      }
    )

    # Export ---------------------------------------------------------------

    observeEvent(

      input$generate_rules,

      {

        rules_tbl <- resolutions()

        conflict_snapshot <- conflicts %>%

          dplyr::select(
            dplyr::all_of(VariableCol),
            Column
          )

        rules_code <- paste0(
          "# Auto-generated harmonization rules\n\n",
          "MergeRules <- list(\n\n",
          "  CodebookNames = c(\n",
          paste0(
            '    "',
            names(codebooks),
            '"',
            collapse = ",\n"
          ),
          "\n  ),\n\n",
          "  ConflictSnapshot = tibble::tribble(\n",
          paste0(
            "    ~",
            VariableCol,
            ", ~Column,\n\n"
          ),
          paste0(
            apply(
              conflict_snapshot,
              1,
              function(x) {
                paste0(
                  '    "',
                  x[1],
                  '", "',
                  x[2],
                  '"'
                )
              }
            ),
            collapse = ",\n"
          ),
          "\n\n  ),\n\n",
          "  VariableRules = tibble::tribble(\n",
          "    ~Variable, ~Column, ~PreferredValue,\n\n",
          if (nrow(rules_tbl) > 0) {
            paste0(
              apply(
                rules_tbl,
                1,
                function(x) {
                  paste0(
                    '    "',
                    x[1],
                    '", "',
                    x[2],
                    '", "',
                    x[3],
                    '"'
                  )
                }
              ),
              collapse = ",\n"
            )
          } else {
            ""
          },
          "\n\n  )\n",
          ")\n"
        )

        shiny::updateTextAreaInput(
          session,
          "rules_text",
          value = rules_code
        )
      }
    )

    observeEvent(

      input$copy_rules,

      {

        shinyjs::runjs("

          var copyText = document.getElementById('rules_text');

          copyText.select();

          copyText.setSelectionRange(0, 99999);

          document.execCommand('copy');

        ")

        shiny::showNotification(
          "Merge rules copied to clipboard.",
          type = "message"
        )
      }
    )

    # Preview --------------------------------------------------------------

    merged_preview <- shiny::reactive({

      merged_long %>%

        dplyr::left_join(
          original_order,
          by = c(
            VariableCol,
            "SourceCodebook"
          )
        ) %>%

        dplyr::arrange(
          SourceCodebook,
          OriginalOrder
        ) %>%

        dplyr::select(
          -OriginalOrder
        )
    })

    output$merged_preview <- DT::renderDT({

      DT::datatable(

        merged_preview(),

        options = list(
          scrollY = "650px",
          paging = FALSE,
          scrollX = TRUE
        ),

        rownames = FALSE
      )
    })

    output$download_merged <- shiny::downloadHandler(

      filename = function() {
        "MergedCodebook.csv"
      },

      content = function(file) {

        utils::write.csv(
          merged_preview(),
          file,
          row.names = FALSE
        )
      }
    )
  }

  shiny::shinyApp(ui, server)
}


# =============================================================================
# MergeCodebooks
# =============================================================================

#' Merge multiple codebooks using harmonization rules
#'
#' Deterministically merge multiple codebooks using optional harmonization
#' rules generated from `CodebookMergeApp()`.
#'
#' If the supplied rules match the current codebooks and conflict structure,
#' warnings are suppressed because the harmonization has already been reviewed.
#'
#' Warnings are only emitted when:
#' - new conflicts appear
#' - codebooks differ from those used to generate rules
#' - conflict structure changes
#'
#' @param codebooks Named list of codebook data frames.
#' @param Rules Optional harmonization rules generated from
#'   `CodebookMergeApp()`.
#' @param VariableCol Name of variable identifier column.
#' @param warn Logical; emit warnings.
#' @param strict Logical; stop on unresolved conflicts.
#'
#' @return A list containing:
#' \describe{
#'   \item{Codebook}{Merged harmonized codebook}
#'   \item{ConflictReport}{Detected conflicts}
#'   \item{AppliedRules}{Applied rules}
#' }
#'
#' @examples
#' \dontrun{
#'
#' MergedCB <- MergeCodebooks(
#'
#'   codebooks = list(
#'     Study1 = cb1,
#'     Study2 = cb2
#'   ),
#'
#'   Rules = MergeRules
#' )
#'
#' }
#'
#' @export

MergeCodebooks <- function(
    codebooks,
    Rules = NULL,
    VariableCol = "Variable",
    warn = TRUE,
    strict = FALSE
) {

  # Validate inputs --------------------------------------------------------

  if (!is.list(codebooks)) {
    stop("codebooks must be a named list.")
  }

  if (is.null(names(codebooks))) {
    stop("codebooks must be a named list.")
  }

  if (any(names(codebooks) == "")) {
    stop("All codebooks must have names.")
  }

  # Preserve original order -----------------------------------------------

  original_order <- purrr::imap_dfr(

    codebooks,

    function(cb, source_name) {

      tibble::tibble(
        SourceCodebook = source_name,
        !!VariableCol := cb[[VariableCol]],
        OriginalOrder = seq_len(nrow(cb))
      )
    }
  )

  # Standardize structure --------------------------------------------------

  all_columns <- unique(
    unlist(
      purrr::map(
        codebooks,
        colnames
      )
    )
  )

  codebooks_std <- purrr::imap(

    codebooks,

    function(cb, source_name) {

      missing_cols <- setdiff(
        all_columns,
        colnames(cb)
      )

      if (length(missing_cols) > 0) {
        cb[missing_cols] <- NA
      }

      cb <- cb[, all_columns]

      cb$SourceCodebook <- source_name

      cb
    }
  )

  merged_long <- dplyr::bind_rows(
    codebooks_std
  )

  metadata_cols <- setdiff(
    colnames(merged_long),
    c(
      VariableCol,
      "SourceCodebook"
    )
  )

  # Long comparison --------------------------------------------------------

  long_compare <- merged_long %>%

    dplyr::mutate(
      dplyr::across(
        everything(),
        as.character
      )
    ) %>%

    tidyr::pivot_longer(
      cols = dplyr::all_of(metadata_cols),
      names_to = "Column",
      values_to = "Value"
    ) %>%

    dplyr::filter(
      !is.na(Value),
      Value != ""
    )

  # Conflict report --------------------------------------------------------

  conflict_report <- long_compare %>%

    dplyr::group_by(
      .data[[VariableCol]],
      Column
    ) %>%

    dplyr::summarise(
      Values = paste(
        unique(Value),
        collapse = " | "
      ),
      n_unique = dplyr::n_distinct(Value),
      .groups = "drop"
    ) %>%

    dplyr::filter(
      n_unique > 1
    )

  # Prepare rules ----------------------------------------------------------

  if (
    !is.null(Rules) &&
    "VariableRules" %in% names(Rules)
  ) {

    rule_table <- Rules$VariableRules

  } else {

    rule_table <- tibble::tibble(
      Variable = character(),
      Column = character(),
      PreferredValue = character()
    )
  }

  # Rename variable column if needed --------------------------------------

  rule_table_renamed <- rule_table %>%

    dplyr::rename(
      !!VariableCol := Variable
    )

  # Merge metadata ---------------------------------------------------------

  merged_codebook <- long_compare %>%

    dplyr::left_join(
      original_order,
      by = c(
        VariableCol,
        "SourceCodebook"
      )
    ) %>%

    dplyr::group_by(
      .data[[VariableCol]],
      Column
    ) %>%

    dplyr::summarise(

      Value = {

        current_variable <- dplyr::first(
          .data[[VariableCol]]
        )

        current_column <- dplyr::first(
          Column
        )

        matching_rule <- rule_table_renamed %>%

          dplyr::filter(
            .data[[VariableCol]] ==
              current_variable,
            Column == current_column
          )

        if (nrow(matching_rule) > 0) {

          matching_rule$PreferredValue[1]

        } else {

          Value[which.min(OriginalOrder)]
        }
      },

      MinOrder = min(OriginalOrder),

      .groups = "drop"
    ) %>%

    dplyr::arrange(MinOrder) %>%

    tidyr::pivot_wider(
      names_from = Column,
      values_from = Value
    ) %>%

    dplyr::select(
      -MinOrder
    )

  # Unresolved conflicts ---------------------------------------------------

  unresolved_conflicts <- conflict_report %>%

    dplyr::anti_join(
      rule_table_renamed,
      by = c(
        VariableCol,
        "Column"
      )
    )

  # Determine whether warnings are needed ---------------------------------

  warnings_needed <- TRUE

  if (!is.null(Rules)) {

    if (
      "ConflictSnapshot" %in% names(Rules) &&
      "CodebookNames" %in% names(Rules)
    ) {

      same_codebooks <- identical(
        sort(Rules$CodebookNames),
        sort(names(codebooks))
      )

      current_snapshot <- conflict_report %>%

        dplyr::select(
          dplyr::all_of(VariableCol),
          Column
        ) %>%

        dplyr::arrange(
          .data[[VariableCol]],
          Column
        )

      saved_snapshot <- Rules$ConflictSnapshot %>%

        dplyr::arrange(
          .data[[VariableCol]],
          Column
        )

      same_conflicts <- identical(
        current_snapshot,
        saved_snapshot
      )

      if (
        same_codebooks &&
        same_conflicts
      ) {

        warnings_needed <- FALSE
      }
    }
  }

  # Warnings ---------------------------------------------------------------

  if (
    warn &&
    warnings_needed &&
    nrow(unresolved_conflicts) > 0
  ) {

    warning(
      paste0(
        nrow(unresolved_conflicts),
        " unresolved conflicts detected."
      )
    )
  }

  if (
    strict &&
    warnings_needed &&
    nrow(unresolved_conflicts) > 0
  ) {

    stop(
      "Unresolved conflicts detected in strict mode."
    )
  }

  # Return result ----------------------------------------------------------

  list(

    Codebook = merged_codebook,

    ConflictReport = conflict_report,

    AppliedRules = rule_table
  )
}
