
#' Create Summary Table
#'
#' Generate a descriptive summary table for specified variables in a dataset.
#'
#' @param Data The dataset containing the variables of interest.
#' @param Variables A character vector specifying the variables for which summary statistics will be calculated.
#' @param numdecimals Number of decimal places to round the summary statistics.
#' @param Relabel Logical, indicating whether to use variable labels as column headers.
#' @return A formatted HTML table displaying summary statistics.
#' @importFrom magrittr %>%
#' @importFrom dplyr select rownames_to_column mutate if_else
#' @importFrom sjlabelled get_label set_label
#' @importFrom summarytools descr
#' @export
CreateSummaryTable <- function(Data, Variables = NULL, numdecimals = 2, Relabel = TRUE, Ordinal = FALSE, ScrollBoxHeight = "700px") {

    if (is.null(Variables)) {
      Variables = colnames(Data)
    }
    suppressWarnings({
      Data <- Data %>% select(all_of(Variables))
      l <- get_label(Data)
      if(Ordinal){
        Data <- ConvertOrdinalToNumeric(Data)
        Data <- lapply(Data, as.numeric) %>% as.data.frame()
        if(Relabel){
          sjlabelled::set_label(Data) <- l
        }
      }

      d <- summarytools::descr(Data)
      statVars <- c("Mean", "Std.Dev", "Median", "IQR", "Min",
                    "Max", "Skewness", "Kurtosis", "N.Valid", "Pct.Valid")
      d2 <- as.data.frame(t(as.data.frame(d)))
      d2 <- d2 %>% select(all_of(statVars))
      d2 <- d2[Variables, ]
      d2[statVars] <- lapply(d2[statVars], round, numdecimals)
      if (Relabel) {
        Data <- ReplaceMissingLabels(Data)
        labels <- sjlabelled::get_label(Data, def.value = colnames(Data)) %>%
          as.data.frame() %>% tibble::rownames_to_column()
        colnames(labels) <- c("Variable", "label")
        d2$label <- labels$label
        d2 <- d2 %>% select(label, all_of(statVars))
      }
      SummaryTable <- d2 %>% tibble::rownames_to_column("Variable") %>%
        mutate(Skewness = kableExtra::cell_spec(round(Skewness, numdecimals), "html",
                                    background = if_else(abs(Skewness) > 10, "yellow",
                                                         "", missing = "grey")), Kurtosis = kableExtra::cell_spec(round(Kurtosis,
                                                                                                            numdecimals), "html", background = if_else(abs(Kurtosis) >
                                                                                                                                               10, "yellow", "", missing = "grey")), IQR = kableExtra::cell_spec(round(IQR,
                                                                                                                                                                                                           numdecimals), "html", background = if_else(abs((IQR/Std.Dev/1.34) -
                                                                                                                                                                                                                                                  1) > 0.5, "yellow", "", missing = "grey")), Pct.Valid = kableExtra::cell_spec(round(Pct.Valid,
                                                                                                                                                                                                                                                                                                                          numdecimals), "html", background = if_else(Pct.Valid <
                                                                                                                                                                                                                                                                                                                                                             70, "red", "", missing = "grey"))) %>% kableExtra::kable(format = "html",
                                                                                                                                                                                                                                                                                                                                                                                                          escape = FALSE, digits = numdecimals, row.names = TRUE, caption = "Descriptive Summary Table. IQR, Skewness, and Kurtosis are highlighted in yellow if they are indicative of a non-normal distribution. Pct.Valid is highlighted in red if over 30% of data is missing") %>%
        kableExtra::kable_styling(bootstrap_options = c("striped", "hover",
                                            "condensed", "responsive")) %>% kableExtra::scroll_box(width = "100%",
                                                                                       height = ScrollBoxHeight)
    })
    return(SummaryTable)
  }
