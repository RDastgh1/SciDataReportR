ReviewQueueFixture <- function(with_id = FALSE) {
  df_Data <- data.frame(ID = sprintf("%03d", 1:40), Visit = as.Date("2020-01-01") + 1:40,
    Site = rep(c("A", "B"), each = 20))
  df_Data$Visit[1:2] <- as.Date(c("2099-01-01", "2098-01-01"))
  attr(df_Data$Visit, "label") <- "Assessment date"
  attr(df_Data$ID, "label") <- "Study identifier"
  attr(df_Data$Site, "label") <- "Assessment site"
  DetectDataAnomalies(df_Data, names(df_Data), id_var = if (with_id) "ID" else NULL, methods = "dates")
}

EditReviewWorkbook <- function(file, transform, sheet = "Review") {
  wb <- openxlsx::loadWorkbook(file)
  df_Sheet <- openxlsx::read.xlsx(wb, sheet, check.names = FALSE, sep.names = " ",
    skipEmptyRows = FALSE, skipEmptyCols = FALSE, na.strings = character())
  df_Sheet <- transform(df_Sheet)
  openxlsx::removeWorksheet(wb, sheet)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df_Sheet)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

