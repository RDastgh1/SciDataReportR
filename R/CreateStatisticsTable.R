
#' Create Statistics Table
#'
#' Generate a table of statistics including means, standard deviations, counts, and p-values.
#'
#' @param Data The data frame containing the variables of interest.
#' @param TargetVar The target variable for which statistics will be calculated.
#' @return A formatted HTML table displaying statistics.
#' @importFrom arsenal tableby tests
#' @importFrom kableExtra kable kable_styling scroll_box cell_spec
#' @importFrom data.table as.data.table
#' @importFrom stats as.formula
#' @export

CreateStatisticsTable <- function(Data,  TargetVar){

  tab1<-tableby(as.formula(paste(TargetVar, "~ .")), data = Data, digits.pct = 1, digits.count = 1, numeric.stats = c("Nmiss2", "meansd"), cat.test = "chisq")
  sd<-as.data.frame(summary(tab1, text = FALSE))
  colnames(sd)[length(colnames(sd))]<-'p_value'
  pvalscript<- sd$p_value
  pvalscript[pvalscript == "< 0.001"]<- "0"
  pvals<-as.numeric(pvalscript)
  pvals[is.na(pvals)]<-1
  colnames(sd)[1]<-'Variable'
  t<- tests(tab1)


  options(kableExtra.auto_format = T)

  sd<-as.data.table(sd)
  StatTable <- sd%>% mutate(p_value = cell_spec(p_value , "html", background = ifelse(pvals<0.05, "yellow", "")))%>%dplyr::rename("p-value" = p_value)%>%
    kable(format = "html", escape = F)%>% kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))%>%
    scroll_box(width = "100%", height = "500px")
  return(StatTable)
}
