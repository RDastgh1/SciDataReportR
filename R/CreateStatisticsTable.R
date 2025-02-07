
#' Create Statistics Table
#'
#' Generate a table of statistics including means, standard deviations, counts, and p-values.
#'
#' @param Data The data frame containing the variables of interest.
#' @param TargetVar The target variable for which statistics will be calculated.
#' @return A formatted HTML table displaying statistics.
#' @importFrom dplyr mutate rename
#' @importFrom kableExtra cell_spec kable kable_styling scroll_box
#' @suggests kableExtra, arsenal, data.table
#' @export

CreateStatisticsTable <- function(Data,  TargetVar){

  tab1<-arsenal::tableby(as.formula(paste(TargetVar, "~ .")), data = Data, digits.pct = 1, digits.count = 1, numeric.stats = c("Nmiss2", "meansd"), cat.test = "chisq")
  sd<-as.data.frame(summary(tab1, text = FALSE))
  colnames(sd)[length(colnames(sd))]<-'p_value'
  pvalscript<- sd$p_value
  pvalscript[pvalscript == "< 0.001"]<- "0"
  pvals<-as.numeric(pvalscript)
  pvals[is.na(pvals)]<-1
  colnames(sd)[1]<-'Variable'
  t<- arsenal::tests(tab1)


  options(kableExtra.auto_format = T)

  sd<-data.table::as.data.table(sd)
  StatTable <- sd%>% mutate(p_value = kableExtra::cell_spec(p_value , "html", background = ifelse(pvals<0.05, "yellow", "")))%>%dplyr::rename("p-value" = p_value)%>%
    kableExtra:: kable(format = "html", escape = F)%>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))%>%
    kableExtra::scroll_box(width = "100%", height = "500px")
  return(StatTable)
}
