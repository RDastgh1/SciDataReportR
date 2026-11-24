ExplorerTables <- function(node) {
  if (inherits(node, "datatables")) return(list(node))
  if (inherits(node, "shiny.tag")) return(unlist(lapply(node$children, ExplorerTables), recursive = FALSE))
  if (is.list(node)) return(unlist(lapply(node, ExplorerTables), recursive = FALSE))
  list()
}

