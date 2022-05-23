#' Main function
#' Contains all the functions
#'
#' @param table.data raw counts matrix
#' @param file.species matrix with species
#' @param design design decided by the user
#'
#' @return obj.dds, list of tables, parameters and datas
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom psych geometric.mean
#' @example
#' obj <- main_ecodif(expr.table, annotation.species, design)
#' @export

main_ecodif <- function(table.data, file.species, design){

  obj_ecodif <- ecodif_loading(table.data, file.species, design)

  obj_ecodif <- ecodif_sizefactor(obj_ecodif)

  obj_ecodif <- ecodif_estidispe(obj_ecodif)

  obj_ecodif <- ecodif_wald(obj_ecodif)
  return(obj_ecodif)
}
