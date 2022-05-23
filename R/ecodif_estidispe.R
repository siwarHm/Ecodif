#' ecodif_estidispe function
#' calculate the dispersion factors
#'
#' @param obj.dds obj.dds, list of tables, parameters and datas
#' @return obj.dds, list of tables, parameters and datas
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom psych geometric.mean
#' @example
#' data(datalist)
#' obj <- ecodif_loading(expr.table, annotation.species, design)
#' obj <- ecodif_sizefactor(obj)
#' obj <- ecodif_estidispe(obj)
#' @export

ecodif_estidispe <- function(obj.dds){
  stopifnot(class(obj.dds)=="ObjDDS")

  dds.global <- obj.dds$dds.global
  dds.global <- estimateDispersions(dds.global)

  obj.dds$dds.global <- dds.global
  return(obj.dds)
}
