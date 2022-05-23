#' ecodif_wald function
#' calculate Wald statistics
#'
#' @param obj.dds obj.dds, list of tables, parameters and datas
#' @param ... arguments to be passed to DESeq2::results
#' @return obj.dds, list of tables, parameters and datas
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom psych geometric.mean
#' @importFrom stats var
#' @example
#' data(datalist)
#' obj <- ecodif_wald(dds)
#' @export

ecodif_wald <- function(obj.dds, ...){
  stopifnot(class(obj.dds)=="ObjDDS")

  dds.global <- obj.dds$dds.global
  dds.global <- nbinomWaldTest(dds.global)

  obj.dds$dds.global <- dds.global

  resu <- DESeq2::results(dds.global, ...)

  list.return <- list(obj.dds=obj.dds, resu=resu)
  class(list.return) = "DESeqResu"
  return(list.return)
}
