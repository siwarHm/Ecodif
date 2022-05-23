#' ecodif_sizefactor function
#' extract the size factors
#'
#' @param obj.dds obj.dds, list of tables, parameters and datas
#' @return obj.dds, list of tables, parameters and datas
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom psych geometric.mean
#' @importFrom igraph %>%
#' @importFrom tibble rownames_to_column
#'
#' @example
#' data(datalist)
#' obj <- ecodif_loading(expr.table, annotation.species, design)
#' obj <- ecodif_sizefactor(obj)
#' @export
utils::globalVariables("Geneid")
ecodif_sizefactor <- function(obj.dds){

  stopifnot(class(obj.dds)=="ObjDDS")
  matrixFactors <- NULL
  target.modif <- obj.dds$target.modif
  All.sizeFactor <- obj.dds$All.sizeFactor
  table.data <- obj.dds$table.data
  matrixFactors = t(do.call(cbind,All.sizeFactor))


  names.col <- colnames(matrixFactors)
  rownames(matrixFactors) <- names(All.sizeFactor)


  matrixFactors <- matrixFactors %>% as.data.frame %>% rownames_to_column("Spe")


  table.data.modif <- table.data %>% as.data.frame %>% dplyr::select(-Geneid)
  #table.data.modif <- table.data.modif %>% rename(Species="Spe")

  colnames(table.data.modif)[ncol(table.data.modif)]<- "Spe"

  nrow(table.data.modif)
  p <- ncol(table.data.modif)

  ## function for the fusion of the two tables

  list.results <- .f_fusion(matrixFactors, table.data.modif)


  expr.table.final <- list.results$expr.table.final
  print(nrow(expr.table.final))
  selected.fusion <- list.results$selected.fusion
  print(nrow(selected.fusion))

  ncol(expr.table.final)
  nrow(target.modif)

  rownames(expr.table.final) <- table.data$Geneid

  ## creation of the global dds

  dds.global <- DESeqDataSetFromMatrix(
    countData = expr.table.final,
    colData = target.modif,
    design = ~ Condition + Time_point + Condition:Time_point)

  normalizationFactors(dds.global) <- selected.fusion

  obj.dds$dds.global <- dds.global

  return(obj.dds)

}
