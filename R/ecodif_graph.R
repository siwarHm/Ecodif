#' ecodif_graph function
#' Make all the graphics
#'
#' @param object obj.dds, list of tables, parameters and datas
#' @return all the plots
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom grDevices rainbow
#' @importFrom psych geometric.mean
#' @example
#' data(datalist)
#' obj <- ecodif_loading(expr.table, annotation.species, design)
#' obj <- ecodif_sizefactor(obj)
#' obj <- ecodif_estidispe(obj)
#' obj <- ecodif_wald(obj)
#' graph <- ecodif_graph(obj)
#' @export
utils::globalVariables(c("species_short" , "file.species" , "objet" , "annotation.species"))
ecodif_graph <- function(object){

  stopifnot(class(object)=="ObjDDS" || class(object)=="DESeqResu")
  ##stopifnot(xor(class(object)=="ObjDDS",class(object)=="DESeqResu"))

  ## Association of a color to each Sample depending on the Condition

  #attach(object)
  classes <- unique(object$target.design$label)
  #classes <- unique(target.design$Label)
  classes.colors <- rainbow(n=length(classes))
  names(classes.colors) <- classes

  ## Add a color column to the target.design table
  ## The color depends of the Label
  object$target.design$color <- classes.colors[object$target.design$Label]
  #target.design$color <- classes.colors[target.design$Label]

  nb.species_short <- c(1:length(object$species_short))
  names(nb.species_short) <- species_short
  #table.data$Spe <- nb.species_short[annotation.species$Species]
  object$table.data$Spe <- nb.species_short[annotation.species$Species]

  ## Boxplot
  print(.f_boxplot(object))

  ## boxplot normalized data
  print(.f_boxplotnorm(object))

  ## Dispersion plot
  ## plots the per-gene dispersion estimates together with the fitted mean-dispersion relationship
  .f_dispe(object)

  #graph.two <- plotDispEsts(dds.global, main="Species normalisation-Ecosytem dispersion")

  ## MA plot
  #.f_maplot(object)

  ## we check if f_wald was loaded

  #if (class(object)=="DESeqResu"){

   # dds.resu <- object$resu

    ## Signification plot
    #print(.f_signiplot(dds.resu))

    ## Volcano plot
    #print(.f_volcano(dds.resu))

  #}
  #if (class(object)=="ObjDDS"){
  #  message("You have to run the function f_wald to obtain these graphics")
 # }

  detach(object)

}
