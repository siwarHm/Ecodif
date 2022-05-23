## Function .f_filtering

## calculate variance and counts the number of zero values per gene


#' Title

#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 scale_y_log10
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme_minimal
#'
#'
utils::globalVariables(c("target.design", "Sample" , "value" , "Label" , "parameters" ,"baseMean" ,"log2FoldChange","Significativity" , "pvalue" ))
.f_filtering <-function(raw.counts,table.data){

  ## Raw counts filter
  ## Zero-variance filter

  gene.var <- apply(raw.counts, 1, var,na.rm=TRUE)

  ## Define a Boolean vector (TRUE/FALSE)
  ## indicating whether each gene has a null variance.
  zero.var.genes <- gene.var == 0

  ## Discard genes having zeros in at least 90% of samples
  ## message("Applying threshold on the percent of non-zero counts per gene: ", parameters$gene.filter.percent.zeros, "%")

  ## A Boolean table indicating whether each count is null (TRUE) or not (FALSE)
  zero.counts <- raw.counts == 0
  max.per.gene <- apply(raw.counts, 1, max,na.rm=TRUE)

  ## Count number of zero values per gene (=row of the expression table)
  percent.zeros <- 100*apply(zero.counts, 1, sum,na.rm=T) / ncol(zero.counts)
  discarded.genes <- data.frame(
    zero.var = zero.var.genes,
    too.many.zeros = (percent.zeros > parameters$gene.filter.percent.zeros),
    too.small.counts = max.per.gene < parameters$gene.filter.min.count)

  ## Genes passing the filters are those for which
  ## all the discarding criteria are FALSE,
  ## i.e. the sum of the row is 0
  filtered.genes <- apply(discarded.genes, 1, FUN=function(x){!any(x)})
  filtered.genes[is.na(filtered.genes)] <- FALSE

  ## Select a matrix with the filtered genes,
  ## i.e. those not discarded by any criterion
  filtered.counts <- table.data[filtered.genes, ]

  test <- table(filtered.counts$Species) > 1 ## A modifier
  if (sum(test)>=2){
    filtered.counts[filtered.counts$Species %in% names(test)[test],]
  } else stop("We need at least 2 species with enough genes")


  return(filtered.counts)
}

## Function f_fusion
## Make the fusion between the matrix factors and the expression table

.f_fusion <- function(matrixFactors, table.data.modif){

  ## innerjoin merge the two tables
  fusion <- dplyr::inner_join(matrixFactors, table.data.modif, by="Spe")

  ## Number of discarded columns
  nb <- ncol(matrixFactors)

  selected.fusion <- fusion[,c(2:nb)]

  ## change col names
  colnames(selected.fusion) <- colnames(matrixFactors)[-1]
  ## transform the table in matrix
  selected.fusion <- as.matrix(selected.fusion)

  ## change NA values in 0
  table.data.modif[is.na(table.data.modif)] <- 0
  table.data.modif$Spe = NULL
  # expr.table.final <- expr.table.modif[,-p]

  listresults <- list(expr.table.final=table.data.modif, selected.fusion=selected.fusion)

  return(listresults)
}

## Function BoxPlot for Raw Counts

.f_boxplot <- function(object){

  transi <- object$obj.dds
  raw.counts <- as.data.frame(transi$raw.counts)
  counts.modif <- reshape2::melt(raw.counts, id.vars = NULL)
  ## Merge target with counts.modif
  colnames(counts.modif)[1] <- "Sample"

  merge.counts <- dplyr::inner_join(counts.modif, target.design, by="Sample")

  lab = unique(target.design$Label)

  ## BOXPLOT raw datas
  plot <- ggplot(merge.counts, aes(x= Sample, y= value, fill=Label )) + ggtitle("Boxplot of Raw Counts") + xlab("Species") + ylab("RawCounts") + geom_boxplot() + scale_y_log10() + scale_fill_igv(name = "Conditions")

  return(plot)
}

## Function BoxPlot for Normalized counts

.f_boxplotnorm <- function(object){

  normalized <- DESeq2::counts(object$obj.dds$dds.global, normalized = TRUE)
  normalized <- reshape2::melt(normalized)
  colnames(normalized)[2] <- "Sample"

  merge.normalized <- dplyr::inner_join(normalized, target.design, by="Sample")

  lab = unique(target.design$Label) ## Déjà fait plus haut à retirer peut être

  plot <- ggplot(merge.normalized, aes(x= Sample, y= value, fill=Label )) + ggtitle("Boxplot of Normalized Counts") +  xlab("Species") + ylab("Normalized Counts") + geom_boxplot() + scale_y_log10() + scale_fill_igv(name = "Conditions", labels = lab)

  return(plot)
}

## function for dispersion boxplot

.f_dispe <- function(object){


  graph.dispe <- plotDispEsts(object$obj.dds$dds.global, main="Species normalisation dispersion", genecol = "black", fitcol = "red", finalcol ="dodgerblue", legend = TRUE)

}

## Function for the MA plot

.f_maplot <- function(object){

  ## Note: we explicitly indicate the package name to avoid confusion with limma::plotMA()

  DESeq2::plotMA(object = object$obj.dds$dds.global)
}

## Function for the Significativity plot

.f_signiplot <- function(dds.resu){

  vecteur <- dds.resu$pvalue > 0.05
  vecteur[is.na(vecteur)] <- FALSE

  dds.resu$Significativity <- vecteur

  dds.resu <- as.data.frame(dds.resu)
  plot <- ggplot(dds.resu, aes(x= baseMean, y= log2FoldChange)) + ggtitle("Significativity Plot") + xlab("") + ylab("") + geom_point(aes(color=Significativity)) + geom_hline(yintercept=0,linetype = "dashed",color="red") + scale_color_uchicago("dark",alpha = 0.9) + theme_minimal()

  return(plot)
  }

## Function for the volcano plot

.f_volcano <- function(dds.resu){

  ## Volcano plot
  dds.resu <- as.data.frame(dds.resu)
  plot <- ggplot(data=dds.resu, aes(x=log2FoldChange, y=-log10(pvalue))) + ggtitle("Volcano Plot") + xlab("Log2 Fold Change") + ylab("-log10(pvalue)") + geom_point() + theme_minimal() + geom_hline(yintercept=-log10(0.05), col="red")

  return(plot)
 }
