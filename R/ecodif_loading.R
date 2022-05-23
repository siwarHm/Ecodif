#' ecodif_loading function
#' Load the data with inputs, tests, filtering and creation of the dds object
#'
#' @param table.data raw counts matrix
#' @param file.species matrix with species
#' @param design design decided by the user
#'
#' @return obj.dds, list of tables, parameters and datas
#' @import DESeq2 tidyverse gdata contrast ggsci
#' @importFrom data.table melt
#' @importFrom psych geometric.mean
#' @importFrom dplyr mutate
#' @importFrom  dplyr across
#' @importFrom dplyr inner_join
#' @importFrom igraph %>%
#' @importFrom tibble column_to_rownames
#' @example
#' data(datalist)
#' obj <- ecodif_loading(expr.table, annotation.species, design)
#' @export


utils::globalVariables("Species")
utils::globalVariables("where")
ecodif_loading <- function(table.data, file.species, design){

  ##-------------------Beginning of the tests------------

  stopifnot(is.character(table.data[,1])) # stop the program if the data is not annotated
  stopifnot(sum(table.data[,1]%in%file.species[,1])>2) # stop the program if the data does not enough correspond to the target
  species_short <- unique(file.species[,2]) # creation of the short list of species

  ## Rename the columns in the user file

  colnames(file.species) <- c("Geneid","Species")
  colnames(table.data)[1]<- "Geneid"

  raw.counts <- table.data[,-1] # remove sample names
  rownames(raw.counts) <- table.data$Geneid
  nblines <- nrow(raw.counts)

  ## Check if there are enough samples

  if (nblines < 1)
  {
    stop("Not enough genes")
  } ## 1 à remplacer


  ## Selection of the design

  target.modif <- design %>% column_to_rownames(colnames(design)[1]) # The column samples becomes row names

  target.modif <- target.modif %>% mutate(across(where(is.character),as.factor)) # transform the data in factors

  ## Check that sample IDs are corresponding between sample descriptions (rows) and count table (columns).

  differing.IDs <- colnames(raw.counts) %in% rownames(target.modif)
  stopifnot(sum(differing.IDs) > 4)
  raw.counts = raw.counts[,differing.IDs]

  differing.IDs.table <- colnames(table.data) %in% rownames(target.modif)
  table.data = table.data[,c(1,which(differing.IDs.table))]
  ## Check in the other way

  differing.IDs <- rownames(target.modif) %in% colnames(raw.counts)
  if (any(!differing.IDs))
  {
    stop("The sample IDs differ between count table (column names) and sample descriptions (row names).")
  }

  ##-------------------End of the tests------------

  ## Creation of the table data

  table.data = table.data %>% inner_join(file.species)

  ## Function for filtering the raw counts

  filtered.counts <- .f_filtering(raw.counts,table.data) ## internal function
  print(nrow(filtered.counts))

  file.species <- filtered.counts[,c(1,ncol(filtered.counts))]
  print(nrow(file.species))


  ## Initialization of the DDS list
  All.dds <- list()
  All.sizeFactor <- list()


  for(i in unique(filtered.counts$Species))
  {

    ##select lines of one species
    filtered.counts.spe <- dplyr::filter(filtered.counts,Species==i)
    genenames <- filtered.counts.spe$Geneid

    ## We have to convert in dataframe
    filtered.counts.spe <- as.data.frame(filtered.counts.spe[,-c(1,ncol(filtered.counts.spe))])
    rownames(filtered.counts.spe)<- genenames

    CT <- as.matrix(filtered.counts.spe)
    CT[which(CT==0)] <- NA
    geomean <- apply(CT,1,geometric.mean,na.rm=TRUE)
    #print(geomean)

    filtered.counts.spe[is.na(filtered.counts.spe)] = 0
    #filtered.counts.spe[(filtered.counts.spe==0)] <- 1

    ncol(filtered.counts.spe)
    nrow(target.modif)

    ddsTemp <- DESeqDataSetFromMatrix(
      countData = filtered.counts.spe,
      colData = target.modif,
      design = ~ Condition + Time_point + Condition:Time_point)

    All.dds[[i]] <- ddsTemp

    ddsTemp <- estimateSizeFactors(ddsTemp,geoMeans=geomean)
    TableSizeFactor <- sizeFactors(ddsTemp)
    TableSizeFactor <- as.data.frame(TableSizeFactor)

    All.sizeFactor[[i]] <- TableSizeFactor
  }

  obj.dds = list(All.sizeFactor=All.sizeFactor, All.dds=All.dds,table.data=filtered.counts, file.species=file.species, design=design, target.modif=target.modif, species_short=species_short, raw.counts=raw.counts)
  class(obj.dds) = "ObjDDS"
  return(obj.dds)
}
