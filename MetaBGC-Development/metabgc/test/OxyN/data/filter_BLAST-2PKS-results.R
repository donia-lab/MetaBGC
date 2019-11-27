#load necessary packages 
require(tidyverse) # 1.2.1

###################################################
#Load BLAST tabular files into a single dataframe. 
###################################################
t2pksBlastfiles<- list.files(path= "/tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/T2PKS-TP-reads",
                           pattern="*tabular-blast", full.names=T, recursive=TRUE)

parseBlastFiles<-function(files){

  df<-data.frame()
  for (x in 1:length(files)) {
    if (!file.size(files[x]) == 0) {
      t <- read.delim(files[x], header=F, stringsAsFactors = FALSE) # load file
      names(t)<-(c("sseqid","slen", "sstart", "send", "qseqid", "qlen", "qstart", "qend", "pident", "evalue"))
      t$Sample<-strsplit(files[x], "/")[[1]][9]
      t$sampleType<-strsplit(files[x], "/")[[1]][8]
      df<-rbind(t, df)
    } else {
      paste("Error file is empty",files[x], sep = ": ")
    }
  }
  return (df)
}
#concat all BLAST data files into a data frame 
t2pksBlastDF<-parseBlastFiles(t2pksBlastfiles)
###################################################
#Filter BLAST dataframe based on a read's 
#positional information
###################################################
filterBlastReads<-function(blastDF){

  filterDF<-data.frame()
  for (i in 1:nrow(blastDF)){
    internal <- NULL
    row <- blastDF[i,]
    smin <- min(row$sstart, row$send)
    smax <- max(row$sstart, row$send)
    qmin <- min(row$qstart, row$qend)
    qmax <- max(row$qstart, row$qend)
    percentReadCoveredByCluster <- ((smax - smin + 1) / (row$slen))* 100
    if(qmin - 1 >= smin - 1  && row$qlen - qmax >= row$slen - smax){
      internal <- TRUE
    }else{
      internal <- FALSE}

    if (internal == TRUE){
      if ((row$pident >= 95.0) && (percentReadCoveredByCluster >= 90)){
        row$readPos <- "internal"
        row$readCov <- percentReadCoveredByCluster
        filterDF<-rbind(row,filterDF )
      }} else if (internal == FALSE){
            if ((row$pident >= 95.0) && (percentReadCoveredByCluster >= 50)){
                row$readPos <- "edge"
                row$readCov <- percentReadCoveredByCluster
                filterDF<-rbind(row,filterDF)}}
  }
  return(filterDF)
}

t2pksBlastDF_filter <- filterBlastReads(t2pksBlastDF)
write_delim(t2pksBlastDF_filter, "t2pksBlastDF_filter-results.txt", col_names=T, delim = "\t")
