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
write_tsv(t2pksBlastDF, "t2pksBlastDF_all-results.txt", col_names=T) 
