require(tidyverse) # 1.2.1

###################################################
#Load BLAST tabular files into a single dataframe.
###################################################
lancBlastfiles<- list.files(path= "/tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/other-BGCs-poc/lantibiotics/hmm-unique-analysis/blast-results",
                             pattern="*tabular-blast", full.names=T, recursive=TRUE)


parseBlastFiles<-function(files){

  df<-data.frame()
  for (x in 1:length(files)) {
    if (!file.size(files[x]) == 0) {
      t <- read.delim(files[x], header=F, stringsAsFactors = FALSE) # load file
      names(t)<-(c("sseqid","slen", "sstart", "send", "qseqid", "qlen", "qstart", "qend", "pident", "evalue"))
      filename <-strsplit(files[x], "/")[[1]][12]
      t$Sample <-strsplit(filename, "-")[[1]][4]
      t$cutoff<-strsplit(files[x], "/")[[1]][11]
      df<-rbind(t, df)
    } else {
      paste("Error file is empty",files[x], sep = ": ")
    }
  }
  return (df)
}
#concat all BLAST data files into a data frame
lancBlastDF<-parseBlastFiles(lancBlastfiles)

write_tsv(lancBlastDF, "lanc_hmm-unique_blast-results.txt", col_names = T )
