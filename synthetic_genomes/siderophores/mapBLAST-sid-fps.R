require(tidyverse)
#R Script to map false positives using BLAST 

setwd("/Users/francinecamacho/Google Drive/T2PKS-paper/computational_analysis/synthetic_genomes/revision-synthetic-genomes/other-BGCs-poc/siderophores/HMM-analysis")

blastx <- read_tsv("false-positives-subfive/combined-siderophore-detected_reads-fps-dereplicated_blastx_VS_BLAST_NR-tabular", col_names =F) 
names(blastx)<-c("sseqid","stitle", "readID","slen","qlen","pident","evalue","qcovs","qstart","qend")

#BLASTx de-replicated cluster info 
clusterDF<- read_tsv("false-positives-subfive/combined-siderophore-detected_reads-fps-dereplicated-cluster_data.txt",col_names = T)
#representative reads per cluster using cd-hit-est
rep_clusters<- clusterDF %>% filter(clstr_rep == 1)  %>% select(c(1:2))
names(rep_clusters)[1]<-"readID"

#Function to map the FPs and TP using text 
checkBLAST <- function(blastdf, clusters){
  TP_blastdf<- blastdf %>% filter( grepl('siderophore|Siderophore', stitle)) %>%  group_by(readID) %>% filter(row_number()==1)%>% select(2:3) %>% ungroup()
  TP_blastdf$readCheck<- "TP"
  FP_blastdf<- blastdf %>% filter(!readID %in% TP_blastdf$readID) %>% group_by(readID) %>% filter(row_number()==1) %>% select(2:3) %>% ungroup()
  FP_blastdf$readCheck<- "FP"
  blastdf_check<-NULL
  if(nrow(FP_blastdf) > 0){
    blastdf_check <-rbind(TP_blastdf, FP_blastdf)
  }else{
    blastdf_check<-TP_blastdf
  }
  final_df <- blastdf_check %>% full_join(.,clusters) 
  final_df[is.na(final_df)] <- "FP" # recode nohits to FPs 
  return(final_df)
  }



#text analysis for TP/FP using BLASTx stitle 
blastx_DF<-checkBLAST(blastx, rep_clusters) 

# Map back the blastx mapping to all reads 

reads_mapped <- clusterDF %>% inner_join(., blastx_DF %>% select(c(1,3,4)), by = "clstr") %>% select(c(1,2,8,9))%>% 
  separate(., id, c("Sample", "readID"), sep = "__", remove = T)

#write_tsv(reads_mapped %>% select(-c(clstr)), "blastx-fp-subfive.txt", col_names = T)

