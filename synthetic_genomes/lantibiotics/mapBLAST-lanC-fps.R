require(tidyverse)
#R Script to map false positives using BLAST 

setwd("/Users/francinecamacho/Google Drive/T2PKS-paper/computational_analysis/synthetic_genomes/revision-synthetic-genomes/other-BGCs-poc/lantibiotics/HMM-analysis")

blastx <- read_tsv("false-positives-subfive/combined-lantibiotics-detected_reads-fps-dereplicated_blastx_VS_BLAST_NR-tabular", col_names =F) 
names(blastx)<-c("sseqid","stitle", "readID","slen","qlen","pident","evalue","qcovs","qstart","qend")

#BLASTx de-replicated cluster info 
clusterDF<- read_tsv("false-positives-subfive/combined-lantibiotics-detected_reads-fps-dereplicated-cluster_data.txt",col_names = T)
#representative reads per cluster using cd-hit-est
rep_clusters<- clusterDF %>% filter(clstr_rep == 1)  %>% select(c(1:2))
names(rep_clusters)[1]<-"readID"

#Function to map the FPs and TP using text 
checkBLAST <- function(blastdf, clusters){
  TP_blastdf<- blastdf %>% filter( grepl('lantibiotic|Lantibiotic|cyclase|Cyclase|lantipeptide|lanthionine', stitle)) %>%  group_by(readID) %>% filter(row_number()==1)%>% select(2:3) %>% ungroup()
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
#blastx_DF<-checkBLAST(blastx, rep_clusters)
blastx_DF<-checkBLAST(blastx, rep_clusters) %>% separate(.,readID, c("Sample", "readID"), sep = "__", remove = T) %>% separate(.,readID,c("scaffold_id", "rest_id"), sep = "-", remove = F)

blastx_DF$lanc_genome <-"other"

lanc_genomes <-c("CP000419.1", "CP001834.1", "CP002365.1", "HE613569.1", "CP008815.1", "gi|1098135721|emb|FOXX01000007.1|",
                 "gi|1098137419|emb|FOXX01000002.1|", "AZXI01000001.1", "AZXI01000044.1", "AZXI01000068.1", "gi|529222834|ref|NC_021985.1|",
                 "gi|640933032|gb|JJOH01000019.1|", "JOFH01000021.1")

blastx_DF$lanc_genome[blastx_DF$scaffold_id %in% lanc_genomes] <-"lantibiotic_genome"

#blastx_DF %>% select(c(3,2,1,6,8)) %>% write_tsv(., "blastx-fp-subfive-dereplicated.txt", col_names = T)

# Map back the blastx mapping to all reads 
reads_mapped <- clusterDF %>% inner_join(., blastx_DF %>% select(c(1,6,7,8)), by = "clstr") %>% select(c(1,2,8,9,10))%>% 
  separate(., id, c("Sample", "readID"), sep = "__", remove = T)



#reads_mapped %>% separate(readID, c("scaffold_id", "rest_id"), sep = "-", remove = F) %>% distinct(scaffold_id) %>% filter(scaffold_id %in% lanc_genomes)
write_tsv(reads_mapped %>% select(-c(clstr)), "blastx-fp-subfive.txt", col_names = T)

