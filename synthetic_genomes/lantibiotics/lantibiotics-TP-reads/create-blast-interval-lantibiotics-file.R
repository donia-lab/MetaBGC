require(tidyverse) 

setwd("/tigress/DONIA/fcamacho/TypeII-PKS-Manuscript/synthetic-genomes/other-BGCs-poc/lantibiotics/lantibiotics-TP-reads")

#load terpene HMM data and recode sampleType for the complexity of the synthetic sample (#of genomes in samples)
lantibiotics <- read_delim("../combined-lantibiotics-spHMM-synthetic_genomes-results.txt",col_names = F, delim = "\t") 
names(lantibiotics) <- c("readID", "sampleType", "sampleID", "cyclaseType", "HMMScore", "window","interval")

lantibiotics$sampleType[str_detect(lantibiotics$sampleID, "high") ==TRUE] <-"high"
lantibiotics$sampleType[str_detect(lantibiotics$sampleID, "low") ==TRUE] <-"low"

lantibiotics$interval <- factor(lantibiotics$interval,levels = c(
  "0_30","10_40","20_50","30_60","40_70","50_80","60_90","70_100","80_110",
  "90_120","100_130","110_140","120_150","130_160","140_170","150_180", "160_190",
  "170_200", "180_210","190_220", "200_230", "210_240", "220_250", "230_260",
  "240_270", "250_280", "260_290", "270_300", "280_310", "290_320", "300_330",
  "310_340", "320_350", "330_360", "340_370", "350_380", "360_390", "370_400",
  "380_410", "390_420", "400_430", "410_440", "420_450"))

# Function to aggregate identitical sample reads located at different reading frames 
# and take the frame with the highest HMM score
formatHMM<-function(hmmdf){
  hmmdfRecoded <- separate(hmmdf, readID, into = c("readIDOnly","F_R_read_frame"), sep = "/", extra = "merge")
  hmmdfRecoded_FR <- separate(hmmdfRecoded, F_R_read_frame, into = c("F_R","frameNumb"), sep = "_", extra = "merge")
  hmmdfRecodedDF<- within(hmmdfRecoded_FR, readID <- paste(readIDOnly,F_R, sep='/'))
  hmmdfRecodedDFUnique<-aggregate(HMMScore ~ readID + sampleID + sampleType + cyclaseType + window + interval, hmmdfRecodedDF, max)
  colnames(hmmdfRecodedDFUnique)<-c("readID","Sample", "sampleType", "cyclaseType", "window", "interval","HMMScore")
  return(hmmdfRecodedDFUnique)
}

#Keep duplicated reads if they are in different reads
lantibiotics_recoded<-formatHMM(lantibiotics)


# BLAST unfiltered reads at 95% pident no readCoverage filter
all_lancBlastDF <- read_delim("lantibioticBlastDF_all-results.txt", col_names = T, delim = "\t")


#load siderophore interval positions data 
lanc_positions <- read_tsv("LanC_gene-interval_position-data.txt",col_names = T)  %>% filter(start != 0)

filter_blast <- function(df, pos_df){
  results <- data.frame()
  for (i in 1:nrow(pos_df)){
    gene_data <- pos_df[i,]
    print(gene_data$gene_name)
    #filters datatframe for edges and internal reads compared to model interval
    interval_df_1 <- df %>% filter(qseqid ==gene_data$gene_name) %>%
      filter(qstart %in% gene_data$start:gene_data$end | qend %in% gene_data$start:gene_data$end)
    # need to get reads that are bigger than the interval 
    interval_df_2 <- df %>% filter(qseqid ==gene_data$gene_name) %>% 
      filter(qstart < gene_data$start & qend > gene_data$end)
    interval_df <- rbind(interval_df_1,interval_df_2)
    if (nrow(interval_df) > 0){
      for (j in 1:nrow(interval_df)){
        in_interval_count <- sum(interval_df[j,]$qstart:interval_df[j,]$qend %in% gene_data$start:gene_data$end)
        interval_count <- length(gene_data$start:gene_data$end)
        interval_cov <- (in_interval_count/interval_count) * 100
        if (interval_cov >=90){
          res_df<-interval_df[j,]
          res_df$model_cov <- interval_cov
          res_df$interval <- gene_data$interval
          results <- rbind(results,res_df)
        }
      }
      
    }
    
  }
  return(results)
}
lanc_blast_intervals<- filter_blast(all_lancBlastDF, lanc_positions)

lanc_blast_intervals %>% write_tsv(.,"lantibiotics_blast_intervals.txt", col_names = T)

print("Completed!")
