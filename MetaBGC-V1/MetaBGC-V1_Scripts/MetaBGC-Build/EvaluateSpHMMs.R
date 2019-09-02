
EvaluateSpHMM <- function(InputFiles.HMMRun,InputFiles.BLAST_TP_NoCov,InputFiles.GeneIntervalPos,InputParam.HMM_Model_Name,InputParam.F1_Threshold,OutputFiles.HMMOutDir,OutputFiles.HMMHighPerfOutDir){

	### Load segmented profiled HMMs for synthetic genomes
	#load HMM data and recode sampleType for the complexity of the synthetic sample (#of genomes in samples)
	hmm_df <- read_delim(InputFiles.HMMRun,col_names = F, delim = "\t") 
	names(hmm_df) <- c("readID", "sampleType", "sampleID", "cyclaseType", "HMMScore", "window","interval")
	hmm_df$sampleType[str_detect(hmm_df$sampleID, "high") ==TRUE] <-"high"
	hmm_df$sampleType[str_detect(hmm_df$sampleID, "low") ==TRUE] <-"low"

	##### Remove the reading frame information from the readID. 
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
	hmm_df_recoded<-formatHMM(hmm_df)

	##### Load the BLAST data for genes against our synthetic dataset. 
	# BLAST unfiltered reads at 95% pident no readCoverage filter
	all_blast_df <- read_delim(InputFiles.BLAST_TP_NoCov, col_names = T, delim = "\t")

	##### Positional information about domains and their locations in respect to the 30_10 spHMM models 
	gene_positions <- read_tsv(InputFiles.GeneIntervalPos,col_names = T)  %>% filter(start != 0)


	#### Positional read analysis in respect to location mapped to siderophore domain 
	##### Keep reads that map to the interval of a given model and covers the model 90% 

	filter_blast <- function(df, pos_df){
	  results <- data.frame()
	  for (i in 1:nrow(pos_df)){
		gene_data <- pos_df[i,]
		#print(gene_data$gene_name)
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

	#Filter BLAST reads that are within the genes intervals and cover 90% of the model interval
	blast_intervals<- filter_blast(all_blast_df, gene_positions)


	##### Determine manual cutoffs for each interval within the model
	compare_reads <- function(hmm_df, blast_df){
	  names(blast_df)[names(blast_df)=="sseqid"] <-"readID"
	  names(blast_df)[names(blast_df)=="qseqid"] <-"bgcName"
	  
	  # remove columns to compare the two dataframe 
	  blastDF <- blast_df %>% select(-c(model_cov, interval))
	  common_reads <- hmm_df %>% semi_join(.,blastDF) 
	  common_reads$readCheck<-"common-read"
	  hmm_unique_reads <- hmm_df %>% anti_join(.,blastDF)
	  hmm_unique_reads$readCheck<-"hmm-unique-read"
	  compared_data <- rbind(common_reads, hmm_unique_reads)
	  return(compared_data)
	}

	blast_bin <- compare_reads(hmm_df_recoded, blast_intervals)

	df_hmm_cutoff_scores <- blast_bin %>% filter(readCheck == "common-read") %>% group_by(interval) %>% mutate(medianScore = round(median(HMMScore))) %>% distinct(interval, readCheck, medianScore) %>% ungroup()
	names(df_hmm_cutoff_scores) <- c("interval", "read_check", "cutoff")

	#Filter lantibiotics data with cutoffs to compare to BLAST interval reads 
	filtered_median <- hmm_df_recoded
	filtered_median <- filtered_median[0,]
	for (i in 1:nrow(df_hmm_cutoff_scores)) {
	  intervalStr <- toString(df_hmm_cutoff_scores[i, "interval"])
	  cutoffScore <- as.numeric(df_hmm_cutoff_scores[i, "cutoff"])
	  filtered_median <- rbind(filtered_median,hmm_df_recoded %>% filter((interval == intervalStr & HMMScore>=cutoffScore)))
	}

	filtered_subfive <- hmm_df_recoded
	filtered_subfive <- filtered_subfive[0,]
	for (i in 1:nrow(df_hmm_cutoff_scores)) {
	  intervalStr <- toString(df_hmm_cutoff_scores[i, "interval"])
	  cutoffScore <- as.numeric(df_hmm_cutoff_scores[i, "cutoff"]) - 5.0
	  filtered_subfive <- rbind(filtered_subfive,hmm_df_recoded %>% filter((interval == intervalStr & HMMScore>=cutoffScore)))
	}

	filtered_plusfive <- hmm_df_recoded
	filtered_plusfive <- filtered_plusfive[0,]
	for (i in 1:nrow(df_hmm_cutoff_scores)) {
	  intervalStr <- toString(df_hmm_cutoff_scores[i, "interval"])
	  cutoffScore <- as.numeric(df_hmm_cutoff_scores[i, "cutoff"]) + 5.0
	  filtered_plusfive <- rbind(filtered_plusfive,hmm_df_recoded %>% filter((interval == intervalStr & HMMScore>=cutoffScore)))
	}

	#return the hmm-unique reads 

	compare_hmm_unique <- function(hmm_df, blast_df, pos_df){
      hmm_unique_df <- hmm_df  %>% inner_join(.,blast_df, by = c("readID"="sseqid", "Sample", "sampleType", "cyclaseType"))
	  intervals <- unique(hmm_unique_df$interval)
	  results <- data.frame()
	  #check that reads are in the same interval to throw out 
	  for (i in 1:length(intervals)){
		gene_interval_data <- pos_df %>% filter(interval == intervals[i])
		for (k in 1:nrow(gene_interval_data)){
			gene_data <- gene_interval_data[k,]
			#filters datatframe for edges and internal reads compared to model interval
			interval_df  <- hmm_unique_df %>% filter(qseqid ==gene_data$gene_name) %>%
			  filter(qstart %in% gene_data$start:gene_data$end | qend %in% gene_data$start:gene_data$end)
			if (nrow(interval_df) > 0 ){
			   results <- rbind(results,interval_df)
			}
		
	  }
	 }  
	  return(results%>% distinct())
	}

	return_hmm_unique <- function(hmm_df, blast_df){
	  oxyn_hmm_df <- hmm_df  %>% select(-c(window))
	  oxyn_hmm_df$interval <- as.character(oxyn_hmm_df$interval)
	  names(blast_df)[names(blast_df)=="sseqid"] <-"readID"
	  OxyN_hmm_unique <- oxyn_hmm_df %>% anti_join(.,blast_df, by= c("readID", "Sample", "sampleType", "interval"))
	  return(OxyN_hmm_unique)
	}

	gene_positions$cyclaseType <- InputParam.HMM_Model_Name
	all_blast_df$cyclaseType <- InputParam.HMM_Model_Name

	median_hmmunique <- return_hmm_unique(filtered_median, blast_intervals)
	median_hmmunique_less_model_cov <- compare_hmm_unique(median_hmmunique, all_blast_df, gene_positions) 
	median_remaining_hmm <- median_hmmunique %>% anti_join(.,median_hmmunique_less_model_cov)
	table(median_remaining_hmm$interval)

	subfive_hmmunique <- return_hmm_unique(filtered_subfive, blast_intervals)
	subfive_hmmunique_less_model_cov<-compare_hmm_unique(subfive_hmmunique,all_blast_df,gene_positions) 
	subfive_remaining_hmm <- subfive_hmmunique %>% anti_join(.,subfive_hmmunique_less_model_cov)
	table(subfive_remaining_hmm$interval)

	plusfive_hmmunique <- return_hmm_unique(filtered_plusfive, blast_intervals)
	plusfive_hmmunique_less_model_cov<-compare_hmm_unique(plusfive_hmmunique,all_blast_df,gene_positions) 
	plusfive_remaining_hmm <- plusfive_hmmunique %>% anti_join(.,plusfive_hmmunique_less_model_cov)
	table(plusfive_remaining_hmm$interval)



	calculate_F1 <- function(hmm_df, hmm_fp, blast_df, intervals){
	  #added this because factor vector 
	  hmm_df$interval <- as.character(hmm_df$interval)
	  hmm_df <- hmm_df %>% select(-c(window))
	  names(blast_df)[names(blast_df)=="sseqid"] <-"readID"
	  
	  results<-data.frame(interval=character(), F1=numeric())
	  for (i in 1:length(intervals)){
		model <- intervals[i]
		model_hmm <- hmm_df %>% filter(interval== model) # filter out model HMM results 
		model_blast <- blast_df %>% filter(interval== model) 
		model_hmm_fp <-hmm_fp %>% filter(interval== model) 
		TP <- model_hmm %>% inner_join(.,model_blast) %>% nrow()
		FP <- model_hmm_fp %>% nrow()
		FN <- model_blast %>% anti_join(.,model_hmm) %>% nrow()
		F1_metric <- (2*TP)/((2*TP)+FP+FN)
		row <- data.frame(cbind(model, F1_metric))
		results<- results %>% add_row(interval = model, F1 = F1_metric)
	  }
	  return(results)
	}

	f1_cutoff_median<- calculate_F1(filtered_median, median_remaining_hmm, blast_intervals, unique(hmm_df_recoded$interval))
	f1_cutoff_median$cutoff<- "median"

	f1_cutoff_plusfive<- calculate_F1(filtered_plusfive, plusfive_remaining_hmm, blast_intervals, unique(hmm_df_recoded$interval))
	f1_cutoff_plusfive$cutoff<- "+5"

	f1_cutoff_subfive<- calculate_F1(filtered_subfive, subfive_remaining_hmm, blast_intervals, unique(hmm_df_recoded$interval))
	f1_cutoff_subfive$cutoff<- "-5"

	f1_df <- rbind(f1_cutoff_median, f1_cutoff_plusfive, f1_cutoff_subfive)
	f1_cutoff_df <- f1_df %>% group_by(interval) %>% top_n(1, F1) %>% ungroup() %>%  filter(F1 >=as.numeric(InputParam.F1_Threshold)) %>% arrange(interval) %>% as.data.frame()

	colnames(f1_cutoff_df)[colnames(f1_cutoff_df)=="cutoff"] <- "cutoff_diff"
	colnames(f1_cutoff_df)[colnames(f1_cutoff_df)=="interval"] <- "Interval"

	var <- ("cutoff_diff")
	f1_cutoff_df[,var] <- sapply(f1_cutoff_df[,var],function(x) ifelse(x=='median',0,x))
	f1_cutoff_df <- base::merge(x=f1_cutoff_df,y=df_hmm_cutoff_scores,by="Interval")
	f1_cutoff_df$FinalCutoff <- as.numeric(f1_cutoff_df$cutoff_diff) + as.numeric(f1_cutoff_df$Cutoff)
	f1_cutoff_df <- subset(f1_cutoff_df, select = -c(cutoff_diff, Cutoff))
	cutoffFileName = paste0(InputParam.HMM_Model_Name, "_F1_Cutoff.txt")
	write_tsv(f1_cutoff_df,file.path(OutputFiles.HMMHighPerfOutDir,cutoffFileName),col_names = T)

	CopyHPModel <- function(x) {
	 spHMMInterval <- x[1]
	 spHMMFileName = paste0(InputParam.HMM_Model_Name, "__30_10__", spHMMInterval, ".hmm")
	 file.copy(file.path(OutputFiles.HMMOutDir,spHMMFileName), OutputFiles.HMMHighPerfOutDir)
	}
	apply(f1_cutoff_df, 1, CopyHPModel)

	supp_fig_lanc_f1 <-ggplot(data = f1_df, mapping = aes(x = interval, y = F1, group= cutoff, colour=cutoff )) + geom_point()+
	geom_line() + 
	ylim(0,1)+ 
	geom_hline(yintercept=c(as.numeric(InputParam.F1_Threshold)), linetype="dashed", color = c("red","green"), size=1) +
	theme_pubclean() +  
	scale_color_npg(name="HMM Score Cutoff")

	ggsave(supp_fig_lanc_f1, file=file.path(OutputFiles.HMMHighPerfOutDir,"plot.eps"), device="eps", width = 20, height = 7)

}

