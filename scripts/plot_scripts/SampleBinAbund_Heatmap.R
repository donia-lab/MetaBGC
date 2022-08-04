require(tidyverse)
require(reshape2)
require(pheatmap)
setwd("C:/Users/ab50/Documents/data/MetaBGCRuns/AcbK-homologs/output/quantify/combined_2")

save_pheatmap_pdf <- function(x, filename, width=10, height=75) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

generate_heatmap <- function(data_wide,metadata,ann_colors,cluster_row,cluster_col,cluster_algo='centroid') {
  heatmap <- pheatmap::pheatmap(data_wide, 
                                cluster_row = cluster_row,
                                cluster_cols = cluster_col,
                                show_colnames = TRUE,
                                show_rownames = FALSE,
                                annotation_row = annotations,
                                annotation_colors = ann_colors,
                                clustering_distance_rows = 'correlation',
                                clustering_method = cluster_algo,
                                scale = "none",
                                breaks=bk,
                                color = hmcols,
                                fontsize = 6.5,
                                fontsize_col = 6, 
                                cellwidth = 5,
                                cellheight = 2)
  return(heatmap)
}

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  return(x)
}

#load data 
hmp_df <- read_tsv("SampleAbundanceBodysite.tsv") 
hmp_df$Subject<-NULL
hmp_df$Bodysite<-NULL
hmp_df$Subject_Status<-NULL
hmp_df <- hmp_df %>% arrange(BodyAggSite, Cohort)

# Select Bins
hmp_df <- hmp_df %>% select(Sample, Cohort, BodyAggSite,"11","494","95","948","2", "36", "139", "1", "24", "92", "21", "223", "725", "30", "90")
hmp_nz  <- apply(hmp_df[c(4:ncol(hmp_df))], 1, function(x) all(x==0))
hmp_df <- hmp_df[!hmp_nz,]

metadata <- hmp_df %>% select(Sample, Cohort, BodyAggSite ) %>% distinct() %>% group_by(Cohort) %>% arrange(desc(BodyAggSite)) %>% ungroup()
annotations<-metadata %>% as.data.frame()
rownames(annotations)<-annotations$Sample
annotations$Sample<-NULL

hmp_df$Cohort<-NULL
hmp_df$BodyAggSite<-NULL
row.names(hmp_df)<-hmp_df$Sample
hmp_df$Sample<-NULL

#Colors for heatmap  
bk = unique(c(seq(0,0, length=100), seq(0,1000,length=1000)))
hmcols<- colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(bk)-1)

# Log10 the data
#hmp_df <- log(hmp_df)
#hmp_df[hmp_df == -Inf] <- 0
#bk = seq(0,max(as.numeric(unlist(hmp_df))),length=1000)
#hmcols<- colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(bk)-1)

#Setup annotation colors
ann_colors = list(
  Cohort = c(Chinese="green",Fiji="blue",HMP="brown","HMP-1_2"="darkred",MetaHit="darkorange"),
  BodyAggSite = c(Gut = "cornflowerblue",Oral = "cornsilk4",Skin = "black",Vaginal = "coral")
)
annotations$Cohort <- as.factor(annotations$Cohort)
annotations$BodyAggSite <- as.factor(annotations$BodyAggSite)
annotations <- annotations[order(annotations$Cohort,annotations$BodyAggSite),]

data_wide<-as.matrix(hmp_df)
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,F)
save_pheatmap_pdf(heatmap, "Heatmap_NoCluster.pdf",10,72)

# Clustered Heatmaps
# 1. Single
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'single')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples_Single.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'single')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins_Single.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'single')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_Single.pdf",10,72)


# 2. UPGMA
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'average')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples.pdf_UPGMA",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'average')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins.pdf_UPGMA",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'average')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_UPGMA.pdf",10,72)


# 3. WPGMA
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'mcquitty')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples_WPGMA.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'mcquitty')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins.pdf_WPGMA",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'mcquitty')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_WPGMA.pdf",10,72)


# 4. WPGMC
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'median')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples_WPGMC.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'median')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins_WPGMC.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'median')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_WPGMC.pdf",10,72)


# 5. UPGMC
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'centroid')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples_UPGMC.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'centroid')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins_UPGMC.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'centroid')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_UPGMC.pdf",10,72)


# 6. complete
heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,F,'complete')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamples_Complete.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,F,T,'complete')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterBins_Complete.pdf",10,72)

heatmap <- generate_heatmap(data_wide,annotations,ann_colors,T,T,'complete')
save_pheatmap_pdf(heatmap, "Heatmap_ClusterSamplesBins_Complete.pdf",10,72)


