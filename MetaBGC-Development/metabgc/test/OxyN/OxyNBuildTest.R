InputFiles.HMMRun <- "data/CombinedHmmSearch.txt"
InputFiles.BLAST_TP_NoCov <- "data/CombinedBLASTSearch.txt"
InputFiles.GeneIntervalPos <- "data/Gene_Interval_Pos.txt"
InputParam.HMM_Model_Name = 'Cyclase_OxyN'
InputParam.F1_Threshold = 0.50
OutputFiles.HMMOutDir <- "output/build/spHMMs"
OutputFiles.HMMHighPerfOutDir <- "output/build/HiPerf"
library(tidyverse)
library(ggsci)
library(ggpubr)
library(dplyr)