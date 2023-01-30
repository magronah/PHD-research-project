library(DESeq2)
library(ggplot2)
library(tidyverse)
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions2.R")
#####################################################################
data_list <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))
#####################################################################
n <-2; filter <- 10
data  = data_list[[n]]; metadata = metadata_list[[n]]

##Deseq processes
dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
keep <- rowSums(counts(dds)) >= filter
dds <- dds[keep,]; data <- data[keep,] 
#############################################################################
#' Remove samples with all zero (deseq needs this to run effectively) 
r <- which(colSums(data) == 0)
if (length(r) > 0){
  data <- data[,-r]; metadata <- metadata[-r,]
  dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
}

#############################################################################
dds$Groups <- relevel(dds$Groups, ref = "NT")
ddsDE <- DESeq(dds, sfType ="poscounts") #4 was good minReplicatesForReplace=Inf
res <- results(ddsDE, cooksCutoff=FALSE, independentFiltering=FALSE) 
reslt <- lfcShrink(ddsDE, res = res, coef="Groups_ASD_vs_NT", type="normal")
result <- data.frame(reslt)

plotMA(reslt)
#############################################################################
k<-which(metadata$Groups == "NT")
control <- data[,k]
stopifnot(colnames(control) == metadata$Samples[k] )
logcontrol <- log2(rowMeans(control))

finit <- which(is.finite(logcontrol))
logcontrol <- logcontrol[finit]

foldchanges <- result$log2FoldChange[finit]
min(foldchanges)
max(foldchanges)

data <- data.frame(control=logcontrol, foldchanges)

ggplot(data, aes(control, foldchanges)) +
  geom_point()  + 
  geom_smooth() +
  #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))

#####################################################################
# I will take it one after the other and get it done. 