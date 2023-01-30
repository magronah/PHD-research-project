## for bioconductor packages:
library(BiocManager)
if (FALSE) BiocManager::install(c("DESeq2", "edgeR"))
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions2.R")
#####################################################################
data_lis <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))

#####################################################################
n <-1
data  = data_list[[n]]; metadata = metadata_list[[n]]

##Deseq processes
dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
filter <- 10
keep <- rowSums(counts(dds)) >= filter
dds <- dds[keep,]; data <- data[keep,] 

#' Remove samples with all zero (deseq needs this to run effectively) 
r <- which(colSums(data) == 0)
if (length(r) > 0){
  data <- data[,-r]; metadata <- metadata[-r,]
  dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
}

dds$Groups <- relevel(dds$Groups, ref = "NT")

ddsDE <- DESeq(dds, sfType ="poscounts",minReplicatesForReplace=Inf) #4 was good 
res <- results(ddsDE, cooksCutoff=TRUE, independentFiltering=TRUE) 

reslt <- lfcShrink(ddsDE, res = res, coef="Groups_ASD_vs_NT", type="ashr")
result <- data.frame(reslt)

K<- which(metadata$Samples == "NT")
control <- data[,K]
logcontrol <- log2(rowSums(control))
finit <- which(is.finite(logcontrol))
View(logcontrol)

###############################################################
foldchange <- result$log2FoldChange
foldchange <- foldchange[finit]
control <- logcontrol[finit]


View(foldchange)


############################################
data <- data.frame(control, foldchange = result$log2FoldChange)
View(result)
ggplot(data, aes(control, foldchange)) +
  geom_point()  + 
  geom_smooth() +
  #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))


##############################
#DeseqResults[[index]] <- result

########## Method 1
library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]

all(rownames(coldata) == colnames(cts))


##############################################
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
dds <- DESeq(dds)
rownames(coldata) <- sub("fb", "", rownames(coldata))
######################
##Filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "untreated")
res <- results(dds)

#####################################################################
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
plotMA(res)
summary(res)
DESeq2::plotMA(res)

#####################################################################
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, main="apeglm")
plotMA(resNorm, main="normal")
plotMA(resAsh, main="ashr")

#####################################################################
########## Method 2

Control_Treat_Split <- function(data_list, metadata_list ){

  index <- 1; L =treat = list(); control <- list()
  for(n in names(metadata_list)){

    data <- data_list[[n]];  meta <- metadata_list[[n]]
    m <-  rownames(data)[1]; p <- colnames(data)[1]
    
    if (nchar(m) < nchar(p) ) {
      data <- t(data)
    } else {
      data <- (data)
    }

    asd <- which(meta$Groups == "ASD" )
    treat_data <- data[asd,]; control_data <- data[-asd,]
    stopifnot(rownames(treat_data) == meta$Samples[asd])
    stopifnot(rownames(control_data) == meta$Samples[-asd])

    L[[index]] <-  (list(treatment = treat_data, control =control_data))
    cont_mean[[index]] <- colMeans()
    treat_mean[[index]] <-
    effects <- log2(colMeans(treat_mean)/colMeans(treat_mean))  
    index <- index + 1
    
  }
  names(L) <- names(metadata_list)
  return(L)
}

v <- Control_Treat_Split(data_list, metadata_list )

data.frame( )
View(v[[1]][["control"]])

for(n in names(v)){
  con <- which(colSums(v[[n]][["control"]]) == 0)
  tre <- which(colSums(v[[n]][["treatment"]]) == 0)
  print(c(length(con), length(tre)))
}
# as.data.frame(DIMS(data_lis))

#colnames(metadata_list[[1]])
###################################################
M <- read_datasets(data_lis,metadata_list) ##'@read_datasets: reads dataset as list
RDataFiles(M,path) ##'@RDataFiles: creates filtered data, control, treatment and metadata files

nam <- c("PRJNA168470", "PRJNA355023", "PRJNA453621", "PRJNA578223", 
         "PRJEB45948","PRJNA644763", "PRJNA589343", "PRJNA624252",
         "PRJNA687773","PRJNA642975")

data_list <- mget(load( paste0(path,"data_new.RData"))) # list of filtered dataset
stopifnot(names(data_list[["data_list"]]) == nam)

metadata_list <- mget(load(paste0(path,"groups.RData"))) #  list of metadata
stopifnot(names(metadata_list) == nam)

control_list <- mget(load(paste0(path,"control.RData")))  # list of control dataset
stopifnot(names(control_list[[names(control_list)]]) == nam)
k <- which(colSums(control_list[[names(control_list)]][[1]]) == 0)
           
length(k)

treat_list <- mget(load(paste0(path,"treat.RData") ))     #list of treatment dataset
stopifnot(names(treat_list[[names(treat_list)]]) == nam)

##############################################
control <- control_list[["cont_list"]]
treat <- treat_list[["treat_list"]]
data_list<- data_list[["data_list"]]
#metadata_list <- metadata_list[["meta_list"]]

##Just check how much datasets is been thrown away from filtering
before_filter <- DIMS(data_lis); after_filter <- DIMS(data_list)

before_filter <- data.frame(before_filter); rownames(before_filter) <- c("samples","taxa")
after_filter <- data.frame(after_filter) ; rownames(after_filter) <- c("samples","taxa")

before_filter
after_filter
######################################################################
cont_treat_list<-read_ContTreat(control,treat) ##'@read_ContTreat: creates a list of control and treatment data
LogMeanData(cont_treat_list,path) ##' @LogMeanData: computes the log means of control and treatment

deseq_res(data_list, metadata_list) ##' @deseq_res: computes foldchanges, Pvalues, etc, from deseq2 
##################################################################
### Read data created by LogMeanData and deseq_res functions
Treat_ControlMeans <- mget(load(paste0(path,"Control_N_TreatmentMeans.RData"))) #list of log mean control and treatments
stopifnot(names(Treat_ControlMeans[[names(Treat_ControlMeans)]] ) == nam)

ContMeans <- mget(load(paste0(path,"ControlMeans.RData"))) #list of log mean of control
stopifnot(names(ContMeans[[names(ContMeans)]] ) == nam)

TreatMeans <- mget(load(paste0(path,"TreatmentMeans.RData"))) #list of log mean of treatments
stopifnot(names(TreatMeans[[names(TreatMeans)]] ) == nam)


DeseqRes <- mget(load(paste0(path,"DeseqResults.RData"))) #list of deseq results
stopifnot(names(DeseqRes[[names(DeseqRes)]] ) == nam)


deseq_list <- DeseqRes[[names(DeseqRes)]]   #extra Deseq data from the list
logcontrol <- ContMeans[[names(ContMeans)]]  #extra control mean datasets from the list
logtreatment <- TreatMeans[[names(TreatMeans)]]  #extra control treatment datasets from the list

plot_path <- "Latex_Files/Proposal/Figs/power/"
Plot_ContEffects(deseq_list,logcontrol,plot_path) ##' @Plot_ContEffects: plots and saves plots for control verse effect sizes

###########Fit normal distribution to control abundances
file_path <- path
plotname <- "His_QQ_control_"
Hist_QQ(logcontrol,plot_path,file_path, plotname) 



dds <- DESeqDataSetFromMatrix(t(data), meta_test, ~Groups)
#ddsDE <- DESeq(dds, sfType ="poscounts",minReplicatesForReplace=Inf)  
#res <- results(ddsDE, alpha = 0.05, contrast = c("Groups", "ASD","NT"), 
#               cooksCutoff=FALSE, independentFiltering=FALSE)


ddsDE <- DESeq(dds,sfType ="poscounts")
#resUnShrank <- data.frame(res)
#resultsNames(ddsDE)

resNorm <- lfcShrink(ddsDE, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")

resLFC <- lfcShrink(ddsDE, coef = "Groups_NT_vs_ASD", type="apeglm")

# res1 <- results(ddsDE, contrast=c("Groups","NT","ASD"))
# res1 <- data.frame(res1)

plotMA(resNorm, main="resNorm")

xlim <- c(1,1e5); ylim <- c(-3,100)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")

resShrink <- data.frame(resLFC)
plot(resUnShrank$log2FoldChange)

Res <- data.frame(res)

View(Res)
