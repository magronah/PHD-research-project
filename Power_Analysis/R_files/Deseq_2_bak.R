library(microbiome)  
library(DESeq2)
library(tidyverse)
library("pasilla")
library(metR)  
library(truncdist)
library(mgcv)
library(distr)
library(fitdistrplus)
source("Power_Project/R_files/Deseq_Functions_bak.R")
path <- "Power_Project/Datasets/"
plot_path <- "Latex_Files/Proposal/Figs/power/"
name <- "deseq2" 

########### Preprocessing
pasCts <- system.file("extdata","pasilla_gene_counts.tsv", 
                      package="pasilla", mustWork=TRUE)

pasAnno <- system.file("extdata","pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)

cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

########### Datasets
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ condition)
keep <- rowSums(counts(dds)) >= 10; dds <- dds[keep,]

dds$condition <- relevel(dds$condition, ref = "untreated")

dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
reslt <- lfcShrink(dds, res=res, coef=2, type="normal") 
v <-data.frame(reslt)

####### Control/ Treatment split
cts <- counts(dds)
control <- cts[,4:7]; treatment <- cts[,1:3]

cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
logcontrol = log2(cont_Means) ; dispers = dispersions(dds)

effects<- v$log2FoldChange

data <- data.frame(Control_abundance=logcontrol, Fold_change =effects)

pp<- ggplot(data, aes(Control_abundance, Fold_change)) +
  geom_point()  +
  geom_smooth() +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
 
############################################################################
alpha = 0.1
power <-(!is.na(v$padj) & v$padj < alpha)
#effects <- effect_#v$log2FoldChange#/v$lfcSE
#power = v$lfcSE

comb <- tibble(abs_lfc = abs(effects),lcontrol = logcontrol, power)
gg_2dim <- (ggplot(comb)
            + aes(lcontrol, abs_lfc)
            + geom_point(aes(color = power), alpha = 0.5)
            + labs(x = "log(control abundance)",
                   y = "absolute(log(fold change))")
)

fit_2d <- bam(power ~ te(lcontrol, abs_lfc),
              data = comb)
#              family = binomial) #
#plot(fit_2d)

pp <- with(comb,
           expand.grid(lcontrol = seq(min(lcontrol),
                                      max(lcontrol),
                                      length = 50),
                       abs_lfc = seq(min(abs_lfc),
                                     max(abs_lfc),
                                     length = 50)))

pp$power <- predict(fit_2d, newdata = pp,
                    type = "response")

################
brkvec <- seq(0,1,0.01)
gg_2dimc <- (gg_2dim
             + geom_contour(data = pp,
                            aes(z=power),
                            breaks = brkvec)
             + geom_label_contour(data = pp, aes(z= power),
                                  breaks = brkvec)
)
gg_2dimc

# not filtering doesn't change the pattern
################

# dispersion is the same, effects is the same and control is also the sames. 
Disp_Parameters <- mget(load("Power_Project/Datasets/Dispersion_Parameters_deseq.RData"))
Cont_Parameters <- mget(load("Power_Project/Datasets/Control_Parameters_deseq.RData"))
Effect_Param <- mget(load(paste0(path,"Effect_Parameters_deseq.RData"))) 

param.cont <- Cont_Parameters[[names(Cont_Parameters)]]
param.disp <- Disp_Parameters[[names(Disp_Parameters)]]
param.effect <- Effect_Param[[names(Effect_Param)]]


scal = data.frame(x =  exp(param.effect[1]+param.effect[2]*logcontrol + 
                             param.effect[3]*logcontrol^2))

effect_ <-apply(scal,1,function(x){
  rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})

mu= param.cont$mu; sigma = param.cont$sigma; lambda = param.cont$lambda
myMix <- UnivarMixingDistribution(Norm(mean=mu[1], sd =sigma[1]),
                                  Norm(mean=mu[2], sd =sigma[2]),
                                  mixCoeff=c(lambda[1],lambda[2]))

rmyMix <- r(myMix)
n <- length(effects)
control <- rmyMix(n)


n_otu <- 1000; n_samples <-  100; n_sim = 300
#data_sim(param.effect,param.disp,param.cont,n_sim,n_samples,n_otu,name,filter=10)

data_list <- mget(load(paste0(path,"Simualated_Data_",name,".RData")))
