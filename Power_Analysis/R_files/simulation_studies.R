library(BiocManager)
if (FALSE) BiocManager::install(c("DESeq2", "edgeR"))
library(nlme)
library(truncdist)
library(DESeq2)
library(patchwork)
library(MASS)
library(tidyverse)
library(sn)
#####
set.seed(101)
n_otu = 1000; samples_per_group = 20; dispersion = runif(n_otu,0,1)
scale_param = 3; 
logcontrol = rsn(n=n_otu, xi=4, omega=scale_param)
control  =  2^logcontrol
MyLogFoldChange = rtrunc(n = n_otu, spec = 'cauchy', a = -5, b = 5,
                       location = 0, scale = 0.5)
FoldChange = 2^MyLogFoldChange
treatment  = control*FoldChange

mean_abundances = (control + treatment)/2 #mean abundance across all samples

n = 2*samples_per_group
x= data.frame(mean_abundances,dispersion)
data = apply(x,1, \(s) rnegbin(n=n, mu = s[1], theta = s[2]) )

sample_names = paste0("sample",1:n)
groups = rep(c("control", "treatment"), each= samples_per_group)


metadata = data.frame(samples = sample_names,Groups = as.factor(groups))
colnames(data) =  paste0("otu",1:n_otu)
rownames(data) = sample_names

data  = t(data) # deseq wants samples by taxa

### No filtering nor shrinkage
dds <- DESeqDataSetFromMatrix(data, metadata,~Groups)
dds$Groups <- relevel(dds$Groups, ref = "control")
dds <- DESeq(dds,sfType ="poscounts")
res <- data.frame(results(dds))

## compare foldchanges 
DeseqLogFoldchange = res$log2FoldChange
plot(MyLogFoldChange,DeseqLogFoldchange)
abline(0,1)

library(ashr)
### Include filtering and shrinkage
dds2 <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=7) 
res2 <- results(dds2, cooksCutoff=TRUE, independentFiltering=TRUE)
reslt2 <- lfcShrink(dds2, res=res2, coef=2, type="normal")  
res2 = data.frame(reslt2)

## compare foldchanges 
DeseqLogFoldchange2 = res2$log2FoldChange
plot(MyLogFoldChange,DeseqLogFoldchange2)
abline(0,1)


### Include filtering only
dds3 <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
res3 <- results(dds3, cooksCutoff=TRUE, independentFiltering=TRUE)
res3 = data.frame(res3)

## compare foldchanges 
DeseqLogFoldchange3 = res3$log2FoldChange
plot(MyLogFoldChange,DeseqLogFoldchange3)
abline(0,1)

######################
dat = data.frame(logcontrol, Myfoldchange= MyLogFoldChange, Deseqfoldchange = DeseqLogFoldchange2 )
p1 <- ggplot(dat, aes(x=logcontrol, y=Myfoldchange)) +
  geom_point(alpha=0.5) +
  geom_smooth()

p2<-ggplot(dat, aes(x=logcontrol, y=Deseqfoldchange)) +
  geom_point(alpha=0.5) +
  geom_smooth()

p1|p2


library(edgeR)
y <- DGEList(counts=data,group=metadata$Groups)
y_full <- y
# keep the old one in case we mess up
head(y$counts)
apply(y$counts, 2, sum) # total OTU counts per sample
#keep <- rowSums(cpm(y)>100) >= 2
#y <- y[keep,]
#dim(y)
y1 <- estimateCommonDisp(y, verbose=T)
et <- exactTest(y1,pair = c( "treatment", "control" ))
y1 <- estimateCommonDisp(y, verbose=T)
y1 <- estimateTagwiseDisp(y1)
edgeRRes = et$table
View(edgeRRes)

plot(MyLogFoldChange, edgeRRes$logFC)
abline(0,1)

