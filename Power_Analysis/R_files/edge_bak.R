library(microbiome)  
library(edgeR)
library(tidyverse)
library("pasilla")
library(MASS)   
library(metR)   
library(mgcv)
library(distr)
library(fitdistrplus)
library(statmod) 
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions.R")
data_list <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))
name <- "edgeR"
################################################################################
path <- "Power_Project/Datasets/"
plot_path <- "Latex_Files/Proposal/Figs/power/"
################################################################################
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
data = cts; group  =  coldata$condition
y <- DGEList(counts=data,group=group) 
keep <- rowSums(y$counts) >= 10
y <- y[keep,];  data <- y$counts

y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y$samples$group <- relevel(y$samples$group, ref="untreated")

design <- model.matrix(~group) 
rownames(design) <- colnames(y)
y2 <- estimateDisp(y, design, robust=TRUE)

dat <- data.frame(countmean = log2(rowMeans(y$counts)), dispersion = y2$trended.dispersion)
mm <- nls(dispersion ~ a + b/exp(countmean), start =list(a=1, b= 1), 
          data=dat)

Dispersion_edgeR = data.frame(asymptDisp = coef(mm)[1],extraPois = coef(mm)[2] )

design1 <- model.matrix(~0+group, data=y$samples)
colnames(design1) <- levels(y$samples$group) 

fit <- glmQLFit(y2, design1)

qlf <- exactTest(y2)
FDR <- p.adjust(qlf$table$PValue, method="BH")
qlf$table$FDR <- FDR
result <- qlf$table
#######################################################################
control <- data[,4:7]; treatment <- data[,1:3]
cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
logcontrol = log2(cont_Means) ; dispers = dispersions(dds)
effects<- result$logFC

data <- data.frame(Control_abundance=logcontrol, Fold_change =effects)

pp<- ggplot(data, aes(Control_abundance, Fold_change)) +
  geom_point()  +
  geom_smooth() +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

pp

############################################################################
alpha = 0.1
pad = result$FDR
power <-(!is.na(pad) & pad < alpha)
power = pad
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
plot(fit_2d)

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
brkvec <- seq(0,1,0.1)
gg_2dimc <- (gg_2dim
             + geom_contour(data = pp,
                            aes(z=power),
                            breaks = brkvec)
             + geom_label_contour(data = pp, aes(z= power),
                                  breaks = brkvec)
)
gg_2dimc

library(truncdist)
library(bbmle)
dat <- data.frame(effects = effects,logcontrol = logcontrol)
mm<-mle2(effects ~ dcauchy(0, scale=exp(a+b*logcontrol + c*logcontrol^2)),
         start =list(a=0.1, b=0.1, c= 0.001), data=dat)

param.effect <-coef(mm)
save(param.effect, file = paste0(path, "Effect_Parameters_",name,".RData")) 


Cont_Parameters <- mget(load("Power_Project/Datasets/Control_Parameters_deseq.RData"))
param.cont <- Cont_Parameters[[names(Cont_Parameters)]]

mu= param.cont$mu; sigma = param.cont$sigma; lambda = param.cont$lambda
myMix <- UnivarMixingDistribution(Norm(mean=mu[1], sd =sigma[1]),
                                  Norm(mean=mu[2], sd =sigma[2]),
                                  mixCoeff=c(lambda[1],lambda[2]))

rmyMix <- r(myMix)
n <- length(effects)
control <- rmyMix(n)

################################################################################
library(truncdist)
library(bbmle)
dat <- data.frame(effects = effects,logcontrol = logcontrol)
mm<-mle2(effects ~ dcauchy(0, scale=exp(a+b*logcontrol + c*logcontrol^2)),
         start =list(a=0.1, b=0.1, c= 0.001), data=dat)

param.effect <-coef(mm)

scal = data.frame(x =  exp(param.effect[1]+param.effect[2]*control + 
                             param.effect[3]*control^2))

effect_ <-apply(scal,1,function(x){
  rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})

hist(effects, prob=T)
lines(density(effects),col="red", lwd=2)
lines(density(effect_),col="blue", lwd=2)
#save(param.effect, file = paste0(path, "Effect_Parameters_",name,".RData")) 



good <- c("PRJNA453621","PRJEB45948","PRJNA589343", "PRJNA687773")
names(good) <- good; name <- names(good) 
data_list <- list(data_list[["PRJNA453621"]], data_list[["PRJEB45948"]],
                  data_list[["PRJNA589343"]], data_list[["PRJNA687773"]])

metadata_list <- list(metadata_list[["PRJNA453621"]], metadata_list[["PRJEB45948"]],
                      metadata_list[["PRJNA589343"]], metadata_list[["PRJNA687773"]])

names(data_list) <- names(good)
names(metadata_list) <- names(good)
################################################################################

edgeR_res <-function(data_list,metadata_list, filter=100){
  
  Dispersion_edgeR = list()
  index = 1
  for(n in names(data_list)){
    n= 2
    
    abund <- log2(rowMeans(data)); effects = result$logFC
    dd <- data.frame(mean_abundance = abund, foldchange= effects)
    ggplot(dd, aes(mean_abundance, foldchange)) +
      geom_point()  + 
      geom_smooth() +
      #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
      ggtitle(paste0("Mean abundance against Fold Changes")) +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
    
    
    aa<- fitdistr(abund, "exponential")
    pr=coef(aa); 
    data.frame(rate = coef(aa))
    #################
    n = 1000
    xr <- rexp(n, rate =  pr[1])
    index= index + 1
    }
}


##############3



b<- fitdist(means, "cauchy")
pp=coef(b); pp
n = 1000
xx <- rtrunc(n, 'cauchy', a=min(means), b=max(means), location=pp[1], scale= pp[2] )

a<- fitdist(means, "lnorm")
p=coef(a);   


x=rlnorm(n, meanlog = p[1] , sdlog = p[2])

hist(means, probability = T, breaks = 30)
lines(density(means), col ="blue")
lines(density(x))
#lines(density(xx), col ="red")
#lines(density(xxx), col ="pink")
lines(density(xr), col ="red")


dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, df)
pt_ls <- function(q, df=1, mu=0, sigma=1)  pt((q - mu)/sigma, df)
qt_ls <- function(p, df=1, mu=0, sigma=1)  qt(p, df)*sigma + mu
rt_ls <- function(n, df=1, mu=0, sigma=1)  rt(n,df)*sigma + mu
fit.t<-fitdist(means, 't_ls', start =list(df=1,mu=mean(means),sigma=sd(means))) 
ppp= coef(fit.t)
xxx = rt(n,ppp[1])*ppp[3] + ppp[2]



fit1 <- glmQLFit(y2, design1) 
qlf1 <- glmQLFTest(fit1, contrast=c(-1,1))
v <- (qlf1$table)

#control abundance 

stopifnot(rownames(y$samples) == colnames(y$counts))
k <- which(y$samples$group == "NT")
control <- y$counts[, k]
logcontrol <-log2(rowMeans(control))
k= which(is.finite(logcontrol))
logcontrol = logcontrol[k]
foldchange = v$logFC[k]


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

counts = cts ; group = coldata$condition 

control <- cts[,4:7]; treatment <- cts[,1:3]
cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)

########### Datasets
y <- DGEList(counts=counts,group=group) #
apply(y$counts, 2, sum) # total OTU counts per sample 
keep <- rowSums(cpm(y)>100) >= 2 
y <- y[keep,]
y_full <- y
dim(y)

y$samples$lib.size <- colSums(y$counts) 
y$samples
y <- calcNormFactors(y)
y$samples$lib.size*y$samples$norm.factors

y1 <- estimateCommonDisp(y, verbose=T)
y1 <- estimateTagwiseDisp(y1)

#use a generalized linear model to estimate the dispersion 
design <- model.matrix(~group) 
rownames(design) <- colnames(y)

library(statmod) 
y2 <- estimateDisp(y, design, robust=TRUE) 
mean <-y2$AveLogCPM
dispersion <- (y2$tagwise.dispersion)
# plot(mean,dispersion)
# plotBCV(y2)
mean <- c(mean)
x <- order(mean); y = x[7]
x[7]
which(mean == 1845)
mean1 <- mean[-which.max(mean)]
mean = mean[-c(1:5)]
View(mean)
library(bbmle)
dat <- data.frame(dispersion = dispersion,mean = mean1)
mm <- nls(dispersion ~ a + b/mean, start =list(a=10, b= 100), 
         data=dat, algorithm = "port", lower = c(a = -Inf, b = 0),
         upper = c(a = Inf, b = Inf))
mm



mm<-mle2(dispersion ~ a + b/mean, start =list(a=0.12345, b = 0.54321), data=dat)


#dispersion = asymptDisp - extraPois / mean 
dispersion <- coef(mm)[1] - coef(mm)[2]/mean
plot(mean,dispersion)


et <- exactTest(y1,pair = c( "treated", "untreated" ))



data_sim<-function(param.effect,param.disp,param.cont,n_sim,n_samples,n_otu,name,filter=10){
  
  Sim_data = Sim_metadata = dispersions =  True_data = list()
  true_treatment = true_effect =  true_control = list()
  
  new_otu <- 2*n_otu
  mu= param.cont$mu; sigma = param.cont$sigma; lambda = param.cont$lambda
  for(i in 1:n_sim){
    
    #simulate control and effects 
    myMix <- UnivarMixingDistribution(Norm(mean=mu[1], sd =sigma[1]),
                                      Norm(mean=mu[2], sd =sigma[2]),
                                      mixCoeff=c(lambda[1],lambda[2]))
    
    rmyMix <- r(myMix)
    
    control <- rmyMix(new_otu)
    cont <- unlist(control)
    
    ## scale parameters values for cauchy 
    scal = data.frame(x =  exp(param.effect[1]+param.effect[2]*cont + 
                                 param.effect[3]*cont^2))
    
    effect <-apply(scal,1,function(x){
      rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
    
    ## treatment simulation
    eff <- unlist(effect)
    treatment <-   eff + cont; treat <- unlist(treatment)
    
    mean_abund <- data.frame(treat = 2^treat, cont = 2^cont)
    mean_abund.data <- data.frame(means = rowMeans(mean_abund))
    
    dispers <- function(mean,c0,c1){c0 + c1/mean}
    disp <- drop(1/apply(mean_abund.data, 2, dispers,  c0 = param.disp[["asymptDisp"]],
                         c1 = param.disp[["extraPois"]]))
    
    
    dispersions <- disp
    df <- data.frame(mean_abund.data,disp)
    
    ## Simulate data
    data <-matrix(apply(df,1,function(x,n=n_samples){rnbinom(n=n,mu=x[1],size=x[2])}),
                  nrow = new_otu,ncol = n_samples)
    
    ## Filter data 
    keep <- rowSums(data) >= filter
    data <- data[keep,]
    effect <- effect[keep]; control <- control[keep]
    treatment <- treatment[keep]; dispersions <- disp[keep]
    
    #Select n_otus after filtering
    select <- sample.int(n_otu)
    data <- data[select,]
    true_effect[[i]] <- effect[select]; true_control[[i]] <- control[select]
    true_treatment[[i]] <- treatment[select]; dispersions <- disp[select] 
    
    
    Sim_metadata[[i]] <- data.frame(Groups = sort(rep(c("ASD","NT"),each = n_samples/2)),
                                    row.names=paste0("Sample",1:n_samples))
    
    rownames(data) <- paste0("otu",1:n_otu)
    colnames(data) <- paste0("Sample",1:n_samples)
    
    Sim_data[[i]] <- data
  } 
  names(true_effect) = names(true_control) = names(true_treatment) =  paste0("Sim_data",1:n_sim)
  names(dispersions) =  names(Sim_metadata) = names(Sim_data) = paste0("Sim_data",1:n_sim)
  
  save(dispersions, file = paste0(path,"Dispersions_", name, ".RData")) 
  save(true_control, file = paste0(path,"True_Control_", name, ".RData")) 
  save(true_effect, file = paste0(path,"True_Effect_", name, ".RData")) 
  save(true_treatment, file = paste0(path,"True_Treatment_", name, ".RData")) 
  
  save(Sim_metadata, file = paste0(path,"Simualated_Metadata_", name, ".RData")) 
  save(Sim_data, file = paste0(path,"Simualated_Data_", name, ".RData")) 
}


Deseq <-function(data_list,metadata_list,name,path){
  result = dispersion_coef = list()
  index = 1
  
  for(n in names(data_list)){
    data = data_list[[n]]; metadata <- metadata_list[[n]]
    dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    
    r <- which(colSums(data) == 0)
    if (length(r) > 0){
      data <- data[,-r]; metadata <- metadata[-r,]
      dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
    }
    dds$Groups <- relevel(dds$Groups, ref = "NT")
    dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=7) 
    res <- results(dds, cooksCutoff=TRUE, independentFiltering=TRUE)
    reslt <- lfcShrink(dds, res=res, coef=2, type="normal")  
    result[[index]] <- data.frame(reslt)
    index = index + 1
  }
  names(result) <- paste0("result",300:499)
  save(result, file = paste0(path, "Sim_Deseq_Results2_",name,".RData")) 
} 
