library(fitdistrplus)
## for bioconductor packages:
library(BiocManager)
if (FALSE) BiocManager::install(c("DESeq2", "edgeR"))
library(nlme)
library(truncdist)
library(DESeq2)
library(tidyverse)
library(patchwork)
library(bbmle)
library(sn)
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions.R")
data_list <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))

################################################################################
path <- "Power_Project/Datasets/"
plot_path <- "Latex_Files/Proposal/Figs/power/"
################################################################################
name <- names(metadata_list); stopifnot(names(data_list) == name)
View(data_list[[1]])
################################################################################
#subset
good <- c("PRJNA453621","PRJEB45948","PRJNA589343", "PRJNA687773")
names(good) <- good
name <- names(good)
data_list <- list(data_list[["PRJNA453621"]], data_list[["PRJEB45948"]],
                  data_list[["PRJNA589343"]], data_list[["PRJNA687773"]])

metadata_list <- list(metadata_list[["PRJNA453621"]], metadata_list[["PRJEB45948"]],
                      metadata_list[["PRJNA589343"]], metadata_list[["PRJNA687773"]])
#############################################################################

View(data_list[[1]])

library(vegan)

S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)

otu_tab = (data_list[[1]]) # sample by taxa
otu_tab = t(otu_tab)
otu_tab = as.data.frame(otu_tab)
str(otu_tab)
raremax1 <- min(rowSums(otu_tab))
otu_tab_rarefy <- rarefy(otu_tab, raremax1)
rare_data = tibble(Samples = names(otu_tab_rarefy), value = otu_tab_rarefy) %>%
  pivot_wider(names_from="Samples", values_from="value") %>%
  as.data.frame()

View(otu_tab)

View(rare_data)
class(otu_tab_rarefy)

rand_name = rownames(otu_tab)
rand_df <- otu_tab %>%
  pivot_long(names_from="rand_name1", values_from="n", values_fill = 0) %>%
  as.data.frame()

rownames(otu_tab_rarefy)
##############################################################################
Control_Treat <- mget(load(paste0(path,"Control_Treat_Data.RData"))) #list of deseq results
Control_Treat <- Control_Treat[[names(Control_Treat)]]
stopifnot(names(Control_Treat) == name)
##############################################################################
names(data_list) <- names(good)
names(metadata_list) <- names(good)
name <- names(metadata_list); stopifnot(names(data_list) == name)
################################################################################
plotname <- "box_plot_effects_"
#v <- deseq_res(data_list, metadata_list,plot_path, plotname, filter=10)
################################################################################
DeseqRes <- mget(load(paste0(path,"DeseqResults.RData"))) #list of deseq results
deseq_list <- DeseqRes[[names(DeseqRes)]]   #extra Deseq data from the list
stopifnot(names(deseq_list) == name)
View(deseq_list[[1]])

Control_Treat <- mget(load(paste0(path,"Control_Treat_Data.RData"))) #list of deseq results
Control_Treat <- Control_Treat[[names(Control_Treat)]]
stopifnot(names(Control_Treat) == name)

logcontrol <- read_data(Control_Treat, "LogControl") #extra control mean datasets from the list
logtreatment <- read_data(Control_Treat, "LogTreat") #extra control treatment datasets from the list
############################################################################
vv= Plot_ContEffects(deseq_list,logcontrol, plot_path)
(vv[[1]]|vv[[2]])/(vv[[3]]|vv[[4]])

 

#############################################################################
# DeseqRes1 <- mget(load(paste0(path,"DeseqResults10.RData"))) #list of deseq results
# DeseqRes2 <- mget(load(paste0(path,"DeseqResults30.RData"))) #list of deseq results
# DeseqRes3 <- mget(load(paste0(path,"DeseqResults50.RData"))) #list of deseq results
# deseq_list1 <- DeseqRes1[[names(DeseqRes1)]]   #extra Deseq data from the list
# deseq_list2 <- DeseqRes2[[names(DeseqRes2)]]   #extra Deseq data from the list
# deseq_list3 <- DeseqRes3[[names(DeseqRes3)]]   #extra Deseq data from the list
# 
# Control_Treat1 <- mget(load(paste0(path,"Control_Treat_Data10.RData"))) #list of deseq results
# Control_Treat2 <- mget(load(paste0(path,"Control_Treat_Data30.RData"))) #list of deseq results
# Control_Treat3 <- mget(load(paste0(path,"Control_Treat_Data50.RData"))) #list of deseq results
# 
# Control_Treat1 <- Control_Treat1[[names(Control_Treat1)]]
# Control_Treat2 <- Control_Treat2[[names(Control_Treat2)]]
# Control_Treat3 <- Control_Treat3[[names(Control_Treat3)]]
# 
# logcontrol1 <- read_data(Control_Treat1, "LogControl") #extra control mean datasets from the list
# logcontrol2 <- read_data(Control_Treat2, "LogControl") 
# logcontrol3 <- read_data(Control_Treat3, "LogControl") 

# vp <-Plot_ContEffects(deseq_list1,logcontrol1, plot_path)
# vp[[1]]|v[[2]]|
##############################################################################


# deseq_Filter_list <- list(filtThresh_10 = deseq_list1, filtThresh_30 = deseq_list2,
#                   filtThresh_50 =deseq_list3)
# 
# logcontrol_Filter_list <- list(filtThresh_10 = logcontrol1, filtThresh_30 = logcontrol2,
#                  filtThresh_50 =logcontrol3)
# 
# save(deseq_Filter_list, file = paste0(path, "deseq_Filter_list.RData")) 
# save(logcontrol_Filter_list, file = paste0(path, "logcontrol_Filter_list.RData")) 

############################################################################
effects <- list(data.frame(effect = deseq_list$PRJNA453621$log2FoldChange), 
                           data.frame(effect = deseq_list$PRJEB45948$log2FoldChange),
                           data.frame(effect = deseq_list$PRJNA589343$log2FoldChange),
                           data.frame(effect =deseq_list$PRJNA687773$log2FoldChange)) 

names(effects) <- names(deseq_list)
plotname  <- "effects_dist_"
HIST(effects,plot_path, plotname)

effect_list <- list(); logcontrol_list <- list()
index <- 1

for(n in names(effects)){
  effect <- deseq_list[[n]]$log2FoldChange; logcontrol1 <-logcontrol[[n]]
  k<-which(is.finite(logcontrol1))
  effect_list[[index]] <- effect[k]; logcontrol_list[[index]] <- logcontrol1[k]
  #data <- data.frame(logcont,effect ) 
  index <- index+1
}
names(effect_list)= names(logcontrol_list)= names(effects)
##########################################################################

Parameters<-Fit_Effects(effect_list,logcontrol_list)

x <- effects[[1]]$effect

ca <- fitdistr(x, "cauchy")
p<-rcauchy(100,-0.014676381,0.148778427)

#do it in a log scale

Plot_Effect_Fit<-function(effect_list,logcontrol_list,Parameters){
  scale_param <- list(); index<- 1
  for(n in names(effect_list)){
    n<- 1
  param <- Parameters[[n]]; logcont <- logcontrol_list[[n]]
  effect <- effect_list[[n]]

  pp<-(exp(param[1,]+param[2,]*logcont + param[3,]*logcont^2))
  #p<-rtrunc(pp, spec="cauchy", a=-5, b=5)
  l <- length(pp)
  p <- rcauchy(n = l, location = 0, scale = pp)
  scale_param[[index]] <- pp
  filename <- paste0(plot_path,"effect_fit_", n, ".pdf", sep="")
  pdf(filename)
  hist(p, breaks = 30, prob = TRUE, cex.main=1) 
  lines(density(p), col="blue", lwd=2) # add a density estimate with defaults
  dev.off()
  
  index <- index + 1
  }
  names(scale_param) = names(effect_list)
  return(scale_param)
}

v<-Plot_Effect_Fit(effect_list,logcontrol_list,Parameters)


###########################################################################
k<-which.max(v[[4]])
vv <- v[[4]][-k]
a <- data.frame('PRJNA453621' = v[[1]]) 
b <- data.frame('PRJEB45948' = v[[2]]) 
c <- data.frame('PRJNA589343' = v[[3]])
d <- data.frame('PRJNA687773' = vv)

k<-c(median(v[[1]]),median(v[[2]]),median(v[[3]]),median(v[[4]]))
plot(k)

filename <- paste0(plot_path,"box_plot_scale_para_", ".pdf", sep="")
pdf(filename)
boxplot(dplyr::bind_rows(a, b,c))
#print(p)
dev.off()

#simulate 
plotname<- "His_QQ_control_"
Hist_QQ(logcontrol,plot_path,path, plotname)
Control_Param <- mget(load(paste0(path,"Control_Parameters.RData")))  
Control_Param <- Control_Param[[names(Control_Param)]]

plotname <- "Effect-Control_Simulated_"
Control_Effect_Simulate(control_param, effect_param, nsim,plot_path, plotname)
  
param <-(Control_Param[[1]])

l <- length(v[[1]])
control <- rnorm(n=l, mean = param$mean, sd = param$sd)
effects <- rtrunc(v[[1]], spec="cauchy", a=-5, b=5) 
data <- data.frame(control,effects)
plts <-  ggplot(data, aes(control, effects)) +
  geom_point()  + 
  geom_smooth() +
  ggtitle(paste0("Control against Fold Changes")) +
  xlab("control abundance (log2 scale)") + ylab("Log2 fold changes") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))

plts
#####################################################################################
##' 1. mu = mean of control abundances 
##' disperson = from deseq
##' 

#####################################################################################
data <- effec[[4]]
c <- fitdist(data, "cauchy",method="qme")
which(data == 0)


control<- (Control_Treat$PRJNA453621$LogControl) 
treatment<- (Control_Treat$PRJNA453621$LogTreat) 
#plot_cnt_treat <- function(){
  
  k1 <-which(is.finite(control))
  control <- control[k1]; treatment <- treatment[k1]
  k2 <-which(is.finite(treatment))
  control <- control[k2]; treatment <- treatment[k2]
  
  plotname <- "cont_treat_plot5_"
  file_name = paste0(plot_path, plotname, "PRJNA453621" , ".pdf", sep="")
  
  pdf(file_name)
  plot(control,treatment)
  abline(0,1)
  dev.off()
  
#}
x <- seq( 0, 3, .1 )
pdf <- dtrunc(x, spec="cauchy", a=1, b=2 )
(pdf)

#######################################
ca$estimate$location
data <- effec[[1]]
fitdistr(data, "truncnorm", list(location=0,scale=1))
v<-fitdistr(data, "cauchy")

summary(ca)
coef(ca)$location
  #######T
a<-rtrunc(100, "cauchy", location=-0.014676381,scale=  0.148778427 ,a = -1, b = 2)

hist(a, breaks=30, prob = TRUE,
     main='Effect size',
     cex.main=1)
lines(density(a), col="blue", lwd=2) # 
library(truncdist)

x <- rtrunc( 500, spec="norm", a=1, b=2 )




# 
# 
# x = seq(0,2,0.02)
# x =1
# 
# unlist(res)
# 

