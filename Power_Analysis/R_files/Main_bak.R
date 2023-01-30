library(fitdistrplus)
## for bioconductor packages:
library(BiocManager)
if (FALSE) BiocManager::install(c("DESeq2", "edgeR"))
library(DESeq2)
library(tidyverse)
library(truncdist)
library(plyr)
library(metR)  
library(sn)
library(mgcv)
library(ggplot2)
library(bbmle)
library(patchwork)
library(extraDistr)
library(truncnorm)
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions.R")
data_list <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))

# p= DIMS(data_list)
# n =data.frame(p, row.names = c("samples", "taxa"))
# View(t(n)) # Is this un
################################################################################
path <- "Power_Project/Datasets/"
plot_path <- "Latex_Files/Proposal/Figs/power/"
################################################################################
name <- names(metadata_list); stopifnot(names(data_list) == name)
################################################################################
#subset
good <- c("PRJNA453621","PRJEB45948","PRJNA589343", "PRJNA687773")
names(good) <- good
data_list <- list(data_list[["PRJNA453621"]], data_list[["PRJEB45948"]],
                  data_list[["PRJNA589343"]], data_list[["PRJNA687773"]])

metadata_list <- list(metadata_list[["PRJNA453621"]], metadata_list[["PRJEB45948"]],
                      metadata_list[["PRJNA589343"]], metadata_list[["PRJNA687773"]])

names(data_list) <- names(good)
names(metadata_list) <- names(good)
name = names(good)
################################################################################
plotname <- "box_plot_effects_"
#v <- deseq_res(data_list, metadata_list,plot_path, plotname, filter=10)

DeseqRes <- mget(load(paste0(path,"DeseqResults.RData"))) #list of deseq results
deseq_list <- DeseqRes[[names(DeseqRes)]]   #extra Deseq data from the list
stopifnot(names(deseq_list) == name)

Control_Treat <- mget(load(paste0(path,"Control_Treat_Data.RData"))) #list of deseq results
Control_Treat <- Control_Treat[[names(Control_Treat)]]
stopifnot(names(Control_Treat) == name)

param.disp <- mget(load(paste0(path,"Dispersion_Parameters.RData"))) 
param.disp <- param.disp[[names(param.disp)]]

logcontrol_list <- read_data(Control_Treat, "LogControl") #extract control mean datasets from the list
logtreatment_list <- read_data(Control_Treat, "LogTreat") #extract control mean datasets from the list

##############################################
##############################################
##Now simulate and run the anlysis
#####################################
file_path <- "control_histograms_"
Fit_control(logcontrol_list,plot_path,file_path, plotname)

param.control <- mget(load(paste0(path,"Control_Parameters.RData"))) #lis lts
param.control <- param.control[[names(param.control)]]

effects_list <- mget(load(paste0(path,"DeseqFoldChanges.RData"))) 
effects_list <- effects_list[[names(effects_list)]]

plotname  <- "effects_dist_"
Fit_Effects(effects_list,logcontrol_list,plot_path,plotname)
param.effects <- mget(load(paste0(path,"Effect_Parameters.RData"))) 
param.effects <- param.effects[[names(param.effects)]]

#################################################################
n_otu <- 300; n_samples <-  20; n_sim = 100
#data_sim(param.effects,param.disp,param.cont=param.control,n_sim,n_samples,n_otu,path,filter=10)
sim_data_list <- mget(load(paste0(path,"Simulated_Data.RData")))
sim_metadata_list <- mget(load(paste0(path,"Simulated_Metadata.RData"))) 
true_control_list <- mget(load(paste0(path,"True_Control.RData"))) 
true_effect_list <- mget(load(paste0(path,"True_Effect.RData"))) 
standard_err_delta_list <- mget(load(paste0(path,"Standard_Err_Delta.RData"))) 

standard_err_delta_list = standard_err_delta_list[[names(standard_err_delta_list)]]
sim_data_list = sim_data_list[[names(sim_data_list)]]
sim_metadata_list = sim_metadata_list[[names(sim_metadata_list)]]
true_control_list = true_control_list[[names(true_control_list)]]
true_effect_list = true_effect_list[[names(true_effect_list)]]
###############################################################
#Deseq1(sim_data_list,sim_metadata_list,path)
deseq_data_list <- mget(load(paste0(path,"Sim_Deseq_RESULT.RData")))
deseq_data_list <- deseq_data_list[[names(deseq_data_list)]]

#########################################################
combined_data =combined_data_list(true_control_list, true_effect_list,
                                 deseq_data_list,standard_err_delta_list, alpha=0.1)
r = Power_Heatmap(combined_data, standardised = TRUE)
plts1 = r[["full_plot"]]; plts2 = r[["reduced_plot2"]] 

###########Show plots
b = ((plts1[[1]]|plts1[[2]])/(plts1[[3]]|plts1[[4]])) & theme(legend.position = "bottom")
b + plot_layout(guides = "collect")

b = (plts1[[3]]|plts1[[4]]) & theme(legend.position = "bottom")
b + plot_layout(guides = "collect")

b2 = (plts2[[1]]|plts2[[2]])/(plts2[[3]]|plts2[[4]]) & theme(legend.position = "bottom")
b2 + plot_layout(guides = "collect")

pr =contour_plots(combined_data,
                  interval=c(0.001,0.0005,0.001,0.001))

pr =contour_plots(combined_data, standardised = TRUE,
                  interval=c(0.006,0.0005,0.005,0.005))

(pr[[1]][["contour_plot"]]|pr[[2]][["contour_plot"]])
(pr[[3]][["contour_plot"]]|pr[[4]][["contour_plot"]]) 
###########################################################################
plts2[[1]]|pr[["PRJNA453621"]][["contour_plot"]]
plts2[[2]]|pr[["PRJEB45948"]][["contour_plot"]]
plts2[[3]]|pr[["PRJNA589343"]][["contour_plot"]]
plts2[[4]]|pr[["PRJNA687773"]][["contour_plot"]]

pr = contour_plots(combined_data, standardised = TRUE)
r=(pr[[1]][["contour_plot"]]|pr[[2]][["contour_plot"]])/(pr[[3]][["contour_plot"]]|pr[[4]][["contour_plot"]])  & theme(legend.position = "bottom")
r + plot_layout(guides = "collect")

#On avearge, there is 30% power for detecting
#standard_lfc = true_effect/sd_error
################################################################################
combined <- tibble(lcontrol = true_control, abs_lfc = abs(true_effect), 
               power = as.numeric(power),pval, count = rep(1,length(true_control)))




p|gg_2dimc

lf1 = true_effect/deseq_err
comb1 <- tibble(abs_lfc = abs(lf1),lcontrol = true_control, power, pval)
n1 = plts(comb1,0.005)

comb2 <- tibble(abs_lfc = abs(standard_lfc),lcontrol = true_control, power, pval)
n2 = plts(comb2,0.01)
################################################################################
deseq_standard_lfc = lfc/deseq_err
deseq_standard_lfc = unlist(deseq_standard_lfc)
comb3 <- tibble(abs_lfc = abs(deseq_standard_lfc),
                lcontrol = true_control, power, pval)

n3 = plts(comb3,0.1)

b = list(n1[["contour_plot"]],n2[["contour_plot"]],n3[["contour_plot"]])

combined <- b[[1]] + b[[2]] + b[[3]] & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")






#' Add the standard error to and write up
# do 100 for each.. 
###################################################################
# pp <-effects_list[[4]]
# 
# 
# hist(p, breaks = 30, probability = TRUE)
# lines(density(p), col="blue",lwd = 2)
# 
# logcont = logcontrol[[4]]
# dd = data.frame(logcont)
# p=apply(dd,1,function(x){ifelse(is.infinite(x), log2(runif(1,0.001,1)), x)})
# 
# k <- which(is.infinite(logcont)) 
# logcon <- logcont[k]; effect <- p[k]

# dat <- data.frame(logcont= logcont ,effect=effect)
# mm<-mle2(effect ~ dcauchy(0, scale=exp(a+b*logcont + c*logcont^2)),
#          start =list(a=1, b=-1, c= 0.01), data=dat)
# length(k)/length(p)
# 
# n = length(k)
# log2(runif(n,0,1))
# 
hist(effect, breaks = 30, probability = TRUE)
lines(density(effect), col="blue",lwd = 2)

Effect.Distribution <- function(effects_list,param.effects,param.control){
  
  for(n in names(param.effects)){
    param.effect = param.effects[[n]]
    control.param = param.control[[n]]
    effect_size   = effects_list[[n]]
    
    control = rnorm(n=length(effect_size),mean=control.param[,1], sd = control.param[,2])
    scal = data.frame(x =  exp(param.effect[,1]+param.effect[,2]*control + 
                             param.effect[,3]*control^2))
    effect <-apply(scal,1,function(x){
      rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
    
    file_name = paste0(plot_path,"effect_dist_",n,".png")
    reso <- 80
    length <- 3.25*reso/72
    png(file_name,res=reso)
    #png(file_name, res = 900)
    hist(effect_size, breaks = 30, probability = TRUE)
    lines(density(effect_size), col="blue",lwd = 2)
    lines(density(effect), col="red",lwd = 2)
    m = max(density(effect_size)$y); n = max(effect_size)
    legend(n-1.2, m- 0.5, legend=c("Observed effect", "Simulated effect"),
           col=c("blue", "red"), lty=1:2, cex=0.8,bty = "n")
    dev.off()
    
  }

}

Effect.Distribution(effects_list,param.effects,param.control)

View(true_effect)
############################################################################





############################################################################
v<-Plot_ContEffects(deseq_list,logcontrol, plot_path)

#v<-Plot_ContEffects(deseq_list,logcontrol,logtreatment, plot_path)
  
View(v$control_treat[[4]])
vv <- v$control_treat
k<-Transform(vv)

# Imputed_data <- k$imputed_data; Transform_data <- k$transform_data
# plotname <- "control_effects_impute_"
# Impute_plot(Imputed_data,plot_path,plotname)
#function(des)
#R <- control_treat_split(data_lis,metadata_list)

# control_treat_split<- function(data_list,metadata_list){
#   
#   L <- list(); index <- 1
#   
#   for (n in names(data_list)){ 
#     
#     ##Check data is taxa  by samples
#     m<-  rownames(data)[1]; p <- colnames(data)[1]
#     if (nchar(m) < nchar(p) ) {
#       data <- data
#     } else {
#       data <- t(data)
#     }
#     
#     data <- data_list[[n]]; metadata<- metadata_list[[n]]
#     
#     k <- which(metadata$Groups == "ASD")
#     control <- data[,-k]; treatment <- data[,k] 
#     stopifnot(colnames(control) == metadata$Samples[-k])
#     stopifnot(colnames(treatment) == metadata$Samples[k])
#     
#     cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
#     logcontrol = log2(cont_Means);   logtreat = log2(treat_Means)
#     
#     L[[index]] <- list(Control_Data = control, Control_Means = cont_Means,
#                        LogControl  =  logcontrol, Treat_Data =  treatment,
#                        Treat_Means  = treat_Means,
#                        LogTreat  = logtreat)
#     index <- index + 1
#   } 
#   
#   names(L) <- names(data_list) 
#   return(L)
#   
# }  
# 
# k <- which(metadata$Groups == "ASD")
# control <- data1[,-k]; treatment <- data1[,k] 
# stopifnot(colnames(control) == metadata$Samples[-k])
# stopifnot(colnames(treatment) == metadata$Samples[k])
# 
# cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
# Control_abundance = log2(cont_Means); Fold_change =result$log2FoldChange
# 
# k1 <-(which(is.finite(log2(cont_Means))))
# length(k1)
# Control_Abundance <- Control_abundance[k1]; Fold_Change = Fold_change[k1]
# 




##########Create threshold file for filtering###########
nam <- c("PRJNA168470", "PRJNA355023", "PRJNA453621", "PRJNA578223", 
          "PRJEB45948","PRJNA644763", "PRJNA589343", "PRJNA624252",
          "PRJNA687773","PRJNA642975")
stopifnot(names(data_lis) == nam)
# No need to filter yet. 


#########################################################
#thresh <- read.csv(paste0(path,"thresholds.csv")) ## data of values used for filtering
M <- read_datasets(data_lis,metadata_list) ##'@read_datasets: reads dataset as list 
RDataFiles(M,path) ##'@RDataFiles: creates filtered data, control, treatment and metadata files

###################################################################
#######Read datasets created by RDataFiles function################
data_list <- mget(load( paste0(path,"data_filtered.RData"))) # list of filtered dataset
stopifnot(names(data_list[["data_list"]]) == nam)

metadata_list <- mget(load(paste0(path,"metadata.RData"))) #  list of metadata
stopifnot(names(metadata_list[[names(metadata_list)]]) == nam)

control_list <- mget(load(paste0(path,"control.RData")))  # list of control dataset
stopifnot(names(control_list[[names(control_list)]]) == nam)

treat_list <- mget(load(paste0(path,"treat.RData") ))     #list of treatment dataset
stopifnot(names(treat_list[[names(treat_list)]]) == nam)

##############################################
control <- control_list[["cont_list"]]
treat <- treat_list[["treat_list"]]
data_list<- data_list[["data_list"]]
metadata_list <- metadata_list[["meta_list"]]

##Just check how much datasets is been thrown away from filtering
before_filter <- DIMS(data_lis); after_filter <- DIMS(data_list)

before_filter <- data.frame(before_filter); rownames(before_filter) <- c("samples","taxa")
after_filter <- data.frame(after_filter) ; rownames(after_filter) <- c("samples","taxa")

before_filter
after_filter
######################################################################
cont_treat_list<-read_ContTreat(control,treat) ##'@read_ContTreat: creates a list of control and treatment data
LogMeanData(cont_treat_list,path) ##' @LogMeanData: computes the log means of control and treatment


source("Power_Project/R_files/Functions.R")
plotname <-"outlier_detect"
deseq_res(data_list, metadata_list,plot_path, plotname) ##' @deseq_res: computes foldchanges, Pvalues, etc, from deseq2 
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

v<-Plot_ContEffects(deseq_list,logcontrol,logtreatment,plot_path) ##' @Plot_ContEffects: plots and saves plots for control verse effect sizes
###################################################################
vv <- v$control_treat; r <-Transform(vv)
Imputed_data <- r$imputed_data; Transform_data <- r$transform_data
###################################################################

des<- deseq_list[[1]]
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")

plot_path
#Imputed_data <- r$imputed_data; Transform_data <- r$transform_data
plotname <- "control_effects_impute_"
k <-Impute_plot(Imputed_data,plot_path,plotname)

plotname <- "control_effects_transform_"
difference_plot(Transform_data,plot_path,plotname)

###########Fit normal distribution to control abundances
effects_list<-v$effects

HIS_Effects <- function(effects_list,plot_path, plotname){
  
  for(n in names(effects_list)){
    dd <- effects_list[[n]]
    
    file_name = paste0(plot_path, plotname, n , ".pdf", sep="")
    pdf(file_name, width=6, height=4)
    hist(dd, breaks=30, prob = TRUE, main='Effect size',cex.main=1)
    lines(density(dd), col="blue", lwd=2) # add a density estimate with defaults
    dev.off()
    
  }
}

plotname <- "HIS_Effect_sizes_"
HIS_Effects(effects_list,plot_path, plotname)

plotname <- "HIS_Effect_sizes_imputed_"
HIS_Effects(k,plot_path, plotname)


file_path <- path
plotname <- "His_QQ_control_"
Hist_QQ(logcontrol,plot_path,file_path, plotname) 

# Fit the 2 and  and see. with limited control and without limited control 
plotname <- "Effect_Simulated_Effect_"
HIST_EFFECT_SIMULATED(deseq_list,logcontrol,file_path, plotname)


logcont<-logcontrol[[n]]; effect_size <-deseq_list[[n]]$log2FoldChange  
a<- which(is.finite(logcont))
logcont <- logcont[a]; effect_size2 <- effect_size[a]



control_param <- mget(load(paste0(path,"Parameters.RData")))#list of parameters for control
control_parameters <- control_param[[names(control_param)]]

effect_param <- mget(load(paste0(path,"EffectSize_Parameters.RData")))#list of parameters for control
effect_parameters <- effect_param[[names(effect_param)]]

n <- 1000

DATA_control <- rnorm(n,mean = control_parameters[[1]]$mean, sd=control_parameters[[1]]$sd)
hist(DATA_control,prob=TRUE, breaks = 30)
lines(density(DATA_control),col="blue", lwd=2)##' Then use it to simulate effect sizes 


plot_path
plotname <- "Effect-Control_Simulated_"
Control_Effect_Simulate(control_parameters, effect_parameters, nsim=n,plot_path, plotname)

to_long <- function(x) {
  pivot_longer(x, everything()) %>% 
    mutate(name = forcats::fct_inorder(name))
}


Box_Plot_Effects  <- function(deseq_list,plot_path, plotname){
  
  df <- list(); index <- 1
  mid <- c()
  for(n in names(deseq_list)){
    Effect_size = abs(deseq_list[[n]]$log2FoldChange)
    df[[index]] <- data.frame(Effect_size)
    mid[index] <- median(Effect_size)
    index <- index + 1
  }
  names(df) <-names(deseq_list)

  df <- lapply(df, to_long)
  df <- bind_rows(df, .id = "id")
  
  file_name = paste0(plot_path, plotname, ".pdf", sep="")
  pdf(file_name)  
  plt <- ggplot(df, aes(name, value)) +
    geom_boxplot() +
    facet_wrap(~id, scales = "free_x")
  dev.off()
  
  #print(mid) 
  plot(mid)
  lines(mid,col = "gray")
  
  print(plt)
  return(df)
}  

plotname <- "box_plot_effects_"
n <-Box_Plot_Effects(deseq_list,plot_path, plotname)
View(n)
n$id <- as.factor(n$id)
unique(n$id)

m <- lapply(logcontrol, to_long)

class(logcontrol[[1]])


plot(m)
lines(m,col = "gray")

summary(n)



#I will use simulation to verify the typical effect size I get I get from
#the actual data 
# We stack up all the effect sizes and  with the minimum length and then do the 
# analysis on it on we simulate effect sizes of the same length and then 
# do the analysis on it 
# What analysis are we going to do? Study their range and the median 
# 2. model the effect sizes with a mixed effect model where control abundance is a 
# mixed effect, maybe? 
# 

# Model the overall effect size against control and treatement abundances 


n <- 1; filter <- 10
data  = data_list[[n]]; metadata = metadata_list[[n]]
##Deseq processes
dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)

keep <- rowSums(counts(dds)) >= filter
dds <- dds[keep,]; data <- data[keep,] 

#' Remove samples with all zero (deseq needs this to run effectively) 
r <- which(colSums(data) == 0)
if (length(r) > 0){
  data <- data[,-r]; metadata <- metadata[-r,]
  dds <- DESeqDataSetFromMatrix(data, metadata, ~Groups)
}
dds$Groups <- relevel(dds$Groups, ref = "NT")

ddsDE <- DESeq(dds, sfType ="poscounts") 
res <- results(ddsDE, alpha = 0.05,
               cooksCutoff=TRUE, independentFiltering=TRUE) 

reslt <- lfcShrink(ddsDE, res = res, coef="Groups_ASD_vs_NT", type="normal")
#file_name = paste0(plot_path, "deseqMAplot_", n , ".pdf", sep="")
#pdf(file_name)
DESeq2::plotMA(reslt, main="normal")
########################
resu <- data.frame(reslt)

dd <- data.frame(control = log2(resu$baseMean), foldchange =resu$log2FoldChange)
ggplot(dd, aes(control, foldchange)) +
  geom_point()  + 
  geom_smooth() +
  #ggtitle(paste0("Control against Fold Changes for"," " , n)) +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))


#######################
compare_models <- function(nIndiv=10, nTaxa=30,Beta= c(1, 0.2, 0.1, 0.2), 
                           nTime = 10,
                           nSim = 2){
  
  Theta = Theta_fun(nTaxa)
  Ncol = length(Beta)
  Other_Effect_Est <- matrix(NA, nrow = nTaxa, ncol = Ncol)
  other_err =  glmmtmb_err =  matrix(NA,  nrow = nSim,ncol = Ncol)
  
  for(i in 1:nSim){
    TrueEffects <- matrix(rep(Beta, each = nTaxa), 
                          nrow = nTaxa,ncol = Ncol)
    
    metadata = metadat(nIndiv,nTime,nTaxa)
    dd = sim(dispersion=1,metadata,Beta,Theta,nIndiv)
    fit_glmmtmb_rr = glmmTMB(count ~ group*time+ rr(-1 + taxa|subject:time,2),
                             data = dd, family = nbinom2)
    glmmtmb_err[i,]  = (Beta - unlist(fixef(fit_glmmtmb_rr)$cond))^2
    
    for(n in 1:nTaxa){
      Taxa <- dd[dd$taxa == n,]
      fit_nbmm<-glmm.nb(count~group + time + group*time,
                        random = ~ 1 | subject/time, data = Taxa,
                        verbose = FALSE #method = "ML", control =list(msMaxIter = 1000, msMaxEval = 1000
                        )
      
      Other_Effect_Est[n,] = fixef(fit_nbmm)
    }
    
    otherErr = (TrueEffects - Other_Effect_Est)^2
    other_err[i,] = colMeans(otherErr)
  }
  
  other_err_dat = data.frame(other_err)
  glmmtmb_err_dat = data.frame(glmmtmb_err) 
  cols = c("intercept", "group", "time", "group_time")
  colnames(other_err_dat) = colnames(glmmtmb_err_dat) = cols
  
  v = list(other_err_data = other_err_dat, glmmtmb_err_data = glmmtmb_err_dat)
  
  nsim = 1:nSim
  inter <- (v
            %>% setNames(c("other", "glmmtmb"))
            %>% purrr::map_dfr(pull, intercept) 
  )
  inter$nsim  = nsim
  grp <- (v
          %>% setNames(c("other", "glmmtmb"))
          %>% purrr::map_dfr(pull, group) 
  )
  
  grp$nsim  = nsim
  tme <- (v
          %>% setNames(c("other", "glmmtmb"))
          %>% purrr::map_dfr(pull, time) 
  )
  tme$nsim  = nsim  
  
  grptime <- (v
              %>% setNames(c("other", "glmmtmb"))
              %>% purrr::map_dfr(pull, group_time) 
  )
  grptime$nsim  = nsim  
  v = list(Intercept = inter, Group  = grp, Time = tme, Group_Time = grptime)
}

v = compare_models(nIndiv=10, nTaxa=10, nTime = 10, Beta=  c(1, 0.2, 0.1, 0.2),  nSim = 100)
## Save these datasets


