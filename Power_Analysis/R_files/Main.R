library(ggplot2)
library(fitdistrplus)
library(truncnorm)
library(BiocManager)
if (FALSE) BiocManager::install(c("DESeq2", "edgeR"))
library(nlme)
library(truncdist)
library(DESeq2)
library(tidyverse)
library(bbmle)
library(patchwork)
library(truncdist)
library(plyr)
library(metR)  
library(mgcv)
library(ggplot2)
library(patchwork)
library(extraDistr)
library(sn)
library(here)
setwd(here())
theme_set(theme_bw())

source("Power_Analysis/R_files/Functions.R")
data_list <- mget(load("Power_Analysis/Datasets/data.RData"))
metadata_list <- mget(load("Power_Analysis/Datasets/groups.RData"))

path <- "Power_Analysis/Datasets/"
plot_path <- "Power_Analysis/Figures/"
############################################################
subset <- c("PRJNA453621","PRJEB45948","PRJNA589343", "PRJNA687773")
names(subset) <- subset
data_list <- list(data_list[["PRJNA453621"]], data_list[["PRJEB45948"]],
                  data_list[["PRJNA589343"]], data_list[["PRJNA687773"]])

metadata_list <- list(metadata_list[["PRJNA453621"]],
                      metadata_list[["PRJEB45948"]],
                      metadata_list[["PRJNA589343"]],
                      metadata_list[["PRJNA687773"]])

names(data_list) <- names(subset)
names(metadata_list) <- names(subset)
name = names(subset)
######################################################################
DeseqRes <- mget(load(paste0(path,"DeseqResults.RData"))) #list of deseq results
deseq_list <- DeseqRes[[names(DeseqRes)]]    
stopifnot(names(deseq_list) == name)

Control_Treat <- mget(load(paste0(path,"Control_Treat_Data.RData")))  
Control_Treat <- Control_Treat[[names(Control_Treat)]]
stopifnot(names(Control_Treat) == name)

param.disp <- mget(load(paste0(path,"Dispersion_Parameters.RData")))
param.disp <- param.disp[[names(param.disp)]]
logcontrol <- read_data(Control_Treat, "LogControl")  

param.control <- mget(load(paste0(path,"Control_Parameters.RData"))) 
param.control <- param.control[[names(param.control)]]

effects_list <- mget(load(paste0(path,"DeseqFoldChanges.RData")))
effects_list <- effects_list[[names(effects_list)]]

#################################################################
n_otu <- 300; n_samples <-  20; n_sim = 100
param.effects <- mget(load(paste0(path,"Effect_Parameters.RData")))
param.effects <- param.effects[[names(param.effects)]]
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

#Generate heatmaps
r = Power_Heatmap(combined_data); plts1 = r[["full_plot"]]; plts2 = r[["reduced_plot2"]]

plot1 = (plts1[[1]]|plts1[[2]]) & theme(legend.position = "bottom")
plot1 + plot_layout(guides = "collect")

plot2 =(plts1[[3]]|plts1[[4]])  & theme(legend.position = "bottom")
plot2 + plot_layout(guides = "collect")





# plts=Plot_ContEffects(deseq_list,logcontrol, plot_path)
# m = plts[[1]]; n = plts[[2]]
# (m[[1]]|m[[2]])/(m[[3]]|m[[4]])
# 
# (n[[1]]|n[[2]])/(n[[3]]|n[[4]]) #& theme(legend.position = "bottom")
# #n + plot_layout(guides = "collect")
# 
# file_path <- "control_histograms_"
# v1 = Fit_control(logcontrol,plot_path,file_path, plotname)
# v1 = (v1[[1]]|v1[[2]])/(v1[[3]]|v1[[4]]) & theme(legend.position = "bottom")
# v1 + plot_layout(guides = "collect")
# 
# plotname  <- "effects_dist_"
# v = Fit_Effects(effects_list,logcontrol,plot_path,plotname)
# vv = (v[[1]]|v[[2]])/(v[[3]]|v[[4]]) & theme(legend.position = "bottom")
# vv + plot_layout(guides = "collect")
