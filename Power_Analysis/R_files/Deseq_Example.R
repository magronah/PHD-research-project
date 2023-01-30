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
name <- "deseq" 
###########################################################

###############################################################################
data_list <- mget(load(paste0(path,"Simualated_Data_",name,".RData")))
metadata_list <- mget(load(paste0(path,"Simualated_Metadata_",name,".RData")))
True_Control <- mget(load(paste0(path,"True_Control_",name,".RData")))
True_Effect <- mget(load(paste0(path,"True_Effect_",name,".RData")))
True_Treatment <- mget(load(paste0(path,"True_Treatment_",name,".RData")))


sim_deseq_list <- mget(load(paste0(path,"Sim_Deseq_Results_",name,".RData")))
res_data <- sim_deseq_list[[names(sim_deseq_list)]]
###############################################################################
lfc <- (res_data
        %>% setNames(paste0("padjust", 1:n_sim))
        %>% purrr::map_dfr(pull, log2FoldChange) 
)

p_error <- (res_data
            %>% setNames(paste0("padjust", 1:n_sim))
            %>% purrr::map_dfr(pull, lfcSE) 
)

lfc = unlist(lfc)/ unlist(p_error)

#store all the padjusted
p_adjust <- (res_data
             ## 'map' won't work without names
             %>% setNames(paste0("padjust", 1:n_sim))
             ## map_dfr = run function on each element, combine the results into a data frame
             %>% purrr::map_dfr(pull, padj) ## 'pull' = tidyverse equivalent of $/[[
)

a<- unlist(p_adjust); alpha = 0.1

power <-(!is.na(a) & a < alpha)

df <- function(lfc, power, control){

  comb <- tibble(abs_lfc = abs(lfc),lcontrol = control, power)
  fit_2d <- bam(power ~ te(lcontrol, abs_lfc),
                data = comb, family = binomial) 
  pp <- with(comb,
             expand.grid(lcontrol = seq(min(lcontrol),
                                        max(lcontrol),
                                        length = 50),
                         abs_lfc = seq(min(abs_lfc),
                                       max(abs_lfc),
                                       length = 50)))
  pp$power <- predict(fit_2d, newdata = pp,
                      type = "response")
  
  result = list(data = comb, predict = pp, model = fit_2d)
  
  return(result)
}

dat1 = df(lfc, power, true_control)
dat2 = df(true_effect, power, true_control)

summary(dat1$data)
summary(dat2$data)
par(mfrow=c(1,2))
plot(abs_lfc ~lcontrol, dat1$data, pch=".")
plot(abs_lfc ~lcontrol, dat2$data, pch=".")


plot(true_control, lfc)


plot()

par(mfrow=c(1,2))
plot(dat1$model)
plot(dat2$model)

summary(dat1$model)
summary(dat2$data)

summary(dat1$predict)
summary(dat2$predict)


