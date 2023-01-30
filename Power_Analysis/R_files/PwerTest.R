library(fitdistrplus)
library(DESeq2)
library(tidyverse)
library(MASS)
library(mgcv)
library(metR)  ## for contour labels (geom_label_contour)
library(purrr)
library(here)
setwd(here())
theme_set(theme_bw())
source("Power_Project/R_files/Functions.R")
path <- "Power_Project/Datasets/"

############################
Control_Treat <- mget(load(paste0(path,"Control_Treat_Data.RData"))) #list of deseq results
Disp_Parameters <- mget(load("Power_Project/Datasets/Dispersion_Parameters.RData"))
Cont_Parameters <- mget(load("Power_Project/Datasets/Control_Parameters.RData"))
DeseqRes <- mget(load(paste0(path,"DeseqResults.RData"))) #list of deseq results

################################################################################
Control_Treat <- Control_Treat[[names(Control_Treat)]]
param.cont <- Cont_Parameters[[names(Cont_Parameters)]]
param.disp <- Disp_Parameters[[names(Disp_Parameters)]]
deseq_list <- DeseqRes[[names(DeseqRes)]]   #extra Deseq data from the list

deseqres = deseq_list[[1]]
param.cont <- param.cont[[1]]
param.disp <- param.disp[[1]]

#Simulate log controls
n_otu <- sum((!is.na(deseqres$padj))); n_samples <-  40; n_sim = 50
#data_sim(param.disp,param.cont,n_sim,n_samples,n_otu)
  
p <- expand.grid(gender = c("male", "female"),
                 Groups = c("NT", "ASD"), rep = 1:floor(n_samples/2)) [1:n_samples,]

metadata <- data.frame(Samples = paste0("Sample",1:n_samples), Groups = p$Groups)

data_list <- mget(load("Power_Project/Datasets/Simualated_Data.RData"))
data_list <- data_list[[names(data_list)]]
metadata_list <-met(metadata)

#k<-Deseq1(data_list,metadata_list,path)
sim_deseq_list <- mget(load("Power_Project/Datasets/Sim_Deseq_Results.RData"))
res_data <- sim_deseq_list[[names(sim_deseq_list)]]
# data_list <- list(data_list[[1]],data_list[[2]],data_list[[3]],
#                   data_list[[4]],data_list[[5]])
# class(res_data)
# res_data <- list(res_data$result1,res_data$result2,res_data$result3,
#                  res_data$result4,res_data$result5)
#store all the padjusted
  p_adjust <- (res_data
               ## 'map' won't work without names
               %>% setNames(paste0("padjust", 1:n_sim))
               ## map_dfr = run function on each element, combine the results into a data frame
               %>% purrr::map_dfr(pull, padj) ## 'pull' = tidyverse equivalent of $/[[
  )
  
  ## observed log fold changes
  obs_lfcmat <- (res_data
                 %>% setNames(paste0("log2FoldChange",1:n_sim))
                 %>% purrr::map_dfr(pull,log2FoldChange)
  )
  obs_abslfc <- rowMeans(abs(obs_lfcmat))
  
  #stopifnot(!any(is.na(p_adjust)))

  # k<-which(rowSums(is.na(p_adjust))>0)
  # p_adjust <-drop_na(p_adjust)
  # obs_abslfc <- obs_abslfc[-k]

## OR
alpha = 0.1
power <- rowMeans(!is.na(p_adjust) & p_adjust < alpha)
#power <- rowMeans(!is.na(p_adjust) < alpha)
#length(power)
##
lfold1 <- deseq_list$PRJNA453621$log2FoldChange
lcont1 <-Control_Treat$PRJNA453621$LogControl


p<- which(!is.na(deseqres$padj) == TRUE) 
lfold2 <- lfold1[p]; lcont2 <- lcont1[p]
lfold3 <- lfold2[-k]; lcont3 <- lcont2[-k]

# remove infinit values in control
pp <-which(is.finite(lcont3)== TRUE)
lfoldchange <- lfold3[pp]; lcontrol <- lcont3[pp]
power <- power[pp]; obs_abslfc <- obs_abslfc[pp]; p_adjust <- p_adjust[pp,]

##Create a table 
comb <- tibble(abs_lfc = abs(lfoldchange), lcontrol, power,
               obs_abslfc)

gg_2dim <- (ggplot(comb)
            + aes(lcontrol, abs_lfc)
            + geom_point(aes(size = power), alpha = 0.5)
            + labs(x = "log(control abundance)",
                   y = "log(fold change)")
)

gg_2dim

fit_2d <- gam(power ~ te(lcontrol, abs_lfc),
              data = comb,
              family = quasibinomial) #
plot(fit_2d)

pp <- with(comb,
           expand.grid(lcontrol = seq(min(lcontrol),
                                      max(lcontrol),
                                      length = 51),
                       abs_lfc = seq(min(abs_lfc),
                                     max(abs_lfc),
                                     length = 51)))

pp$power <- predict(fit_2d, newdata = pp,
                    type = "response")

# gg0 <- ggplot(comb, aes(lcontrol, abs_lfc)) +
#   geom_point() +
#   theme_classic()
# gg0 + geom_bin_2d()

brkvec <- c(0.01, 0.02, 0.05, (1:9)/10)

gg_2dimc <- (gg_2dim
             + geom_contour(data = pp,
                            aes(z=power),
                            breaks = brkvec)
             + geom_label_contour(data = pp, aes(z= power),
                                  breaks = brkvec)
)
gg_2dimc 


## power vs log fold change
gg_abslfc <- (ggplot(comb)
              + aes(abs_lfc, power)
              + geom_point()
              + geom_smooth()
)

gg_abslfc + geom_bin_2d()

gg_lcontrol <- (ggplot(comb)
                + aes(lcontrol, power)
                + geom_point()
                + geom_smooth()
)
gg_lcontrol
#
###############################################################################
##' I will finish and send him the results today and send it




