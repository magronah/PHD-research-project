library(microbiome)  
library(DESeq2)
library(tidyverse)
library("pasilla")
library(metR)  
library(truncdist)
library(mgcv)
library(plyr)
library(patchwork)
library(distr)
library(fitdistrplus)
library(ggplot2)
source("Power_Project/R_files/Deseq_Functions_bak.R")
path <- "Power_Project/Datasets/"
plot_path <- "Latex_Files/Proposal/Figs/power/"
name <- "deseq" 
############################Preprocessing#################################
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
########################Datasets#################################
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ condition)
keep <- rowSums(counts(dds)) >= 50; dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "untreated")
dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
reslt <- lfcShrink(dds, res=res, coef=2, type="normal") 
v <-data.frame(reslt)
#################Control/Treatment split######################
cts <- counts(dds)
control <- cts[,4:7]; treatment <- cts[,1:3]
cont_Means <- rowMeans(control); treat_Means <- rowMeans(treatment)
logcontrol = log2(cont_Means); dispers = dispersions(dds)
logtreatment = log2(treat_Means)
################Fit Dispersions###############################
dat <- data.frame(logcontrol =logcontrol, dispersion = dispers)
mm <- nls(dispersion ~ a + b/(2^logcontrol), start =list(a=1, b= 1), 
          data=dat)
param.disp = data.frame(asymptDisp = coef(mm)[1], extraPois =  coef(mm)[2]) 
#save(param.disp, file = paste0(path, "Dispersion_Parameters_",name,".RData")) 
################Fit Control###############################
library(mixtools)
mixmdl = normalmixEM(logcontrol)
##########Plot
plotname <- "control_fit_"
file_name = paste0(plot_path, plotname, name , ".pdf", sep="")
pdf(file_name, width=6, height=5)
plot(mixmdl,which=2)
lines(density(logcontrol), lty=2, lwd=2)
legend(-2, 0.14, legend=c("fitted gaussian mixtures", " Observed control"),
       col = c("red","black"), lty=1:2, cex=0.8, box.lwd=0 )
dev.off()
########Save fitted values
lambda<-mixmdl$lambda; mu <-mixmdl$mu; sigma <- mixmdl$sigma
dat = data.frame(lambda  = mixmdl$lambda, mu = mixmdl$mu, sigma =mixmdl$sigma)
#save(dat, file = paste0(path, "Control_Parameters_",name,".RData")) 
################Fit effects###############################
library(truncdist)
library(bbmle)
dat <- data.frame(effects = v$log2FoldChange,logcontrol = logcontrol)
mm<-mle2(effects ~ dcauchy(0, scale=exp(a+b*logcontrol + c*logcontrol^2)),
         start =list(a=0.1, b=0.1, c= 0.001), data=dat)
param.effect <-coef(mm)
#save(param.effect, file = paste0(path, "Effect_Parameters_",name,".RData")) 
scal = data.frame(x =  exp(param.effect[1]+param.effect[2]*logcontrol + 
                             param.effect[3]*logcontrol^2))

effect_ <-apply(scal,1,function(x){
  rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})

########Plot
effects = v$log2FoldChange
plotname <- "observed_sim_effect_"
file_name = paste0(plot_path, plotname, name , ".pdf", sep="")
pdf(file_name, width=6, height=5)
hist(effects, prob=T)
lines(density(effects),col="red", lwd=2)
lines(density(effect_),col="blue", lwd=2)
legend(-4, 0.95, legend=c("Observed effect", "Simulated effect"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()
####################Reading data############################################
Disp_Parameters <- mget(load("Power_Project/Datasets/Dispersion_Parameters_deseq.RData"))
Cont_Parameters <- mget(load("Power_Project/Datasets/Control_Parameters_deseq.RData"))
Effect_Param <- mget(load(paste0(path,"Effect_Parameters_deseq.RData"))) 

param.cont <- Cont_Parameters[[names(Cont_Parameters)]]
param.disp <- Disp_Parameters[[names(Disp_Parameters)]]
param.effect <- Effect_Param[[names(Effect_Param)]]

##############################################################
n_otu <- 1000; n_samples <-  100; n_sim = 300
#data_sim(param.effect,param.disp,param.cont,n_sim,n_samples,n_otu,name,filter=10)

####################Reading data############################################
Simulated_Data <- mget(load("Power_Project/Datasets/Simulated_Data_deseq.RData"))
Simulated_Metadata <- mget(load("Power_Project/Datasets/Simulated_Metadata_deseq.RData"))
True_Control <- mget(load("Power_Project/Datasets/True_Control_deseq.RData"))
True_Treatment <- mget(load("Power_Project/Datasets/True_Treatment_deseq.RData"))
True_Effect <- mget(load("Power_Project/Datasets/True_Effect_deseq.RData"))
Standard_error <- mget(load("Power_Project/Datasets/Standard_error_deseq.RData"))
##############################################################
simulated_data_list <- Simulated_Data[[names(Simulated_Data)]]
simulated_metadata_list <- Simulated_Metadata[[names(Simulated_Metadata)]]

true_effect <- True_Effect[[names(True_Effect)]]
true_effect <- unlist(true_effect)

true_control <- True_Control[[names(True_Control)]]
true_control <- unlist(true_control)

true_treatment <- True_Treatment[[names(True_Treatment)]]
true_treatment <- unlist(true_treatment)

standard_error = Standard_error[[names(Standard_error)]]
vv =standard_error[[1]]; v = list()
for(i in 1:length(standard_error)){ 
  v[[i]] = vv[[1]]
}
names(v) = paste0("Sim_data",1:n_sim)
sd_error = unlist(v)

#plot(v[[1]],v[[200]]); abline(0,1)
###############################################################################
#Deseq(simulated_data_list,simulated_metadata_list,name,path)
sim_deseq_list <- mget(load(paste0(path,"Sim_Deseq_Results_deseq.RData")))
res_data <- sim_deseq_list[[names(sim_deseq_list)]]

##############################################################################
lfc <- (res_data
        %>% setNames(paste0("padjust", 1:n_sim))
        %>% purrr::map_dfr(pull, log2FoldChange) 
)

lfc = unlist(lfc)
standard_lfc = true_effect/sd_error

deseq_err <- (res_data
        %>% setNames(paste0("padjust", 1:n_sim))
        %>% purrr::map_dfr(pull, lfcSE) 
)
deseq_err =  unlist(deseq_err)
p_adjust <- (res_data
             ## 'map' won't work without names
             %>% setNames(paste0("padjust", 1:n_sim))
             ## map_dfr = run function on each element, combine the results into a data frame
             %>% purrr::map_dfr(pull, padj) ## 'pull' = tidyverse equivalent of $/[[
)

basemean <- (res_data
             ## 'map' won't work without names
             %>% setNames(paste0("padjust", 1:n_sim))
             ## map_dfr = run function on each element, combine the results into a data frame
             %>% purrr::map_dfr(pull, baseMean) ## 'pull' = tidyverse equivalent of $/[[
)
basemean =  unlist(basemean) 

pval<- unlist(p_adjust)
alpha = 0.1
power <-(!is.na(pval) & pval < alpha)

################################################################################
comb <- tibble(lcontrol = true_control, abs_lfc = abs(true_effect), 
               power = as.numeric(power),pval, count = rep(1,length(true_control)))

## GRAPH COLOR AND TEXT CONTROLS
showText<-TRUE
txtSize<-3
heatmap.low<-"lightgreen"
heatmap.high<-"orangered"

## SUMMARIZE TO NEW X/Y GRID
xblocks<-10
yblocks<-10

# so the best places are those with high total number of otus and high number of power.
# Can you color by both of them?

## CALL ddply to roll-up the data and calculate summary means, SDs,ec
dfe.plot<-ddply(comb,
                .(logControl=cut(comb$lcontrol,xblocks),
                  logFoldChange=cut(comb$abs_lfc,yblocks)),
                summarize,
                Power=mean(power),
                Sum = sum(count) )

## BUILD THE SUMMARY CHART
g<-ggplot(dfe.plot) +
  geom_raster(aes(logControl,logFoldChange,fill=Sum),alpha=0.75) +
  scale_fill_gradient(low=heatmap.low, high=heatmap.high) +
  theme_bw() + theme(axis.text.x=element_text(angle=-90)) +
  ggtitle(paste(xblocks,
                " X ",
                yblocks,
                " grid of Data\nbetween ( ",
                min(comb$lcontrol),
                " : ",
                min(comb$abs_lfc),
                " ) and ( ",
                max(comb$lcontrol),
                " : ",
                max(comb$abs_lfc),
                " )\n\n",
                sep=""))

g

if(showText)g<-g+geom_text(aes(logControl,logFoldChange,label=paste("SUM=",round(Sum,0),
                                                    "\nPower=",round(Power,1)),
                               fontface=c("italic")),
                           color="black",size=txtSize)
g

dfe.plot  = dfe.plot[dfe.plot$Sum >=100, ]
#On avearge, there is 30% power for detecting


############################################################################
################################################################################
plts <-function(comb,interval){ 
  gg_2dim <- (ggplot(comb)
            + aes(lcontrol, abs_lfc)
            + geom_point(aes(color = power), alpha = 0.5)
            + labs(x = "log(control abundance)",
                   y = "absolute(log(fold change))")
            )

  fit_2d <- bam(power ~ te(lcontrol, abs_lfc),
              data = comb,
              family = binomial) #
  
  pp <- with(comb,
           expand.grid(lcontrol = seq(min(lcontrol),
                                      max(lcontrol),
                                      length = 50),
                       abs_lfc = seq(min(abs_lfc),
                                     max(abs_lfc),
                                     length = 50)))

  pp$power <- predict(fit_2d, newdata = pp,
                    type = "response")
################################################################################
  brkvec <- seq(0,1,interval)
  gg_2dimc <- (gg_2dim
             + geom_contour(data = pp,
                            aes(z=power),
                            breaks = brkvec)
             + geom_label_contour(data = pp, aes(z= power),
                                  breaks = brkvec)
        )
  v = list(gg_2dim,pp,fit_2d,gg_2dimc)
  names(v) = c("point_plot","tibble","bam_predict_plot","contour_plot")
  v
}
lf1 = true_effect/deseq_err
comb1 <- tibble(abs_lfc = abs(lf1),lcontrol = true_control, power, pval)
n1 = plts(comb1,0.005)
n1
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

################################################################################
# plotname <- "power_"
# file_name = paste0(plot_path, plotname, name , ".pdf", sep="")
# pdf(file_name, width=6, height=5)
# print(gg_2dimc) 
# dev.off()

MM <- comb %>% plyr::mutate(xbin  = cut(lcontrol, breaks = 10), ybin = cut(abs_lfc, breaks=10))%>%
  group_by(xbin,ybin) %>% count(cut_x, cut_y)


ggplot(data=MM, aes(x=(xbin), y=(ybin), z = power)) +
  geom_tile(aes(fill=power))

ggplot(data=MM, aes(x=(xbin), y=(ybin), z = power)) +
  geom_contour()

ggplot(data=MM, aes(x=(xbin), y=(ybin), z = power)) +
  geom_contour(aes(group = 1))

summary(MM)
colnames(comb)
## power vs log fold change
gg_abslfc <- (ggplot(comb)
              + aes(abs_lfc, power)
              + geom_point()
              + geom_smooth()
)

gg_abslfc

gg_lcontrol <- (ggplot(comb)
                + aes(lcontrol, power)
                + geom_point()
                + geom_smooth()
)

gg_lcontrol

mean <- (res_data
         ## 'map' won't work without names
         %>% setNames(paste0("padjust", 1:n_sim))
         ## map_dfr = run function on each element, combine the results into a data frame
         %>% purrr::map_dfr(pull, baseMean) ## 'pull' = tidyverse equivalent of $/[[
)

basemean <- unlist(mean)
basemean <- log2(basemean)
View(comb)# plotname <- "number_of_OTUs_"
# file_name = paste0(plot_path, plotname, name , ".pdf", sep="")
# pdf(file_name, width=6, height=5)
hist(count, probability = T, breaks = 30)
lines(density(count))
# dev.off()

hist(comb$power, probability = T, breaks = 30)
lines(density(comb$power))

k <- which(count == 0 )
length(count)
count <- count[-k]

# plotname <- "number_of_OTUs_nozeros_"
# file_name = paste0(plot_path, plotname, name , ".pdf", sep="")
# pdf(file_name, width=6, height=5)
# hist(count, probability = T, breaks = 30)
# lines(density(count))
# dev.off()



###############################################################################
effects<- v$log2FoldChange
data <- data.frame(Control_abundance=logcontrol, Fold_change =effects)
pp<- ggplot(data, aes(Control_abundance, Fold_change)) +
  geom_point()  +
  geom_smooth() +
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
###############################################################################

library(akima)
 library(reshape2)
 d1 <- with(comb, interp(x = lcontrol, y = abs_lfc, z = power))
# d2 <- melt(d1$z, na.rm = TRUE)
# names(d2) <- c("x", "y", "power")
# d2$lcontrol <- d1$x[d2$x]
# d2$abs_lfc <- d1$y[d2$y]
# ggplot(data = d2, aes(x = lcontrol, y = abs_lfc, z = power)) +
#   geom_contour()
# 
# grid <- akima::interp(comb$lcontrol, comb$abs_lfc, comb$power)
# griddf <- data.frame(x = rep(grid$x, ncol(grid$z)), 
#                      y = rep(grid$y, each = nrow(grid$z)), 
#                      z = as.numeric(grid$z))


# interpolate data to regular grid
d1 <- with(all_merge, interp(x = NMDS1, y = NMDS2, z = shannon))

