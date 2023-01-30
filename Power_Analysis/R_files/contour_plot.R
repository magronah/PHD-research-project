library(DESeq2)
## library(HMP)
library(MASS)
library(mgcv)
library(metR)  ## for contour labels (geom_label_contour)
library(tidyverse) ## includes ggplot2
library(here)
setwd(here())
theme_set(theme_bw())

## set working directory to head of repo
source("R_files/Functions.R")


data_list <- mget(load("Power_Project/Datasets/data.RData"))
metadata_list <- mget(load("Power_Project/Datasets/groups.RData"))
################################################################################
data_list <- list(data_list[["PRJNA453621"]])
names(data_list) <- "PRJNA453621"
dim(data_list[["PRJNA453621"]])
###########################################
get_file <- function(f) {
  r <- read.csv(file.path(path, f), row.names = 1)
  if (ncol(r) == 1) return(r[[1]])
  r
}


#Get result file for  each otu
res_data <- list()
for (i in 1:n_sim) {
  res_data[[i]] <-  get_file(paste0("result",i,".csv"))
}

dims <- sapply(res_data, dim)
stopifnot(
  dims[1,] == n_otu,
  dims[2,] == 6   ## baseMean, l2fc, lfoldchangeSE, stat, pvalue, padj
)

#store all the padjusted
p_adjust <- matrix(NA, nrow = n_otu, ncol = n_sim)
for (i in 1:n_sim) {
  p_adjust[,i] <- (res_data[[i]])$padj
}
## or: sapply(res_data, "[[", "padj")
## 'sapply' = 'lapply + [s]implify'
## do.call("rbind", lapply(res_data, "[[", "padj"))
## OR: library(purrr)

#
tmp <- (res_data
        ## 'map' won't work without names
        %>% setNames(paste0("sim", 1:n_sim))
        ## map_dfr = run function on each element, combine the results into a data frame
        %>% purrr::map_dfr(pull, padj) ## 'pull' = tidyverse equivalent of $/[[
)

## observed log fold changes
obs_lfcmat <- matrix(0,nrow(res_data[[1]]), n_sim)
for (i in 1:n_sim){
  obs_lfcmat[,i] <- (res_data[[i]])$log2FoldChange
}
obs_abslfc <- rowMeans(abs(obs_lfcmat))
stopifnot(!any(is.na(p_adjust)))

## read Chapter 2 'Growing Objects'
## https://www.burns-stat.com/pages/Tutor/R_inferno.pdf
power <- rep(NA, n_otu)
for (i in 1:nrow(p_adjust)){
  power[i] <- sum(p_adjust[i,] < alpha)/ncol(p_adjust)
}
## OR
power <- rowMeans(p_adjust < alpha)

png(filename=file.path(path, "plots/orderpower.png"),width = 500, height = 480) 
plot(power[order(power)])
dev.off()

comb <- tibble(abs_lfc = abs(lfoldchange), lcontrol, power,
               obs_abslfc)

gg_2dim <- (ggplot(comb)
            + aes(lcontrol, abs_lfc)
            + geom_point(aes(size = power), alpha = 0.5)
            + labs(x = "log(control abundance)",
                   y = "log(fold change)")
)


png(filename=file.path(path, "plots/lfc-lca.png"),width = 500, height = 480) 
print(gg_2dim)
dev.off()


fit_2d <- gam(power ~ te(lcontrol, abs_lfc),
              data = comb,
              family = quasibinomial)
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

brkvec <- c(0.01, 0.02, 0.05, (1:9)/10)
gg_2dimc <- (gg_2dim
             + geom_contour(data = pp,
                            aes(z=power),
                            breaks = brkvec)
             + geom_label_contour(data = pp, aes(z= power),
                                  breaks = brkvec)
)

png(filename=file.path(path, "plots/with_contour.png"),width = 500, height = 480) 
print(gg_2dimc)
dev.off()

## power vs log fold change
gg_abslfc <- (ggplot(comb)
              + aes(abs_lfc, power)
              + geom_point()
              + geom_smooth(method = "gam",
                            method.args = list(family = "quasibinomial"))
)

png(filename=file.path(path, "plots/pwr_fold.png"),width = 500, height = 480) 
print(gg_abslfc)
dev.off()

gg_lcontrol <- (ggplot(comb)
                + aes(lcontrol, power)
                + geom_point()
                + geom_smooth(method = "gam",
                              method.args = list(family = "quasibinomial"))
)
png(filename=file.path(path, "plots/pwr_contabund.png"),width = 500, height = 480) 
print(gg_lcontrol)
dev.off()

## BMB: stopped here
######################################

#MA Begins 
# what is the relationshio
df <- tibble()
###########Power verse Effects

eff <- abs(obs_lfcmat)
mean  <- rowMeans(eff)
lower <- c(); upper <- c()
for (i in 1:nrow(eff)){
  lower[i] <- min(eff[i,])
  upper[i] <- max(eff[i,])
}

df <- data.frame(power, mean, lower, upper)
#png(filename= paste(path,"plots/","pwr_vr_meaneffect1.png",sep=""),width = 480, height = 480)
ggplot(df, aes(x = power,y=mean, ymin= lower, ymax= upper)) +
  geom_point() +
  #geom_errorbar() +
  geom_line(color = "red") +
  theme(text = element_text(size = 20))  +
  labs(title="power and effect sizes  for \n individual taxa",
       x ="power", y = "mean effect sizes") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

####################################
x<-factor(1:length(power))
df <- data.frame(samples = x, power)

newdf <- df[order(power),]
#png(filename= paste(path,"plots/","power_vr_samples1.png",sep=""),width = 480, height = 480,)
ggplot(newdf, aes(x = fct_inorder(x), y= power))+
  geom_point() +
  theme(text = element_text(size = 8)) +
  labs(title="Ordered power for individual taxa",
       x ="taxa", y = "power") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()


#############
index <- which(power >= 0.8)
index
mm1 <- abs(LFOLDCHANGE[index,])
n <- 1#nrow(mm1)
l1 <- c(min(mm1)); u1 <- c(max(mm1))
#for (i in 1:n){
# l1[i] <- min(mm1[i,])
# u1[i] <- max(mm1[i,])
#}
m1 <- mean(mm1)

##########
index <- which(power >= 0.7 & power < 0.8)
index
mm2 <- abs(LFOLDCHANGE[index,])
l2 <- c(min(mm2)); u2 <- c(max(mm2))
m2<-mean(mm2)
#########
index <- which(power >= 0.6 & power < 0.7)
index
mm3 <- abs(LFOLDCHANGE[index,])
l3 <- c(min(mm3)); u3 <- c(max(mm3))
m3<-mean(mm3)
#########
index <- which(power >= 0.5 & power < 0.6)
index
mm4 <- abs(LFOLDCHANGE[index,])
l4 <- c(min(mm4)); u4 <- c(max(mm4))
m4<-mean(mm4)
#########
index <- which(power < 0.5)
index
l5 <- c() ; u5 <- c()
mm5 <- abs(LFOLDCHANGE[index,])
k <- rowMeans(mm5)
n <- nrow(mm5)
for (i in 1:n){
  l5[i] <- min(mm5[i,])
  u5[i] <- max(mm5[i,])
}

m5 <- mean(k)
l5m <- mean(l5)
u5m <- mean(u5)

means<-c(m5,m3,m2); lower <- c(l5m,l3,l2); upper <- c(u5m,u3,u2)
df <- data.frame(power= c("<50%", "60 - 69%", "70 - 79%"),effects = means, lower, upper)

png(filename= paste(path,"plots/","groupedpwr_vr_effect1.png",sep=""),width = 480, height = 480)
ggplot(df, aes(x =fct_inorder(power), y=effects, ymin=lower, ymax =upper)) +
  geom_point() +
  geom_errorbar() +
  theme(text = element_text(size = 20))  +
  labs(title="mean and range of effect sizes for \n various power groupings",
       x ="power", y = "effect sizes") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

##################

men<- rowMeans(mm2)
mean_men <- mean(men)
l1 <- c( min(mm2[1,]), min(mm2[2,]))
lm <- mean(l1)
u1 <-c( max(mm2[1,]), max(mm2[2,]))
um <- mean(u1)

k2<-rowMeans(mm2)
lower2 <- c( min(mm2[1,]), min(mm2[2,]) )
upper2 <- - c( max(mm2[1,]), max(mm2[2,]) )

index <- which(power >= 0.5 & power < 0.6)
index
mm3 <- abs(LFOLDCHANGE[index,])
k3<-rowMeans(mm3)

# index <- which(power >= 0.6 & power < 0.7)
# mm4 <- abs(LFOLDCHANGE[index,])
# k4<-rowMeans(mm4)
#
#index <- which(power >= 0.5 & power < 0.6)

# mm5 <- abs(LFOLDCHANGE[index,])
# k5<-rowMeans(mm5)

############
ndf <- data.frame(power=power,means,lower,upper)
base <- ggplot(ndf, aes(x=power, y=means, ymin = lower, ymax = upper))
#png(filename= paste(path,"plots/","effect_vr_power.png",sep=""),width = 480, height = 480)
base + geom_point() +
  geom_errorbar() +
  geom_line(color="red") +
  theme(text = element_text(size = 20))  +
  labs(title="power and  \n size for power groupings",
       x ="power", y = "effect size") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#############
m2
mean1
###################
#plot power against effect sizes and see how it will look
library(tidyverse)
df <- data.frame(power=c("<=50","60%-80%",">80%"), effects=c(ave_mean,m2,mean1),
                 l= c(lmean,l2,l1),u=c(umean,u2,u1))

base <- ggplot(df, aes(x=fct_inorder(power), y=effects, ymin = l, ymax = u ))
png(filename= paste(path,"plots/","mean_vr_effect.png",sep=""),width = 480, height = 480,)
base + geom_point() +
  geom_errorbar() +
  theme(text = element_text(size = 20))  +
  labs(title="Averages of means and ranges of effect \n size for power groupings",
       x ="power", y = "effect size") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
# base + geom_crossbar()
# base + geom_pointrange()
# base + geom_errorbar()

##############
mdf <- data.frame(power=c("60%-80%"), effects=mean_men, lm,um)

base <- ggplot(mdf, aes(power, effects, ymin = lm, ymax = um ))
base + geom_point() +
  geom_errorbar() +
  theme(text = element_text(size = 20))  +
  labs(title="Mean and range of effect size for taxa with power >80%",
       x ="power", y = "effect size") +
  theme(plot.title = element_text(hjust = 0.5))


#################


#################################
mdf2 <- data.frame(power=c("<=50%", "60%-80%"), effects=c(ave_mean,mean_men),
                   l= c(lmean,lm),u =c(umean,um))
png(filename= paste(path,"plots/","mean_vr_effect.png",sep=""),width = 480, height = 480,)
base <- ggplot(mdf2, aes(power, effects, ymin=l, ymax= u ))
base + geom_point() +
  geom_errorbar() +
  theme(text = element_text(size = 20))  +
  labs(title="Mean and range of effect sizes for taxa \n with different power",
       x ="power", y = "effect size") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

###########################
base <- ggplot(mdf2, aes(power, effects, ymin=l, ymax= u ))
base + geom_point() +
  geom_errorbar() +
  theme(text = element_text(size = 20))  +
  labs(title="Mean and range of effect sizes for taxa \n with different power",
       x ="power", y = "effect size") +
  theme(plot.title = element_text(hjust = 0.5))
######################

means <- c(k2, k6)

lowers <- c(lower2,lower)
uppers <- c(upper2,upper)
df <- data.frame(index = 1:133,mean = means,lower=lowers,upper =uppers )

base <- ggplot(df, aes(index, mean, ymin = lowers, ymax = uppers ))
base + geom_point()
base + geom_smooth()
base + geom_crossbar()
base + geom_pointrange()
base + geom_errorbar()

length(means)

View(means)
k6<-rowMeans(mm6)
k6
data.frame(por = c("<50","40-60"), y=k6, z=k2)
k2
#max(k6)
x = seq(1:length(power))
par(mfrow = c(2, 3))
par(cex = 0.6)
par(mar = c(3, 3, 0, 0), oma = c(1, 1, 1, 1))
plot(k1[order(k1)])
plot(k2[order(k2)])
plot(k3[order(k3)])
plot(k4[order(k4)])
plot(k5[order(k5)])
plot(k6[order(k6)])

# set.seed(101); r <- sample(c(-1,1), size = 1000, replace= TRUE)
# *rexp(1000); hist(r, freq=FALSE, breaks = 100);
# curve(dnorm(x, sd = sd(r)), add=TRUE, col = 2, lwd = 3)


max(k2[order(k3)])
min(k1[order(k1)])
##########


n<- nrow(mm)
r<- matrix(0,n ,2)

for (i in 1:n){
  r[i,1] = range(mm[i, ])[1]
  r[i,2] = range(mm[i, ])[2]
}

m <- rowMeans(mm)
x <- seq(1,n)
df <- data.frame(index , y = m, lower = r[,1], upper=r[,2])

base <- ggplot(df, aes(index, y, ymin = lower, ymax = upper ))
base + geom_crossbar()
base + geom_pointrange()
base + geom_crossbar()
base + geom_errorbar()
base + geom_linerange()

#They have a wide range? what does this tell me? Nothing
# 80 -100
# 70 - 80
# 60 - 70
# < 50

##################################################################
####Parameters for control abundance
#shape 24.151835   2.941457
#rate   3.286305   0.404418
shape<- 24.151835
rate <- 3.286305

n = 133 #sample size in control group/per group
control_abund<-rgamma(n,shape= shape, rate = rate)
write.csv(control_abund, file= file.path(path,"Sim_cont.csv"))
lcontrol<-read.csv(file.path(path,"Sim_cont.csv"), row.names = 1)
#######
n <- nrow(lcontrol)
#(Intercept)  cont_abund
#1.1155695  -0.1401233
inter <- 1.1155695
cont_coef  <- -0.1401233
sd <- 0.623
a<- rep(1,n)
means <- inter*a +  cont_coef*lcontrol$x + rnorm(n,mean=0, sd=sd)

effect_sim <- rnorm(n,mean = means, sd = 0.05)
write.csv(effect_sim, file= file.path(path,"True_lfoldchange.csv"))

lfoldchange <-read.csv(file.path(path,"True_lfoldchange.csv"), row.names = 1)
treat_sim <- lfoldchange$x + lcontrol$x

write.csv(treat_sim, file= file.path(path,"Sim_treat.csv"))
ltreat <- read.csv(file= file.path(path,"Sim_treat.csv"),row.names = 1)
##### Dispersion
av_means <- rowMeans(data.frame(2^(lcontrol$x), 2^(ltreat$x)))
av_means <- rowMeans(data.frame(2^(lcontrol), 2^(ltreat)))


#write.csv(av_means, file= paste(path,"average_means.csv",sep=""))

#asymptDisp   extraPois
# 0.6735505 218.4284917
asymptDisp <-0.6735505; extraPois <- 218.4284917

dispersion1 <- c()
for (i in 1:length(av_means)){
  dispersion1[i] =  extraPois/av_means[i] + asymptDisp
}

write.csv(dispersion1, file= file.path(path,"dispersion1.csv"))

################################################################
data <- data.frame(av_means,dispersion)

plot(m,av_means)
plot(av_means,dispersion)
points(m)
av_means
length(DD)

##################################

n_sim <- 11; sample_per_Grp <- 120
disp_treat <- dispersion;  disp_cont <- dispersion;
treat_means <- 2^(lcontrol$x); cont_means <- 2^(control_abund$x)

Data_simulate(path,treat_means,disp_treat,cont_means,disp_cont,n_sim,sample_per_Grp)

data <- read.csv(file= paste(path,"sim1.csv",sep=""),row.names = 1)
#View(data)

my_data <- list()
for (i in 1:n_sim){
  my_data[[i]]<-  read.csv(paste(path,paste0("sim",i),".csv",sep=""),
                           row.names = 1)
}



results_data(path, my_data,conditions)


##########
res_data <- list()
for (i in 1:n_sim){
  res_data[[i]]<-  read.csv(paste(path,paste0("result",i),".csv",sep=""),
                            row.names = 1)
}

View(res_data[[1]])
p_adjust <- matrix(0,nrow(res_data[[1]]), n_sim)
for (i in 1:n_sim){
  p_adjust[,i] <- (res_data[[i]])$padj
}


LFOLDCHANGE <- matrix(0,nrow(res_data[[1]]), n_sim)
for (i in 1:n_sim){
  LFOLDCHANGE[,i] <- (res_data[[i]])$log2FoldChange
}
#p_adjust<-read.csv(paste(path,"p_adjust.csv",sep=""),row.names = 1)

power <- c()
for (i in 1:nrow(p_adjust)){
  power[i] <-sum(p_adjust[i,] < 0.05)/ncol(p_adjust)
}

index <- which(power == 1)
full_pwr <- LFOLDCHANGE[index,]
View(full_pwr)

x<- seq(0, 1, length.out = length(pwr))

pwr <- power[order(power)]
max(pwr)
full_pwrr <- abs(full_pwr)
View(full_pwrr)
plot(full_pwrr)
