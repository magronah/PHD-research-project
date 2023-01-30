library(fitdistrplus)
library(BiocManager)
library(nlme)
library(truncdist)
library(DESeq2)
library(tidyverse)
library(patchwork)
library(bbmle)

#data
x = seq(1,10,0.02)
df = data.frame(x)
y <-apply(df,1,function(x){
  rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x )})
dat = data.frame(x,y)
###########################################
TruncCauchy = function(x, loc, scale, low, high){
  z = (x - loc) / scale;
  A = (atan((high - loc) / scale) - atan((low - loc) / scale))
  if (low <= x && x <= high){
    y = 1 / (pi * scale * (1 + z^2) * A)
  }
  else{
    y = 0
  }
  y
}

mm<-mle2(y ~ TruncCauchy(loc=0,scale=exp(a+b*x + c*x^2),low =l, high = h),
         start =list(a=0.1, b=0.1, c= 0.001, l=-5, h= 5), data=dat)

########################################### 
# truccauchy <- function(mu,s,A,B,x){
#   (1/s)*(1/(1+((x-mu)/s)^2))/(atan((B -mu)/s) - atan((A-mu)/s))
# }
# 
# 
# truccauchy <- function(par,mu=0,data){
#   s=exp(par[1]+par[2]*x + par[3]*x^2)
#   sum((1/s)*(1/(1+((x)/s)^2))/(atan((par[5])/s) - atan((par[4])/s)) - y)^2
# }
# 
# (m1 = optim(par = c(0.5,0.2,0.001,3,16),fn = truccauchy, data = dat, method ="BFGS"))
# s = exp(-0.77188813 -0.69554310*x-0.07589402*x^2 )
# plot(s)

# rtrunc(1, 'cauchy', a=-5, b=5, location=0, scale=x)
# truccauchy2 <- function(par,data, log=FALSE){
#   sum((1/par)*(x) - y)^2
# }




TruncCauchy = function(x, loc=0, par){
  scale=exp(par[1]+par[2]*x + par[3]*x^2);
  high  = par[4]; low = par[5]
  z = (x - loc) / scale;
  A = (atan((high - loc) / scale) - atan((low - loc) / scale))
  if (low <= x && x <= high){
    y = 1 / (pi * scale * (1 + z^2) * A)
  }
  else{
    y = 0
  }
  y
}



x = data.frame(seq(0,10,0.2))
x = seq(0,10,0.02)
res <- apply(eff,1, "TruncCauchy",loc=0, scale=0.7, low=-5, high=5)
res <- apply(eff,1, "truccauchy",mu=0, s=0.4, A=-5, B=5)

max(res)
hist(res, breaks = 30,prob=T)
lines(density(res))
logcontrol <- read_data(Control_Treat, "LogControl") 
logcontrol = unlist(logcontrol[[3]])
logcontr = which(is.finite(logcontrol))
logcont  = logcontrol[logcontr]; eff = (abs(deseq_list[[3]]$log2FoldChange))[logcontr]
dat = data.frame(logcont, eff)

res <- cauchy_sample(n, params){
 repeat {
             bad <- res < lower_bound | res > upper_bound
             if (!any(bad)) break
             res[bad] <- cauchy_sample(sum(bad), params)
             } 
}




mm<-mle2(eff ~ TruncCauchy(loc=b*logcont,0.4,5,2),
         start =list(b=1), data=dat)

dcauchy1 <- function(x, location=0, scale=1, log = FALSE)
  .Call(C_dcauchy, x, location, scale, log)








TruncCauchy2 = function(x, loc, scale, low, high){

  A = scale/(atan( (high - loc) / scale) - atan((low - loc) / scale))
  if (low <= x && x <= high){
    y = A/ (scale^2 + (x - loc)^2)
  }
  else{
    y = 0
  }
  y
}


set.seed(101)
n = 20; notu = 2000
Lfoldchange = rtrunc(notu, 'cauchy', a=-5, b=5, location=0, scale=0.5)
foldchange = 2^Lfoldchange
control_mean = 2^rnorm(notu,5, 2); treat_mean = control_mean*foldchange


#########
dispersion = runif(notu,0,1)
countmeans  =  rowMeans(data.frame(control_mean,treat_mean))
df = data.frame(countmeans ,dispersion)

data <- data.frame(apply(df,1,function(x){
    rnegbin(n=2*n, mu = x[1], theta = x[2])})) 
  rownames(data)  <-  paste0("Sample",1:(2*n))
  colnames(data)  <-  paste0("OTU",1:notu)

metadata =  data.frame(Sample=paste0("Sample",1:(2*n)),
                       Groups=factor(rep(c("control", "treatment"),each=n)))

remove = which(colMeans(data) == 0)
data  =  data[,-remove]

#Filtering
keep <- rowSums(data) >= 50; data <- data[keep,]
dds <- DESeqDataSetFromMatrix(t(data),metadata,~Groups) 
dds$Groups <- relevel(dds$Groups, ref = "control")
dds <- DESeq(dds,sfType ="poscounts", minReplicatesForReplace=Inf) 
res <- results(dds, cooksCutoff=TRUE, independentFiltering=TRUE)
reslt <- lfcShrink(dds, res=res, coef=2, type="apeglm") 
result <- data.frame(reslt)

#### Control and treatment means
ControlData =  data[1:n,]; TreatData = data[-c(1:n),]
ControlMean =  colMeans(ControlData); TreatMean = colMeans(TreatData)

r1 = which(ControlMean == 0); r2= which(TreatMean == 0)
r =  c(r1, r2)
r

if(length(r)>0){
  ControlMean = ControlMean[-r]; TreatMean = TreatMean[-r]
  logfolchange = res$log2FoldChange[-r]
}

logcontrol = log2(ControlMean); logtreat  = log2(TreatMean)
stopifnot(which(is.infinite(logcontrol)) == NULL)
stopifnot(which(is.infinite(logtreat)) == NULL)
stopifnot(which((logfolchange) == 0) == NULL)

dat = data.frame(logcontrol, effect = (logfolchange))
mm<-mle2(effect ~ dcauchy(0, scale=exp(a+b*logcontrol + c*logcontrol^2)),
         start =list(a=0.1, b=0.1, c= 0.001), data=dat)

coef(mm)
#Scale-location plot was used to validate the relationship and we showed that the variance
#seem to follow a quadractic relationship. 
pl


## Notes 
df = data.frame(logcontrol, logfolchangeSq = (logfolchange)^2)
ggplot(df, aes(x=logcontrol , y=logfolchangeSq)) +
  geom_point()  + 
  geom_smooth(method = "gam") +
  geom_smooth(method = "lm", formula = y ~poly(x,2), color="red")+
  ggtitle(paste0("Control against Fold Changes")) +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 16))



# yi^2 is a measure of an estimates of the variance of f(xi) so I will use that to check
# THis is just a check and not for you to be worried about. 

#Model with linear and quadratic and see 


plot(log2(control_mean), res$log2FoldChange)

Splot(log(countmeans), log(res$baseMean))
abline(0,1)

#Validate how to model the relationship. 


