library(ggplot2)
library(dplyr)
## uses denovo_simulation branch
library(glmmTMB)
library(MASS)
## remotes::install_github("chvlyl/ZIBR")
library(ZIBR)
## remotes::install_github("nyiuab/NBZIMM")
library(NBZIMM)
## library(splinectomeR)
library(lme4)
library(patchwork)

## Data Simulation
nIndiv <- 10; nTime <- 10; nTaxa <- 30
metadat <-function(nIndiv,nTime,nTaxa){

dd <- expand.grid(subject = factor(1:nIndiv),
                    taxa =factor(1:nTaxa), time = (1:nTime))
  
  grp = (rep(factor(c("control", "treatment")),length.out=nIndiv))
  dd$group <- grp
  dd
}

dd =metadat(nIndiv,nTime,nTaxa)

rowIndices <- rep(1:nTaxa, 1:nTaxa); colIndices <- sequence(1:nTaxa)
template <- sparseMatrix(rowIndices, colIndices, x =
                           if_else(rowIndices==colIndices,1,0))
theta <- template@x
dispersion =1
bet <- c(1, 0.2, 0.1, 0.2)

sim <- function(dispersion,data=dd,beta= bet, thet = theta,n=nIndiv){
  p  = simulate(~group*time+(-1 + taxa|subject:time),
                newdata = data,
                newparams = list(beta= beta, theta =thet),
                family=negative.binomial(theta=dispersion))
  data$count = unlist(p)
  data
}

dd =sim(dispersion = 1)
training <- dd[dd$time < floor(nTime*0.8),]
test <- anti_join(dd, training, by='time') 

#############Fit data
glmmTMBfit_rr <- glmmTMB(count ~ group*time+ rr(-1 + taxa|subject:time,3),
                         data = training, family = nbinom2)

#' @Question: thinking about the predict function gives an alternating values
predict <- predict(glmmTMBfit_rr,
                   newdata = test, 
                   type = "response",
                   allow.new.levels=TRUE)#re.form = NULL

summary(glmmTMBfit_rr)
str(predict)

#############Investigating if method is unbiased and asymptotically consistent

#' @Objective: showing that method is unbiased, 
#' 1: Need to create the random component to be added to simulated counts
#' to get the true values for calculating mean(true_values -  predict) 
#' @Question1: My understanding for showing a that an estimator is biased is  
#'            to show that for completely different sets true data 
#'            (ie, random + count components), mean(true_values -  predict) 
#'            are relatively the same. But I am asking myself hpw to define
#'            relatively (ie, what "relatively the same" would mean ) since, 
#'            they will surely not be the same 
#'
#' @Objective: Showing that a method is asymptotically consistent
#' @Question2: I have the same question here as in @Question1, as sample increasing, 
#'             I expect the mean(true_values -  predict)  to be similar. I wonder if
#'             I am thinking of this right. 
#'             
nsim = 20; err = c()
for(i in 1:nsim){
  
  dd =sim(dispersion = 1)
  training <- dd[dd$time < floor(nTime*0.8),]
  test <- anti_join(dd, training, by='time') 
  test_x <- test %>% dplyr::select(-c("count"))
  
  glmmTMBfit_rr <- glmmTMB(count ~ group*time+ rr(-1 + taxa|subject:time,3),
                           data = training, family = nbinom2)
  
  predict <- predict(glmmTMBfit_rr,
                     newdata = test[], 
                     type = "response",
                     allow.new.levels=TRUE)#re.form = NULL
  
  summary(glmmTMBfit_rr)
  str(predict)
  test$predicted <- predict
  View(test)
  err[i] = mean((test$count - test$predicted))
}
data = data.frame(mean_squared_errors =err, index = 1:length(err))

ggplot(data, aes(x = index,y = log(mean_squared_errors))) +
  geom_point() +   
  geom_line() 

