library(nlme)
library(bbmle)

simulate.pilot = function(m, mn, sigma,c0, nu0, nu1) {
  
  W=rnorm(m, mean=mn, sd=sigma)
  pilot.dat=data.frame(W)
  scal=(exp(c0+nu0*W+nu1*W^2))
  pilot.dat=transform(pilot.dat, Y=rnorm(m, mean=0,sd=scal))
  
  pilot.dat   
} 

View(effects[[1]])
mn=3.3;sigma=sqrt(0.5); m=200
c0=-5.413;   nu0=-0.2; nu1=-0.3

simu.dat=simulate.pilot(m, mn, sigma, c0, nu0, nu1)

p<-mle2(Y ~ dnorm(0, sd=exp(a+b*W + c*W^2)),
     start =list(a=-5.413, b=-0.2, c= -0.3), method = "Nelder-Mead" ,dat=simu.dat)
summary(p)

####################################################################
fit1= gls(Y~W, data=simu.dat, weights=varComb(varExp(form=~W),
                                                 varExp(form=~I(W^2))))

nu0.hat=2*fit1$modelStruct$varStruct$A[1]
nu1.hat=2*fit1$modelStruct$varStruct$B[1]

nu0.hat
nu1.hat










##########################################
set.seed(123)
n<-1000

z0<-data.frame(days=runif(n,-3,3))

z0_mat <- model.matrix(~ days, z0)
View(z0_mat)
betas<-c(10,.85)
a0<-1.5
b0<-8
b1<- -0.2
x_target<-60
alpha<-(log(b0)-log(a0))/x_target
z0$weight <- z0_mat%*%betas+rnorm(n,0,a0*exp(alpha*z0_mat[,2]+b1*z0_mat[,2]^2))
z0 %>% ggplot(aes(x=days,y=weight))+geom_point()

View(z0)
v1  <- varComb(
  varExp(form=~fitted(.)),
  varExp(form=~I(fitted(.)^2)))
m01 <- gls(weight~ days, weights= v1, data = z0) 
m01

scal <- runif(100,0,1) 
weight <- rcauchy(100, location=0, scale=scal)
days <- rnorm(100,0,0.3)
data <- data.frame(weight,days)
mle2(weight ~ dcauchy(0, scale=exp(a+b*days)),
     start =list(a=0.5, b=0.51),data=data) 

p<- -0.94942590 + 0.05137998*days
mean(p)


#expected error in DaDA2 212


