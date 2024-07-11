rgen.exp <- function(n, alpha, lambda)
{
  if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
    stop("non-numeric argument to mathematical function")
  #if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
  #stop("Invalid arguments")
  return(-(1.0/lambda) * log(1.0 - (runif(n) ^ (1/alpha))))
}
dgen.exp <- function (x, alpha, lambda, log = FALSE)
{
  # if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(x)))
  #   stop("non-numeric argument to mathematical function")
  # if((min(alpha) <= 0) || (min(lambda) <= 0) || (x <= 0))    
  #   stop("Invalid arguments")
  u <- exp(log(lambda) + log(x))
  pdf <- exp(log(alpha) + log(lambda) - u + (alpha - 1.0) * log(1.0 - exp(-u)))   
  if(log)
    pdf <- log(pdf)
  return(pdf)
}

library(MASS)
library(maxLik)
x=c(594.4, 202.75, 168.37, 574.86, 225.65, 76.38, 156.67, 127.81, 813.87, 562.39,
    468.47, 135.09, 72.24, 497.94, 355.56, 569.07, 640.48, 200.76, 550.42, 748.75, 489.66, 678.06, 457.71,
    106.73, 716.3, 42.66, 80.4, 339.22, 70.09, 193.42)
y=c(71.46,419.02,284.64,585.57,456.60,113.85,187.85,688.16,662.66,45.58,578.62, 756.70,594.29,166.49,99.72,707.36,765.14,187.13,145.96,350.70,547.44,116.99,375.81,581.60, 119.86,48.01,200.16,36.75,244.53, 83.55)
m=length(x)
n=length(y)

LL=function(param){
  alpha1=param[1]
  alpha2=param[2]
  lambda1=param[3]
  lambda2=param[4]
  sum(dgen.exp(x,alpha=alpha1,lambda=lambda1,log=TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda2,log=TRUE))
}
mle1=maxLik(logLik=LL,start=c(0.81,0.81,0.01,0.01))
u1=mle1$estimate[1]
u2=mle1$estimate[2]
u3=mle1$estimate[3]
u4=mle1$estimate[4]

LL0=function(param){
  alpha1=param[1]
  alpha2=param[2]
  lambda=param[3]
  sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log=TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log=TRUE))
}
mle2=maxLik(logLik=LL0,start=c(0.81,0.81,0.01))
v1=mle2$estimate[1]
v2=mle2$estimate[2]
v3=mle2$estimate[3]

L=prod(dgen.exp(x,u1,u3))*prod(dgen.exp(y,u2,u4))
L0=prod(dgen.exp(x,v1,v3))*prod(dgen.exp(y,v2,v3))
lL=-2*log(L0/L)
cat("\n",lL)
