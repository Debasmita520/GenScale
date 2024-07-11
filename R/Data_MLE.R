library(maxLik)
set.seed(123)
rgen.exp <- function(n, alpha, lambda)
{
  if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(n)))
    stop("non-numeric argument to mathematical function")
  if((min(alpha) <= 0) || (min(lambda) <= 0) || (n <= 0))
    stop("Invalid arguments")
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
pgen.exp <- function (q, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
  # if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
  #   stop("non-numeric argument to mathematical function")
  # if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
  #   stop("Invalid arguments")
  u <- exp(log(lambda) + log(q))
  cdf <- exp(alpha * log(1.0 - exp(-u)))   
  if(!lower.tail)
    cdf <- 1.0 - cdf
  if(log.p)
    cdf <- log(cdf)
  return(cdf)
}



library(goftest)
library(stats)
library(fitdistrplus)
library(MASS)
x=c(594.4, 202.75, 168.37, 574.86, 225.65, 76.38, 156.67, 127.81, 813.87, 562.39,
    468.47, 135.09, 72.24, 497.94, 355.56, 569.07, 640.48, 200.76, 550.42, 748.75, 489.66, 678.06, 457.71,
    106.73, 716.3, 42.66, 80.4, 339.22, 70.09, 193.42)
y=c(71.46, 419.02, 284.64, 585.57, 456.60, 113.85, 187.85, 688.16, 662.66, 45.58,
    578.62, 756.70, 594.29, 166.49, 99.72, 707.36, 765.14, 187.13, 145.96, 350.70, 547.44, 116.99, 375.81,
    581.60, 119.86, 48.01, 200.16, 36.75, 244.53, 83.55)
alpha1ml=c()
alpha2ml=c()
lamdaml=c()


# 
# fnewton=function(x){
#   y=numeric(2)
#   y[1]=x[1]+(25/(sum(log(1-exp(-x[2]*z[i:(i+24)])))))
#   y[2]=x[2]-1/((sum((z[i:(i+24)]*exp(-x[2]*z[i:(i+24)]))/(1-exp(-x[2]*z[i:(i+24)])))/(sum(log(1-exp(-x[2]*z[i:(i+24)])))))+(1/(length(z[i:(i+24)])))*(sum(z[i:(i+24)]/(1-exp(-x[2]*z[i:(i+24)])))))
#   y
# }
# alphaml[i]=nleqslv(c(3,2),fnewton,control=list(btol=0.0001),method="Newton")$x[1]
# lamdaml[i]=nleqslv(c(3,2),fnewton,control=list(btol=0.0001),method="Newton")$x[2]

logLikFun=function(param){
  alpha1=param[1]
  alpha2=param[2]
  lamda=param[3]
  sum(dgen.exp(x,alpha=alpha1,lambda=lamda,log=TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lamda,log=TRUE))
}

mle=maxLik(logLik = logLikFun, start = c(0.81,0.81,0.01))
alpha1ml=mle$estimate[1]
alpha2ml=mle$estimate[2]
lamdaml=mle$estimate[3]



rule1ml=log(alpha1ml/alpha2ml)+(alpha1ml-alpha2ml)*log(1-exp(-lamdaml*x))
p11=length(which(rule1ml>=0))/length(x)
rule2ml=log(alpha1ml/alpha2ml)+(alpha1ml-alpha2ml)*log(1-exp(-lamdaml*y))
p22=length(which(rule2ml<0))/length(y)
epc=0.5*(p11+p22)
epm=1-epc


cat("\n",format(round(mean(alpha1ml),4),nsmall=4),"\t",format(round(mean(alpha2ml),4),nsmall=4),"\t",format(round(mean(lamdaml),4),nsmall=4))
cat("\n",format(round(epc,4),nsmall = 4))
