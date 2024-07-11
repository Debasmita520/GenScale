library(maxLik)
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

qgen.exp <- function (p, alpha, lambda, lower.tail = TRUE, log.p = FALSE)
{
  # if((!is.numeric(alpha)) || (!is.numeric(lambda)) || (!is.numeric(q)))
  #   stop("non-numeric argument to mathematical function")
  # if((min(alpha) <= 0) || (min(lambda) <= 0) || (q <= 0))    
  #   stop("Invalid arguments")
  u <- exp(log(lambda) + log(p))
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
library(nleqslv)
library(MASS)
z=c(594.4, 202.75, 168.37, 574.86, 225.65, 76.38, 156.67, 127.81, 813.87, 562.39,
    468.47, 135.09, 72.24, 497.94, 355.56, 569.07, 640.48, 200.76, 550.42, 748.75, 489.66, 678.06, 457.71,
    106.73, 716.3, 42.66, 80.4, 339.22, 70.09, 193.42)
zbar=mean(z)
alphaml=numeric(1)
lamdaml=numeric(1)
exml=1/zbar
  
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
    alpha=param[1]
    lamda=param[2]
    sum(dgen.exp(z,alpha=alpha,lambda=lamda,log=TRUE))
  }
  
  mle=maxLik(logLik = logLikFun, start = c(0.81,0.01))
  alphaml=mle$estimate[1]
  lamdaml=mle$estimate[2]
  
  
  
  
  
  print(paste("P-value:", ad.test(z, "pgen.exp", alpha=alphaml, lambda=lamdaml)$p.value))
    cat("\n",ad.test(z,"pexp",rate=exml)$p.value)
    #cat("\n",format(round(mean(alphaml),4),nsmall = 4),format(round(mean(lamdaml),4),nsmall = 4))
    
    # GenMLE<-fitdist(z,"gen.exp",start = list(alpha= alphaml,lambda=lamdaml))
    # GenMLE
    # cdfcomp(GenMLE,main  = "",xlab = "data",ylab = " ")
    #denscomp(GenMLE,main="",xlab = "data",ylab = " ")
    #ppcomp(GenMLE,main="",xlab = "Theoretical Probabilities",ylab = "Empirical Probabilities")
