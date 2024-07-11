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
m=20
n=20
M=20000
alpha1=1.0
alpha2=0.5
lambda=10
x=matrix(0,m,M)
y=matrix(0,n,M)
z1=matrix(0,1,M)
z2=matrix(0,1,M)
x1=array(0,M)
x2=array(0,M)
x3=array(0,M)
x4=array(0,M)
x5=array(0,M)
x6=array(0,M)
y1=array(0,M)
y2=array(0,M)
y3=array(0,M)
y4=array(0,M)
y5=array(0,M)
y6=array(0,M)
g1=array(0,M)
g2=array(0,M)
g3=array(0,M)
g4=array(0,M)
count1=0
count2=0


for (j in 1:M) 
{
  x[,j]=rgen.exp(m,alpha1,lambda) 
  y[,j]=rgen.exp(n,alpha2,lambda)
  z1[,j]=rgen.exp(1,alpha1,lambda)
  z2[,j]=rgen.exp(1,alpha2,lambda)
  
  LA1=function(param){
    alpha1=param[1]
    alpha2=param[2]
    lambda=param[3]
    sum(dgen.exp(z1[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(x[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y[,j],alpha=alpha2,lambda=lambda,log = TRUE))
  }
  mle1=maxLik(logLik=LA1,start=c(1,3,2))
  x1[j]=mle1$estimate[1]
  x2[j]=mle1$estimate[2]
  x3[j]=mle1$estimate[3]
  
  LA2=function(param){
    alpha1=param[1]
    alpha2=param[2]
    lambda=param[3]
    sum(dgen.exp(z1[,j],alpha=alpha2,lambda=lambda,log = TRUE))+sum(dgen.exp(x[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y[,j],alpha=alpha2,lambda=lambda,log = TRUE))
  }
  mle2=maxLik(logLik=LA2,start=c(1,3,2))
  x4[j]=mle2$estimate[1]
  x5[j]=mle2$estimate[2]
  x6[j]=mle2$estimate[3]
  
  LA3=function(param){
    alpha1=param[1]
    alpha2=param[2]
    lambda=param[3]
    sum(dgen.exp(z2[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(x[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y[,j],alpha=alpha2,lambda=lambda,log = TRUE))
  }
  mle3=maxLik(logLik=LA3,start=c(1,3,2))
  y1[j]=mle3$estimate[1]
  y2[j]=mle3$estimate[2]
  y3[j]=mle3$estimate[3]
  
  LA4=function(param){
    alpha1=param[1]
    alpha2=param[2]
    lambda=param[3]
    sum(dgen.exp(z2[,j],alpha=alpha2,lambda=lambda,log = TRUE))+sum(dgen.exp(x[,j],alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y[,j],alpha=alpha2,lambda=lambda,log = TRUE))
  }
  mle4=maxLik(logLik=LA4,start=c(1,3,2))
  y4[j]=mle4$estimate[1]
  y5[j]=mle4$estimate[2]
  y6[j]=mle4$estimate[3]
  
  g1[j]=dgen.exp(z1[,j],x1[j],x3[j])*prod(dgen.exp(x[,j],x1[j],x3[j]))*prod(dgen.exp(y[,j],x2[j],x3[j]))
  
  g2[j]=dgen.exp(z1[,j],x5[j],x6[j])*prod(dgen.exp(x[,j],x4[j],x6[j]))*prod(dgen.exp(y[,j],x5[j],x6[j]))
    
  g3[j]=dgen.exp(z2[,j],y1[j],y3[j])*prod(dgen.exp(x[,j],y1[j],y3[j]))*prod(dgen.exp(y[,j],y2[j],y3[j]))
    
  g4[j]=dgen.exp(z2[,j],y5[j],y6[j])*prod(dgen.exp(x[,j],y4[j],y6[j]))*prod(dgen.exp(y[,j],y5[j],y6[j]))
  
  if(g1[j]>g2[j]){count1=count1+1}
  if(g3[j]<g4[j]){count2=count2+1}
}
p11=count1/M
p22=count2/M
epc=0.5*(p11+p22)
epm=1-epc
cat("\n",epm)



