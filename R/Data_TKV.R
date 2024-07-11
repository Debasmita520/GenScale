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
library(goftest)
library(stats)
library(fitdistrplus)
library(MASS)
c=1
x=c(594.4, 202.75, 168.37, 574.86, 225.65, 76.38, 156.67, 127.81, 813.87, 562.39,
    468.47, 135.09, 72.24, 497.94, 355.56, 569.07, 640.48, 200.76, 550.42, 748.75, 489.66, 678.06, 457.71,
    106.73, 716.3, 42.66, 80.4, 339.22, 70.09, 193.42)
y=c(71.46, 419.02, 284.64, 585.57, 456.60, 113.85, 187.85, 688.16, 662.66, 45.58,
    578.62, 756.70, 594.29, 166.49, 99.72, 707.36, 765.14, 187.13, 145.96, 350.70, 547.44, 116.99, 375.81,
    581.60, 119.86, 48.01, 200.16, 36.75, 244.53, 83.55)
m=length(x)
n=length(y)

LL=function(thetaa,x,y){
  alpha1=thetaa[1]
  alpha2=thetaa[2]
  lambda=thetaa[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))
  -loglik
}
LL1=function(theta,x,y,c){
  alpha1=theta[1]
  alpha2=theta[2]
  lambda=theta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)
  -loglik
}
LL2=function(thet,x,y,c){
  alpha1=thet[1]
  alpha2=thet[2]
  lambda=thet[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(alpha1)
  -loglik
}
LL3=function(thett,x,y,c){
  alpha1=thett[1]
  alpha2=thett[2]
  lambda=thett[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(alpha2)
  -loglik
}
LL4=function(thetta,x,y,c){
  alpha1=thetta[1]
  alpha2=thetta[2]
  lambda=thetta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(lambda)
  -loglik
}
xstart=c(1.84,1.61,0.0039)
z1=optim(xstart,LL,x=x,y=y)$par[1]
z2=optim(xstart,LL,x=x,y=y)$par[2]
z3=optim(xstart,LL,x=x,y=y)$par[3]
x1=optim(xstart,LL1,x=x,y=y,c=c)$par[1]
x2=optim(xstart,LL1,x=x,y=y,c=c)$par[2]
x3=optim(xstart,LL1,x=x,y=y,c=c)$par[3]
x4=optim(xstart,LL2,x=x,y=y,c=c)$par[1]
x5=optim(xstart,LL2,x=x,y=y,c=c)$par[2]
x6=optim(xstart,LL2,x=x,y=y,c=c)$par[3]
y1=optim(xstart,LL3,x=x,y=y,c=c)$par[1]
y2=optim(xstart,LL3,x=x,y=y,c=c)$par[2]
y3=optim(xstart,LL3,x=x,y=y,c=c)$par[3]
y4=optim(xstart,LL4,x=x,y=y,c=c)$par[1]
y5=optim(xstart,LL4,x=x,y=y,c=c)$par[2]
y6=optim(xstart,LL4,x=x,y=y,c=c)$par[3]

H=-LL2(c(x4,x5,x6),x=x,y=y,c=c)+LL1(c(x1,x2,x3),x=x,y=y,c=c)

I=-LL3(c(y1,y2,y3),x=x,y=y,c=c)+LL1(c(x1,x2,x3),x=x,y=y,c=c)

J=-LL4(c(y4,y5,y6),x=x,y=y,c=c)+LL1(c(x1,x2,x3),x=x,y=y,c=c)

MLtko=optim(c(x1,x2,x3),LL1,x=x,y=y,c=c)$par
Ho=optim(c(x1,x2,x3),LL1,x=x,y=y,c=c,hessian = TRUE)$hessian
MLtkalpha1=optim(c(x4,x5,x6),LL2,x=x,y=y,c=c)$par
Halpha1=optim(c(x4,x5,x6),LL2,x=x,y=y,c=c,hessian = TRUE)$hessian
MLtkalpha2=optim(c(y1,y2,y3),LL3,x=x,y=y,c=c)$par
Halpha2=optim(c(y1,y2,y3),LL3,x=x,y=y,c=c,hessian = TRUE)$hessian
MLtklambda=optim(c(y4,y5,y6),LL4,x=x,y=y,c=c)$par
Hlambda=optim(c(y4,y5,y6),LL4,x=x,y=y,c=c,hessian = TRUE)$hessian
sigma0=det(solve(Ho))
sigmaalpha1=det(solve(Halpha1))
sigmaalpha2=det(solve(Halpha2))
sigmalambda=det(solve(Hlambda))
k1=sigmaalpha1/sigma0
TK1=sqrt(k1)*exp(H)
k2=sigmaalpha2/sigma0
TK2=sqrt(k2)*exp(I)
k3=sigmalambda/sigma0
TK3=sqrt(k3)*exp(J)

Tk1=mean(TK1)
Tk2=mean(TK2)
Tk3=mean(TK3)

rule1tk=log(Tk1/Tk2)+(Tk1-Tk2)*log(1-exp(-Tk3*x))
p11=length(which(rule1tk>=0))/length(x)
rule2tk=log(Tk1/Tk2)+(Tk1-Tk2)*log(1-exp(-Tk3*y))
p22=length(which(rule2tk<0))/length(y)
epc=0.5*(p11+p22)
epm=1-epc

cat("\n",format(round(mean(Tk1),4),nsmall=4),"\t",format(round(mean(Tk2),4),nsmall=4),"\t",format(round(mean(Tk3),4),nsmall=4))
cat("\n",format(round(epc,4),nsmall = 4))







