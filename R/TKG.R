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
m=20
n=20
M=20000
a1=1
b1=1
a2=1
b2=1
a3=1
b3=1
alpha1=1
alpha2=0.5
lambda=8
x=matrix(0,m,M)
y=matrix(0,n,M)
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
z1=array(0,M)
z2=array(0,M)
z3=array(0,M)
biasc1=array(0,M)
biasc2=array(0,M)
biasc3=array(0,M)
msec1=array(0,M)
msec2=array(0,M)
msec3=array(0,M)
H=array(0,M)
I=array(0,M)
J=array(0,M)
MLtko=matrix(0,3,M)
MLtkalpha1=matrix(0,3,M)
MLtkalpha2=matrix(0,3,M)
MLtklambda=matrix(0,3,M)
sigma0=array(0,M)
sigmaalpha1=array(0,M)
sigmaalpha2=array(0,M)
sigmalambda=array(0,M)
k1=array(0,M)
k2=array(0,M)
k3=array(0,M)
TK1=array(0,M)
TK2=array(0,M)
TK3=array(0,M)
biasTk1=array(0,M)
biasTk2=array(0,M)
biasTk3=array(0,M)
mseTk1=array(0,M)
mseTk2=array(0,M)
mseTk3=array(0,M)
x1test=array(0,M)
x2test=array(0,M)

rule1tk=array(0,M)
rule2tk=array(0,M)

count1tk=0
count2tk=0
LL=function(thetaa,x,y){
  alpha1=thetaa[1]
  alpha2=thetaa[2]
  lambda=thetaa[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))
  -loglik
}
LL1=function(theta,x,y,m,n,a1,b1,a2,b2,a3,b3){
  alpha1=theta[1]
  alpha2=theta[2]
  lambda=theta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+dgamma(alpha1,a1,b1,log = TRUE)+dgamma(alpha2,a2,b2,log = TRUE)+dgamma(lambda,a3,b3,log = TRUE)
  -loglik
}
LL2=function(thet,x,y,m,n,a1,b1,a2,b2,a3,b3){
  alpha1=thet[1]
  alpha2=thet[2]
  lambda=thet[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+dgamma(alpha1,a1,b1,log = TRUE)+dgamma(alpha2,a2,b2,log = TRUE)+dgamma(lambda,a3,b3,log = TRUE)+log(alpha1)
  -loglik
}
LL3=function(thett,x,y,m,n,a1,b1,a2,b2,a3,b3){
  alpha1=thett[1]
  alpha2=thett[2]
  lambda=thett[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+dgamma(alpha1,a1,b1,log = TRUE)+dgamma(alpha2,a2,b2,log = TRUE)+dgamma(lambda,a3,b3,log = TRUE)+log(alpha2)
  -loglik
}
LL4=function(thetta,x,y,m,n,a1,b1,a2,b2,a3,b3){
  alpha1=thetta[1]
  alpha2=thetta[2]
  lambda=thetta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+dgamma(alpha1,a1,b1,log = TRUE)+dgamma(alpha2,a2,b2,log = TRUE)+dgamma(lambda,a3,b3,log = TRUE)+log(lambda)
  -loglik
}
for (j in 1:M) 
{
  x[,j]=rgen.exp(m,alpha1,lambda)
  y[,j]=rgen.exp(n,alpha2,lambda)
  xstart=c(4,3,2)
  z1[j]=optim(xstart,LL,x=x[,j],y=y[,j])$par[1]
  z2[j]=optim(xstart,LL,x=x[,j],y=y[,j])$par[2]
  z3[j]=optim(xstart,LL,x=x[,j],y=y[,j])$par[3]
  x1[j]=optim(xstart,LL1,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[1]
  x2[j]=optim(xstart,LL1,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[2]
  x3[j]=optim(xstart,LL1,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[3]
  x4[j]=optim(xstart,LL2,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[1]
  x5[j]=optim(xstart,LL2,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[2]
  x6[j]=optim(xstart,LL2,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[3]
  y1[j]=optim(xstart,LL3,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[1]
  y2[j]=optim(xstart,LL3,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[2]
  y3[j]=optim(xstart,LL3,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[3]
  y4[j]=optim(xstart,LL4,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[1]
  y5[j]=optim(xstart,LL4,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[2]
  y6[j]=optim(xstart,LL4,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par[3]
  biasc1[j]=(z1[j]-alpha1)
  biasc2[j]=(z2[j]-alpha2)
  biasc3[j]=(z3[j]-lambda)
  msec1[j]=(z1[j]-alpha1)^2
  msec2[j]=(z2[j]-alpha2)^2
  msec3[j]=(z3[j]-lambda)^2
  
  #.............................TK approximation
  H[j]=-LL2(c(x4[j],x5[j],x6[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)
    
  I[j]=-LL3(c(y1[j],y2[j],y3[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)
 
  J[j]=-LL4(c(y4[j],y5[j],y6[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)
  
  
  MLtko[,j]=optim(c(x1[j],x2[j],x3[j]),LL1,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par
  Ho=optim(c(x1[j],x2[j],x3[j]),LL1,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3,hessian = TRUE)$hessian
  MLtkalpha1[,j]=optim(c(x4[j],x5[j],x6[j]),LL2,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par
  Halpha1=optim(c(x4[j],x5[j],x6[j]),LL2,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3,hessian = TRUE)$hessian
  MLtkalpha2[,j]=optim(c(y1[j],y2[j],y3[j]),LL3,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par
  Halpha2=optim(c(y1[j],y2[j],y3[j]),LL3,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3,hessian = TRUE)$hessian
  MLtklambda[,j]=optim(c(y4[j],y5[j],y6[j]),LL4,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3)$par
  Hlambda=optim(c(y4[j],y5[j],y6[j]),LL4,x=x[,j],y=y[,j],m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3,hessian = TRUE)$hessian
  sigma0[j]=det(solve(Ho))
  sigmaalpha1[j]=det(solve(Halpha1))
  sigmaalpha2[j]=det(solve(Halpha2))
  sigmalambda[j]=det(solve(Hlambda))
  k1[j]=sigmaalpha1[j]/sigma0[j]
  TK1[j]=sqrt(k1[j])*exp(H[j])
  k2[j]=sigmaalpha2[j]/sigma0[j]
  TK2[j]=sqrt(k2[j])*exp(I[j])
  k3[j]=sigmalambda[j]/sigma0[j]
  TK3[j]=sqrt(k3[j])*exp(J[j])
  biasTk1[j]=(TK1[j]-alpha1)
  biasTk2[j]=(TK2[j]-alpha2)
  biasTk3[j]=(TK3[j]-lambda)
  mseTk1=(TK1[j]-alpha1)^2
  mseTk2=(TK2[j]-alpha2)^2
  mseTk3=(TK3[j]-lambda)^2
  
  x1test[j]=rgen.exp(1,alpha1,lambda)
  x2test[j]=rgen.exp(1,alpha2,lambda)
  
  rule1tk[j]=log(TK1[j]/TK2[j]) +(TK1[j]-TK2[j])*log(1-exp(-TK3[j]*x1test[j]))
  if(rule1tk[j]>=0){count1tk=count1tk+1}
  rule2tk[j]=log(TK1[j]/TK2[j]) +(TK1[j]-TK2[j])*log(1-exp(-TK3[j]*x2test[j]))
  if(rule2tk[j]<0){count2tk=count2tk+1} 
  
}
p11tk=count1tk/M
p22tk=count2tk/M

epctk=(1/2)*(p11tk+p22tk)
epmtk=1-epctk

cat("\n",format(round(mean(z1),3),nsmall=3),"\t",format(round(mean(z2),3),nsmall=3),"\t",format(round(mean(z3),3),nsmall=3))
cat("\n",format(round(mean(biasc1),3),nsmall=3),"\t",format(round(mean(biasc2),3),nsmall=3),"\t",format(round(mean(biasc3),3),nsmall=3))
cat("\n",format(round(mean(msec1),3),nsmall=3),"\t",format(round(mean(msec2),3),nsmall=3),"\t",format(round(mean(msec3),3),nsmall=3))

cat("\n",format(round(mean(TK1),3),nsmall=3),"\t",format(round(mean(TK2),3),nsmall=3),"\t",format(round(mean(TK3),3),nsmall=3))
cat("\n",format(round(mean(biasTk1),3),nsmall=3),"\t",format(round(mean(biasTk2),3),nsmall=3),"\t",format(round(mean(biasTk3),3),nsmall=3))
cat("\n",format(round(mean(mseTk1),3),nsmall=3),"\t",format(round(mean(mseTk2),3),nsmall=3),"\t",format(round(mean(mseTk3),3),nsmall=3))


cat("\n",format(round(epmtk,3),nsmall = 3))





