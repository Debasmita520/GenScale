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
for (alpha2 in c(0.5,1,2,3,4,5)) {
m=20
n=20
M=20000
alpha1=1
#alpha2=0.5
lambda=1.5
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
LL1=function(theta,x,y){
  alpha1=theta[1]
  alpha2=theta[2]
  lambda=theta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(alpha2)-log(lambda)
  -loglik
}
LL2=function(thet,x,y){
  alpha1=thet[1]
  alpha2=thet[2]
  lambda=thet[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha2)-log(lambda)
  -loglik
}
LL3=function(thett,x,y){
  alpha1=thett[1]
  alpha2=thett[2]
  lambda=thett[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(lambda)
  -loglik
}
LL4=function(thetta,x,y){
  alpha1=thetta[1]
  alpha2=thetta[2]
  lambda=thetta[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(alpha2)
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
  x1[j]=optim(xstart,LL1,x=x[,j],y=y[,j])$par[1]
  x2[j]=optim(xstart,LL1,x=x[,j],y=y[,j])$par[2]
  x3[j]=optim(xstart,LL1,x=x[,j],y=y[,j])$par[3]
  x4[j]=optim(xstart,LL2,x=x[,j],y=y[,j])$par[1]
  x5[j]=optim(xstart,LL2,x=x[,j],y=y[,j])$par[2]
  x6[j]=optim(xstart,LL2,x=x[,j],y=y[,j])$par[3]
  y1[j]=optim(xstart,LL3,x=x[,j],y=y[,j])$par[1]
  y2[j]=optim(xstart,LL3,x=x[,j],y=y[,j])$par[2]
  y3[j]=optim(xstart,LL3,x=x[,j],y=y[,j])$par[3]
  y4[j]=optim(xstart,LL4,x=x[,j],y=y[,j])$par[1]
  y5[j]=optim(xstart,LL4,x=x[,j],y=y[,j])$par[2]
  y6[j]=optim(xstart,LL4,x=x[,j],y=y[,j])$par[3]
  biasc1[j]=(z1[j]-alpha1)
  biasc2[j]=(z2[j]-alpha2)
  biasc3[j]=(z3[j]-lambda)
  msec1[j]=(z1[j]-alpha1)^2
  msec2[j]=(z2[j]-alpha2)^2
  msec3[j]=(z3[j]-lambda)^2
  
  #.............................TK approximation
  H[j]=-LL2(c(x4[j],x5[j],x6[j]),x=x[,j],y=y[,j])+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j])
  
  I[j]=-LL3(c(y1[j],y2[j],y3[j]),x=x[,j],y=y[,j])+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j] )
  
  J[j]=-LL4(c(y4[j],y5[j],y6[j]),x=x[,j],y=y[,j])+LL1(c(x1[j],x2[j],x3[j]),x=x[,j],y=y[,j] )
  
  
  MLtko[,j]=optim(c(x1[j],x2[j],x3[j]),LL1,x=x[,j],y=y[,j] )$par
  Ho=optim(c(x1[j],x2[j],x3[j]),LL1,x=x[,j],y=y[,j] ,hessian = TRUE)$hessian
  MLtkalpha1[,j]=optim(c(x4[j],x5[j],x6[j]),LL2,x=x[,j],y=y[,j] )$par
  Halpha1=optim(c(x4[j],x5[j],x6[j]),LL2,x=x[,j],y=y[,j] ,hessian = TRUE)$hessian
  MLtkalpha2[,j]=optim(c(y1[j],y2[j],y3[j]),LL3,x=x[,j],y=y[,j] )$par
  Halpha2=optim(c(y1[j],y2[j],y3[j]),LL3,x=x[,j],y=y[,j] ,hessian = TRUE)$hessian
  MLtklambda[,j]=optim(c(y4[j],y5[j],y6[j]),LL4,x=x[,j],y=y[,j] )$par
  Hlambda=optim(c(y4[j],y5[j],y6[j]),LL4,x=x[,j],y=y[,j] ,hessian = TRUE)$hessian
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


LG1=round(mean(TK1),2)
LG2=round(mean(TK2),2)
LG3=round(mean(TK3),2)
biasLG1=round(mean(biasTk1),3)
biasLG2=round(mean(biasTk2),3)
biasLG3=round(mean(biasTk3),3)
mseLG1=round(mean(mseTk1),3)
mseLG2=round(mean(mseTk2),3)
mseLG3=round(mean(mseTk3),3)
cat("\n",alpha2,"\t",LG1,"\t",LG2,"\t",LG3,"\t",biasLG1,"\t",biasLG2,"\t",biasLG3,"\t",mseLG1,"\t",mseLG2,"\t",mseLG3,"\t",epmtk)

# cat("\n",format(round(mean(TK1),3),nsmall=3),"\t",format(round(mean(TK2),3),nsmall=3),"\t",format(round(mean(TK3),3),nsmall=3))
# cat("\n",format(round(mean(biasTk1),3),nsmall=3),"\t",format(round(mean(biasTk2),3),nsmall=3),"\t",format(round(mean(biasTk3),3),nsmall=3))
# cat("\n",format(round(mean(mseTk1),3),nsmall=3),"\t",format(round(mean(mseTk2),3),nsmall=3),"\t",format(round(mean(mseTk3),3),nsmall=3))
#cat("\n",format(round(epmtk,3),nsmall = 3))
}





