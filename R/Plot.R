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
N=10000
a1=1
b1=1
a2=1
b2=1
a3=1
b3=1
c=1
alpha1=1
alpha2=0.5
lambda=1.5
x=matrix(0,m,M)
y=matrix(0,n,M)
z1=array(0,M)
z2=array(0,M)
z3=array(0,M)

U=array(0,M)
V=array(0,M)
W=array(0,M)
S=array(0,M)
M1=array(0,M)
D1=array(0,M)
sig11=array(0,M)
sig12=array(0,M)
sig13=array(0,M)
sig22=array(0,M)
sig23=array(0,M)
sig33=array(0,M)
L111=array(0,M)
L222=array(0,M)
L333=array(0,M)
L331=array(0,M)
L332=array(0,M)
A11=array(0,M)
B11=array(0,M)
C11=array(0,M)
q1=array(0,M)
q2=array(0,M)
q3=array(0,M)
m1=array(0,M)
m2=array(0,M)
m3=array(0,M)
Lg1=array(0,M)
Lg2=array(0,M)
Lg3=array(0,M)

m11=array(0,M)
m21=array(0,M)
m31=array(0,M)
Lg11=array(0,M)
Lg21=array(0,M)
Lg31=array(0,M)

m12=array(0,M)
m22=array(0,M)
m32=array(0,M)
Lg12=array(0,M)
Lg22=array(0,M)
Lg32=array(0,M)

m13=array(0,M)
m23=array(0,M)
m33=array(0,M)
Lg13=array(0,M)
Lg23=array(0,M)
Lg33=array(0,M)


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

x1VJ=array(0,M)
x2VJ=array(0,M)
x3VJ=array(0,M)
x4VJ=array(0,M)
x5VJ=array(0,M)
x6VJ=array(0,M)
y1VJ=array(0,M)
y2VJ=array(0,M)
y3VJ=array(0,M)
y4VJ=array(0,M)
y5VJ=array(0,M)
y6VJ=array(0,M)

x1V=array(0,M)
x2V=array(0,M)
x3V=array(0,M)
x4V=array(0,M)
x5V=array(0,M)
x6V=array(0,M)
y1V=array(0,M)
y2V=array(0,M)
y3V=array(0,M)
y4V=array(0,M)
y5V=array(0,M)
y6V=array(0,M)

x1J=array(0,M)
x2J=array(0,M)
x3J=array(0,M)
x4J=array(0,M)
x5J=array(0,M)
x6J=array(0,M)
y1J=array(0,M)
y2J=array(0,M)
y3J=array(0,M)
y4J=array(0,M)
y5J=array(0,M)
y6J=array(0,M)

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

HVJ=array(0,M)
IVJ=array(0,M)
JVJ=array(0,M)
MLtkoVJ=matrix(0,3,M)
MLtkalpha1VJ=matrix(0,3,M)
MLtkalpha2VJ=matrix(0,3,M)
MLtklambdaVJ=matrix(0,3,M)
sigma0VJ=array(0,M)
sigmaalpha1VJ=array(0,M)
sigmaalpha2VJ=array(0,M)
sigmalambdaVJ=array(0,M)
k1VJ=array(0,M)
k2VJ=array(0,M)
k3VJ=array(0,M)
TK1VJ=array(0,M)
TK2VJ=array(0,M)
TK3VJ=array(0,M)

HV=array(0,M)
IV=array(0,M)
JV=array(0,M)
MLtkoV=matrix(0,3,M)
MLtkalpha1V=matrix(0,3,M)
MLtkalpha2V=matrix(0,3,M)
MLtklambdaV=matrix(0,3,M)
sigma0V=array(0,M)
sigmaalpha1V=array(0,M)
sigmaalpha2V=array(0,M)
sigmalambdaV=array(0,M)
k1V=array(0,M)
k2V=array(0,M)
k3V=array(0,M)
TK1V=array(0,M)
TK2V=array(0,M)
TK3V=array(0,M)

HJ=array(0,M)
IJ=array(0,M)
JJ=array(0,M)
MLtkoJ=matrix(0,3,M)
MLtkalpha1J=matrix(0,3,M)
MLtkalpha2J=matrix(0,3,M)
MLtklambdaJ=matrix(0,3,M)
sigma0J=array(0,M)
sigmaalpha1J=array(0,M)
sigmaalpha2J=array(0,M)
sigmalambdaJ=array(0,M)
k1J=array(0,M)
k2J=array(0,M)
k3J=array(0,M)
TK1J=array(0,M)
TK2J=array(0,M)
TK3J=array(0,M)


alpha1mc=array(0,M)
alpha2mc=array(0,M)
lamdamc=array(0,M)

alpha11mc=array(0,M)
alpha21mc=array(0,M)
lamda1mc=array(0,M)

alpha12mc=array(0,M)
alpha22mc=array(0,M)
lamda2mc=array(0,M)

alpha13mc=array(0,M)
alpha23mc=array(0,M)
lamda3mc=array(0,M)

#........................MLE
LL=function(thetaa,x,y){
  alpha1=thetaa[1]
  alpha2=thetaa[2]
  lambda=thetaa[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))
  -loglik
}
#.........................TK(Gamma)
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

#.................................TK(Vague)
LL9=function(the5,x,y,c){
  alpha1=the5[1]
  alpha2=the5[2]
  lambda=the5[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)
  -loglik
}
LL10=function(the6,x,y,c){
  alpha1=the6[1]
  alpha2=the6[2]
  lambda=the6[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(alpha1)
  -loglik
}
LL11=function(the7,x,y,c){
  alpha1=the7[1]
  alpha2=the7[2]
  lambda=the7[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(alpha2)
  -loglik
}
LL12=function(the8,x,y,c){
  alpha1=the8[1]
  alpha2=the8[2]
  lambda=the8[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-c*log(lambda)+log(lambda)
  -loglik
}

#.................................TK(Vague Jeffreys)
LL5=function(the1,x,y,m,n){
  alpha1=the1[1]
  alpha2=the1[2]
  lambda=the1[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(alpha2)-log(lambda)
  -loglik
}
LL6=function(the2,x,y,m,n){
  alpha1=the2[1]
  alpha2=the2[2]
  lambda=the2[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha2)-log(lambda)
  -loglik
}
LL7=function(the3,x,y,m,n){
  alpha1=the3[1]
  alpha2=the3[2]
  lambda=the3[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(lambda)
  -loglik
}
LL8=function(the4,x,y,m,n){
  alpha1=the4[1]
  alpha2=the4[2]
  lambda=the4[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))-log(alpha1)-log(alpha2)
  -loglik
}
#.................................TK(Jeffreys)
LL13=function(the9,x,y,I1){
  alpha1=the9[1]
  alpha2=the9[2]
  lambda=the9[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+(1/2)*log(I1)
  -loglik
}
LL14=function(the10,x,y,I1){
  alpha1=the10[1]
  alpha2=the10[2]
  lambda=the10[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+(1/2)*log(I1)+log(alpha1)
  -loglik
}
LL15=function(the11,x,y,I1){
  alpha1=the11[1]
  alpha2=the11[2]
  lambda=the11[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+(1/2)*log(I1)+log(alpha2)
  -loglik
}
LL16=function(the12,x,y,I1){
  alpha1=the12[1]
  alpha2=the12[2]
  lambda=the12[3]
  loglik=sum(dgen.exp(x,alpha=alpha1,lambda=lambda,log = TRUE))+sum(dgen.exp(y,alpha=alpha2,lambda=lambda,log = TRUE))+(1/2)*log(I1)+log(lambda)
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
  
  #Lindley's approximation(gamma)..........#
  U[j]=(m/(z1[j])^2)
  V[j]=-sum((x[,j]*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j]))))
  W[j]=(n/(z2[j])^2)
  S[j]=-sum((y[,j]*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j]))))
  M1[j]=(m/(z3[j])^2)+(z1[j]-1)*sum(((x[,j])^2*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j])))^2)+(n/(z3[j])^2)+(z2[j]-1)*sum(((y[,j])^2*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j])))^2)
  D1[j]=(U[j]*((W[j]*M1[j])-(S[j])^2))-((V[j])^2*W[j])
  
  
  sig11[j]=(((W[j]*M1[j])-(S[j])^2)/D1[j])
  sig12[j]=((V[j]*S[j])/D1[j])
  sig13[j]=-((V[j]*W[j])/D1[j])
  sig22[j]=(((U[j]*M1[j])-(V[j])^2)/D1[j])
  sig23[j]=-((U[j]*S[j])/D1[j])
  sig33[j]=((U[j]*W[j])/D1[j])
  
  
  L111[j]=((2*m)/(z1[j])^3)
  L222[j]=((2*n)/(z2[j])^3)
  L333[j]=((2*m)/(z3[j])^3)+(z1[j]-1)*sum(((x[,j])^3*exp(-(z3[j]*x[,j]))*(1+exp(-(z3[j]*x[,j]))))/((1-exp(-(z3[j]*x[,j])))^3))+((2*n)/(z3[j])^3)+(z2[j]-1)*sum(((y[,j])^3*exp(-(z3[j]*y[,j]))*(1+exp(-(z3[j]*y[,j]))))/((1-exp(-(z3[j]*y[,j])))^3))
  L331[j]=-sum(((x[,j])^2*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j])))^2)
  L332[j]=-sum(((y[,j])^2*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j])))^2)
  
  
  A11[j]=((sig11[j]*L111[j])+(sig33[j]*L331[j]))
  B11[j]=((sig22[j]*L222[j])+(sig33[j]*L332[j]))
  C11[j]=(2*sig13[j]*L331[j])+(2*sig23[j]*L332[j])+(sig33[j]*L333[j])
  
  
  q1[j]=(0.5*((A11[j]*sig11[j])+(B11[j]*sig12[j])+(C11[j]*sig13[j])))
  q2[j]=(0.5*((A11[j]*sig12[j])+(B11[j]*sig22[j])+(C11[j]*sig23[j])))
  q3[j]=(0.5*((A11[j]*sig13[j])+(B11[j]*sig23[j])+(C11[j]*sig33[j])))
  
  
  m1[j]=(((b1-1)/z1[j])-a1)
  m2[j]=(((b2-1)/z2[j])-a2)
  m3[j]=(((b3-1)/z3[j])-a3)
  
  Lg1[j]=z1[j]+(m1[j]*sig11[j])+(m2[j]*sig12[j])+(m3[j]*sig13[j])+q1[j]
  Lg2[j]=z2[j]+(m1[j]*sig12[j])+(m2[j]*sig22[j])+(m3[j]*sig23[j])+q2[j]
  Lg3[j]=z3[j]+(m1[j]*sig13[j])+(m2[j]*sig23[j])+(m3[j]*sig33[j])+q3[j]
  
  #Lindley's approximation(vague)..........#
  m11[j]=0
  m21[j]=0
  m31[j]=-(c/z3[j])
  
  Lg11[j]=z1[j]+(m11[j]*sig11[j])+(m21[j]*sig12[j])+(m31[j]*sig13[j])+q1[j]
  Lg21[j]=z2[j]+(m11[j]*sig12[j])+(m21[j]*sig22[j])+(m31[j]*sig23[j])+q2[j]
  Lg31[j]=z3[j]+(m11[j]*sig13[j])+(m21[j]*sig23[j])+(m31[j]*sig33[j])+q3[j]
  
  #Lindley's approximation(vague Jeffreys)..........#
  m12[j]=-(1/z1[j])
  m22[j]=-(1/z2[j])
  m32[j]=-(1/z3[j])
  
  Lg12[j]=z1[j]+(m12[j]*sig11[j])+(m22[j]*sig12[j])+(m32[j]*sig13[j])+q1[j]
  Lg22[j]=z2[j]+(m12[j]*sig12[j])+(m22[j]*sig22[j])+(m32[j]*sig23[j])+q2[j]
  Lg32[j]=z3[j]+(m12[j]*sig13[j])+(m22[j]*sig23[j])+(m32[j]*sig33[j])+q3[j]
  
  #Lindley's approximation(Jeffreys)..........#
  C1=(n/lambda^2)*(1+((alpha2*(alpha2-1))/(alpha2-2))*(trigamma(1)-trigamma(alpha2-1)+(digamma(alpha2-1)-digamma(1))^2))-((n*alpha2)/(lambda^2))*((trigamma(1)-trigamma(alpha2))+(digamma(alpha2)-digamma(1))^2)
  
  p1111=sqrt(m*(((n*C1)/alpha2^2)+((n^2)/(lambda^2))*((alpha2/(alpha2-1))*(digamma(alpha2)-digamma(1))-(digamma(alpha2+1)-digamma(1)))^2)+((n*m^2)/(alpha2^2*lambda^2))*(0.639996))
  
  
  f =expression(p1111)
  p_x = Deriv(f, "alpha1")
  p_y = Deriv(f, "alpha2")
  p_z= Deriv(f, "lambda")
  
  m13[j] = eval(p_x, list(alpha1 = z1[j], alpha2 = z2[j], lambda=z3[j]))
  m23[j] = eval(p_y, list(alpha1 = z1[j], alpha2 = z2[j], lambda=z3[j]))
  m33[j]=eval(p_z, list(alpha1 = z1[j], alpha2 = z2[j], lambda=z3[j]))
  Lg13[j]=z1[j]+(m13[j]*sig11[j])+(m23[j]*sig12[j])+(m33[j]*sig13[j])+q1[j]
  Lg23[j]=z2[j]+(m13[j]*sig12[j])+(m23[j]*sig22[j])+(m33[j]*sig23[j])+q2[j]
  Lg33[j]=z3[j]+(m13[j]*sig13[j])+(m23[j]*sig23[j])+(m33[j]*sig33[j])+q3[j]
  
  #...........................TK(Gamma)
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
  
  #...........................TK(Vague)
  x1V[j]=optim(xstart,LL9,x=x[,j],y=y[,j],c=c)$par[1]
  x2V[j]=optim(xstart,LL9,x=x[,j],y=y[,j],c=c)$par[2]
  x3V[j]=optim(xstart,LL9,x=x[,j],y=y[,j],c=c)$par[3]
  x4V[j]=optim(xstart,LL10,x=x[,j],y=y[,j],c=c)$par[1]
  x5V[j]=optim(xstart,LL10,x=x[,j],y=y[,j],c=c)$par[2]
  x6V[j]=optim(xstart,LL10,x=x[,j],y=y[,j],c=c)$par[3]
  y1V[j]=optim(xstart,LL11,x=x[,j],y=y[,j],c=c)$par[1]
  y2V[j]=optim(xstart,LL11,x=x[,j],y=y[,j],c=c)$par[2]
  y3V[j]=optim(xstart,LL11,x=x[,j],y=y[,j],c=c)$par[3]
  y4V[j]=optim(xstart,LL12,x=x[,j],y=y[,j],c=c)$par[1]
  y5V[j]=optim(xstart,LL12,x=x[,j],y=y[,j],c=c)$par[2]
  y6V[j]=optim(xstart,LL12,x=x[,j],y=y[,j],c=c)$par[3]
  
  #...........................TK(Vague Jeffreys)
  x1VJ[j]=optim(xstart,LL5,x=x[,j],y=y[,j],m=m,n=n)$par[1]
  x2VJ[j]=optim(xstart,LL5,x=x[,j],y=y[,j],m=m,n=n)$par[2]
  x3VJ[j]=optim(xstart,LL5,x=x[,j],y=y[,j],m=m,n=n)$par[3]
  x4VJ[j]=optim(xstart,LL6,x=x[,j],y=y[,j],m=m,n=n)$par[1]
  x5VJ[j]=optim(xstart,LL6,x=x[,j],y=y[,j],m=m,n=n)$par[2]
  x6VJ[j]=optim(xstart,LL6,x=x[,j],y=y[,j],m=m,n=n)$par[3]
  y1VJ[j]=optim(xstart,LL7,x=x[,j],y=y[,j],m=m,n=n)$par[1]
  y2VJ[j]=optim(xstart,LL7,x=x[,j],y=y[,j],m=m,n=n)$par[2]
  y3VJ[j]=optim(xstart,LL7,x=x[,j],y=y[,j],m=m,n=n)$par[3]
  y4VJ[j]=optim(xstart,LL8,x=x[,j],y=y[,j],m=m,n=n)$par[1]
  y5VJ[j]=optim(xstart,LL8,x=x[,j],y=y[,j],m=m,n=n)$par[2]
  y6VJ[j]=optim(xstart,LL8,x=x[,j],y=y[,j],m=m,n=n)$par[3]
  
  #...........................TK(Jeffreys)
  VJ=(n/lambda^2)*(1+((alpha2*(alpha2-1))/(alpha2-2))*(trigamma(1)-trigamma(alpha2-1)+(digamma(alpha2-1)-digamma(1))^2))-((n*alpha2)/(lambda^2))*((trigamma(1)-trigamma(alpha2))+(digamma(alpha2)-digamma(1))^2)
  
  VJ1=sqrt(m*(((n*VJ)/alpha2^2)+((n^2)/(lambda^2))*((alpha2/(alpha2-1))*(digamma(alpha2)-digamma(1))-(digamma(alpha2+1)-digamma(1)))^2)+((n*m^2)/(alpha2^2*lambda^2))*(0.639996))
  
  
  x1J[j]=optim(xstart,LL13,x=x[,j],y=y[,j],I1=VJ1)$par[1]
  x2J[j]=optim(xstart,LL13,x=x[,j],y=y[,j],I1=VJ1)$par[2]
  x3J[j]=optim(xstart,LL13,x=x[,j],y=y[,j],I1=VJ1)$par[3]
  x4J[j]=optim(xstart,LL14,x=x[,j],y=y[,j],I1=VJ1)$par[1]
  x5J[j]=optim(xstart,LL14,x=x[,j],y=y[,j],I1=VJ1)$par[2]
  x6J[j]=optim(xstart,LL14,x=x[,j],y=y[,j],I1=VJ1)$par[3]
  y1J[j]=optim(xstart,LL15,x=x[,j],y=y[,j],I1=VJ1)$par[1]
  y2J[j]=optim(xstart,LL15,x=x[,j],y=y[,j],I1=VJ1)$par[2]
  y3J[j]=optim(xstart,LL15,x=x[,j],y=y[,j],I1=VJ1)$par[3]
  y4J[j]=optim(xstart,LL16,x=x[,j],y=y[,j],I1=VJ1)$par[1]
  y5J[j]=optim(xstart,LL16,x=x[,j],y=y[,j],I1=VJ1)$par[2]
  y6J[j]=optim(xstart,LL16,x=x[,j],y=y[,j],I1=VJ1)$par[3]
  
  
  #.............................TK approximation(Gamma)
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
  
  
  
  #.............................TK approximation(Vague)
  HV[j]=-LL10(c(x4V[j],x5V[j],x6V[j]),x=x[,j],y=y[,j],c=c)+LL9(c(x1V[j],x2V[j],x3V[j]),x=x[,j],y=y[,j],c=c)
  
  IV[j]=-LL11(c(y1V[j],y2V[j],y3V[j]),x=x[,j],y=y[,j],c=c)+LL9(c(x1V[j],x2V[j],x3V[j]),x=x[,j],y=y[,j],c=c)
  
  JV[j]=-LL12(c(y4V[j],y5V[j],y6V[j]),x=x[,j],y=y[,j],c=c)+LL9(c(x1V[j],x2V[j],x3V[j]),x=x[,j],y=y[,j],c=c)
  
  
  MLtkoV[,j]=optim(c(x1V[j],x2V[j],x3V[j]),LL9,x=x[,j],y=y[,j],c=c)$par
  HoV=optim(c(x1V[j],x2V[j],x3V[j]),LL9,x=x[,j],y=y[,j],c=c,hessian = TRUE)$hessian
  MLtkalpha1V[,j]=optim(c(x4V[j],x5V[j],x6V[j]),LL10,x=x[,j],y=y[,j],c=c)$par
  Halpha1V=optim(c(x4V[j],x5V[j],x6V[j]),LL10,x=x[,j],y=y[,j],c=c,hessian = TRUE)$hessian
  MLtkalpha2V[,j]=optim(c(y1V[j],y2V[j],y3V[j]),LL11,x=x[,j],y=y[,j],c=c)$par
  Halpha2V=optim(c(y1V[j],y2V[j],y3V[j]),LL11,x=x[,j],y=y[,j],c=c,hessian = TRUE)$hessian
  MLtklambdaV[,j]=optim(c(y4V[j],y5V[j],y6V[j]),LL12,x=x[,j],y=y[,j],c=c)$par
  HlambdaV=optim(c(y4V[j],y5V[j],y6V[j]),LL12,x=x[,j],y=y[,j],c=c,hessian = TRUE)$hessian
  sigma0V[j]=det(solve(HoV))
  sigmaalpha1V[j]=det(solve(Halpha1V))
  sigmaalpha2V[j]=det(solve(Halpha2V))
  sigmalambdaV[j]=det(solve(HlambdaV))
  k1V[j]=sigmaalpha1V[j]/sigma0V[j]
  TK1V[j]=sqrt(k1V[j])*exp(HV[j])
  k2V[j]=sigmaalpha2V[j]/sigma0V[j]
  TK2V[j]=sqrt(k2V[j])*exp(IV[j])
  k3V[j]=sigmalambdaV[j]/sigma0V[j]
  TK3V[j]=sqrt(k3V[j])*exp(JV[j])
  
  #.............................TK approximation(Vague Jeffreys)
  HVJ[j]=-LL6(c(x4VJ[j],x5VJ[j],x6VJ[j]),x=x[,j],y=y[,j],m=m,n=n)+LL5(c(x1VJ[j],x2VJ[j],x3VJ[j]),x=x[,j],y=y[,j],m=m,n=n)
  
  IVJ[j]=-LL7(c(y1VJ[j],y2VJ[j],y3VJ[j]),x=x[,j],y=y[,j],m=m,n=n)+LL5(c(x1VJ[j],x2VJ[j],x3VJ[j]),x=x[,j],y=y[,j],m=m,n=n)
  
  JVJ[j]=-LL8(c(y4VJ[j],y5VJ[j],y6VJ[j]),x=x[,j],y=y[,j],m=m,n=n)+LL5(c(x1VJ[j],x2VJ[j],x3VJ[j]),x=x[,j],y=y[,j],m=m,n=n)
  
  
  MLtkoVJ[,j]=optim(c(x1VJ[j],x2VJ[j],x3VJ[j]),LL5,x=x[,j],y=y[,j],m=m,n=n)$par
  HoVJ=optim(c(x1VJ[j],x2VJ[j],x3VJ[j]),LL5,x=x[,j],y=y[,j],m=m,n=n,hessian = TRUE)$hessian
  MLtkalpha1VJ[,j]=optim(c(x4VJ[j],x5VJ[j],x6VJ[j]),LL6,x=x[,j],y=y[,j],m=m,n=n)$par
  Halpha1VJ=optim(c(x4VJ[j],x5VJ[j],x6VJ[j]),LL6,x=x[,j],y=y[,j],m=m,n=n,hessian = TRUE)$hessian
  MLtkalpha2VJ[,j]=optim(c(y1VJ[j],y2VJ[j],y3VJ[j]),LL7,x=x[,j],y=y[,j],m=m,n=n)$par
  Halpha2VJ=optim(c(y1VJ[j],y2VJ[j],y3VJ[j]),LL7,x=x[,j],y=y[,j],m=m,n=n,hessian = TRUE)$hessian
  MLtklambdaVJ[,j]=optim(c(y4VJ[j],y5VJ[j],y6VJ[j]),LL8,x=x[,j],y=y[,j],m=m,n=n)$par
  HlambdaVJ=optim(c(y4VJ[j],y5VJ[j],y6VJ[j]),LL8,x=x[,j],y=y[,j],m=m,n=n,hessian = TRUE)$hessian
  sigma0VJ[j]=det(solve(HoVJ))
  sigmaalpha1VJ[j]=det(solve(Halpha1VJ))
  sigmaalpha2VJ[j]=det(solve(Halpha2VJ))
  sigmalambdaVJ[j]=det(solve(HlambdaVJ))
  k1VJ[j]=sigmaalpha1VJ[j]/sigma0VJ[j]
  TK1VJ[j]=sqrt(k1VJ[j])*exp(HVJ[j])
  k2VJ[j]=sigmaalpha2VJ[j]/sigma0VJ[j]
  TK2VJ[j]=sqrt(k2VJ[j])*exp(IVJ[j])
  k3VJ[j]=sigmalambdaVJ[j]/sigma0VJ[j]
  TK3VJ[j]=sqrt(k3VJ[j])*exp(JVJ[j])
  
  #.............................TK approximation(Jeffreys)
  HJ[j]=-LL14(c(x4J[j],x5J[j],x6J[j]),x=x[,j],y=y[,j],I1=VJ1)+LL13(c(x1J[j],x2J[j],x3J[j]),x=x[,j],y=y[,j],I1=VJ1)
  
  IJ[j]=-LL15(c(y1J[j],y2J[j],y3J[j]),x=x[,j],y=y[,j],I1=VJ1)+LL13(c(x1J[j],x2J[j],x3J[j]),x=x[,j],y=y[,j],I1=VJ1)
  
  JJ[j]=-LL16(c(y4J[j],y5J[j],y6J[j]),x=x[,j],y=y[,j],I1=VJ1)+LL13(c(x1J[j],x2J[j],x3J[j]),x=x[,j],y=y[,j],I1=VJ1)
  
  
  MLtkoJ[,j]=optim(c(x1J[j],x2J[j],x3J[j]),LL13,x=x[,j],y=y[,j],I1=VJ1)$par
  HoJ=optim(c(x1J[j],x2J[j],x3J[j]),LL13,x=x[,j],y=y[,j],I1=VJ1,hessian = TRUE)$hessian
  MLtkalpha1J[,j]=optim(c(x4J[j],x5J[j],x6J[j]),LL14,x=x[,j],y=y[,j],I1=VJ1)$par
  Halpha1J=optim(c(x4J[j],x5J[j],x6J[j]),LL14,x=x[,j],y=y[,j],I1=VJ1,hessian = TRUE)$hessian
  MLtkalpha2J[,j]=optim(c(y1J[j],y2J[j],y3J[j]),LL15,x=x[,j],y=y[,j],I1=VJ1)$par
  Halpha2J=optim(c(y1J[j],y2J[j],y3J[j]),LL15,x=x[,j],y=y[,j],I1=VJ1,hessian = TRUE)$hessian
  MLtklambdaJ[,j]=optim(c(y4J[j],y5J[j],y6J[j]),LL16,x=x[,j],y=y[,j],I1=VJ1)$par
  HlambdaJ=optim(c(y4J[j],y5J[j],y6J[j]),LL16,x=x[,j],y=y[,j],I1=VJ1,hessian = TRUE)$hessian
  sigma0J[j]=det(solve(HoJ))
  sigmaalpha1J[j]=det(solve(Halpha1J))
  sigmaalpha2J[j]=det(solve(Halpha2J))
  sigmalambdaJ[j]=det(solve(HlambdaJ))
  k1J[j]=sigmaalpha1J[j]/sigma0J[j]
  TK1J[j]=sqrt(k1J[j])*exp(HJ[j])
  k2J[j]=sigmaalpha2J[j]/sigma0J[j]
  TK2J[j]=sqrt(k2J[j])*exp(IJ[j])
  k3J[j]=sigmalambdaJ[j]/sigma0J[j]
  TK3J[j]=sqrt(k3J[j])*exp(JJ[j])
  #...............................................MCMC(Gamma)
  MCMC=function(alpha1int,alpha2int,lamdaint,sd1,sd2,sd3,X1,X2,alpha1,alpha2,lamda,m,n,a1,b1,a2,b2,a3,b3,N){
    alpha1mc=array(0,N)
    alpha1mc[1]=alpha1int
    alpha2mc=array(0,N)
    alpha2mc[1]=alpha2int
    lamdamc=array(0,N)
    lamdamc[1]=lamdaint
    
    for (i in 2:N) {
      repeat{
        alpha1str=rnorm(1,alpha1int,sd1)
        if(alpha1str>0){
          break
        }
      }
      repeat{
        alpha2str=rnorm(1,alpha2int,sd2)
        if(alpha2str>0){
          break
        }
      }
      repeat{
        lamdastr=rnorm(1,lamdaint,sd3)
        if(lamdastr>0){
          break
        }
      }
      
      w11=array(0,m)
      w12=array(0,m)
      w21=array(0,n)
      w22=array(0,n)
      for (k1 in 1:m) {
        w11[k1]=dgen.exp(X1[k1],alpha1str,lamdastr)
        w12[k1]=dgen.exp(X1[k1],alpha1int,lamdaint)
        
      }
      for (k2 in 1:n) {
        w21[k2]=dgen.exp(X2[k2],alpha2str,lamdastr)
        w22[k2]=dgen.exp(X2[k2],alpha2int,lamdaint)
      }
      p11=dgamma(alpha1str,a1,b1)
      p12=dgamma(alpha1int,a1,b1)
      p21=dgamma(alpha2str,a2,b2)
      p22=dgamma(alpha2int,a2,b2)
      p31=dgamma(lamdastr,a3,b3)
      p32=dgamma(lamdaint,a3,b3)
      
      si1=prod(w11)*prod(w21)*p11*p21*p31
      si2=prod(w12)*prod(w22)*p12*p22*p32
      A=si1/si2
      si0=min(1,A) 
      u=runif(1)
      
      
      if(u<=si0){
        alpha1mc[i]=alpha1str
        alpha2mc[i]=alpha2str
        lamdamc[i]=lamdastr
      }else{
        alpha1mc[i]=alpha1mc[i-1]
        alpha2mc[i]=alpha2mc[i-1]
        lamdamc[i]=lamdamc[i-1]
      }
      
    }
    mylist=list("alpha1mc"=alpha1mc,"alpha2mc"=alpha2mc,"lamdamc"=lamdamc)
    return(mylist)
    
  }
  
  mcmc=MCMC(alpha1int=1.0,alpha2int=0.5,lamdaint=1.5,sd1=0.3,sd2=0.3,sd3=0.3,X1=x[,j],X2=y[,j],alpha1=alpha1,alpha2=alpha2,lamda=lamda,m=m,n=n,a1=a1,b1=b1,a2=a2,b2=b2,a3=a3,b3=b3,N=N)
  alpha1mc[j]=mean(mcmc$alpha1mc[N/2:N])
  alpha2mc[j]=mean(mcmc$alpha2mc[N/2:N])
  lamdamc[j]=mean(mcmc$lamdamc[N/2:N])
  
  
  #...............................................MCMC(Vague)
  MCMC2=function(alpha12int,alpha22int,lamda2int,sd12,sd22,sd32,Z1,Z2,alpha1,alpha2,lamda,m,n,c,N){
    alpha12mc=array(0,N)
    alpha12mc[1]=alpha12int
    alpha22mc=array(0,N)
    alpha22mc[1]=alpha22int
    lamda2mc=array(0,N)
    lamda2mc[1]=lamda2int
    
    for (i in 2:N) {
      repeat{
        alpha12str=rnorm(1,alpha12int,sd12)
        if(alpha12str>0){
          break
        }
      }
      repeat{
        alpha22str=rnorm(1,alpha22int,sd22)
        if(alpha22str>0){
          break
        }
      }
      repeat{
        lamda2str=rnorm(1,lamda2int,sd32)
        if(lamda2str>0){
          break
        }
      }
      
      w112=array(0,m)
      w122=array(0,m)
      w212=array(0,n)
      w222=array(0,n)
      for (k5 in 1:m) {
        w112[k5]=dgen.exp(Z1[k5],alpha12str,lamda2str)
        w122[k5]=dgen.exp(Z1[k5],alpha12int,lamda2int)
        
      }
      for (k6 in 1:n) {
        w212[k6]=dgen.exp(Z2[k6],alpha22str,lamda2str)
        w222[k6]=dgen.exp(Z2[k6],alpha22int,lamda2int)
      }
      p112=1
      p122=1
      p212=1
      p222=1
      p312=1/(lamda2str)^(c)
      p322=1/(lamda2int)^(c)
      
      
      si12=prod(w112)*prod(w212)*p112*p212*p312
      si22=prod(w122)*prod(w222)*p122*p222*p322
      C=si12/si22
      si02=min(1,C) 
      w=runif(1)
      
      
      if(w<=si02){
        alpha12mc[i]=alpha12str
        alpha22mc[i]=alpha22str
        lamda2mc[i]=lamda2str
      }else{
        alpha12mc[i]=alpha12mc[i-1]
        alpha22mc[i]=alpha22mc[i-1]
        lamda2mc[i]=lamda2mc[i-1]
      }
      
    }
    mylist=list("alpha12mc"=alpha12mc,"alpha22mc"=alpha22mc,"lamda2mc"=lamda2mc)
    return(mylist)
    
  }
  
  mcmc=MCMC2(alpha12int=1.0,alpha22int=0.5,lamda2int=1.5,sd12=0.3,sd22=0.3,sd32=0.3,Z1=x[,j],Z2=y[,j],alpha1=alpha1,alpha2=alpha2,lamda=lamda,m=m,n=n,c=c,N=N)
  alpha12mc[j]=mean(mcmc$alpha12mc[N/2:N])
  alpha22mc[j]=mean(mcmc$alpha22mc[N/2:N])
  lamda2mc[j]=mean(mcmc$lamda2mc[N/2:N])
  
  #...............................................MCMC(Vague Jeffreys)  
  MCMC1=function(alpha11int,alpha21int,lamda1int,sd11,sd21,sd31,Y1,Y2,alpha1,alpha2,lamda,m,n,N){
    alpha11mc=array(0,N)
    alpha11mc[1]=alpha11int
    alpha21mc=array(0,N)
    alpha21mc[1]=alpha21int
    lamda1mc=array(0,N)
    lamda1mc[1]=lamda1int
    
    for (i in 2:N) {
      repeat{
        alpha11str=rnorm(1,alpha11int,sd11)
        if(alpha11str>0){
          break
        }
      }
      repeat{
        alpha21str=rnorm(1,alpha21int,sd21)
        if(alpha21str>0){
          break
        }
      }
      repeat{
        lamda1str=rnorm(1,lamda1int,sd31)
        if(lamda1str>0){
          break
        }
      }
      
      w111=array(0,m)
      w121=array(0,m)
      w211=array(0,n)
      w221=array(0,n)
      for (k3 in 1:m) {
        w111[k3]=dgen.exp(Y1[k3],alpha11str,lamda1str)
        w121[k3]=dgen.exp(Y1[k3],alpha11int,lamda1int)
        
      }
      for (k4 in 1:n) {
        w211[k4]=dgen.exp(Y2[k4],alpha21str,lamda1str)
        w221[k4]=dgen.exp(Y2[k4],alpha21int,lamda1int)
      }
      p111=1/(alpha11str)
      p121=1/(alpha11int)
      p211=1/(alpha21str)
      p221=1/(alpha21int)
      p311=1/(lamda1str)
      p321=1/(lamda1int)
      
      si11=prod(w111)*prod(w211)*p111*p211*p311
      si21=prod(w121)*prod(w221)*p121*p221*p321
      B=si11/si21
      si01=min(1,B)
      v=runif(1)
      
      
      if(v<=si01){
        alpha11mc[i]=alpha11str
        alpha21mc[i]=alpha21str
        lamda1mc[i]=lamda1str
      }else{
        alpha11mc[i]=alpha11mc[i-1]
        alpha21mc[i]=alpha21mc[i-1]
        lamda1mc[i]=lamda1mc[i-1]
      }
      
    }
    mylist=list("alpha11mc"=alpha11mc,"alpha21mc"=alpha21mc,"lamda1mc"=lamda1mc)
    return(mylist)
    
  }
  mcmc=MCMC1(alpha11int=1.0,alpha21int=0.5,lamda1int=1.5,sd11=0.3,sd21=0.3,sd31=0.3,Y1=x[,j],Y2=y[,j],alpha1=alpha1,alpha2=alpha2,lamda=lamda,m=m,n=n,N=N)
  alpha11mc[j]=mean(mcmc$alpha11mc[N/2:N])
  alpha21mc[j]=mean(mcmc$alpha21mc[N/2:N])
  lamda1mc[j]=mean(mcmc$lamda1mc[N/2:N])
  #...............................................MCMC(Jeffreys) 
  
  MCMC3=function(alpha13int,alpha23int,lamda3int,sd13,sd23,sd33,Z3,Z4,alpha1,alpha2,lamda,m,n,N){
    alpha13mc=array(0,N)
    alpha13mc[1]=alpha13int
    alpha23mc=array(0,N)
    alpha23mc[1]=alpha23int
    lamda3mc=array(0,N)
    lamda3mc[1]=lamda3int
    
    for (i in 2:N) {
      repeat{
        alpha13str=rnorm(1,alpha13int,sd13)
        if(alpha13str>0){
          break
        }
      }
      repeat{
        alpha23str=rnorm(1,alpha23int,sd23)
        if(alpha23str>0){
          break
        }
      }
      repeat{
        lamda3str=rnorm(1,lamda3int,sd33)
        if(lamda3str>0){
          break
        }
      }
      w113=array(0,m)
      w123=array(0,m)
      w213=array(0,n)
      w223=array(0,n)
      for (k7 in 1:m) {
        w113[k7]=dgen.exp(Z3[k7],alpha13str,lamda3str)
        w123[k7]=dgen.exp(Z3[k7],alpha13int,lamda3int)
        
      }
      for (k8 in 1:n) {
        w213[k8]=dgen.exp(Z4[k8],alpha23str,lamda3str)
        w223[k8]=dgen.exp(Z4[k8],alpha23int,lamda3int)
      }
      
      A1=(n/lamda3str^2)*(1+((alpha23str*(alpha23str-1))/(alpha23str-2))*(trigamma(1)-trigamma(alpha23str-1)+(digamma(alpha23str-1)-digamma(1))^2))-((n*alpha23str)/(lamda3str^2))*((trigamma(1)-trigamma(alpha23str))+(digamma(alpha23str)-digamma(1))^2)
      B1=(n/lamda3int^2)*(1+((alpha23int*(alpha23int-1))/(alpha23int-2))*(trigamma(1)-trigamma(alpha23int-1)+(digamma(alpha23int-1)-digamma(1))^2))-((n*alpha23int)/(lamda3int^2))*((trigamma(1)-trigamma(alpha23int))+(digamma(alpha23int)-digamma(1))^2)
      p113=sqrt(m*(((n*A1)/alpha23str^2)+((n^2)/(lamda3str^2))*((alpha23str/(alpha23str-1))*(digamma(alpha23str)-digamma(1))-(digamma(alpha23str+1)-digamma(1)))^2)+((n*m^2)/(alpha23str^2*lamda3str^2))*(0.639996))
      p123=sqrt(m*(((n*B1)/alpha23int^2)+((n^2)/(lamda3int^2))*((alpha23int/(alpha23int-1))*(digamma(alpha23int)-digamma(1))-(digamma(alpha23int+1)-digamma(1)))^2)+((n*m^2)/(alpha23int^2*lamda3int^2))*(0.639996))
      
      si13=prod(w113)*prod(w213)*p113
      si23=prod(w123)*prod(w223)*p123
      D=si13/si23
      si03=min(1,D)
      w1=runif(1)
      
      
      if(w1<=si03){
        alpha13mc[i]=alpha13str
        alpha23mc[i]=alpha23str
        lamda3mc[i]=lamda3str
      }else{
        alpha13mc[i]=alpha13mc[i-1]
        alpha23mc[i]=alpha23mc[i-1]
        lamda3mc[i]=lamda3mc[i-1]
      }
      
    }
    
    mylist=list("alpha13mc"=alpha13mc,"alpha23mc"=alpha23mc,"lamda3mc"=lamda3mc)
    return(mylist)
    
  }
  
  mcmc=MCMC3(alpha13int=1.0,alpha23int=0.5,lamda3int=1.5,sd13=0.3,sd23=0.3,sd33=0.3,Z3=x[,j],Z4=y[,j],alpha1=alpha1,alpha2=alpha2,lamda=lamda,m=m,n=n,N=N)
  alpha13mc[j]=mean(mcmc$alpha13mc[N/2:N])
  alpha23mc[j]=mean(mcmc$alpha23mc[N/2:N])
  lamda3mc[j]=mean(mcmc$lamda3mc[N/2:N])
  
}

plot(cumsum(z1)/(1:M),type='l',col='darkred',xlab="Iteration",ylab="Values of the Estimators",xlim= c(1,M),ylim=c(0.82,1.8))
lines(cumsum(Lg1)/(1:M),type='l', col='skyblue',)
lines(cumsum(Lg11)/(1:M),type='l', col='black',)
lines(cumsum(Lg12)/(1:M),type='l', col='orange',)
lines(cumsum(Lg13)/(1:M),type='l', col='navy',)
lines(cumsum(TK1)/(1:M),type='l', col='brown')
lines(cumsum(TK1V)/(1:M),type='l', col='green')
lines(cumsum(TK1VJ)/(1:M),type='l', col='darkgrey')
lines(cumsum(TK1J)/(1:M),type='l', col='blue')
lines(cumsum(alpha1mc)/(1:M),type='l', col='red')
lines(cumsum(alpha11mc)/(1:M),type='l', col='purple')
lines(cumsum(alpha12mc)/(1:M),type='l', col='darkgreen')
lines(cumsum(alpha13mc)/(1:M),type='l', col='deeppink')
# # # # legend(60,0.95,legend=c("Line 1", "Line 2","Line 3","Line 4","Line 5","Line 6","Line 7"),
# # # #        fill = c("red", "blue","purple","orange","green","black","yellow"),lty=1:1)
legend("topright",text.width = 0.6,title="(a)",legend=c(expression(alpha[1][ML]), expression(alpha[1][LG]), expression(alpha[1][LV]), expression(paste(alpha[1][LJ]^"*")), expression(alpha[1][LJ]),expression(alpha[1][TKG]),expression(alpha[1][TKV]),expression(paste(alpha[1][TKJ]^"*")),expression(alpha[1][TKJ]),expression(alpha[1][MCG]),expression(alpha[1][MCV]),expression(paste(alpha[1][MCJ]^"*")),expression(alpha[1][MCJ])),
       fill = c("darkred", "skyblue","black","orange","navy","brown","green","darkgrey","blue","red","purple","darkgreen","deeppink"),
       lty=0.5:0.5,ncol=4,cex=0.8,adj=0.93)


#legend("topright",legend=c(as.expression(bquote('alpha'['1ML']), expression(alpha[BG]^TK),expression(alpha[BJ]^TK),expression(alpha[BV]^TK),expression(alpha[BG]^MG),expression(alpha[BJ]^MJ),expression(alpha[BV]^MV)),
#      fill = c("red", "black","green","blue","deeppink","purple","darkgreen"),lty=1:1)
# legend("topright",legend=c(expression(alpha[1][ML]), expression(alpha[1][LG]), expression(alpha[1][LV]), expression(paste(alpha[1][LJ]^"*")), expression(alpha[1][LJ]),expression(alpha[1][TKG]),expression(alpha[1][TKV]),expression(paste(alpha[1][TKJ]^"*")),expression(alpha[1][TKJ]),expression(alpha[1][MCG]),expression(alpha[1][MCV]),expression(paste(alpha[1][MCJ]^"*")),expression(alpha[1][MCJ])),
#        fill = c("red", "black","green","blue","deeppink","purple","darkgreen"),lty=1:1)

# plot(cumsum(z2)/(1:M),type='l',col='darkred',xlab="Iteration",ylab="Values of the Estimators",xlim= c(1,M),ylim=c(0.1,1.5))
# lines(cumsum(Lg2)/(1:M),type='l', col='skyblue',)
# lines(cumsum(Lg21)/(1:M),type='l', col='black',)
# lines(cumsum(Lg22)/(1:M),type='l', col='orange',)
# lines(cumsum(Lg23)/(1:M),type='l', col='navy',)
# lines(cumsum(TK2)/(1:M),type='l', col='brown',)
# lines(cumsum(TK2V)/(1:M),type='l', col='green')
# lines(cumsum(TK2VJ)/(1:M),type='l', col='darkgrey')
# lines(cumsum(TK2J)/(1:M),type='l', col='blue')
# lines(cumsum(alpha2mc)/(1:M),type='l', col='red')
# lines(cumsum(alpha21mc)/(1:M),type='l', col='purple')
# lines(cumsum(alpha22mc)/(1:M),type='l', col='darkgreen')
# lines(cumsum(alpha23mc)/(1:M),type='l', col='deeppink')
# legend("topright",text.width = 0.6,title="(b)",legend=c(expression(alpha[2][ML]), expression(alpha[2][LG]), expression(alpha[2][LV]), expression(paste(alpha[2][LJ]^"*")), expression(alpha[2][LJ]),expression(alpha[2][TKG]),expression(alpha[2][TKV]),expression(paste(alpha[2][TKJ]^"*")),expression(alpha[2][TKJ]),expression(alpha[2][MCG]),expression(alpha[2][MCV]),expression(paste(alpha[2][MCJ]^"*")),expression(alpha[2][MCJ])),
#        fill = c("darkred", "skyblue","black","orange","navy","brown","green","darkgrey","blue","red","purple","darkgreen","deeppink"),
#        lty=0.5:0.5,ncol=4,cex=0.8,adj=0.9)

# plot(cumsum(z3)/(1:M),type='l',col='darkred',xlab="Iteration",ylab="Values of the Estimators",xlim= c(1,M),ylim=c(0.7,3))
# lines(cumsum(Lg3)/(1:M),type='l', col='skyblue',)
# lines(cumsum(Lg31)/(1:M),type='l', col='black',)
# lines(cumsum(Lg32)/(1:M),type='l', col='orange',)
# lines(cumsum(Lg33)/(1:M),type='l', col='navy',)
# lines(cumsum(TK3)/(1:M),type='l', col='brown',)
# lines(cumsum(TK3V)/(1:M),type='l', col='green')
# lines(cumsum(TK3VJ)/(1:M),type='l', col='darkgrey')
# lines(cumsum(TK3J)/(1:M),type='l', col='blue')
# lines(cumsum(lamdamc)/(1:M),type='l', col='red')
# lines(cumsum(lamda1mc)/(1:M),type='l', col='purple')
# lines(cumsum(lamda2mc)/(1:M),type='l', col='darkgreen')
# lines(cumsum(lamda3mc)/(1:M),type='l', col='deeppink')
# legend("topright",text.width = 0.6,title="(c)", legend=c(expression(lambda[ML]), expression(lambda[LG]), expression(lambda[LV]), expression(paste(lambda[LJ]^"*")), expression(lambda[LJ]), expression(lambda[TKG]),expression(lambda[TKV]),expression(paste(lambda[TKJ]^"*")),expression(lambda[TKJ]),expression(lambda[MCG]),expression(lambda[MCV]),expression(paste(lambda[MCJ]^"*")),expression(lambda[MCJ])),
#        fill = c("darkred", "skyblue","black","orange","navy","brown","green","darkgrey","blue","red","purple","darkgreen","deeppink"),
#        lty=0.5:0.5,ncol=4,cex=0.8,adj=1 )



