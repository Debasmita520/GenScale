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


MCMC=function(alpha1int,alpha2int,lamdaint,sd1,sd2,sd3,X1,X2,alpha1,alpha2,lamda,n1,n2,c,N){
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
    
    w11=array(0,n1)
    w12=array(0,n1)
    w21=array(0,n2)
    w22=array(0,n2)
    for (k1 in 1:n1) {
      w11[k1]=dgen.exp(X1[k1],alpha1str,lamdastr)
      w12[k1]=dgen.exp(X1[k1],alpha1int,lamdaint)
      
    }
    for (k2 in 1:n2) {
      w21[k2]=dgen.exp(X2[k2],alpha2str,lamdastr)
      w22[k2]=dgen.exp(X2[k2],alpha2int,lamdaint)
    }
    p11=1/(alpha1str)
    p12=1/(alpha1int)
    p21=1/(alpha2str)
    p22=1/(alpha2int)
    p31=1/(lamdastr)
    p32=1/(lamdaint)
    
    
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

M=20000
n1=20
n2=20

alpha1=1.0
alpha2=0.5
lamda=1.0

bias1=array(0,M)
bias2=array(0,M)
bias3=array(0,M)
MSE1=array(0,M)
MSE2=array(0,M)
MSE3=array(0,M)
N=5000
x1=matrix(0,n1,M)
x2=matrix(0,n2,M)

alpha1mc=c()
alpha2mc=c()
lamdamc=c()
x1test=array(0,M)
x2test=array(0,M)

rule1mc=array(0,M)
rule2mc=array(0,M)

count1mc=0
count2mc=0

for (j in 1:M) {
  x1[,j]=rgen.exp(n1,alpha1,lamda)
  x2[,j]=rgen.exp(n2,alpha2,lamda)
  
  mcmc=MCMC(alpha1int=1.0,alpha2int=0.5,lamdaint=1.0,sd1=0.3,sd2=0.3,sd3=0.3,X1=x1[,j],X2=x2[,j],alpha1=alpha1,alpha2=alpha2,lamda=lamda,n1=n1,n2=n2,c=c,N=N)
  alpha1mc[j]=mean(mcmc$alpha1mc[N/2:N])
  bias1[j]=alpha1mc[j]-alpha1
  MSE1[j]=(alpha1mc[j]-alpha1)^2
  
  alpha2mc[j]=mean(mcmc$alpha2mc[N/2:N])
  bias2[j]=alpha2mc[j]-alpha2
  MSE2[j]=(alpha2mc[j]-alpha2)^2
  
  lamdamc[j]=mean(mcmc$lamdamc[N/2:N])
  bias3[j]=lamdamc[j]-lamda
  MSE3[j]=(lamdamc[j]-lamda)^2
  
  x1test[j]=rgen.exp(1,alpha1,lamda)
  x2test[j]=rgen.exp(1,alpha2,lamda)
  
  rule1mc[j]=log(alpha1mc[j]/alpha2mc[j]) +(alpha1mc[j]-alpha2mc[j])*log(1-exp(-lamdamc[j]*x1test[j]))
  if(rule1mc[j]>=0){count1mc=count1mc+1}
  rule2mc[j]=log(alpha1mc[j]/alpha2mc[j]) +(alpha1mc[j]-alpha2mc[j])*log(1-exp(-lamdamc[j]*x2test[j]))
  if(rule2mc[j]<0){count2mc=count2mc+1} 
  
}

p11mc=count1mc/M
p22mc=count2mc/M

epcmc=(1/2)*(p11mc+p22mc)
epmmc=1-epcmc


cat("\n",format(round(mean(alpha1mc),4),nsmall=4),"\t",format(round(mean(alpha2mc),4),nsmall=4),"\t",format(round(mean(lamdamc),4),nsmall=4))
cat("\n",format(round(mean(bias1),4),nsmall=4),"\t",format(round(mean(bias2),4),nsmall=4),"\t",format(round(mean(bias3),4),nsmall=4))
cat("\n",format(round(mean(MSE1),4),nsmall=4),"\t",format(round(mean(MSE2),4),nsmall=4),"\t",format(round(mean(MSE3),4),nsmall=4))
cat("\n",format(round(epmmc,4),nsmall = 4))