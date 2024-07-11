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
    p11=1
    p12=1
    p21=1
    p22=1
    p31=1/(lamdastr)^(c)
    p32=1/(lamdaint)^(c)
    
    
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
c=1
N=10000
x1=c(594.4, 202.75, 168.37, 574.86, 225.65, 76.38, 156.67, 127.81, 813.87, 562.39,
    468.47, 135.09, 72.24, 497.94, 355.56, 569.07, 640.48, 200.76, 550.42, 748.75, 489.66, 678.06, 457.71,
    106.73, 716.3, 42.66, 80.4, 339.22, 70.09, 193.42)
x2=c(71.46, 419.02, 284.64, 585.57, 456.60, 113.85, 187.85, 688.16, 662.66, 45.58,
    578.62, 756.70, 594.29, 166.49, 99.72, 707.36, 765.14, 187.13, 145.96, 350.70, 547.44, 116.99, 375.81,
    581.60, 119.86, 48.01, 200.16, 36.75, 244.53, 83.55)
n1=length(x1)
n2=length(x2)

mcmc=MCMC(alpha1int=1.84,alpha2int=1.61,lamdaint=0.004,sd1=0.3,sd2=0.3,sd3=0.3,X1=x1,X2=x2,n1=n1,n2=n2,c=c,N=N)
alpha1mc=mean(mcmc$alpha1mc[N/2:N])
alpha2mc=mean(mcmc$alpha2mc[N/2:N])
lamdamc=mean(mcmc$lamdamc[N/2:N])


rule1mc=log(alpha1mc/alpha2mc)+(alpha1mc-alpha2mc)*log(1-exp(-lamdamc*x1))
p11=length(which(rule1mc>=0))/length(x1)
rule2mc=log(alpha1mc/alpha2mc)+(alpha1mc-alpha2mc)*log(1-exp(-lamdamc*x2))
p22=length(which(rule2mc<0))/length(x2)
epc=0.5*(p11+p22)
epm=1-epc


cat("\n",format(round(mean(alpha1mc),4),nsmall=4),"\t",format(round(mean(alpha2mc),4),nsmall=4),"\t",format(round(mean(lamdamc),4),nsmall=4))
cat("\n",format(round(epc,4),nsmall = 4))

  
  
  
  
  
  
  
  
  
  
  
  
  