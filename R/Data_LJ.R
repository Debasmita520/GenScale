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
y=c(71.46, 419.02, 284.64, 585.57, 456.60, 113.85, 187.85, 688.16, 662.66, 45.58,
    578.62, 756.70, 594.29, 166.49, 99.72, 707.36, 765.14, 187.13, 145.96, 350.70, 547.44, 116.99, 375.81,
    581.60, 119.86, 48.01, 200.16, 36.75, 244.53, 83.55)
  m=length(x)
  n=length(y)
  
    
    logLikfun=function(param){
      alp1=param[1]
      alp2=param[2]
      lamda=param[3]
      sum(dgen.exp(x,alpha=alp1,lambda=lamda,log=TRUE))+sum(dgen.exp(y,alpha=alp2,lambda=lamda,log=TRUE))
    }
    mle=maxLik(logLik=logLikfun,start=c(4,3,0.005),method = "BFGS")
    z1=coef(mle)[1]
    z2=coef(mle)[2]
    z3=coef(mle)[3]
    
    # z1=optim(c(alp1,alp2,lamda),LL,sx=sum(x),sy=sum(y),sp=sum(log(1-exp(-(lamda*x)))),sq=sum(log(1-exp(-(lamda*y)))),m=m,n=n)$par[1]  
    # z2=optim(c(alp1,alp2,lamda),LL,sx=sum(x),sy=sum(y),sp=sum(log(1-exp(-(lamda*x)))),sq=sum(log(1-exp(-(lamda*y)))),m=m,n=n)$par[2]
    # z3=optim(c(alp1,alp2,lamda),LL,sx=sum(x),sy=sum(y),sp=sum(log(1-exp(-(lamda*x)))),sq=sum(log(1-exp(-(lamda*y)))),m=m,n=n)$par[3]  
    
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Lindley's approximation..........#
    U=(m/(z1)^2)
    V=-sum((x*exp(-(z3*x)))/(1-exp(-(z3*x))))
    W=(n/(z2)^2)
    S=-sum((y*exp(-(z3*y)))/(1-exp(-(z3*y))))
    M1=(m/(z3)^2)+(z1-1)*sum(((x)^2*exp(-(z3*x)))/(1-exp(-(z3*x)))^2)+(n/(z3)^2)+(z2-1)*sum(((y)^2*exp(-(z3*y)))/(1-exp(-(z3*y)))^2)
    D=(U*((W*M1)-(S)^2))-((V)^2*W)
    
    
    sig11=(((W*M1)-(S)^2)/D)
    sig12=((V*S)/D)
    sig13=-((V*W)/D)
    sig22=(((U*M1)-(V)^2)/D)
    sig23=-((U*S)/D)
    sig33=((U*W)/D)
    
    
    L111=((2*m)/(z1)^3)
    L222=((2*n)/(z2)^3)
    L333=((2*m)/(z3)^3)+(z1-1)*sum(((x)^3*exp(-(z3*x))*(1+exp(-(z3*x))))/((1-exp(-(z3*x)))^3))+((2*n)/(z3)^3)+(z2-1)*sum(((y)^3*exp(-(z3*y))*(1+exp(-(z3*y))))/((1-exp(-(z3*y)))^3))
    L331=-sum(((x)^2*exp(-(z3*x)))/(1-exp(-(z3*x)))^2)
    L332=-sum(((y)^2*exp(-(z3*y)))/(1-exp(-(z3*y)))^2)
    
    
    A=((sig11*L111)+(sig33*L331))
    B=((sig22*L222)+(sig33*L332))
    C=(2*sig13*L331)+(2*sig23*L332)+(sig33*L333)
    
    
    q1=(0.5*((A*sig11)+(B*sig12)+(C*sig13)))
    q2=(0.5*((A*sig12)+(B*sig22)+(C*sig23)))
    q3=(0.5*((A*sig13)+(B*sig23)+(C*sig33)))
    
    A=(m/z3^2)*(1+((z1*(z1-1))/(z1-2))*(trigamma(1)-trigamma(z1-1)+(digamma(z1-1)-digamma(1))^2))-((m*z1)/(z3^2))*((trigamma(1)-trigamma(z1))+(digamma(z1)-digamma(1))^2)+(n/z3^2)*(1+((z2*(z2-1))/(z2-2))*(trigamma(1)-trigamma(z2-1)+(digamma(z2-1)-digamma(1))^2))-((n*z2)/(z3^2))*((trigamma(1)-trigamma(z2))+(digamma(z2)-digamma(1))^2)
    
    p11=sqrt((m/z1^2)*(((n*A)/z2^2)+((n^2)/(z3^2))*((z2/(z2-1))*(digamma(z2)-digamma(1))-(digamma(z2+1)-digamma(1)))^2)+((n*m^2)/(z2^2*z3^2))*(((z1/(z1-1))*(digamma(z1)-digamma(1))-(digamma(z1+1)-digamma(1)))^2))
    
    
    f =expression(p11)
    p_x = Deriv(f, "alp1")
    p_y = Deriv(f, "alp2")
    p_z= Deriv(f, "lamda")
    
    m1 = eval(p_x, list(alp1 = z1, alp2 = z2, lamda=z3))
    m2 = eval(p_y, list(alp1 = z1, alp2 = z2, lamda=z3))
    m3=eval(p_z, list(alp1 = z1, alp2 = z2, lamda=z3))
    
    #m1=0
    #m2=0
    #m3=-(c/z3)
    Lg1=z1+(m1*sig11)+(m2*sig12)+(m3*sig13)+q1
    Lg2=z2+(m1*sig12)+(m2*sig22)+(m3*sig23)+q2
    Lg3=z3+(m1*sig13)+(m2*sig23)+(m3*sig33)+q3
    LG1=mean(Lg1)
    LG2=mean(Lg2)
    LG3=mean(Lg3)
    #..........................................................................
    
    rule1L=log(LG1/LG2) +(LG1-LG2)*log(1-exp(-LG3*x))
    p11=length(which(rule1L>=0))/length(x)
    rule2L=log(LG1/LG2) +(LG1-LG2)*log(1-exp(-LG3*y))
    p22=length(which(rule2L<0))/length(y)
    epc=0.5*(p11+p22)
  
    cat("\n",format(round(epc,3),nsmall = 3))