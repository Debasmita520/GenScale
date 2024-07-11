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
for (lamda in c(0.5,1,1.5,2,3,3.5,4,5,8,10)) {
  m=20
  n=20
  M=20000

  alp1=1
  alp2=0.5
  #lamda=1.5
  x=matrix(0,m,M)
  y=matrix(0,n,M)
  z1=array(0,M)
  z2=array(0,M)
  z3=array(0,M)
  biasc1=array(0,M)
  biasc2=array(0,M)
  biasc3=array(0,M)
  msec1=array(0,M)
  msec2=array(0,M)
  msec3=array(0,M)
  U=array(0,M)
  V=array(0,M)
  W=array(0,M)
  S=array(0,M)
  M1=array(0,M)
  D=array(0,M)
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
  A=array(0,M)
  B=array(0,M)
  C=array(0,M)
  q1=array(0,M)
  q2=array(0,M)
  q3=array(0,M)
  m1=array(0,M)
  m2=array(0,M)
  m3=array(0,M)
  Lg1=array(0,M)
  Lg2=array(0,M)
  Lg3=array(0,M)
  biasLg1=array(0,M)
  biasLg2=array(0,M)
  biasLg3=array(0,M)
  mseLg1=array(0,M)
  mseLg2=array(0,M)
  mseLg3=array(0,M)

  x1test=array(0,M)
  x2test=array(0,M)

  rule1L=array(0,M)
  rule2L=array(0,M)

  count1L=0
  count2L=0

  for (j in 1:M)
  {
    x[,j]=rgen.exp(m,alp1,lamda)
    y[,j]=rgen.exp(n,alp2,lamda)


    logLikfun=function(param){
      alp1=param[1]
      alp2=param[2]
      lamda=param[3]
      sum(dgen.exp(x[,j],alp1,lamda,log=TRUE))+sum(dgen.exp(y[,j],alp2,lamda,log=TRUE))
    }
    mle=maxLik(logLik=logLikfun,start=c(alp1,alp2,lamda),method = "BFGS")
    z1[j]=coef(mle)[1]
    z2[j]=coef(mle)[2]
    z3[j]=coef(mle)[3]

    biasc1[j]=(z1[j]-alp1)
    biasc2[j]=(z2[j]-alp2)
    biasc3[j]=(z3[j]-lamda)
    msec1[j]=(z1[j]-alp1)^2
    msec2[j]=(z2[j]-alp2)^2
    msec3[j]=(z3[j]-lamda)^2

    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Lindley's approximation..........#
    U[j]=(m/(z1[j])^2)
    V[j]=-sum((x[,j]*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j]))))
    W[j]=(n/(z2[j])^2)
    S[j]=-sum((y[,j]*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j]))))
    M1[j]=(m/(z3[j])^2)+(z1[j]-1)*sum(((x[,j])^2*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j])))^2)+(n/(z3[j])^2)+(z2[j]-1)*sum(((y[,j])^2*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j])))^2)
    D[j]=(U[j]*((W[j]*M1[j])-(S[j])^2))-((V[j])^2*W[j])


    sig11[j]=(((W[j]*M1[j])-(S[j])^2)/D[j])
    sig12[j]=((V[j]*S[j])/D[j])
    sig13[j]=-((V[j]*W[j])/D[j])
    sig22[j]=(((U[j]*M1[j])-(V[j])^2)/D[j])
    sig23[j]=-((U[j]*S[j])/D[j])
    sig33[j]=((U[j]*W[j])/D[j])


    L111[j]=((2*m)/(z1[j])^3)
    L222[j]=((2*n)/(z2[j])^3)
    L333[j]=((2*m)/(z3[j])^3)+(z1[j]-1)*sum(((x[,j])^3*exp(-(z3[j]*x[,j]))*(1+exp(-(z3[j]*x[,j]))))/((1-exp(-(z3[j]*x[,j])))^3))+((2*n)/(z3[j])^3)+(z2[j]-1)*sum(((y[,j])^3*exp(-(z3[j]*y[,j]))*(1+exp(-(z3[j]*y[,j]))))/((1-exp(-(z3[j]*y[,j])))^3))
    L331[j]=-sum(((x[,j])^2*exp(-(z3[j]*x[,j])))/(1-exp(-(z3[j]*x[,j])))^2)
    L332[j]=-sum(((y[,j])^2*exp(-(z3[j]*y[,j])))/(1-exp(-(z3[j]*y[,j])))^2)


    A[j]=((sig11[j]*L111[j])+(sig33[j]*L331[j]))
    B[j]=((sig22[j]*L222[j])+(sig33[j]*L332[j]))
    C[j]=(2*sig13[j]*L331[j])+(2*sig23[j]*L332[j])+(sig33[j]*L333[j])


    q1[j]=(0.5*((A[j]*sig11[j])+(B[j]*sig12[j])+(C[j]*sig13[j])))
    q2[j]=(0.5*((A[j]*sig12[j])+(B[j]*sig22[j])+(C[j]*sig23[j])))
    q3[j]=(0.5*((A[j]*sig13[j])+(B[j]*sig23[j])+(C[j]*sig33[j])))

    A=(n/lamda^2)*(1+((alp2*(alp2-1))/(alp2-2))*(trigamma(1)-trigamma(alp2-1)+(digamma(alp2-1)-digamma(1))^2))-((n*alp2)/(lamda^2))*((trigamma(1)-trigamma(alp2))+(digamma(alp2)-digamma(1))^2)

    p11=sqrt(m*(((n*A)/alp2)+((n^2)/(lamda^2))*((alp2/(alp2-1))*(digamma(alp2)-digamma(1))-(digamma(alp2+1)-digamma(1)))^2)+((n*m^2)/(alp2^2*lamda^2))*(0.639996))


    f =expression(p11)
    p_x = Deriv(f, "alp1")
    p_y = Deriv(f, "alp2")
    p_z= Deriv(f, "lamda")

    m1[j] = eval(p_x, list(alp1 = z1[j], alp2 = z2[j], lamda=z3[j]))
    m2[j] = eval(p_y, list(alp1 = z1[j], alp2 = z2[j], lamda=z3[j]))
    m3[j]=eval(p_z, list(alp1 = z1[j], alp2 = z2[j], lamda=z3[j]))

    Lg1[j]=z1[j]+(m1[j]*sig11[j])+(m2[j]*sig12[j])+(m3[j]*sig13[j])+q1[j]
    Lg2[j]=z2[j]+(m1[j]*sig12[j])+(m2[j]*sig22[j])+(m3[j]*sig23[j])+q2[j]
    Lg3[j]=z3[j]+(m1[j]*sig13[j])+(m2[j]*sig23[j])+(m3[j]*sig33[j])+q3[j]
    biasLg1[j]=Lg1[j]-alp1
    biasLg2[j]=Lg2[j]-alp2
    biasLg3[j]=Lg3[j]-lamda
    mseLg1[j]=(Lg1[j]-alp1)^2
    mseLg2[j]=(Lg2[j]-alp2)^2
    mseLg3[j]=(Lg3[j]-lamda)^2
    #..........................................................................
    x1test[j]=rgen.exp(1,alp1,lamda)
    x2test[j]=rgen.exp(1,alp2,lamda)

    rule1L[j]=log(Lg1[j]/Lg2[j]) +(Lg1[j]-Lg2[j])*log(1-exp(-Lg3[j]*x1test[j]))
    if(rule1L[j]>=0){count1L=count1L+1}
    rule2L[j]=log(Lg1[j]/Lg2[j]) +(Lg1[j]-Lg2[j])*log(1-exp(-Lg3[j]*x2test[j]))
    if(rule2L[j]<0){count2L=count2L+1}

  }

  p11L=count1L/M
  p22L=count2L/M

  epcL=(1/2)*(p11L+p22L)
  epmL=1-epcL

  ML1=round(mean(z1),2)
  ML2=round(mean(z2),2)
  ML3=round(mean(z3),2)
  biasML1=round(mean(biasc1),2)
  biasML2=round(mean(biasc2),2)
  biasML3=round(mean(biasc3),3)
  mseML1=round(mean(msec1),2)
  mseML2=round(mean(msec2),2)
  mseML3=round(mean(msec3),3)
  LG1=round(mean(Lg1),2)
  LG2=round(mean(Lg2),2)
  LG3=round(mean(Lg3),2)
  biasLG1=round(mean(biasLg1),3)
  biasLG2=round(mean(biasLg2),3)
  biasLG3=round(mean(biasLg3),3)
  mseLG1=round(mean(mseLg1),3)
  mseLG2=round(mean(mseLg2),3)
  mseLG3=round(mean(mseLg3),3)
  cat("\n",lamda,"\t",LG1,"\t",LG2,"\t",LG3,"\t",biasLG1,"\t",biasLG2,"\t",biasLG3,"\t",mseLG1,"\t",mseLG2,"\t",mseLG3,"\t",epmL)
  #cat("\n",ML1,"\t",ML2,"\t",ML3)
  #cat("\n",format(round(epmL,3),nsmall = 3))
}
