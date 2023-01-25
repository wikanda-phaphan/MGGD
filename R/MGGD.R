#' @title The Mixture Generalized Gamma Distribution
#'
#' @description Density, distribution function, random generation, and parameter estimation function for the mixture generalized gamma distribution as Suksaengrakcharoen and Bodhisuwan  (2014) introduced.
#'
#' @param X vector of data
#' @param n a number of observations
#' @param lambda a value of the parameter lambda
#' @param beta a value of the parameter beta
#' @param alpha a value of the parameter alpha
#' @param p a value of the parameter p
#'
#' @return ?
#'
#' @examples
#' lambda <-0.7
#' beta <-1
#' alpha <-3
#' p <- 0.5
#' n <- 30
#' X <- rMGG(n,lambda,beta,alpha,p)
#' X
#' dMGG(X,lambda,beta,alpha,p)
#' pMGG(X,lambda,beta,alpha,p)
#' EM(X,n,lambda,beta,alpha,p)

#'
###############The  Mixture Generalized Gamma Distribution#####################
fGG=function(X,lambda,beta,alpha)
{
  ((lambda*beta)/(gamma(alpha)))*(lambda*X)^((alpha*beta)-1)*
    (exp(-(lambda*X)^(beta)))
}
fLBGG=function(X,lambda,beta,alpha)
{
  (lambda*beta)/(gamma(alpha+(1/beta)))*((lambda*X)^(alpha*beta))*
    (exp(-(lambda*X)^(beta)))
}
fMGG=function(X,lambda,beta,alpha,p){
  p*fGG(X,lambda,beta,alpha)+(1-p)*fLBGG(X,lambda,beta,alpha)
}
###########pdf
dMGG=function(X,lambda,beta,alpha,p){
  p*fGG(X,lambda,beta,alpha)+(1-p)*fLBGG(X,lambda,beta,alpha)
}
#dMGG(X,lambda,beta,alpha,p)
#####cdf
#install.packages('zipfR')
#library(zipfR)
pMGG=function(X,lambda,beta,alpha,p){
  1-(p*Rgamma(alpha, (lambda*X)^beta)/gamma(alpha))-((1-p)*Rgamma(alpha+(1/beta),(lambda*X)^beta)/gamma(alpha+(1/beta)))
}
#pMGG(X,lambda,beta,alpha,p)
######### MGG-random numbers generation procedure ### 1 time
fM=function(X,lambda,beta,alpha,p){
  fMGG(X,lambda,beta,alpha,p)/dlnorm(X, meanlog = lambda, sdlog = beta, log = FALSE)
}

rMGG=function(n,lambda,beta,alpha,p){
  M <-optimise(f=function(X){fM(X,lambda,beta,alpha,p)},interval=c(0,45),maximum=T)$objective
  X=NULL
  while(length(X)<n){
    y=rlnorm(n*M, meanlog = lambda, sdlog = beta)
    u=runif(n*M,0,M)
    X=c(X,y[u<fM(y,lambda,beta,alpha,p)])
  }
  X=X[1:n]
  MGG<-matrix(X, ncol=1)
  return(MGG)
}
#rMGG(n,lambda,beta,alpha,p)
#X<-rMGG(n,lambda,beta,alpha,p)

EM=function(X,n,lambda,beta,alpha,p){
  flcomplete1<-function(beta)
  {
    
    (n*(ak-1))*log(gamma(beta[3]+(1/beta[2])))+(n*((beta[3]*beta[2])-1-ak)*log(beta[1]))+(n*log(beta[2]))+sum(n*((beta[3]*beta[2])-1)*log(X))-sum((beta[1]*X)^beta[2])-(n*ak*log(gamma(beta[3]))) 
    
  } 
  ####################### find initial value#################################
  Xbar <-apply(X,2,mean)
  alpha.old <- 0.5/(log(Xbar)-mean(log(X)))
  beta.old <- Xbar/alpha.old
  p.old <- 0.5
  #lambda.old <- 0.5
  lambda.old <- Xbar/(apply(X,2,var))
  
  cond=1
  cond2=1
  cond3=1
  cond4=1
  #while(cond>0.000001){
  while((cond>0.01)&(cond2>0.00001)&(cond2>0.01)&(cond2>0.0001)){
    Ez=(p.old)*fGG(X,lambda.old,beta.old,alpha.old)/((1-p.old)*fLBGG(X,lambda.old,beta.old,alpha.old)+p.old*fGG(X,lambda.old,beta.old,alpha.old))
    ak=mean(Ez)
    p.new <- ak
    jib <-nlminb(c(lambda.old,beta.old,alpha.old), flcomplete1,scale = 100)
    #jib <-nlminb(c(lambda.old,beta.old,alpha.old,p.old), flcomplete2,scale = 10)
    lambda.new<-jib$par[1]
    beta.new<-jib$par[2]
    alpha.new<-jib$par[3]
    #p.new <-jib$par[4]
    
    #cond=abs(p.new-p.old)+abs(lambda.new-lambda.old)+abs(beta.new-beta.old)+abs(alpha.new-alpha.old)
    cond=abs(p.new-p.old)
    cond2=abs(lambda.new-lambda.old)
    cond3=abs(beta.new-beta.old)
    cond4=abs(alpha.new-alpha.old)
    
    p.old=p.new
    lambda.old=lambda.new
    beta.old=beta.new
    alpha.old=alpha.new
  }
  return(list(p.hat=p.new,lambda.hat=lambda.new,beta.hat=beta.new,alpha.hat=alpha.new))
}
#EM(X,n,lambda,beta,alpha,p)

