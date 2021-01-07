#' @title UMAU regression coefficient intervals
#' 
#' @description Compute the usual t-intervals for the coefficients of a 
#' regression model
#'
#' @details This function computes the 'usual' uniformly most 
#' accurate unbiased confidence interval for each coefficient in 
#' a linear regression model. 
#' 
#' @param y a numeric vector of data 
#' @param X a design matrix
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' 
#' @return A matrix where each row corresponds to the interval and OLS 
#' estimate of a coefficient.  
#' 
#' @author Peter Hoff
#'
#' @export
umauregCI<-function(y,X,alpha=.05)
{
  fit<-summary(lm(y~ -1 + X))$coef
  bh<-fit[,1]
  se<-fit[,2]
  sweep( outer(  se, c(-1,0,1)*qt(1-alpha/2,nrow(X)-ncol(X)) ) ,1,bh,"+")
}



#' @title Empirical Bayes estimation of hyperparameters
#' 
#' @description Compute emprirical Bayes estimates of the 
#' error variance and distribution of the regression coefficients. 
#'
#' @details This function computes the adaptive FAB 
#' confidence interval for each coefficient in 
#' a linear regression model. 
#' 
#' @param y a numeric vector of data 
#' @param X a design matrix
#' @param emu (logical) estimate mean of coefficient (TRUE) or assume it is 
#' zero (FALSE)? 
#' @param dof degrees of freedom to use for the t-quantiles (the remainder 
#' go to adaptive estimation of the prior)
#' 
#' @return A list (s,sigma2,tau2,mu) where 
#' \enumerate{
#' \item s an estimate of the error standard deviation 
#' \item sigma2 an estimate of the error variance, independent of s
#' \item tau2 an estimate of the coefficient variance, independent of s
#' \item mu an estimate of the coefficient mean, independent of s
#' } 
#' 
#' @author Peter Hoff
#'
#' @export
ebayes_est<-function(y,X,emu=FALSE,dof=min(50,round(.5*(dim(X)[1]-dim(X)[2]))))
{
  sX<-svd(X) ; p<-sum(zapsmall(sX$d)>0)
  U<-sX$u[,1:p]
  H<-U%*%t(U) ; G<-MASS::Null(svd(X)$u)

  yh<-H%*%y
  e<-t(G)%*%y

  idx1<-seq(1,dof,length=dof)
  idx2<-seq(dof+1,length(e),length=length(e)-dof)
  s<-sqrt(mean(e[idx1]^2))
  sigma2<-mean(e[idx2]^2)
  mu<-emu*sum(yh)/sum(X)

  tau2<- (sum( (yh-mu)^2) - p*sigma2)/( sum( diag(crossprod(X)) ) - emu)
  tau2<-max(tau2,0)
  list(s=s,sigma2=sigma2,tau2=tau2,mu=mu)
}




#' @title FAB regression coefficient intervals
#' 
#' @description Compute the adaptive FAB t-intervals for the 
#' coefficients of a regression model. 
#'
#' @details This function computes the adaptive FAB 
#' confidence interval for each coefficient in 
#' a linear regression model. 
#' 
#' @param y a numeric vector of data 
#' @param X a design matrix
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' @param dof degrees of freedom to use for the t-quantiles (the remainder 
#' go to adaptive estimation of the prior) 
#' @param verbose logical, print progress or not
#' 
#' @return A matrix where each row corresponds to the interval and OLS 
#' estimate of a coefficient.
#' 
#' @author Peter Hoff
#'
#' @export
fabregCI<-function(y,X,alpha=.05,dof=min(50,round(.5*(dim(X)[1]-dim(X)[2]))),
          verbose=TRUE)
{
  p<-dim(X)[2]
  sX<-svd(X) ; H<-sX$u%*%diag(1/sX$d)%*%t(sX$v)
  beta_ols<-t(H)%*%y
  CI<-NULL
  for(j in 1:p)
  {
    G<-MASS::Null(H[,j])
    spsi<-ebayes_est(t(G)%*%y,t(G)%*%X,dof=dof)
    h<-sqrt(sum(H[,j]^2))
    CI<-rbind(CI, fabtzCI(beta_ols[j],spsi$s*h,dof,alpha=alpha,
        psi=list(mu=spsi$mu,tau2=spsi$tau2,sigma2=spsi$sigma2*h^2))) 
    if(verbose){ cat(j," of ",p," intervals computed\n") }
  }
  cbind(CI[,1],beta_ols,CI[,2])
}



#' @title z-optimal FAB t-interval 
#' 
#' @description Computation of a 1-alpha FAB t-interval using 
#' z-optimal spending function
#' 
#' @param y a numeric scalar, a normally distributed statistic
#' @param s a numeric scalar, the standard error of y
#' @param dof positive integer, degrees of freedom for s
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' @param psi a list of parameters for the spending function, including 
#' \enumerate{
#' \item mu, the prior expectation of E[y]
#' \item tau2, the prior variance of E[y]
#' \item sigma2 the variance of y
#' }
#'
#' @examples
#' n<-10 
#' y<-rnorm(n) 
#' fabtzCI(mean(y),sqrt(var(y)/n),n-1)  
#' t.test(y)$conf.int 
#' @export
fabtzCI<-function(y,s,dof,alpha=.05,psi=list(mu=0,tau2=1e5,sigma2=1)) 
{ 
  mu<-psi$mu ; tau2<-psi$tau2 ; sigma2<-psi$sigma2
  if(tau2<=0)
  {
    thetaL<-min(mu,y+s*qt(alpha,dof)) 
    thetaU<-max(mu,y+s*qt(1-alpha,dof)) 
  }

  if(tau2>0)
  {
    
    root<-function(theta)
    { 
      sfabz(theta,alpha=alpha,psi=psi) - pt( (y-theta)/s,dof )/alpha 
    }
    a<-b<-y+s*qt(1-alpha,dof) 
    #while(root(a)>0){ a<- a + s*qnorm(alpha)*.25 }
    #while(root(b)<0){ b<- b + s*qnorm(1-alpha)*.25 }  
    while(root(a)>0){ a<- a - s } 
    while(root(b)<0){ b<- b + s } 
    thetaU<-uniroot(root,c(a,b))$root 

    root<-function(theta)
    { 
      sfabz(theta,alpha=alpha,psi=psi) - (1- pt( (theta-y)/s,dof )/alpha) 
    }
    a<-b<-y+s*qt(alpha,dof) 
    #while(root(a)>0){ a<- a + s*qnorm(alpha)*.25 }
    #while(root(b)<0){ b<- b + s*qnorm(1-alpha)*.25 }
    while(root(a)>0){ a<- a - s } 
    while(root(b)<0){ b<- b + s } 
    thetaL<-uniroot(root,c(a,b))$root  
  }

c(thetaL,thetaU) 
}

#' @title Bayes-optimal spending function
#'
#' @description Compute Bayes optimal spending function
#'
#' @details This function computes the 
#' value of s that minimizes the acceptance probability of a 
#' biased level-alpha test for a normal population with 
#' known variance, under a specified  prior
#' predictive distribution.
#' 
#' @param theta value of theta being tested
#' @param psi a list of parameters for the spending function, including 
#' \enumerate{
#' \item mu, the prior expectation of E[y]
#' \item tau2, the prior variance of E[y]
#' \item sigma2 the variance of y
#' }
#' @param alpha level of test
#'
#' @author Peter Hoff 
#' 
#' @export
sfabz<-function(theta,psi,alpha=.05)
{ 
  mu<-psi$mu ; tau2<-psi$tau2 ; sigma2<-psi$sigma2 

  s<-1*(theta>mu)  
  if(tau2>0) 
  {
    igfun<-function(x,alpha)
    {
    gsmx <-function(s){ qnorm(alpha*s) - qnorm(alpha*(1-s)) - x }
    uniroot(gsmx,interval=c(0,1),maxiter=2000,tol=.Machine$double.eps^0.5)$root
    }
    s<-igfun( 2*sqrt(sigma2)*(theta-mu)/tau2,alpha) 
  } 
  s
}




