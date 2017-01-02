#' @title FAB z-interval
#' 
#' @description Computation of a 1-alpha FAB z-interval. 
#'  
#' @details A FAB interval is the "frequentist" interval procedure
#' that is Bayes optimal: It  minimizes the prior expected
#' interval width among all interval procedures with 
#' exact 1-alpha frequentist coverage. This function computes 
#' the FAB z-interval for the mean of a normal population with an 
#' known variance, given a user-specified prior distribution 
#' determined by \code{psi}. The prior is that the population mean 
#' is normally distributed. 
#' Referring to the elements of \code{psi}
#' as mu, t2, s2, the prior and population variance are
#' determined as follows:
#' \enumerate{
#' \item mu is the prior expectation of the mean
#' \item t2 is the prior variance of the  mean
#' \item s2 is the population variance
#' }
#' 
#' @param y a numeric scalar
#' @param mu a numeric scalar
#' @param t2 a positive numeric scalar
#' @param s2 a positive numeric scalar
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' 
#' @author Peter Hoff  
#' 
#' @keywords htest 
#' 
#' @examples
#' y<-0
#' fabzCI(y,0,10,1)
#' fabzCI(y,0,1/10,1)
#' fabzCI(y,2,10,1)
#' fabzCI(y,0,1/10,1) 
#' 
#' @export
fabzCI<-function(y,mu,t2,s2,alpha=.05)
{ 
  if(!is.numeric(y) | length(y)>1){ stop("y must be a numeric scalar") }
  if(!is.numeric(mu) | length(mu)>1){ stop("mu must be a numeric scalar") }
  if(!is.numeric(t2) | length(t2)>1 | t2<0){stop("t2 must be a positive numeric scalar") }
  if(!is.numeric(s2) | length(s2)>1 | s2<0){stop("s2 must be a positive numeric scalar") }
  if(!is.numeric(alpha) | length(alpha)>1){ stop("alpha must be a numeric scalar") }
  if(alpha<=0 | alpha>=1 ){ stop("alpha must be between 0 and 1") } 

  s<-sqrt(s2) 

  ubroot<-function(theta)
  {  
    (y+s*qnorm(min(1,1-alpha+pnorm((y-theta)/s)))+ mu*2*s2/t2)/(1+2*s2/t2) - 
    theta
  }
  a<-b<- y + s*qnorm(1-alpha)  
  while(ubroot(a)<0){ a<- a - 1e-12  }
  while(ubroot(b)>0){ b<- b + s*qnorm(1-alpha)*.25 }
  thetaU<-uniroot(ubroot,c(a,b))$root


  lbroot<-function(theta)
  {  
    (y+s*qnorm(max(0,alpha-pnorm((theta-y)/s)))+ mu*2*s2/t2)/(1+2*s2/t2) - 
     theta
  }
  a<-b<- y + s*qnorm(alpha) 
  while(lbroot(a)<0){ a<- a + s*qnorm(alpha)*.25 }
  while(lbroot(b)>0){ b<- b + 1e-12  }
  thetaL<-uniroot(lbroot,c(a,b))$root

  c(thetaL,thetaU)
}

#' @title Bayes-optimal w-function
#'
#' @description Compute Bayes optimal w-function
#'
#' @details This function computes the 
#' value of w that minimizes the acceptance probability of a 
#' biased level-alpha test for a normal population with 
#' known variance, under a specified  prior
#' predictive distribution.
#' 
#' @param theta value of theta being tested
#' @param mu prior mean of theta
#' @param t2 prior variance of theta
#' @param s2 population variance
#' @param alpha level of test
#'
#' @author Peter Hoff 
#' 
#' @keywords internal 
#' 
#' @export
wfabz<-function(theta,mu,t2,s2,alpha=.05)
{
  igfun<-function(x,alpha)
  {
    gwmx <-function(w){ qnorm(alpha*w) - qnorm(alpha*(1-w)) - x }
    uniroot(gwmx,interval=c(0,1),maxiter=2000,tol=.Machine$double.eps^0.5)$root
  }
  igfun( 2*sqrt(s2)*(theta-mu)/t2,alpha)
}



