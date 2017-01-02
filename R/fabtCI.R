#' @title FAB t-interval
#' 
#' @description Computation of a 1-alpha FAB t-interval. 
#'  
#' @details A FAB interval is the "frequentist" interval procedure
#' that is Bayes optimal: It  minimizes the prior expected
#' interval width among all interval procedures with 
#' exact 1-alpha frequentist coverage. This function computes 
#' the FAB t-interval for the mean of a normal population with an 
#' unknown variance, given a user-specified prior distribution 
#' determined by \code{psi}. The prior is that the population mean 
#' and variance are independently distributed as normal and 
#' inverse-gamma random variables. 
#' Referring to the elements of \code{psi}
#' as mu, t2, s20, nu0, the prior is determined as follows:
#' \enumerate{
#' \item mu is the prior expectation of the mean
#' \item t2 is the prior variance of the  mean
#' \item the population variance is inverse-gamma(nu0/2,nu0 s20/2)
#' }
#' 
#' @param y a numeric vector with at least two non-missing values
#' @param psi a length-four vector of hyperparameters for the prior
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' 
#' @author Peter Hoff  
#' 
#' @keywords htest 
#' 
#' @examples
#' y<-rnorm(10)
#' fabtCI(y,c(0,10,1,5)) 
#' fabtCI(y,c(0,1/10,1,5))
#' fabtCI(y,c(2,10,1,5))
#' fabtCI(y,c(0,1/10,1,5)) 
#' 
#' @export
fabtCI<-function(y,psi=c(0,100,1,2),alpha=.05)
{
  if(!is.numeric(y)){ stop("y must be numeric") }     
  if(!is.numeric(psi)){ stop("psi must be numeric") } 
  if(!is.numeric(alpha)){ stop("alpha must be numeric") }  
  if(alpha<=0 | alpha>=1 ){ stop("alpha must be between 0 and 1") }  
  if(length(psi)!=4 | any(psi[2:4]<0)) { stop("psi not specified correctly") }

  ybar<-mean(y,na.rm=TRUE) ; s<-sd(y,na.rm=TRUE) ; n<-sum(!is.na(y)) 
  mu<-psi[1] ; t2<-psi[2] ; s20<-psi[3] ; nu0<-psi[4] 

  thetaL<-thetaU<-NA  
  if(n>1)
  {
    ubroot<-function(theta)
    {  
      w<-wfabt(theta,mu,t2,s20,nu0,n,alpha)
      ybar + s*qt(1-alpha*w,n-1)/sqrt(n) - theta
    } 
    a<-b<-ybar + .99*(s/sqrt(n))*qt(1-alpha,n-1)
    while(ubroot(b)>0){ b<- b + (s/sqrt(n))*qnorm(1-alpha)*n/(n+4) }
    thetaU<-uniroot(ubroot,c(a,b))$root

    lbroot<-function(theta)
    {
      w<-wfabt(theta,mu,t2,s20,nu0,n,alpha)
      ybar + s*qt(alpha*(1-w),n-1)/sqrt(n) - theta
    }
    a<-b<-ybar + .99*(s/sqrt(n))*qt(alpha,n-1)
    while(lbroot(a)<0){ a<- a + (s/sqrt(n))*qnorm(alpha)*n/(n+4) }
    thetaL<-uniroot(lbroot,c(a,b))$root
  } 
 
  c(thetaL,thetaU)
}


#' @title Prior predictive acceptance probability 
#' 
#' @description Compute prior probability of falling into an acceptance region 
#'
#' @details Internal function for fabtCI. This function computes the prior 
#' probability that the value theta will fall in the acceptance 
#' region of a biased level-alpha test, under a specified  prior 
#' predictive distribution.  
#'  
#' @param w the weight in [0,1] determining acceptance region endpoints
#' @param theta value of theta being tested
#' @param mu prior mean of theta 
#' @param t2 prior variance of theta
#' @param s20 prior parameter for population variance
#' @param nu0 prior parameter for population variance 
#' @param n sample size 
#' @param alpha level of test  
#'
#' @author Peter Hoff
#' 
#' @keywords internal 
#'
#' @export
pppaccept<-function(w,theta,mu,t2,s20,nu0,n,alpha)
{
  a<-nu0/2 ; b<-nu0*s20/2

  f<-function(s2)
  {
    c<-sqrt(s2/n)/sqrt(s2/n+t2)
    ncp<- c*(mu - theta)/(sqrt(s2/n))
    suppressWarnings(
    pacc<- pt( c*qt(1-alpha*(1-w),n-1),n-1,ncp ) -
           pt( c*qt(alpha*w,n-1),n-1,ncp )
                    )

    ds2<-exp( a*log(b) - lgamma(a) - (a+1)*log(s2) - b/s2 )

    pacc*ds2
  }

  int<-try(integrate(f,lower=0,upper=Inf)$val,silent=TRUE)
  if(is.character(int)){int<-1}
  int
}

#' @title Bayes-optimal w-function
#'
#' @description Compute Bayes optimal w-function
#'
#' @details Internal function for fabtCI. This function computes the 
#' value of w that minimizes the acceptance probability of a 
#' biased level-alpha test under a specified  prior
#' predictive distribution.
#' 
#' @param theta value of theta being tested
#' @param mu prior mean of theta
#' @param t2 prior variance of theta
#' @param s20 prior parameter for population variance
#' @param nu0 prior parameter for population variance
#' @param n sample size
#' @param alpha level of test
#'
#' @author Peter Hoff 
#' 
#' @keywords internal 
#' 
#' @export
wfabt<-function(theta,mu,t2,s20,nu0,n,alpha)
{
  optimize( pppaccept, lower=0,upper=1,
            theta=theta,mu=mu,t2=t2,s20=s20,nu0=nu0,n=n,alpha=alpha)$min
}



