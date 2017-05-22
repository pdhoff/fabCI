#' @title Hierarchical heteroscedastic model estimates
#' 
#' @description Estimate across-group heterogeneity of means and variances. 
#'
#' @details This function estimates 
#' parameters in a hierarchical model for normally distributed
#' groups, where the across-group model for means is normal 
#' and the across group model for variances is inverse-gamma.
#'  
#' @param y a numeric vector of data 
#' @param g a group membership vector, of the same length as y 
#' 
#' @return A vector (mu,t2,s20,nu0), where 
#' \enumerate{
#' \item mu is the mean of the group means 
#' \item t2 is the variance of the group means
#' \item the the distribution of group variances is inverse-gamma(nu0/2,nu0 s20/2)
#' }
#'
#' @author Peter Hoff
#'
#' @keywords htest 
#' 
#' @export
hhetmodel<-function(y,g)
{
  n<-c(table(g))
  ss<-(n-1)*c(tapply(y,g,var))
  ybar<-tapply(y,g,mean)

  n<-n[!is.na(ss)]
  ss<-ss[!is.na(ss)]
  ybar<-ybar[!is.na(ss)]

  ## -- get mmle for prior over s2
  mllab<-function(lab)
  {
    a<-exp(lab[1]) ; b<-exp(lab[2])
    ap<-a+(n-1)/2
    bp<-b+ss/2
    -sum( ( a*log(b)-lgamma(a) ) - ( ap*log(bp)-lgamma(ap) ) )
  }
  l<-(n-1)/ss ; lab0<-log(c( mean(l)^2/var(l) ,mean(l)/var(l) ))
  ab<-exp(optim(lab0, mllab)$par)

  ## -- get posterior modes for each s2 
  a<-ab[1] ; b<-ab[2]
  ap<-a+(n-1)/2
  bp<-b+ss/2
  s2<-bp/(ap+1)

  ## -- mmle of mu, t2 with plugin s2
  mllmut2<-function(mut2)
  {
    -sum(dnorm(ybar,mut2[1],sqrt(s2/n+exp(mut2[2])),log=TRUE) )
  }
  mut20<-c(mean(ybar),log(var(ybar)) )
  mut2<-optim(mut20,mllmut2)$par
  mut2[2]<-exp(mut2[2])

  psi<-c(mut2,ab[2]/ab[1],ab[1]*2) ; names(psi)<-c("mu","t2","s20","nu0") 
  psi
}


#' @title Multigroup FAB t-intervals
#' 
#' @description Computation of 1-alpha FAB t-intervals for heteroscedastic
#' multigroup data.
#'
#' @details For each group j, this function computes an estimate of 
#' the parameters in a hierarchical model for means and variances
#' from data other than group j, and uses this information to 
#' construct a FAB t-interval for group j. These intervals have 
#' 1-alpha frequentist coverage, assuming within-group normality. 
#' 
#' @param y a numeric vector of data
#' @param g a group membership vector, of the same length as y
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#'
#' @author Peter Hoff  
#' 
#' @keywords htest
#' 
#' @examples 
#' ## -- simulated data
#' p<-10 ; n<-10
#' y<-rnorm(n*p) ; g<-rep(1:p,n) 
#'
#' ## -- more interesting data takes longer 
#' # data(radon) ; y<-radon[,2] ; g<-radon[,1] 
#'
#' ## -- FAB t-intervals
#' FCI<-multifabCI(y,g) 
#'
#' ## -- UMAU t-intervals 
#' ybar<-tapply(y,g,mean) ; ssd<-tapply(y,g,sd) ; n<-table(g) 
#' qtn<-cbind( qt(.025,n-1),  qt(.975,n-1) ) 
#' UCI<-sweep(sweep(qtn,1,ssd/sqrt(n),"*"),1,ybar,"+") 
#' 
#' mean( (UCI[,2]-UCI[,1])/(FCI[,2]-FCI[,1]) , na.rm=TRUE)
#' 
#' @export
multifabCI<-function(y,g,alpha=.05)
{
  FCI<-NULL
  for(j in sort(unique(g)))
  {
    fci<-c(-Inf,Inf)
    yj<-y[g==j] 
    if(length(yj)>1)
    {
      psij<-hhetmodel(y[g!=j],g[g!=j])     
      fci<-fabtCI(y[g==j],psij,alpha=alpha)  
    }
    FCI<-rbind(FCI,fci)
  }
  rownames(FCI)<-sort(unique(g)) 
FCI
}


