#' @title FAB t-interval with z-optimal w-function
#' 
#' @description Computation of a 1-alpha FAB t-interval using the w-function of
#' the optimal FAB z-interval
#'  
#' @details Internal function for multifabCIhom. It's similar to fabtCI, the 
#' difference being that fabtwzCI uses the w-function that leads to the 
#' optimal CI if the sampling variance were known. 
#' 
#' @param y a numeric vector with at least two non-missing values
#' @param s2 a variance estimate
#' @param df degrees of freedom corresponding to \code{s2} 
#' @param muw prior expectation of the mean
#' @param t2w prior variance of the  mean
#' @param s2w assumed population variance
#' @param alpha the type I error rate, so 1-alpha is the coverage rate  
#' 
#' @author Chaoyu Yu  
#' 
#' @keywords internal 
#' 
#' @examples
#' y<-rnorm(10)
#' fabtwzCI(y,10,9,1,5,10)
#' fabtwzCI(y,1/10,9,1,5,1/10)
#' fabtwzCI(y,10,1,9,5,10)
#' fabtwzCI(y,1/10,9,1,5,1/10)
#' 
#' @export
fabtwzCI<-function(y,s2,df,muw,t2w,s2w,alpha=.05)
{ 
  y = y[!is.na(y)]
  ybar = mean(y)
  n = length(y)
  s<-sqrt(s2)
  g<-function(w,alpha){ qnorm(alpha*w) - qnorm(alpha*(1-w)) }
  
  ubroot<-function(theta)
  {
    w<-max(min(  pt((ybar-theta)/(s/sqrt(n)),df)/alpha, 1) , 0)
    g(w , alpha) - 2*(theta-muw)*sqrt(s2w)/t2w
  }
  
  a<-b<- ybar + (s/sqrt(n))*qt(1-alpha,df)
  while(ubroot(a)<0){ a<- a - 1e-12  }
  while(ubroot(b)>0){ b<- b + s*qnorm(1-alpha)*.25 }
  thetaU<-uniroot(ubroot,c(a,b))$root
  
  
  lbroot<-function(theta)
  {
    w<-max(min(  1-pt((theta-ybar)/(s/sqrt(n)),df)/alpha, 1), 0)
    g(w, alpha) - 2*(theta-muw)*sqrt(s2w)/t2w
  }
  
  a<-b<- ybar + (s/sqrt(n))*qt(alpha,df)
  while(lbroot(a)<0){ a<- a + s*qnorm(alpha)*.25 }
  while(lbroot(b)>0){ b<- b + 1e-12  }
  thetaL<-uniroot(lbroot,c(a,b))$root
  
  c(thetaL,thetaU)
}



#' @title Hierarchical homoscedastic model estimates 
#' 
#' @description Estimate across-group heterogeneity of means
#'
#' @details This function estimates 
#' parameters in a hierarchical model for normally distributed
#' groups, where the across-group model for means is normal 
#' and the variance is the same across groups.
#'  
#' @param y a numeric vector of data 
#' @param g a group membership vector, of the same length as y 
#' @param group the index of the group 
#' @param p1 number of groups used to pool sample variance
#' 
#' @return A vector (s2,df,muw,t2w,s2w), where 
#' \enumerate{
#' \item s2 is the pooled variance 
#' \item df is the degree of freedom of the t-quantiles
#' \item muw is the estimate mean of the group means 
#' \item t2w is the estimate variance of the group means
#' \item s2w is the estimate within-group variance
#' }
#'
#' @author Chaoyu Yu
#'
#' @keywords htest 
#' 
#' @export
hhommodel<-function(y,g,group,p1)
{
  ## -- estimate for tilde S2
  hmmodel1<-function(y,g){
    n<-c(table(g))   
    no<-unique(g)
    s2<-sum((n-1)*c(tapply(y,g,var)))/sum(n-1)
    df = sum(n-1)
    return(c(s2,df))
  }
  
  ## -- estimate for gamma
  hmmodel2<-function(y,g)
  {
    n<-c(table(g))
    ss<-(n-1)*c(tapply(y,g,var))
    ybar<-tapply(y,g,mean)
    
    n<-n[!is.na(ss)]
    ss<-ss[!is.na(ss)]
    ybar<-ybar[!is.na(ss)]
    
    s2<-mean(tapply(y,g,var))
    
    ## -- mmle of mu, t2 with plugin s2
    mllmut2<-function(mut2)
    {
      -sum(dnorm(ybar,mut2[1],sqrt(s2/n+exp(mut2[2])),log=TRUE) )
    }
    mut20<-c(mean(ybar),log(var(ybar)) )
    mut2<-optim(mut20,mllmut2)$par
    mut2[2]<-exp(mut2[2])
    
    psi<-c(mut2,s2)  
    psi
  } 
  
  nos = unique(g)
  for(i in nos){
    if(length(g[g==i])<2){
      y=y[!(g==i)]
      g=g[!(g==i)]
    }
  }
  
  nos = unique(g)
  nos = nos[order(table(g))]
  nos1 = nos[nos!=group]
  pp = length(nos1)
  nos1 = nos1[c(seq(1,pp,2),seq(2,pp,2))]
  ind = nos1[1:(p1-1)]
  ind = c(ind,group)
  
  y1 = y[g %in% ind]
  g1 = g[g %in% ind]
  
  
  y2 = y[!(g %in% ind)]
  g2 = g[!(g %in% ind)]
  
  return(c(hmmodel1(y1,g1),hmmodel2(y2,g2)))
  
}


#' @title Multigroup FAB t-intervals for the homoscedastic model
#' 
#' @description Computation of 1-alpha FAB t-intervals for homoscedastic 
#' multigroup data. 
#'
#' @details For each group j, this function computes an estimate of 
#' the parameters in a hierarchical model for means 
#' using data from other groups, and uses this information to 
#' construct a FAB t-interval for group j. These intervals have 
#' 1-alpha frequentist coverage, assuming within-group normality and that
#' the within group variance is the same across groups. 
#' 
#' @param y a numeric vector of data
#' @param g a group membership vector, of the same length as y
#' @param alpha the type I error rate, so 1-alpha is the coverage rate 
#' @param prop the proportion of groups to obtain the sample variance estimate
#'
#' @author Chaoyu Yu 
#' 
#' @keywords htest
#' 
#' @examples 
#' ## -- simulate the data
#' mu = 0; sigma2 = 10; tau2 = 1; p =100; 
#' theta = rnorm(p,mu,sqrt(tau2))
#' ns = round(runif(p,2,18))
#' Y=c()
#' for(i in 1:p){
#'  d2 = rnorm(ns[i],theta[i],sqrt(sigma2))
#'  d1 = rep(i,ns[i])
#'  d = cbind(d1,d2)
#'  Y = rbind(Y,d)}
#' y = Y[,2]
#' g = Y[,1]
#' 
#' ## -- FAB t-intervals
#' FCI = multifabCIhom(y,g)  
#' 
#' ## -- UMAU t-intervals 
#' ybar<-tapply(y,g,mean) ; ssd<-tapply(y,g,sd) ; n<-table(g) 
#' qtn<-cbind( qt(.025,n-1),  qt(.975,n-1) ) 
#' UCI<-sweep(sweep(qtn,1,ssd/sqrt(n),"*"),1,ybar,"+") 
#' 
#' mean( (UCI[,2]-UCI[,1])/(FCI[,2]-FCI[,1]) , na.rm=TRUE)
#' 
#' 
#' @export
multifabCIhom<-function(y,g,alpha=.05,prop=0.5)
{
  if(prop<=0 | prop>=1){ stop("prop must be between 0 to 1") }
  ns<-c(table(g))
  qq = sum(ns>1)
  if(qq<4){ stop("too few groups") }
  p1 = max(round(prop*qq),2)
  if(qq-p1<2){p1 = p1-1}
  FCI<-NULL
  nos<-sort(unique(g))
  
  for(j in nos)
  {
    fci<-c(-Inf,Inf)
    yj<-y[g==j] 
    if(length(yj)>1)
    {
      psij<-hhommodel(y,g,j,p1)     
      fci<-fabtwzCI(y[g==j],psij[1],psij[2],psij[3],psij[4],psij[5],alpha=alpha)  
    }
    FCI<-rbind(FCI,fci)
  }
  rownames(FCI)<-sort(unique(g)) 
  FCI
}

