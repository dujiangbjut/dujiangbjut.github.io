  ESD<-function(x,y,tau,B=500){
    n=length(y)
    p=ncol(x)
    beta.hat=rq(y~x+0,tau)$coef
    BETA=matrix(0,B,p)
    for(i in 1:B){
      W=rexp(n)
      BETA[i,]=rq(y~x+0,tau,weight=W)$coef-beta.hat
    }
    res=list(beta.sd=sqrt(diag(cov(BETA) )))
  return(res) 
}