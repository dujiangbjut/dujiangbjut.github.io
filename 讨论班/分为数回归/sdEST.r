sdEST<-function(X,y,tau){
  require('quantreg')
  p=ncol(X)
  n=length(y)
  h<- bandwidth.rq(tau, n, hs = TRUE)
  beta1<-IPRQ(X,y,tau+h)$beta.hat
  beta2<-IPRQ(X,y,tau-h)$beta.hat
  f.hat=2*h/(X%*%beta1-X%*%beta2 )
  Sigma1=1/n*t(X)%*%diag(as.vector(f.hat))%*%(X)
  Sigma2=1/n*tau*(1-tau)*t(X)%*%X
  sd=solve(Sigma1)%*%Sigma2%*%solve(Sigma1)
  res=list(beta.sd=diag(sd))
  return (res)
}