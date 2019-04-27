VSRQ<- function(X,y,tau,lambda){
  require(SparseM)
  require(quantreg)	 
  X <- as.matrix(x)
  p <- ncol(X)	 
  n <- length(y) 
  pn=2*p+2*n
  if(length(lambda)==1) lambda=rep(lemabda,p)
  e.n=rep(1,n)
  I.n=diag(n)
  f.obj=c(lambda,lambda,tau*e.n,(1-tau)*e.n)
  XE=cbind(X,-X,I.n,-I.n)
  f.con=rbind(XE,diag(pn))
  f.dir= c(rep("=", n), rep(">=", pn)) 
  f.rhs=c(y,rep(0,pn))
  out=lp ("min", f.obj, f.con, f.dir, f.rhs)
  theta.hat=out$solution
  beta.hat=theta.hat[1:p]-theta.hat[(p+1):(2*p)]
  result = list(beta.hat= beta.hat)
  return(result)
}