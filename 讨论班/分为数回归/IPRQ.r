IPRQ<-function(X,Y,tau){
  X=as.matrix(X)
  n=length(Y)
  p=ncol(X)
  library(lpSolve)
  e.n=rep(1,n)
  I.n=diag(n)
  n.t=2*p+2*n
  AE=cbind(X,-X,I.n,-I.n)
  AI=diag(n.t)
  A=rbind(AE,AI)
  f.obj <- c(rep(0,2*p),tau*e.n,(1-tau)*e.n)
  f.con <-A
  ce=rep("=",n)
  ci=rep(">=",n.t)
  f.dir=c(ce,ci)
  f.rhs <- c(Y,rep(0,n.t))
  out=lp ("min", f.obj, f.con, f.dir, f.rhs)
  theta.hat=out$solution
  beta.hat= theta.hat[1:(p)]-theta.hat[(p+1):(2*p)]
  result=list(beta.hat=beta.hat)
  return(result)
}