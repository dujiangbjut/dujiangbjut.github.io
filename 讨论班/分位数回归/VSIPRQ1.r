VSIPRQ1<-function(X,Y,tau,lambda){
  X=as.matrix(X)
  n=length(Y)
  p=ncol(X) 
  library(lpSolve) 
  e.n=rep(1,n)
  I.n=diag(n)
  if(length(lambda)==1) lambda=rep(lambda,p)
  AE1=cbind(X,-X,I.n,-I.n,matrix(0,n,2*p))
  AE2=cbind(diag(lambda),-diag(lambda), matrix(0,p,2*n),diag(p),-diag(p))
  AE=rbind(AE1,AE2) 
  AI=diag(2*n+4*p)
  A=rbind(AE,AI)
  f.obj <- c(rep(0,2*p),tau*e.n,(1-tau)*e.n,rep(1,2*p)) 
  f.con <-A
  ce=rep("=",n+p)
  ci=rep(">=",2*n+4*p)
  f.dir=c(ce,ci)
  f.rhs <- c(Y,rep(0,5*p+2*n))
  out=lp ("min", f.obj, f.con, f.dir, f.rhs)
  theta.hat=out$solution
  beta.hat= theta.hat[1:p]-theta.hat[(p+1):(2*p)]
  result=list(beta.hat=beta.hat)
  return(result)
}




