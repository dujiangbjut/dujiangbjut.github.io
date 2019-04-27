DCRQ <- function(X,y,tau){
  require(SparseM)
  require(quantreg)
  n=length(y)
  p=ncol(X)
  M=length(tau)
  library(Matrix)
  e.n=t(rep(1,n))
  e.M=t(rep(1,M))
  D.M=diag(M)%x%e.n
  temp=cbind(matrix(1,M,1)%x%X,diag(M)%x%t(e.n))
  X <- as.matrix(temp)
  D=as.matrix.csr(X) 
  Y=rep(y,M)
  a <- t(temp)%*%t((e.M-tau)%*%(D.M) )
  fit=rq.fit.sfn(D,Y,rhs=a)
  beta.hat=fit$coefficients
  return(list(beta.hat=beta.hat))
}