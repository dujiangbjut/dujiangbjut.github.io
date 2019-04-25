VSDCRQ2 <- function(X,y,tau,lambda){
  require(SparseM)
  require(quantreg)
  n=length(y)
  p=ncol(X)
  M=length(tau)
  library(Matrix)
  e.n=t(rep(1,n))
  e.M=t(rep(1,M))
  e.pM=t(rep(1,p+M))
  D.M=diag(M)%x%e.n
  temp=cbind(matrix(1,M,1)%x%X,diag(M)%x%t(e.n))
  X <- as.matrix(rbind(temp,2*diag(lambda)))
  D=as.matrix.csr(X) 
  Y=c(rep(y,M),rep(0,p+M))
  a <- t(temp)%*%t((e.M-tau)%*%(D.M) )+diag(lambda)%*%t(e.pM)
  fit=rq.fit.sfn(D,Y,rhs=a)
  beta.hat=fit$coefficients
  return(list(beta.hat=beta.hat))
}