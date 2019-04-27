DRQ <- function(X,y,tau){
  require(SparseM)
  require(quantreg)	 
  X <- as.matrix(x)
  p <- ncol(X)	 
  n <- length(y) 
  D=as.matrix.csr(X) 
  a <- (1-tau)*(t(X)%*%rep(1,n)) 
  fit=rq.fit.sfn(D,y,rhs=a)
  result = list(beta.hat=fit$coef)
  return(result)
}