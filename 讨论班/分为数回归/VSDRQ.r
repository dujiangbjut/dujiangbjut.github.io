VSDRQ <- function(X,y,tau,lambda,acc=1e-3){
  require(SparseM)
  require(quantreg)
  if(length(lambda)==1) lambda=rep(lemabda,p)
  temp=X
  X=rbind(X,2*diag(lambda))
  X <- as.matrix(X)
  p <- ncol(X)
  n <- length(y) 
  D <- as.matrix.csr(X)
  Y <- c(y,rep(0,p))
  a <- (1-tau)*(t(temp)%*%matrix(1,n,1)) +lambda
  fit=rq.fit.sfn(D,Y,rhs=a)
  
  beta.hat=fit$coef
  beta.hat[abs(beta.hat)<=acc]=0
  result=list(beta.hat=beta.hat)
  
  return(result)
}
