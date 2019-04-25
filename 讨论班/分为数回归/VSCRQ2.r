VSCRQ2<-function(X,Y,tau,lambda){
  n=length(Y)
  p=ncol(X)
  M=length(tau)
  library(Matrix)
  e.n=t(rep(1,n))
  e.M=t(rep(1,M))
  D.M=diag(M)%x%e.n
  f.obj <- c(rep(0,2*(p+M)),tau%*%D.M,(e.M-tau)%*%D.M, rep(1,2*(p+M)) )
  temp=cbind(matrix(1,M,1)%x%X,diag(M)%x%t(e.n))
  AE1=cbind(temp,-temp,diag(M*n),-diag(M*n),matrix(0,M*n,2*(p+M)) )
  AE2=cbind(diag(lambda),-diag(lambda),matrix(0,p+M,2*M*n),diag(p+M),-diag(p+M))
  AI=diag(2*M*n+4*(p+M))
  A=rbind(AE1,AE2,AI)
  f.con <-A
  ce=rep("=",M*n+p+M)
  ci=rep(">=",2*M*n+4*(p+M))
  f.dir=c(ce,ci)
  f.rhs <- c(rep(Y,M),rep(0,(p+M)),rep(0,2*M*n+4*(p+M)))
  out=lp ("min", f.obj, f.con, f.dir, f.rhs)
  theta.hat=out$solution
  beta.hat=theta.hat[1:(2*(p+M))]
  beta=beta.hat[1:(p+M)] - beta.hat[(p+M+1):(2*(p+M))]
  result=list(beta.hat=beta)
  return(result) 
}