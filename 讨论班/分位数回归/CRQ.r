CRQ<-function(X,Y,tau){
    n=length(Y)
    p=ncol(X)
    M=length(tau)
    library(Matrix)
    e.n=t(rep(1,n))
    e.M=t(rep(1,M))
    D.M=diag(M)%x%e.n
    f.obj <- c(rep(0,2*(p+M)),tau%*%D.M,(e.M-tau)%*%D.M)
    temp=cbind(matrix(1,M,1)%x%X,diag(M)%x%t(e.n))
    AE=cbind(temp,-temp,diag(M*n),-diag(M*n))
    AI=diag(2*M*n+2*(p+M))
    A=rbind(AE,AI)
    f.con <-A
    ce=rep("=",M*n)
    ci=rep(">=",2*M*n+2*(p+M))
    f.dir=c(ce,ci)
    f.rhs <- c(rep(Y,M),rep(0,(2*(p+M)+2*M*n)))
    out=lp ("min", f.obj, f.con, f.dir, f.rhs)
    theta.hat=out$solution
    beta.hat=theta.hat[1:(2*(p+M))]
    beta=beta.hat[1:(p+M)] - beta.hat[(p+M+1):(2*(p+M))]
    result=list(beta.hat=beta)
  return(result) 
}