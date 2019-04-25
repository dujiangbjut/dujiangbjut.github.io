MMCRQ=function(x,y,tau,delta=1e-3){
  n=length(y)
  K=length(tau)
  x=as.matrix(x)
  p=ncol(x)
  e.n=t(rep(1,n))
  x=cbind(matrix(1,K,1)%x%x,diag(K)%x%t(e.n))
  y=rep(y,K)
  e.new=e.old=rep(0,n*K)
  tau=matrix(tau%x%rep(1,n))
  beta.new=beta.old=rep(0,p+K)
  loop=0
  Loop=1000
  esp=1e-3
  error=1
  while(error>esp&&loop<Loop){
    loop=loop+1
    beta.old=beta.new     
    e.hat=y-x%*%beta.old      
    W=diag(as.vector(1/(delta+abs(e.hat))))
    beta.new=beta.old+solve(t(x)%*%W%*%x)%*%t(x)%*%(2*tau-1+e.hat/(abs(e.hat)+delta))
    #beta.new=solve(t(x)%*%W%*%x)%*%t(x)%*%(2*tau-1+y/(abs(e.hat)+delta))
    error=sum(abs(beta.new-beta.old))
  }
  ind=sum(abs(beta.new-beta.old))#1*(loop<Loop)
  ind=loop
  loss=0 
  temp=y-x%*%beta.new
  loss=sum(temp*(tau-1*(temp<0)))
  tmp.out = list(beta.hat=beta.new, loss=loss,ind=ind)    
  return(tmp.out)
}


