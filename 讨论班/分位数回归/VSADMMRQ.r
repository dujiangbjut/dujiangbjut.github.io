VSADMMRQ <- function(X,y,tau,lambda,iter=200,esp=1e-03){
  n=length(y)
  e.new=e.old=rep(0,n)
  X=as.matrix(X)
  p=ncol(X)
  beta.new=beta.old=rep(0,p)
  u.new=u.old=rep(0,n)
  step=0
  error=1
  while(step<iter&error>esp){ 
    beta.old=beta.new
    e.old=e.new
    u.old=u.new
    step=step+1
    c=y-X%*%beta.old+u.old/1.2
    c1=c-(2*tau-1)/1.2
    e.new=pmax(c1-1/1.2,0)-pmax(-c1-1/1.2,0)
    beta.new=solve(t(X)%*%X)%*%t(X)%*%(y-e.new+u.old/1.2)+lambda*(pmax(beta.old,0)-pmax(-beta.old,0))
    temp=y-e.new-X%*%beta.new
    u.new=u.old+temp*1.2
    temp1=y-e.new-X%*%beta.new
    temp2=1.2*t(X)%*%(e.new-e.old)
    error1=max(sqrt(sum((X%*%beta.new)^2)),sqrt(sum(e.new^2)),sqrt(sum(y^2)))
    error2=sqrt(p)*0.01+esp*sqrt(sum((X%*%beta.new)^2))
    if(sqrt(sum(temp1^2))<sqrt(n)*0.01+esp*error1&sqrt(sum(temp2^2))<error2) step=iter+1
    error=sqrt(sum((beta.new-beta.old)^2))
  }
  e.hat=y-X%*%beta.new
  loss=sum(e.hat*(tau-1*(e.hat<0)))
  tmp.out = list(beta.hat=beta.new, loss=loss,ind=step)    
  return(tmp.out)
}