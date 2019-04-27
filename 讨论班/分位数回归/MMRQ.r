  MMRQ=function(x,y,tau,delta=1e-3){
    p=ncol(x)
    n=length(y)    
    beta.hat=rnorm(p)  
    beta.new=beta.hat
    beta.old=rep(n,length(beta.new))    
    loop=0
    Loop=1000
    esp=p*0.01
    while(sum(abs(beta.new-beta.old))>esp&&loop<Loop){
      loop=loop+1
      beta.old=beta.new     
      e.hat=y-x%*%beta.old      
      W=diag(as.vector(1/(delta+abs(e.hat))))
      beta.new=beta.old+solve(t(x)%*%W%*%x)%*%t(x)%*%(2*tau-1+e.hat/(abs(e.hat)+delta))
     #beta.new=solve(t(x)%*%W%*%x)%*%t(x)%*%(2*tau-1+y/(abs(e.hat)+delta))
    }
    ind=sum(abs(beta.new-beta.old))#1*(loop<Loop)
    ind=loop
    loss=0 
    temp=y-x%*%beta.new
    loss=sum(temp*(tau-1*(temp<0)))
    tmp.out = list(beta.hat=beta.new, loss=loss,ind=ind)    
  return(tmp.out)
}
 

