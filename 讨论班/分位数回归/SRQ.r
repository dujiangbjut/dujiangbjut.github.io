SRQ=function(x,y,tau){
    require(MASS)
    p=ncol(x)
    n=length(y)
    h=log(n)/sqrt(n)
    beta.new=rq(y~x+0,tau)$coef
    beta.old=rep(0,p)
    Loop=100
    loop=0
    error=1
    while(loop<Loop&error>1e-3){
      loop=loop+1
      beta.old=beta.new 
      e.hat=(y-x%*%beta.old)/h
      L1=t(x)%*%(tau-pnorm(-e.hat))
      L2=-t(x)%*%diag(as.vector(dnorm(e.hat)))%*%x/h
      beta.new=beta.old-ginv(L2)%*%L1
      e=(y-x%*%beta.old)     
      error=sum(abs(beta.old-beta.new))
    }
     D=tau*(1-tau)*t(x)%*%x/n 
     H=ginv(L2/n)%*%D%*%ginv(L2/n)/n 
     #H=ginv(L2/n)%*%L1%*%t(L1)%*%ginv(L2/n) 
     beta.sd=sqrt(diag(H))
     tmp.out = list(beta.hat=beta.new, beta.sd= beta.sd)    
     return(tmp.out)
}

