ISRQ=function(x,y,tau){
    p=ncol(x)
    n=length(y)
    beta.new=rq(y~x+0,tau)$coef
    beta.old=rep(0,p)
    Loop=100
    loop=0
    error=1
    H=diag(p)/n
    while(loop<Loop&error>1e-3){
      loop=loop+1
      beta.old=beta.new
      h=sqrt(diag((x)%*%H%*%t(x)))
      e.hat=-(y-x%*%beta.old)/h
      L1=t(x)%*%(pnorm(e.hat)-tau)
      L2=t(x)%*%diag(as.vector(dnorm(e.hat)/h))%*%x
      beta.new=beta.old-ginv(L2)%*%L1
      e=(y-x%*%beta.old)
      Temp= t(x)%*%e
      D=tau*(1-tau)*t(x)%*%x/n 
      #H=ginv(L2)%*%L1%*%t(L1)%*%ginv(L2)
      H=ginv(L2/n)%*%D%*%ginv(L2/n)/n 
      error=sum(abs(beta.old-beta.new))
    }
    beta.sd=sqrt(diag(H))
    tmp.out = list(beta.hat=beta.new, beta.sd= beta.sd)    
    return(tmp.out)
}

