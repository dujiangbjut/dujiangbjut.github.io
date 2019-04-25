setwd("F:\\latex\\files\\导师\\分位数回归")
source('CRQ.r')
source('DCRQ.r')
source('VSCRQ.r')
source('VSCRQ2.r')
source('VSDCRQ2.r')
source('VSDCRQ.r')
source('IPRQ.r')
source('VSRQ.r')
source('DRQ.r')
source('VSIPRQ1.r')
source('VSDRQ.r')
source("ADMMRQ.r")
source("ADMMCRQ.r")
source("MMCRQ.r")
source("MMRQ.r")
source('sdEST.r')
source('ESD.r')
source('SRQ.r')
source('ISRQ.r')
source('VSADMMRQ.r')
library(lpSolve)
library(Matrix)
require(xtable)
##分位数  不含惩罚 table 1
n=c(200,400,600,800,1000,2000)
p=5
M=50
tau=0.3
numofmethod=5
options(digits = 20)
result=matrix(0,length(n),2*numofmethod)
for (i in 1:length(n)) {
  for (j in 1:M) {
  set.seed(201801)
  x=rnorm(n[i]*p)
  x=matrix(x,n[i],p)
  e=rnorm(n[i])
  beta.true=runif(p,-1,1)
  y=x%*%beta.true+e 
  #IP quantreg
  ptm=proc.time()
  beta.hat=rq(y~x+0,tau)$coef
  Time=proc.time()-ptm
  result[i,1]=sum(abs(beta.hat-beta.true))+result[i,1]
  result[i,2]=Time[1]+result[i,2]
  #IPRQ
  ptm=proc.time()
  beta.hat=IPRQ(x,y,tau)$beta.hat
  Time=proc.time()-ptm
  result[i,3]=sum(abs(beta.hat-beta.true))+result[i,3]
  result[i,4]=Time[1]+result[i,4]
  #DRQ
  ptm=proc.time()
  beta.hat=DRQ(x,y,tau)$beta.hat
  Time=proc.time()-ptm
  result[i,5]=sum(abs(beta.hat-beta.true))+result[i,5]
  result[i,6]=Time[1]+result[i,6]
  #ADMMRQ
  ptm=proc.time()
  beta.hat=ADMMRQ(x,y,tau)$beta.hat
  Time=proc.time()-ptm
  result[i,7]=sum(abs(beta.hat-beta.true))+result[i,7]
  result[i,8]=Time[1]+result[i,8]
  #MM
  ptm=proc.time()
  beta.hat=MMRQ(x,y,tau)$beta.hat
  Time=proc.time()-ptm
  result[i,9]=sum(abs(beta.hat-beta.true))+result[i,9]
  result[i,10]=Time[1]+result[i,10]
  }
}
result=result/M
Rname=NULL
for (i in 1:length(n)) {
  Rname[i]= paste('(',n[i],',',p,')')
}
name=c("error","time")
Cname=rep(name,numofmethod)
rownames(result)<-Rname
colnames(result)<-Cname
xtable((result),digits=3)
#adpLasso






##复合分位数  #无惩罚 table 2 3
n=c(200,400,600,800,1000)
p=5
M=50
tau=c(2,5,8)
tau=tau/10
numofmethod=4
options(digits = 20)
result=matrix(0,length(n),2*numofmethod)
cqr=matrix(0,length(n),6)
for (i in 1:length(n)) {
  for (j in 1:M) {
    set.seed(201801)
    x=rnorm(n[i]*p)
    x=matrix(x,n[i],p)
    e=rnorm(n[i])
    beta.true=runif(p,-1,1)
    y=x%*%beta.true+e 
    #CRQ
    ptm=proc.time()
    beta.hat=CRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,1]=sum(abs(beta.hat[1:p]-beta.true))+result[i,1]
    result[i,2]=Time[1]+result[i,2]
    #DCRQ
    ptm=proc.time()
    beta.hat=DCRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,3]=sum(abs(beta.hat[1:p]-beta.true))+result[i,3]
    result[i,4]=Time[1]+result[i,4]
    #ADMMCRQ
    ptm=proc.time()
    beta.hat=ADMMCRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,5]=sum(abs(beta.hat[1:p]-beta.true))+result[i,5]
    result[i,6]=Time[1]+result[i,6]
    #MMCRQ
    ptm=proc.time()
    beta.hat=MMCRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,7]=sum(abs(beta.hat[1:p]-beta.true))+result[i,7]
    result[i,8]=Time[1]+result[i,8]
    ####cqr.fit
    #ip
    ptm=proc.time()
    beta.hat=cqr.fit(x,y,tau,method = "ip")$beta
    Time=proc.time()-ptm
    cqr[i,1]=sum(abs(beta.hat-beta.true))+result[i,1]
    cqr[i,2]=Time[1]+result[i,2]
    #mm
    ptm=proc.time()
    beta.hat=cqr.fit(x,y,tau,method = "mm")$beta
    Time=proc.time()-ptm
    cqr[i,3]=sum(abs(beta.hat-beta.true))+result[i,3]
    cqr[i,4]=Time[1]+result[i,4]
    #admm
    ptm=proc.time()
    beta.hat=cqr.fit(x,y,tau,method = "admm")$beta
    Time=proc.time()-ptm
    cqr[i,5]=sum(abs(beta.hat-beta.true))+result[i,5]
    cqr[i,6]=Time[1]+result[i,6]
  }
}
result=result/M
cqr=cqr/M
Rname=NULL
for (i in 1:length(n)) {
  Rname[i]= paste('(',n[i],',',p,')')
}
name=c("error","time")
Cname=rep(name,numofmethod)
rownames(result)<-Rname
colnames(result)<-Cname
xtable((result),digits=3)

rownames(cqr)=Rname
colnames(cqr)=rep(name,3)
xtable(cqr)

##lasso惩罚 单一分位数
n=c(200,400,600,800,1000,2000)
p=5
M=50
tau=0.3
numofmethod=5
options(digits = 20)
result=matrix(0,length(n),2*numofmethod)
for (i in 1:length(n)) {
  for (j in 1:M) {
    set.seed(201801)
    x=rnorm(n[i]*p)
    x=matrix(x,n[i],p)
    e=rnorm(n[i])
    beta.true=runif(p,-1,1)
    y=x%*%beta.true+e 
    #IP quantreg
    ptm=proc.time()
    beta.hat=rq(y~x+0,tau)$coef
    Time=proc.time()-ptm
    result[i,1]=sum(abs(beta.hat-beta.true))+result[i,1]
    result[i,2]=Time[1]+result[i,2]
    #IPRQ
    ptm=proc.time()
    beta.hat=IPRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,3]=sum(abs(beta.hat-beta.true))+result[i,3]
    result[i,4]=Time[1]+result[i,4]
    #DRQ
    ptm=proc.time()
    beta.hat=DRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,5]=sum(abs(beta.hat-beta.true))+result[i,5]
    result[i,6]=Time[1]+result[i,6]
    #ADMMRQ
    ptm=proc.time()
    beta.hat=ADMMRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,7]=sum(abs(beta.hat-beta.true))+result[i,7]
    result[i,8]=Time[1]+result[i,8]
    #MM
    ptm=proc.time()
    beta.hat=MMRQ(x,y,tau)$beta.hat
    Time=proc.time()-ptm
    result[i,9]=sum(abs(beta.hat-beta.true))+result[i,9]
    result[i,10]=Time[1]+result[i,10]
  }
}
result=result/M
Rname=NULL
for (i in 1:length(n)) {
  Rname[i]= paste('(',n[i],',',p,')')
}
name=c("error","time")
Cname=rep(name,numofmethod)
rownames(result)<-Rname
colnames(result)<-Cname
xtable((result),digits=3)




##方差估计
cmpInterval<-function(beta,beta.sd,alpha,n){
  l=qnorm(alpha/2)
  beta.l=beta + l*beta.sd/sqrt(n)
  r=qnorm(1-alpha/2)
  beta.r=beta+r*beta.sd/sqrt(n)
  return(cbind(beta.l,beta.r))
}
library(quantreg)
options(digits=5)
n=100
M=5
pn=2
p1=2
p2=2
beta=c(rep(-1,pn),rep(1,p1),rep(0,p2))
p=p1+p2+pn
alpha=0.05
Resinterval=matrix(0,p,8)
for (i in 1:M) {
x=rnorm(n*p)
x=matrix(x,n,p)
e=rnorm(n)
tau=0.3
e=e-quantile(e,0.5)
y=x%*%beta+e
beta.hat=rq(y~x+0,tau)$coef
fit1=ESD(x,y,tau)
Resinterval[,1:2]=cmpInterval(beta.hat,fit1$beta.sd,alpha ,n)+Resinterval[,1:2]
fit2=sdEST(x,y,tau)
Resinterval[,3:4]=cmpInterval(beta.hat,fit2$beta.sd,alpha,n)+Resinterval[,3:4]
fit3=SRQ(x,y,tau)
Resinterval[,5:6]=cmpInterval(beta.hat,fit3$beta.sd,alpha,n)+Resinterval[,5:6]
fit4=ISRQ(x,y,tau)
Resinterval[,7:8]=cmpInterval(beta.hat,fit4$beta.sd,alpha,n)+Resinterval[,7:8]

}
Resinterval=Resinterval/M
format(Resinterval,digits=2)
xtable(Resinterval)

##power




fit2=DCRQ(x,y,tau)

CRQ=fit1$beta.hat
DCRQ=fit2$beta.hat
lambda=log(n)/(abs(CRQ)+1/n) 


fit3=VSCRQ(x,y,tau,lambda)
VSCRQ=fit3$beta.hat

fit4=VSDCRQ(x,y,tau,lambda)
VSDCRQ=fit4$beta.hat

fit5=VSCRQ2(x,y,tau,lambda)
VSCRQ2=fit5$beta.hat

fit6=VSDCRQ2(x,y,tau,lambda)
VSDCRQ2=fit6$beta.hat
result=cbind(CRQ,DCRQ,VSCRQ,VSDCRQ,VSCRQ2,VSDCRQ2)


fit0=CRQ(x,y,tau)
fit0$beta.hat
lambda=log(n)/(abs(fit0$beta.hat)+1/n) 

fit00=VSDRQ2(x,y,tau,lambda)
fit00$beta.hat


N=100
rho=seq(from=0.1,to=10,by=0.1)

temp1=temp2=matrix(0,1,N)
result1=result2=matrix(0,length(rho),1)
for (i in 1:length(rho)) {
  
  for (j in 1:N) {
    x=rnorm(n*p)
    x=matrix(x,n,p)
    e=rnorm(n)
    e=e-quantile(e,0.5)
    y=x%*%beta+e
    fit1=ADMMRQ(x,y,tau)
    temp1[j]=fit1$ind
    fit2=ADMM(x,y,tau,rho=rho[i])
    temp2[j]=fit2$ind
  }
  
  result1[i]=mean(temp1)
  result2[i]=mean(temp2)
}


