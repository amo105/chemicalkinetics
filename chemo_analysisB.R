library(tseries)
library(deSolve)
library(mvtnorm)
library(MCMCpack)
library(tmvtnorm)
library(lhs)
library(fields)
library(msm)
library(mgcv)
library(chemopac)
source("kron.prod.R")

source("data.R")
source("antony_ode.R")
source("estimate_r.R")
source("adaption.R")

ppp<-5

Y<-log(resp)
y<-as.vector(Y)
missing<-(1:length(y))[y==(-Inf)]

lower<-c(22.5,91.4,26.42,25,31.28,0)
upper<-c(45,91.59,26.47,40,32.56,3000)

# fty<-function(u){
# sum((qnorm(p=c(0.025,0.975),mean=u[1],sd=sqrt(u[2]))-log(c(10^(-8),10^(-4))))^2)}
# opt1<-optim(fn=fty,par=c(-13,1),method="L-BFGS-B",lower=c(-Inf,0))

# fty<-function(u){
# sum((qnorm(p=c(0.025,0.975),mean=u[1],sd=sqrt(u[2]))-log(c(100,1000000)))^2)}
# opt2<-optim(fn=fty,par=c(10,1),method="L-BFGS-B",lower=c(-Inf,0))

prior.loc<-c(rep(-13.8,ppp-2),rep(9.21,2))
prior.var<-c(rep(5.52,ppp-2),rep(5.52,2))
## Limits on the design parameters and prior hyperparameters.

design<-0*diag(5)
for(i in 1:5){
design[i,]<-as.vector(DATASET[min((1:dim(DATASET)[1])[marker==i]),-c(3,7)])}
times<-c(dataset$time[marker<6],dataset$time[(length(dataset$time)-1):length(dataset$time)])
what.times<-c(marker[marker<5],rep(5,20))

uni.times<-sort(unique(times))

MU<-function(inputs,times){
mu5(thet=inputs[1:5],x=inputs[6:10],tim=times)}

M<-function(theta){
output<-mu5(thet=theta,x=design[1,],tim=times[what.times==1])
for(j in 2:4){
output<-rbind(output,mu5(thet=theta,x=design[j,],tim=times[what.times==j]))}
output<-rbind(output,mu5(thet=theta,x=design[5,],tim=sort(times[what.times==5]))[c(1:16,18,20,17,19),])
output}

m<-function(theta){
as.vector(M(theta))}

set.seed(1)
INI<-matrix(0,nrow=0,ncol=ppp)
INIt<-matrix(0,nrow=0,ncol=ppp)
Z<-matrix(0,nrow=0,ncol=3)
while(dim(INI)[1]<50){
if(dim(INI)[1]==0){
ini.u<-runif(ppp)
ini.theta<-qnorm(p=ini.u,mean=prior.loc,sd=sqrt(prior.var))} else{
ini.u<-augmentLHS(lhs=INIt,m=1)
ini.u<-ini.u[dim(ini.u)[1],]
ini.theta<-qnorm(p=ini.u,mean=prior.loc,sd=sqrt(prior.var))}
try(z<-MU(inputs=c(ini.theta,design[1,]),times=uni.times))
try(z<-rbind(z,MU(inputs=c(ini.theta,design[2,]),times=uni.times)))
try(z<-rbind(z,MU(inputs=c(ini.theta,design[3,]),times=uni.times)))
try(z<-rbind(z,MU(inputs=c(ini.theta,design[4,]),times=uni.times)))
try(z<-rbind(z,MU(inputs=c(ini.theta,design[5,]),times=uni.times)))
if(!any(is.na(z)) & dim(z)[1]==155){
Z<-rbind(Z,z)
INIt<-rbind(INIt,ini.u)
INI<-rbind(INI,ini.theta)}}

n_MO<-dim(Z)[1]
kk<-dim(Z)[2]
m_MO<-1

X.arrayA<-array(0,dim=c(dim(INI)[1],ppp,dim(INI)[1]))
for(i in 1:dim(INI)[1]){
for(j in 1:ppp){
X.arrayA[i,j,]<-abs(INI[i,j]-INI[,j])}}

X.arrayB<-array(0,dim=c(dim(design)[1],5,dim(design)[1]))
for(i in 1:dim(design)[1]){
for(j in 1:5){
X.arrayB[i,j,]<-abs(design[i,j]-design[,j])}}

X.arrayC<-array(0,dim=c(length(uni.times),length(uni.times)))
for(i in 1:length(uni.times)){
X.arrayC[i,]<-abs(uni.times[i]-uni.times)}

q<-2

X.arrayA2<-X.arrayA^q
X.arrayB2<-X.arrayB^q
X.arrayC2<-X.arrayC^q

marg1<-estimate_r5(X.arrayA2=X.arrayA2,X.arrayB2=X.arrayB2,X.arrayC2=X.arrayC2,Z=Z,ini=c(rep(0,ppp),rep(-4,5),0),lower=rep(c(-5,-20,-5),c(5,5,1)))
logr<-marg1$par

DESIGN<-INI
Z<-Z[1:(dim(INI)[1]*length(uni.times)*5),]

GPF<-GP_functions5(X.arrayA2=X.arrayA2,X.arrayB2=X.arrayB2,X.arrayC2=X.arrayC2,Z=Z,opt.logr=logr,DESIGN=DESIGN,design=design,uni.times=uni.times)

diffs<-c()
for(ii in 1:dim(DESIGN)[1]){
diffs[ii]<-sum(((y-as.vector(X%*%matrix(GPF$mhat(DESIGN[ii,]),ncol=1)))[-missing])^2)}

#ini.theta<-DESIGN[(1:dim(DESIGN)[1])[diffs<=sort(diffs)[5]],1:ppp]
ini.theta<-matrix(rep(DESIGN[which.min(diffs),],5),nrow=5,byrow=TRUE)

KT<-5
AC<-matrix(0,ncol=KT,nrow=0)
AC.rho<-matrix(0,ncol=KT,nrow=0)
AC.psi<-matrix(0,ncol=KT,nrow=0)

error<-c()

eps<-rep(0.9,KT)
eps2<-eps^2
eps.rho<-rep(3000,KT)
eps.rho2<-eps.rho^2
eps.psi<-rep(1.1,KT)
eps.psi2<-eps.psi^2

while(dim(DESIGN)[1]<100){

if(dim(INI)[1]==dim(DESIGN)[1]){
curr.theta<-ini.theta
curr.missing<-matrix(rep(rep(c,length(missing)),KT),ncol=length(missing),byrow=TRUE)
curr.ome<-array(0,dim=c(KT,kk,kk))
icurr.ome<-array(0,dim=c(KT,kk,kk))
curr.sig<-array(0,dim=c(KT,kk,kk))
icurr.sig<-array(0,dim=c(KT,kk,kk))
for(kt in 1:KT){
curr.ome[kt,,]<-diag(kk)
icurr.ome[kt,,]<-solve.cpp(curr.ome[kt,,])
curr.sig[kt,,]<-diag(kk)
icurr.sig[kt,,]<-solve.cpp(curr.sig[kt,,])}
curr.rho<-rep(5,KT)
curr.psi<-cbind(rep(5,KT),rep(5,KT))

curr.T<-array(0,dim=c(KT,sNi,sNi))
curr.iT<-curr.T
curr.dT<-curr.T
curr.dT2<-curr.T
curr.detT<-rep(0,KT)
for(kt in 1:KT){
for(i in 1:6){
curr.T[kt,marker==i,marker==i]<-exp(-curr.rho[kt]*T.array[marker==i,marker==i])
curr.iT[kt,marker==i,marker==i]<-chol2inv(chol(curr.T[kt,marker==i,marker==i]))
curr.dT[kt,marker==i,marker==i]<--curr.rho[kt]*T.array[marker==i,marker==i]*curr.T[kt,marker==i,marker==i]
curr.dT2[kt,marker==i,marker==i]<-(curr.rho[kt]*T.array[marker==i,marker==i]-1)*curr.rho[kt]*T.array[marker==i,marker==i]*curr.T[kt,marker==i,marker==i]
curr.detT[kt]<-curr.detT[kt]+determinant(curr.T[kt,marker==i,marker==i])$modulus[1]}}

curr.A<-array(0,dim=c(KT,mm,mm))
curr.iA<-curr.A
curr.dA<-array(0,dim=c(KT,2,mm,mm))
curr.dA2<-array(0,dim=c(KT,3,mm,mm))
curr.detA<-rep(0,KT)
for(kt in 1:KT){
curr.A[kt,,]<-exp(-curr.psi[kt,1]*A.array1-curr.psi[kt,2]*A.array2)
curr.iA[kt,,]<-chol2inv(chol(curr.A[kt,,]))
curr.dA[kt,1,,]<--curr.psi[kt,1]*A.array1*curr.A[kt,,]
curr.dA[kt,2,,]<--curr.psi[kt,2]*A.array2*curr.A[kt,,]
curr.dA2[kt,1,,]<-(curr.psi[kt,1]*A.array1-1)*curr.psi[kt,1]*A.array1*curr.A[kt,,]
curr.dA2[kt,2,,]<-(curr.psi[kt,2]*A.array2-1)*curr.psi[kt,2]*A.array2*curr.A[kt,,]
curr.dA2[kt,3,,]<-curr.psi[kt,1]*curr.psi[kt,2]*A.array1*A.array2*curr.A[kt,,]
curr.detA[kt]<-determinant(curr.A[kt,,])$modulus[1]}

curr.R<-array(0,dim=c(KT,n,n))
curr.iR<-curr.R
curr.S<-array(0,dim=c(KT,mm*kk,mm*kk))
curr.iS<-curr.S
for(kt in 1:KT){
curr.R[kt,,]<-kronecker(curr.ome[kt,,],curr.T[kt,,])
curr.iR[kt,,]<-kronecker(icurr.ome[kt,,],curr.iT[kt,,])
curr.S[kt,,]<-kronecker(curr.sig[kt,,],curr.A[kt,,])
curr.iS[kt,,]<-kronecker(icurr.sig[kt,,],curr.iA[kt,,])}

ini.temp.rhos<-rep(-5,KT-1)

curr.dstar<-matrix(0,ncol=mm*kk,nrow=KT)
curr.m<-matrix(0,ncol=mm*kk,nrow=KT)
curr.eta<-matrix(0,ncol=n,nrow=KT)
curr.dm<-array(0,dim=c(KT,mm*kk,ppp))
curr.dm2<-array(0,dim=c(KT,ppp,mm*kk,ppp))
for(kt in 1:KT){
curr.dstar[kt,]<-rep(0,mm*kk)
curr.dm[kt,,]<-GPF$dmhat(curr.theta[kt,])
curr.dm2[kt,,,]<-GPF$d2mhat(curr.theta[kt,])
curr.m[kt,]<-GPF$mhat(curr.theta[kt,])
curr.eta[kt,]<-as.vector(X%*%matrix(curr.dstar[kt,],ncol=1))}

} else{
curr.m<-matrix(0,ncol=mm*kk,nrow=KT)
curr.eta<-matrix(0,ncol=n,nrow=KT)
curr.dm<-array(0,dim=c(KT,mm*kk,ppp))
curr.dm2<-array(0,dim=c(KT,ppp,mm*kk,ppp))
for(kt in 1:KT){
curr.dm[kt,,]<-GPF$dmhat(curr.theta[kt,])
curr.dm2[kt,,,]<-GPF$d2mhat(curr.theta[kt,])
curr.m[kt,]<-GPF$mhat(curr.theta[kt,])
curr.eta[kt,]<-as.vector(X%*%matrix(curr.dstar[kt,],ncol=1))}

}

###########################################################################################################################################################
###########################################################################################################################################################

fm<-function(theta,kt){
GPF$mhat(theta)}

LEN<-50
acc<-matrix(0,ncol=KT,nrow=0)
acc.rho<-matrix(0,ncol=KT,nrow=0)
acc.psi<-matrix(0,ncol=KT,nrow=0)
for(len in 1:LEN){

TEMP<-1
for(jj in 2:KT){
TEMP[jj]<-prod(exp(exp(ini.temp.rhos[1:(jj-1)])))}

whi<-sample(x=1:2,size=1)

if(whi==1){
inner<-rep(0,KT)
inner.rho<-rep(0,KT)
inner.psi<-rep(0,KT)

source("MCMC_B.R")

acc<-rbind(acc,inner)
acc.rho<-rbind(acc.rho,inner.rho)
acc.psi<-rbind(acc.psi,inner.psi)}

if(whi==2){

source("SWAP_B.R")}

source("UUU_B.R")

ini.temp.rhos<-adaption(rhos=ini.temp.rhos, u = UUU, iterations = len*(dim(DESIGN)[1]-dim(INI)[1])/KT+len)

}

###########################################################################################################################################################
###########################################################################################################################################################

z<-matrix(0,nrow=0,ncol=kk)
zhat<-z
for(kt in 1:KT){
for(ii in 1:5){
z<-rbind(z,mu5(thet=curr.theta[kt,],x=design[ii,],tim=uni.times))
zhat<-rbind(zhat,GPF$mu2hat(thet=curr.theta[kt,],x=design[ii,],tim=uni.times))}}

error<-c(error,log(tapply(X=as.vector((z-zhat)^2),INDEX=rep(rep(1:KT,each=length(uni.times)*5),kk),FUN=sum)))

AC<-rbind(AC,apply(acc,2,mean))
AC.rho<-rbind(AC.rho,apply(acc.rho,2,mean))
AC.psi<-rbind(AC.psi,apply(acc.psi,2,mean))

cat("accept",AC[dim(AC)[1],],"\n")
cat("accept rho",AC.rho[dim(AC.rho)[1],],"\n")
cat("accept psi",AC.psi[dim(AC.psi)[1],],"\n")

DESIGN<-rbind(DESIGN,curr.theta)
Z<-rbind(Z,z)

par(mfrow=c(3,2),mai=c(0.5,0.5,0.1,0.1),mgp=c(1,0.5,0))
ts.plot(DESIGN[,1])
ts.plot(DESIGN[,2])
ts.plot(DESIGN[,3])
ts.plot(DESIGN[,4])
ts.plot(DESIGN[,5])
ts.plot(error)
par(mfrow=c(1,1))

X.arrayA<-array(0,dim=c(dim(DESIGN)[1],ppp,dim(DESIGN)[1]))
for(i in 1:dim(DESIGN)[1]){
for(j in 1:ppp){
X.arrayA[i,j,]<-abs(DESIGN[i,j]-DESIGN[,j])}}
X.arrayA2<-X.arrayA^q

GPF<-GP_functions5(X.arrayA2=X.arrayA2,X.arrayB2=X.arrayB2,X.arrayC2=X.arrayC2,Z=Z,opt.logr=logr,DESIGN=DESIGN,design=design,uni.times=uni.times)

}

fm<-function(theta,kt){
if(kt==1){
ret<-m(theta=theta)}
if(kt>1){
ret<-GPF$mhat(theta=theta)}
ret}

TEMP<-c(1,TEMP)
KT<-length(TEMP)

curr.theta<-rbind(curr.theta[1,],curr.theta)
curr.missing<-rbind(curr.missing[1,],curr.missing)
curr.rho<-c(curr.rho[1],curr.rho)
curr.psi<-rbind(curr.psi[1,],curr.psi)
inter.ome<-curr.ome
inter.iome<-icurr.ome
inter.sig<-curr.sig
inter.isig<-icurr.sig
curr.ome<-array(0,dim=c(KT,kk,kk))
icurr.ome<-array(0,dim=c(KT,kk,kk))
curr.sig<-array(0,dim=c(KT,kk,kk))
icurr.sig<-array(0,dim=c(KT,kk,kk))
curr.ome[1,,]<-inter.ome[1,,]
icurr.ome[1,,]<-inter.iome[1,,]
curr.sig[1,,]<-inter.sig[1,,]
icurr.sig[1,,]<-inter.isig[1,,]
for(kt in 2:KT){
curr.ome[kt,,]<-inter.ome[kt-1,,]
icurr.ome[kt,,]<-inter.iome[kt-1,,]
curr.sig[kt,,]<-inter.sig[kt-1,,]
icurr.sig[kt,,]<-inter.isig[kt-1,,]}
rm(inter.ome);rm(inter.iome);rm(inter.sig);rm(inter.isig)

inter.T<-curr.T
inter.iT<-curr.iT
inter.dT<-curr.dT
inter.dT2<-curr.dT2
inter.detT<-curr.detT
curr.T<-array(0,dim=c(KT,sNi,sNi))
curr.iT<-curr.T
curr.dT<-curr.T
curr.dT2<-curr.T
curr.detT<-rep(0,KT)
curr.T[1,,]<-inter.T[1,,]
curr.iT[1,,]<-inter.iT[1,,]
curr.dT[1,,]<-inter.dT[1,,]
curr.dT2[1,,]<-inter.dT2[1,,]
curr.detT[1]<-inter.detT[1]
for(kt in 2:KT){
curr.T[kt,,]<-inter.T[kt-1,,]
curr.iT[kt,,]<-inter.iT[kt-1,,]
curr.dT[kt,,]<-inter.dT[kt-1,,]
curr.dT2[kt,,]<-inter.dT2[kt-1,,]
curr.detT[kt]<-inter.detT[kt-1]}

inter.A<-curr.A
inter.iA<-curr.iA
inter.dA<-curr.dA
inter.dA2<-curr.dA2
inter.detA<-curr.detA
curr.A<-array(0,dim=c(KT,mm,mm))
curr.iA<-curr.A
curr.dA<-array(0,dim=c(KT,2,mm,mm))
curr.dA2<-array(0,dim=c(KT,3,mm,mm))
curr.detA<-rep(0,KT)
curr.A[1,,]<-inter.A[1,,]
curr.iA[1,,]<-inter.iA[1,,]
curr.dA[1,,,]<-inter.dA[1,,,]
curr.dA2[1,,,]<-inter.dA2[1,,,]
curr.detA[1]<-inter.detA[1]
for(kt in 2:KT){
curr.A[kt,,]<-inter.A[kt-1,,]
curr.iA[kt,,]<-inter.iA[kt-1,,]
curr.dA[kt,,,]<-inter.dA[kt-1,,,]
curr.dA2[kt,,,]<-inter.dA2[kt-1,,,]
curr.detA[kt]<-inter.detA[kt-1]}

rm(inter.T);rm(inter.iT);rm(inter.dT);rm(inter.dT2);rm(inter.detT)
rm(inter.A);rm(inter.iA);rm(inter.dA);rm(inter.dA2);rm(inter.detA)

inter.R<-curr.R
inter.iR<-curr.iR
inter.S<-curr.S
inter.iS<-curr.iS
curr.R<-array(0,dim=c(KT,n,n))
curr.iR<-curr.R
curr.S<-array(0,dim=c(KT,mm*kk,mm*kk))
curr.iS<-curr.S
curr.R[1,,]<-inter.R[1,,]
curr.iR[1,,]<-inter.iR[1,,]
curr.S[1,,]<-inter.S[1,,]
curr.iS[1,,]<-inter.iS[1,,]
for(kt in 2:KT){
curr.R[kt,,]<-inter.R[kt-1,,]
curr.iR[kt,,]<-inter.iR[kt-1,,]
curr.S[kt,,]<-inter.S[kt-1,,]
curr.iS[kt,,]<-inter.iS[kt-1,,]}
rm(inter.R);rm(inter.iR);rm(inter.S);rm(inter.iS)

curr.dstar<-rbind(curr.dstar[1,],curr.dstar)
curr.m<-rbind(curr.m[1,],curr.m)
curr.eta<-rbind(curr.eta[1,],curr.eta)

inter.dm<-curr.dm
inter.dm2<-curr.dm2
curr.dm<-array(0,dim=c(KT,mm*kk,ppp))
curr.dm2<-array(0,dim=c(KT,ppp,mm*kk,ppp))
curr.dm[1,,]<-inter.dm[1,,]
curr.dm2[1,,,]<-inter.dm2[1,,,]
for(kt in 2:KT){
curr.dm[kt,,]<-inter.dm[kt-1,,]
curr.dm2[kt,,,]<-inter.dm2[kt-1,,,]}
rm(inter.dm);rm(inter.dm2)

THETA<-matrix(curr.theta[1,],nrow=1)
MISSING<-matrix(curr.missing[1,],nrow=1)
DSTAR<-matrix(curr.dstar[1,],nrow=1)
OME<-curr.ome[1,,]
SIG<-curr.sig[1,,]
RHO<-curr.rho[1]
PSI<-matrix(curr.psi[1,],nrow=1)
MMM<-matrix(curr.m[1,],nrow=1)
AC<-matrix(0,nrow=1,ncol=KT)
AC.rho<-matrix(0,nrow=1,ncol=KT)
AC.psi<-matrix(0,nrow=1,ncol=KT)

eps<-c(eps[1],eps)
eps2<-eps^2
eps.rho<-c(eps.rho[1],eps.rho)
eps.rho2<-eps.rho^2
eps.psi<-c(eps.psi[1],eps.psi)
eps.psi2<-eps.psi^2

counter<-1
while(counter<50000){

whi<-sample(x=1:2,size=1)

if(whi==1){

inner<-rep(0,KT)
inner.rho<-rep(0,KT)
inner.psi<-rep(0,KT)

source("MCMC_B.R")

AC<-rbind(AC,inner)
AC.rho<-rbind(AC.rho,inner.rho)
AC.psi<-rbind(AC.psi,inner.psi)}

if(whi==2){

source("SWAP_B.R")

}

THETA<-rbind(THETA,curr.theta[1,])
MISSING<-rbind(MISSING,curr.missing[1,])
DSTAR<-rbind(DSTAR,curr.dstar[1,])
OME<-rbind(OME,curr.ome[1,,])
SIG<-rbind(SIG,curr.sig[1,,])
RHO<-c(RHO,curr.rho[1])
PSI<-rbind(PSI,curr.psi[1,])
MMM<-rbind(MMM,curr.m[1,])

counter<-counter+1

if(ceiling(counter/100)==(counter/100)){
cat("Iteration", counter, "\n")
write.table(file="./mcmc/B_THETA.txt",x=THETA,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_MISSING.txt",x=MISSING,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_DSTAR.txt",x=DSTAR,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_OME.txt",x=OME,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_SIG.txt",x=SIG,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_RHO.txt",x=RHO,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_PSI.txt",x=PSI,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_MMM.txt",x=MMM,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_AC.txt",x=AC,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_ACrho.txt",x=AC.rho,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file="./mcmc/B_ACpsi.txt",x=AC.psi,row.names=FALSE,col.names=FALSE,append=TRUE)
THETA<-matrix(0,nrow=0,ncol=dim(THETA)[2])
MISSING<-matrix(0,nrow=0,ncol=dim(MISSING)[2])
DSTAR<-matrix(0,nrow=0,ncol=dim(DSTAR)[2])
OME<-matrix(0,nrow=0,ncol=dim(OME)[2])
SIG<-matrix(0,nrow=0,ncol=dim(SIG)[2])
RHO<-c()
AC<-matrix(0,nrow=0,ncol=dim(AC)[2])
AC.rho<-matrix(0,nrow=0,ncol=dim(AC.rho)[2])
AC.psi<-matrix(0,nrow=0,ncol=dim(AC.psi)[2])
PSI<-matrix(0,nrow=0,ncol=dim(PSI)[2])
MMM<-matrix(0,nrow=0,ncol=dim(MMM)[2])}

}

