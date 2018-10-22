
estimate_r5<-function(X.arrayA2,X.arrayB2,X.arrayC2,Z,ini,lower=rep(c(-5,-16,-5),c(5,5,1)),upper=rep(8,length(ini))){

drt<-apply(Z,2,mean)
Zr<-Z-matrix(rep(drt,dim(Z)[1]),ncol=length(drt),byrow=TRUE)

SIG<-var(Zr)

EIG<-eigen(x=SIG,symmetric=TRUE)
UUU<-EIG$vectors
LAM<-EIG$values
CCC<-UUU%*%diag(1/sqrt(LAM))
Zr<-Zr%*%CCC

n_MO<-dim(Z)[1]
lut<-dim(X.arrayC2)[1]
m_MO<-1
kk<-dim(Z)[2]
ddd<-dim(X.arrayA2)[1]
eee<-dim(X.arrayB2)[1]

marg.post<-function(paras){
#fff<<-rbind(fff,paras)
r<-exp(paras)
aaA<-exp(makeA5.cpp(D1=X.arrayA2[,1,],D2=X.arrayA2[,2,],D3=X.arrayA2[,3,],D4=X.arrayA2[,4,],D5=X.arrayA2[,5,],R=r[1:5]))
aaB<-exp(makeA5.cpp(D1=X.arrayB2[,1,],D2=X.arrayB2[,2,],D3=X.arrayB2[,3,],D4=X.arrayB2[,4,],D5=X.arrayB2[,5,],R=r[6:10]))
aaC<-exp(-r[11]*(X.arrayC2))
# iaaA<-solve.cpp(aaA)
# iaaB<-solve.cpp(aaB)
# iaaC<-solve.cpp(aaC)
iaaA<-solve(aaA)
iaaB<-solve(aaB)
iaaC<-solve(aaC)
iaaZ<-cbind(kron.prod(y=Zr[,1],matrices=list(iaaC,iaaB,iaaA)),kron.prod(y=Zr[,2],matrices=list(iaaC,iaaB,iaaA)),kron.prod(y=Zr[,3],matrices=list(iaaC,iaaB,iaaA)))
detaa<-n_MO*((1/ddd)*determinant(aaA)$modulus[1]+(1/eee)*determinant(aaB)$modulus[1]+(1/lut)*determinant(aaC)$modulus[1])
ret<--0.5*kk*detaa-0.5*sum(iaaZ*Zr)
ret}

dmarg.post<-function(paras){
#ggg<<-rbind(ggg,paras)
	r<-exp(paras)
aaA<-exp(makeA5.cpp(D1=X.arrayA2[,1,],D2=X.arrayA2[,2,],D3=X.arrayA2[,3,],D4=X.arrayA2[,4,],D5=X.arrayA2[,5,],R=r[1:5]))
aaB<-exp(makeA5.cpp(D1=X.arrayB2[,1,],D2=X.arrayB2[,2,],D3=X.arrayB2[,3,],D4=X.arrayB2[,4,],D5=X.arrayB2[,5,],R=r[6:10]))
aaC<-exp(-r[11]*(X.arrayC2))
# iaaA<-solve.cpp(aaA)
# iaaB<-solve.cpp(aaB)
# iaaC<-solve.cpp(aaC)
iaaA<-solve(aaA)
iaaB<-solve(aaB)
iaaC<-solve(aaC)
der<-c()
for(i in 1:5){
dA<--r[i]*X.arrayA2[,i,]*aaA
iAdA<-iaaA%*%dA
iAdAiA<-iAdA%*%iaaA
iaaZ<-cbind(kron.prod(y=Zr[,1],matrices=list(iaaC,iaaB,iAdAiA)),kron.prod(y=Zr[,2],matrices=list(iaaC,iaaB,iAdAiA)),kron.prod(y=Zr[,3],matrices=list(iaaC,iaaB,iAdAiA)))
der[i]<-0.5*sum(iaaZ*Zr)-0.5*kk*lut*eee*sum(diag(iAdA))}
for(i in 1:5){
dB<--r[i+5]*X.arrayB2[,i,]*aaB
iBdB<-iaaB%*%dB
iBdBiB<-iBdB%*%iaaB
iaaZ<-cbind(kron.prod(y=Zr[,1],matrices=list(iaaC,iBdBiB,iaaA)),kron.prod(y=Zr[,2],matrices=list(iaaC,iBdBiB,iaaA)),kron.prod(y=Zr[,3],matrices=list(iaaC,iBdBiB,iaaA)))
der[i+5]<-0.5*sum(iaaZ*Zr)-0.5*kk*lut*ddd*sum(diag(iBdB))}
dC<--r[11]*X.arrayC2*aaC
iCdC<-iaaC%*%dC
iCdCiC<-iCdC%*%iaaC
iaaZ<-cbind(kron.prod(y=Zr[,1],matrices=list(iCdCiC,iaaB,iaaA)),kron.prod(y=Zr[,2],matrices=list(iCdCiC,iaaB,iaaA)),kron.prod(y=Zr[,3],matrices=list(iCdCiC,iaaB,iaaA)))
der[11]<-0.5*sum(iaaZ*Zr)-0.5*kk*ddd*eee*sum(diag(iCdC))
der}

opt<-optim(fn=marg.post,gr=dmarg.post,par=ini,control=list(fnscale=-1,maxit=200),method="L-BFGS-B",lower=lower,upper=upper)

return(opt)}

GP_functions5<-function(X.arrayA2,X.arrayB2,X.arrayC2,Z,opt.logr,DESIGN,design,uni.times){

drt<-apply(Z,2,mean)
Zr<-Z-matrix(rep(drt,dim(Z)[1]),ncol=length(drt),byrow=TRUE)

n_MO<-dim(Z)[1]
lut<-dim(X.arrayC2)[1]
m_MO<-1
kk<-dim(Z)[2]
ddd<-dim(X.arrayA2)[1]
eee<-dim(X.arrayB2)[1]

r<-exp(opt.logr)
aaA<-exp(makeA5.cpp(D1=X.arrayA2[,1,],D2=X.arrayA2[,2,],D3=X.arrayA2[,3,],D4=X.arrayA2[,4,],D5=X.arrayA2[,5,],R=r[1:5]))
aaB<-exp(makeA5.cpp(D1=X.arrayB2[,1,],D2=X.arrayB2[,2,],D3=X.arrayB2[,3,],D4=X.arrayB2[,4,],D5=X.arrayB2[,5,],R=r[6:10]))
aaC<-exp(-r[11]*(X.arrayC2))
iaaA<-solve.cpp(aaA)
iaaB<-solve.cpp(aaB)
iaaC<-solve.cpp(aaC)

iaaZ<-cbind(kron.prod(y=Zr[,1],matrices=list(iaaC,iaaB,iaaA)),kron.prod(y=Zr[,2],matrices=list(iaaC,iaaB,iaaA)),kron.prod(y=Zr[,3],matrices=list(iaaC,iaaB,iaaA)))
residualv<-as.vector(Zr)
s.hat<-t(Zr)%*%iaaZ

srt1<-c();srt2<-c();srt3<-c();srt4<-c();srt5<-c()
for(bob in times[what.times==1]){
srt1<-c(srt1,(1:length(uni.times))[uni.times==bob])}
for(bob in times[what.times==2]){
srt2<-c(srt2,(1:length(uni.times))[uni.times==bob])}
for(bob in times[what.times==3]){
srt3<-c(srt3,(1:length(uni.times))[uni.times==bob])}
for(bob in times[what.times==4]){
srt4<-c(srt4,(1:length(uni.times))[uni.times==bob])}
for(bob in times[what.times==5]){
srt5<-c(srt5,(1:length(uni.times))[uni.times==bob])}
srt<-list(srt1,srt2,srt3,srt4,srt5)

SRT<-c(srt[[1]],31+srt[[2]],62+srt[[3]],93+srt[[4]],124+srt[[5]])
SRT<-c(SRT,155+SRT,2*155+SRT)

#######################################################################################################

sigtilde<-function(THETA,INI_X,TIM){
aA<-makeaa.cpp(DESIGN=DESIGN,THETA=THETA,R=r[1:5])
aB<-makeaa.cpp(DESIGN=design,THETA=INI_X,R=r[6:10])
aC<-makeaa.cpp(DESIGN=matrix(uni.times,ncol=1),THETA=matrix(TIM,ncol=1),R=r[11])
ata<-aA%*%iaaA
btb<-aB%*%iaaB
ctc<-aC%*%iaaC
tata<-t(ata)
1-diags.cpp(X=aA,Y=tata)*as.vector(aB%*%t(btb))*as.vector(aC%*%t(ctc))}

mu2hat<-function(thet,x,tim){
aA<-matrix(exp(-rowSums.cpp(matrix(rep(r[1:5],ddd),ncol=5,byrow=TRUE)*(matrix(rep(thet,ddd),ncol=5,byrow=TRUE)-DESIGN)^2)),nrow=1)
aB<-matrix(exp(-rowSums.cpp(matrix(rep(r[6:10],eee),ncol=5,byrow=TRUE)*(matrix(rep(x,eee),ncol=5,byrow=TRUE)-design)^2)),nrow=1)
aC<-matrix(exp(-r[11]*(rep(tim,each=length(uni.times))-rep(uni.times,length(tim)))^2),nrow=length(tim),byrow=TRUE)
ata<-aA%*%iaaA
btb<-aB%*%iaaB
ctc<-aC%*%iaaC
out<-cbind(kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,1]),kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,2]),kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,3]))
matrix(rep(drt,length(tim)),ncol=kk,byrow=TRUE)+out}

mu2hatB<-function(THETA,INI_X,TIM){
aA<-makeaa.cpp(DESIGN=DESIGN,THETA=THETA,R=r[1:5])
aB<-makeaa.cpp(DESIGN=design,THETA=INI_X,R=r[6:10])
aC<-makeaa.cpp(DESIGN=matrix(uni.times,ncol=1),THETA=matrix(TIM,ncol=1),R=r[11])
ata<-aA%*%iaaA
btb<-aB%*%iaaB
ctc<-aC%*%iaaC
out<-cbind(kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,1]),kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,2]),kron3cpp(aa0=ata,aa1=btb,aa2=ctc,yy=Zr[,3]))
matrix(rep(drt,dim(out)[1]),ncol=kk,byrow=TRUE)+out}

Mhat<-function(theta){
nnn<-length(uni.times)
ddd<-dim(DESIGN)
diffs<-matrix(rep(theta,ddd[1]),ncol=ddd[2],byrow=TRUE)-DESIGN
AAA<-Mhat5.cpp(RESIDUAL=Zr,D2=diffs^2,R=r[1:5],iaaA=iaaA)
matrix(rep(drt,mm),ncol=kk,byrow=TRUE)+AAA[SRT[1:mm],]}

mhat<-function(theta){
as.vector(Mhat(theta))}

#########################################################################################################

dMhat<-function(theta){
nnn<-length(uni.times)
ddd<-dim(DESIGN)
diffs<-matrix(rep(theta,dim(DESIGN)[1]),ncol=dim(DESIGN)[2],byrow=TRUE)-DESIGN
out<-dMhat5.cpp(RESIDUAL=Zr,iaaA=iaaA,D=diffs,D2=diffs^2,R=r[1:5])
out2<-array(0,dim=c(nnn*5,5,kk))
for(j in 1:5){
out2[,j,]<-out[((j-1)*nnn*5+1):(j*nnn*5),]}
out2[SRT[1:mm],,]}

dmhat<-function(theta){
out<-dMhat(theta)
out2<-cbind(as.vector(out[,1,]),as.vector(out[,2,]),as.vector(out[,3,]),as.vector(out[,4,]),as.vector(out[,5,]))
out2}

############################################################################################################

d2Mhat<-function(theta){
nnn<-length(uni.times)
ddd<-dim(DESIGN)
diffs<-matrix(rep(theta,ddd[1]),ncol=ddd[2],byrow=TRUE)-DESIGN
out<-d2Mhat5.cpp(RESIDUAL=Zr,iaaA=iaaA,D=diffs,D2=diffs^2,R=r[1:5])
output<-array(0,dim=c(nnn*5,5,5,kk))
for(ii in 1:5){
for(jj in 1:5){
output[,ii,jj,]<-out[((ii-1)*5+jj-1)*155+(1:155),]}}
output<-output[SRT[1:mm],,,]
output}

d2mhat<-function(theta){
out<-d2Mhat(theta)
out2<-array(0,dim=c(5,mm*kk,5))
for(j in 1:5){
out2[j,,]<-cbind(as.vector(out[,j,1,]),as.vector(out[,j,2,]),as.vector(out[,j,3,]),as.vector(out[,j,4,]),as.vector(out[,j,5,]))}
out2}

###############################################################################################################

list(s.hat=s.hat,sigtilde=sigtilde,Mhat=Mhat,mhat=mhat,dMhat=dMhat,dmhat=dmhat,d2Mhat=d2Mhat,d2mhat=d2mhat,mu2hat=mu2hat,mu2hatB=mu2hatB)

}

