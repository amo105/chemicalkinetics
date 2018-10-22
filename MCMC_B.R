nu<-4

for(kt in 1:KT){

curr.y<-y
curr.y[missing]<-curr.missing[kt,]

############################ dstar ################################################

GT<-t(G)%*%curr.iT[kt,,]
HRH<-icurr.ome[kt,,]%x%(GT%*%G)
iVd<-HRH+icurr.sig[kt,,]%x%curr.iA[kt,,]
Vd<-chol2inv(chol(iVd))

mud<-as.vector(Vd%*%matrix(kron2cpp(aa0=icurr.ome[kt,,],aa1=GT,yy=curr.y)+kron2cpp(aa0=icurr.sig[kt,,],aa1=curr.iA[kt,,],yy=curr.m[kt,]),ncol=1))
curr.dstar[kt,]<-as.vector(rmvnorm(n=1,mean=mud,sigma=TEMP[kt]*Vd))
curr.eta[kt,]<-kron2cpp(aa0=diag(kk),aa1=G,yy=curr.dstar[kt,])

############################## theta ##############################################

curr.jack<-as.vector(matrix(curr.dstar[kt,]-curr.m[kt,],nrow=1)%*%curr.iS[kt,,]%*%curr.dm[kt,,])-(curr.theta[kt,]-prior.loc)/prior.var
curr.G<-t(curr.dm[kt,,])%*%curr.iS[kt,,]%*%curr.dm[kt,,]+diag(ppp)/prior.var
curr.iG<-chol2inv(chol(curr.G))

curr.Gj<-array(0,dim=c(ppp,ppp,ppp))
for(j in 1:ppp){
curr.Gj[j,,]<-2*t(curr.dm2[kt,j,,])%*%curr.iS[kt,,]%*%curr.dm[kt,,]}

curr.bit1<-c()
curr.bit2<-c()
for(i in 1:ppp){
inner1<-c()
inner2<-c()
for(j in 1:ppp){
GLOB<-curr.iG%*%curr.Gj[j,,]
inner1[j]<-(GLOB%*%curr.iG)[i,j]
inner2[j]<-curr.iG[i,j]*sum(diag(GLOB))}
curr.bit1[i]<-sum(inner1)
curr.bit2[i]<-sum(inner2)}

prop.mean<-curr.theta[kt,]+0.5*eps2[kt]*as.vector(curr.iG%*%matrix(curr.jack,ncol=1))-eps2[kt]*curr.bit1+0.5*eps2[kt]*curr.bit2
prop.var<-TEMP[kt]*eps2[kt]*curr.iG
iprop.var<-curr.G/(TEMP[kt]*eps2[kt])

prop.theta<-as.vector(rmvnorm(n=1,mean=prop.mean,sigma=prop.var))
prop.m<-fm(theta=prop.theta,kt=kt)

prop.dm<-GPF$dmhat(prop.theta)
prop.dm2<-GPF$d2mhat(prop.theta)

prop.jack<-as.vector(matrix(curr.dstar[kt,]-prop.m,nrow=1)%*%curr.iS[kt,,]%*%prop.dm)-(prop.theta-prior.loc)/prior.var
prop.G<-t(prop.dm)%*%curr.iS[kt,,]%*%prop.dm+diag(ppp)/prior.var
prop.iG<-chol2inv(chol(prop.G))

prop.Gj<-array(0,dim=c(ppp,ppp,ppp))
for(j in 1:ppp){
prop.Gj[j,,]<-2*t(prop.dm2[j,,])%*%curr.iS[kt,,]%*%prop.dm}

prop.bit1<-c()
prop.bit2<-c()
for(i in 1:ppp){
inner1<-c()
inner2<-c()
for(j in 1:ppp){
GLOB<-prop.iG%*%prop.Gj[j,,]
inner1[j]<-(GLOB%*%prop.iG)[i,j]
inner2[j]<-prop.iG[i,j]*sum(diag(GLOB))}
prop.bit1[i]<-sum(inner1)
prop.bit2[i]<-sum(inner2)}

curr.mean<-prop.theta+0.5*eps2[kt]*as.vector(prop.iG%*%matrix(prop.jack,ncol=1))-eps2[kt]*prop.bit1+0.5*eps2[kt]*prop.bit2
curr.var<-TEMP[kt]*eps2[kt]*prop.iG
icurr.var<-prop.G/(TEMP[kt]*eps2[kt])

pY<-matrix(curr.dstar[kt,]-prop.m,ncol=kk)
cY<-matrix(curr.dstar[kt,]-curr.m[kt,],ncol=kk)
top<-(1/TEMP[kt])*(-0.5*sum(diag(icurr.sig[kt,,]%*%t(pY)%*%curr.iA[kt,,]%*%pY))-0.5*sum(((prop.theta-prior.loc)^2)/prior.var))-0.5*determinant(curr.var)$modulus[1]-0.5*as.vector(matrix(curr.theta[kt,]-curr.mean,nrow=1)%*%icurr.var%*%matrix(curr.theta[kt,]-curr.mean,ncol=1))
bot<-(1/TEMP[kt])*(-0.5*sum(diag(icurr.sig[kt,,]%*%t(cY)%*%curr.iA[kt,,]%*%cY))-0.5*sum(((curr.theta[kt,]-prior.loc)^2)/prior.var))-0.5*determinant(prop.var)$modulus[1]-0.5*as.vector(matrix(prop.theta-prop.mean,nrow=1)%*%iprop.var%*%matrix(prop.theta-prop.mean,ncol=1))
prob.theta<-exp(top-bot)
if(prob.theta>=runif(1)){
inner[kt]<-1
curr.theta[kt,]<-prop.theta
curr.m[kt,]<-prop.m
curr.dm[kt,,]<-prop.dm
curr.dm2[kt,,,]<-prop.dm2}

######################### rho ############################################################################

curr.residual1<-curr.y-curr.eta[kt,]

curr.iTdT<-curr.iT[kt,,]%*%curr.dT[kt,,]
curr.iTdTiT<-curr.iTdT%*%curr.iT[kt,,]

curr.jack<-1-curr.rho[kt]+0.5*as.vector(matrix(curr.residual1,nrow=1)%*%(icurr.ome[kt,,]%x%curr.iTdTiT)%*%matrix(curr.residual1,ncol=1))-0.5*kk*sum(diag(curr.iTdT))
curr.G<-curr.rho[kt]+0.5*kk*sum(diag(curr.iTdTiT%*%curr.iT[kt,,]))
curr.iG<-1/curr.G
curr.3G<-curr.rho[kt]+kk*sum(diag(curr.iTdTiT%*%(curr.dT2[kt,,]-curr.dT[kt,,]%*%curr.iTdT)))

prop.mean<-log(curr.rho[kt])+0.5*eps.rho2[kt]*curr.jack*curr.iG-eps.rho2[kt]*(curr.iG^2)*curr.3G+0.5*eps.rho2[kt]*(curr.iG^2)*curr.3G
prop.var<-TEMP[kt]*eps.rho2[kt]*curr.iG

prop.rho<-exp(rnorm(n=1,mean=prop.mean,sd=sqrt(prop.var)))

if(prop.rho>0.001){
prop.T<-0*diag(sNi)
prop.iT<-prop.T
prop.dT<-prop.T
prop.dT2<-prop.T
prop.detT<-0
for(i in 1:6){
prop.T[marker==i,marker==i]<-exp(-prop.rho*T.array[marker==i,marker==i])
prop.iT[marker==i,marker==i]<-chol2inv(chol(prop.T[marker==i,marker==i]))
prop.dT[marker==i,marker==i]<--prop.rho*T.array[marker==i,marker==i]*prop.T[marker==i,marker==i]
prop.dT2[marker==i,marker==i]<-(prop.rho*T.array[marker==i,marker==i]-1)*prop.rho*T.array[marker==i,marker==i]*prop.T[marker==i,marker==i]
prop.detT<-prop.detT+determinant(prop.T[marker==i,marker==i])$modulus[1]}

prop.iTdT<-prop.iT%*%prop.dT
prop.iTdTiT<-prop.iTdT%*%prop.iT

prop.jack<-1-prop.rho+0.5*as.vector(matrix(curr.residual1,nrow=1)%*%(icurr.ome[kt,,]%x%prop.iTdTiT)%*%matrix(curr.residual1,ncol=1))-0.5*kk*sum(diag(prop.iTdT))
prop.G<-prop.rho+0.5*kk*sum(diag(prop.iTdTiT%*%prop.iT))
prop.iG<-1/prop.G
prop.3G<-prop.rho+kk*sum(diag(prop.iTdTiT%*%(prop.dT2-prop.dT%*%prop.iTdT)))

curr.mean<-log(prop.rho)+0.5*eps.rho2[kt]*prop.jack*prop.iG-eps.rho2[kt]*(prop.iG^2)*prop.3G+0.5*eps.rho2[kt]*(prop.iG^2)*prop.3G
curr.var<-TEMP[kt]*eps.rho2[kt]*prop.iG

top<-TEMP[kt]*(-0.5*sum(diag(icurr.ome[kt,,]%*%matrix(curr.residual1,nrow=kk,byrow=TRUE)%*%prop.iT%*%matrix(curr.residual1,ncol=kk)))-0.5*kk*prop.detT+log(prop.rho)-prop.rho)+dnorm(x=log(curr.rho[kt]),mean=curr.mean,sd=sqrt(curr.var),log=TRUE)
bot<-TEMP[kt]*(-0.5*sum(diag(icurr.ome[kt,,]%*%matrix(curr.residual1,nrow=kk,byrow=TRUE)%*%curr.iT[kt,,]%*%matrix(curr.residual1,ncol=kk)))-0.5*kk*curr.detT[kt]+log(curr.rho[kt])-curr.rho[kt])+dnorm(x=log(prop.rho),mean=prop.mean,sd=sqrt(prop.var),log=TRUE)
prob.rho<-exp(top-bot)
if(prob.rho>=runif(1)){
inner.rho[kt]<-1
curr.R[kt,,]<-curr.ome[kt,,]%x%prop.T
curr.iR[kt,,]<-icurr.ome[kt,,]%x%prop.iT
curr.T[kt,,]<-prop.T
curr.iT[kt,,]<-prop.iT
curr.dT[kt,,]<-prop.dT
curr.dT2[kt,,]<-prop.dT2
curr.detT[kt]<-prop.detT
curr.rho[kt]<-prop.rho}}

######################### omega ############################################################################

Shat<-(1/TEMP[kt])*(diag(kk)+matrix(curr.residual1,nrow=kk,byrow=TRUE)%*%curr.iT[kt,,]%*%matrix(curr.residual1,ncol=kk))
nuhat<-(nu+sNi)
nuhat<-((nuhat+kk+1)/TEMP[kt])-kk-1

curr.ome[kt,,]<-riwish(v=nuhat,S=Shat)
icurr.ome[kt,,]<-solve.cpp(curr.ome[kt,,])

curr.R[kt,,]<-curr.ome[kt,,]%x%curr.T[kt,,]
curr.iR[kt,,]<-icurr.ome[kt,,]%x%curr.iT[kt,,]

######################### psi ############################################################################

curr.residual2<-curr.dstar[kt,]-curr.m[kt,]

curr.iAdA1<-curr.iA[kt,,]%*%curr.dA[kt,1,,]
curr.iAdAiA1<-curr.iAdA1%*%curr.iA[kt,,]
curr.iAdA2<-curr.iA[kt,,]%*%curr.dA[kt,2,,]
curr.iAdAiA2<-curr.iAdA2%*%curr.iA[kt,,]
C1<-curr.iAdA1
C2<-curr.iAdA2

curr.jack<-c()
curr.jack[1]<-1-curr.psi[kt,1]-0.5*kk*sum(diag(curr.iAdA1))+0.5*as.vector(matrix(curr.residual2,nrow=1)%*%(icurr.sig[kt,,]%x%curr.iAdAiA1)%*%matrix(curr.residual2,ncol=1))
curr.jack[2]<-1-curr.psi[kt,2]-0.5*kk*sum(diag(curr.iAdA2))+0.5*as.vector(matrix(curr.residual2,nrow=1)%*%(icurr.sig[kt,,]%x%curr.iAdAiA2)%*%matrix(curr.residual2,ncol=1))
curr.G<-0*diag(2)
curr.G[1,1]<-0.5*kk*sum(diag(curr.iAdA1%*%curr.iAdA1))+curr.psi[kt,1]
curr.G[2,2]<-0.5*kk*sum(diag(curr.iAdA2%*%curr.iAdA2))+curr.psi[kt,2]
curr.G[1,2]<-0.5*kk*sum(diag(curr.iAdA1%*%curr.iAdA2))
curr.G[2,1]<-curr.G[1,2]
curr.iG<-solve.cpp(curr.G)

curr.Gj<-array(0,dim=c(2,2,2))
curr.Gj[1,1,1]<-curr.psi[kt,1]+kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,1,,]-C1%*%C1)%*%C1))
curr.Gj[1,1,2]<-0.5*kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,1,,]-C1%*%C1)%*%C2))+0.5*kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,3,,]-C1%*%C2)%*%C1))
curr.Gj[1,2,1]<-curr.Gj[1,1,2]
curr.Gj[1,2,2]<-kk*sum(diag(C2%*%(curr.iA[kt,,]%*%curr.dA2[kt,3,,]-C1%*%C2)))
curr.Gj[2,1,1]<-kk*sum(diag(C1%*%(curr.iA[kt,,]%*%curr.dA2[kt,3,,]-C2%*%C1)))
curr.Gj[2,1,2]<-0.5*kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,2,,]-C2%*%C2)%*%C1))+0.5*kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,3,,]-C2%*%C1)%*%C2))
curr.Gj[2,2,1]<-curr.Gj[2,1,2]
curr.Gj[2,2,2]<-curr.psi[kt,2]+kk*sum(diag((curr.iA[kt,,]%*%curr.dA2[kt,2,,]-C2%*%C2)%*%C2))

curr.bit1<-c()
curr.bit2<-c()
for(i in 1:2){
inner1<-c()
inner2<-c()
for(j in 1:2){
GLOB<-curr.iG%*%curr.Gj[j,,]
inner1[j]<-(GLOB%*%curr.iG)[i,j]
inner2[j]<-curr.iG[i,j]*sum(diag(GLOB))}
curr.bit1[i]<-sum(inner1)
curr.bit2[i]<-sum(inner2)}

prop.mean<-log(curr.psi[kt,])+0.5*eps.psi2[kt]*as.vector(curr.iG%*%matrix(curr.jack,ncol=1))-eps.psi2[kt]*curr.bit1+0.5*eps.psi2[kt]*curr.bit2
prop.var<-TEMP[kt]*eps.psi2[kt]*curr.iG
iprop.var<-curr.G/(TEMP[kt]*eps.psi2[kt])

prop.psi<-as.vector(exp(rmvnorm(n=1,mean=prop.mean,sigma=prop.var)))

prop.dA<-array(0,dim=c(2,mm,mm))
prop.dA2<-array(0,dim=c(3,mm,mm))
prop.A<-exp(-prop.psi[1]*A.array1-prop.psi[2]*A.array2)
prop.iA<-c()
try(prop.iA<-chol2inv(chol(prop.A)),silent=TRUE)
if(!is.null(prop.iA)){
prop.dA[1,,]<--prop.psi[1]*A.array1*prop.A
prop.dA[2,,]<--prop.psi[2]*A.array2*prop.A
prop.dA2[1,,]<-(prop.psi[1]*A.array1-1)*prop.psi[1]*A.array1*prop.A
prop.dA2[2,,]<-(prop.psi[2]*A.array2-1)*prop.psi[2]*A.array2*prop.A
prop.dA2[3,,]<-prop.psi[1]*prop.psi[2]*A.array1*A.array2*prop.A
prop.detA<-determinant(prop.A)$modulus[1]

prop.iAdA1<-prop.iA%*%prop.dA[1,,]
prop.iAdAiA1<-prop.iAdA1%*%prop.iA
prop.iAdA2<-prop.iA%*%prop.dA[2,,]
prop.iAdAiA2<-prop.iAdA2%*%prop.iA
P1<-prop.iAdA1
P2<-prop.iAdA2

prop.jack<-c()
prop.jack[1]<-1-prop.psi[1]-0.5*kk*sum(diag(prop.iAdA1))+0.5*as.vector(matrix(curr.residual2,nrow=1)%*%(icurr.sig[kt,,]%x%prop.iAdAiA1)%*%matrix(curr.residual2,ncol=1))
prop.jack[2]<-1-prop.psi[2]-0.5*kk*sum(diag(prop.iAdA2))+0.5*as.vector(matrix(curr.residual2,nrow=1)%*%(icurr.sig[kt,,]%x%prop.iAdAiA2)%*%matrix(curr.residual2,ncol=1))
prop.G<-0*diag(2)
prop.G[1,1]<-0.5*kk*sum(diag(prop.iAdA1%*%prop.iAdA1))+prop.psi[1]
prop.G[2,2]<-0.5*kk*sum(diag(prop.iAdA2%*%prop.iAdA2))+prop.psi[2]
prop.G[1,2]<-0.5*kk*sum(diag(prop.iAdA1%*%prop.iAdA2))
prop.G[2,1]<-prop.G[1,2]
prop.iG<-solve.cpp(prop.G)

prop.Gj<-array(0,dim=c(2,2,2))
prop.Gj[1,1,1]<-prop.psi[1]+kk*sum(diag((prop.iA%*%prop.dA2[1,,]-P1%*%P1)%*%P1))
prop.Gj[1,1,2]<-0.5*kk*sum(diag((prop.iA%*%prop.dA2[1,,]-P1%*%P1)%*%P2))+0.5*kk*sum(diag((prop.iA%*%prop.dA2[3,,]-P1%*%P2)%*%P1))
prop.Gj[1,2,1]<-prop.Gj[1,1,2]
prop.Gj[1,2,2]<-kk*sum(diag(P2%*%(prop.iA%*%prop.dA2[3,,]-P1%*%P2)))
prop.Gj[2,1,1]<-kk*sum(diag(P1%*%(prop.iA%*%prop.dA2[3,,]-P2%*%P1)))
prop.Gj[2,1,2]<-0.5*kk*sum(diag((prop.iA%*%prop.dA2[2,,]-P2%*%P2)%*%P1))+0.5*kk*sum(diag((prop.iA%*%prop.dA2[3,,]-P2%*%P1)%*%P2))
prop.Gj[2,2,1]<-prop.Gj[2,1,2]
prop.Gj[2,2,2]<-prop.psi[2]+kk*sum(diag((prop.iA%*%prop.dA2[2,,]-P2%*%P2)%*%P2))

prop.bit1<-c()
prop.bit2<-c()
for(i in 1:2){
inner1<-c()
inner2<-c()
for(j in 1:2){
GLOB<-prop.iG%*%prop.Gj[j,,]
inner1[j]<-(GLOB%*%prop.iG)[i,j]
inner2[j]<-prop.iG[i,j]*sum(diag(GLOB))}
prop.bit1[i]<-sum(inner1)
prop.bit2[i]<-sum(inner2)}

curr.mean<-log(prop.psi)+0.5*eps.psi2[kt]*as.vector(prop.iG%*%matrix(prop.jack,ncol=1))-eps.psi2[kt]*prop.bit1+0.5*eps.psi2[kt]*prop.bit2
curr.var<-TEMP[kt]*eps.psi2[kt]*prop.iG
icurr.var<-prop.G/(TEMP[kt]*eps.psi2[kt])

top<-TEMP[kt]*(-0.5*sum(diag(icurr.sig[kt,,]%*%matrix(curr.residual2,nrow=kk,byrow=TRUE)%*%prop.iA%*%matrix(curr.residual2,ncol=kk)))-0.5*kk*prop.detA+sum(log(prop.psi)-prop.psi))+dmvnorm(x=log(curr.psi[kt,]),mean=curr.mean,sigma=curr.var,log=TRUE)
bot<-TEMP[kt]*(-0.5*sum(diag(icurr.sig[kt,,]%*%matrix(curr.residual2,nrow=kk,byrow=TRUE)%*%curr.iA[kt,,]%*%matrix(curr.residual2,ncol=kk)))-0.5*kk*curr.detA[kt]+sum(log(curr.psi[kt,])-curr.psi[kt,]))+dmvnorm(x=log(prop.psi),mean=prop.mean,sigma=prop.var,log=TRUE)
prob.psi<-exp(top-bot)
if(prob.psi>=runif(1)){
inner.psi[kt]<-1
curr.S[kt,,]<-curr.sig[kt,,]%x%prop.A
curr.iS[kt,,]<-icurr.sig[kt,,]%x%prop.iA
curr.A[kt,,]<-prop.A
curr.iA[kt,,]<-prop.iA
curr.dA[kt,,,]<-prop.dA
curr.dA2[kt,,,]<-prop.dA2
curr.detA[kt]<-prop.detA
curr.psi[kt,]<-prop.psi}}

######################## sigma ##########################################################################

Shat<-(1/TEMP[kt])*(diag(kk)+matrix(curr.residual2,nrow=kk,byrow=TRUE)%*%curr.iA[kt,,]%*%matrix(curr.residual2,ncol=kk))
nuhat<-(nu+mm)
nuhat<-((nuhat+kk+1)/TEMP[kt])-kk-1

curr.sig[kt,,]<-riwish(v=nuhat,S=Shat)
icurr.sig[kt,,]<-solve.cpp(curr.sig[kt,,])

curr.S[kt,,]<-curr.sig[kt,,]%x%curr.A[kt,,]
curr.iS[kt,,]<-icurr.sig[kt,,]%x%curr.iA[kt,,]

######################### missing terms ##################################################################

for(ii in 1:6){
curr.y5<-curr.y[rep(marker,kk)==ii]
marker5<-(1:length(curr.y5))[as.vector(Y[marker==ii,])==(-Inf)]
curr.mu5<-curr.eta[kt,rep(marker,kk)==ii]
curr.SIG5<-curr.ome[kt,,]%x%curr.T[kt,marker==ii,marker==ii]
for(jj in 1:length(marker5)){
invhere<-solve.cpp(curr.SIG5[-marker5[jj],-marker5[jj]])
muD<-curr.mu5[marker5[jj]]+as.vector(curr.SIG5[marker5[jj],-marker5[jj]]%*%invhere%*%matrix(curr.y5[-marker5[jj]]-curr.mu5[-marker5[jj]],ncol=1))
sigD<-curr.SIG5[marker5[jj],marker5[jj]]-as.vector(curr.SIG5[marker5[jj],-marker5[jj]]%*%invhere%*%curr.SIG5[-marker5[jj],marker5[jj]])
curr.y5[marker5[jj]]<-rtnorm(n=1,mean=muD,sd=sqrt(TEMP[kt]*sigD),upper=c)}
curr.y[rep(marker,kk)==ii]<-curr.y5}

curr.missing[kt,]<-curr.y[missing]


}














