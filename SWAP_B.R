nu<-3

L1<-sample(x=1:(KT-1),size=1)
L2<-L1+1

UUU<-function(s){
curr.y1<-y
curr.y1[missing]<-curr.missing[s,]

curr.iR1<-kronecker(icurr.ome[s,,],curr.iT[s,,])
curr.iS1<-kronecker(icurr.sig[s,,],curr.iA[s,,])

EVAL_LIK1<--0.5*sNi*determinant(curr.ome[s,,])$modulus[1]-0.5*kk*curr.detT[s]-0.5*as.vector(matrix(curr.y1-curr.eta[s,],nrow=1)%*%curr.iR1%*%matrix(curr.y1-curr.eta[s,],ncol=1))
EVAL_PRIORA<--0.5*mm*determinant(curr.sig[s,,])$modulus[1]-0.5*kk*curr.detA[s]-0.5*as.vector(matrix(curr.dstar[s,]-curr.m[s,],nrow=1)%*%curr.iS1%*%matrix(curr.dstar[s,]-curr.m[s,],ncol=1))
EVAL_PRIOR1<-EVAL_PRIORA+sum(dnorm(x=curr.theta[s,],mean=prior.loc,sd=sqrt(prior.var),log=TRUE))-curr.rho[s]-0.5*(nu+kk+1)*determinant(curr.ome[s,,])$modulus[1]-0.5*sum(diag(icurr.ome[s,,]))-sum(curr.psi[s,])-0.5*(nu+kk+1)*determinant(curr.sig[s,,])$modulus[1]-0.5*sum(diag(icurr.sig[s,,]))

EVAL1<-EVAL_LIK1+EVAL_PRIOR1
EVAL1}

EVAL1<-UUU(L1)
EVAL2<-UUU(L2)

prob.swap<-exp(EVAL1/TEMP[L2]+EVAL2/TEMP[L1]-EVAL1/TEMP[L1]-EVAL2/TEMP[L2])
if(prob.swap>=runif(1)){

inter<-curr.theta[L1,]
curr.theta[L1,]<-curr.theta[L2,]
curr.theta[L2,]<-inter

inter<-curr.dstar[L1,]
curr.dstar[L1,]<-curr.dstar[L2,]
curr.dstar[L2,]<-inter

inter<-curr.missing[L1,]
curr.missing[L1,]<-curr.missing[L2,]
curr.missing[L2,]<-inter

inter<-curr.ome[L1,,]
curr.ome[L1,,]<-curr.ome[L2,,]
curr.ome[L2,,]<-inter
inter<-icurr.ome[L1,,]
icurr.ome[L1,,]<-icurr.ome[L2,,]
icurr.ome[L2,,]<-inter

inter<-curr.sig[L1,,]
curr.sig[L1,,]<-curr.sig[L2,,]
curr.sig[L2,,]<-inter
inter<-icurr.sig[L1,,]
icurr.sig[L1,,]<-icurr.sig[L2,,]
icurr.sig[L2,,]<-inter

inter<-curr.rho[L1]
curr.rho[L1]<-curr.rho[L2]
curr.rho[L2]<-inter

inter<-curr.T[L1,,]
curr.T[L1,,]<-curr.T[L2,,]
curr.T[L2,,]<-inter
inter<-curr.iT[L1,,]
curr.iT[L1,,]<-curr.iT[L2,,]
curr.iT[L2,,]<-inter
inter<-curr.dT[L1,,]
curr.dT[L1,,]<-curr.dT[L2,,]
curr.dT[L2,,]<-inter
inter<-curr.dT2[L1,,]
curr.dT2[L1,,]<-curr.dT2[L2,,]
curr.dT2[L2,,]<-inter
inter<-curr.detT[L1]
curr.detT[L1]<-curr.detT[L2]
curr.detT[L2]<-inter

inter<-curr.R[L1,,]
curr.R[L1,,]<-curr.R[L2,,]
curr.R[L2,,]<-inter
inter<-curr.iR[L1,,]
curr.iR[L1,,]<-curr.iR[L2,,]
curr.iR[L2,,]<-inter

inter<-curr.psi[L1,]
curr.psi[L1,]<-curr.psi[L2,]
curr.psi[L2,]<-inter

inter<-curr.A[L1,,]
curr.A[L1,,]<-curr.A[L2,,]
curr.A[L2,,]<-inter
inter<-curr.iA[L1,,]
curr.iA[L1,,]<-curr.iA[L2,,]
curr.iA[L2,,]<-inter
inter<-curr.dA[L1,,,]
curr.dA[L1,,,]<-curr.dA[L2,,,]
curr.dA[L2,,,]<-inter
inter<-curr.dA2[L1,,,]
curr.dA2[L1,,,]<-curr.dA2[L2,,,]
curr.dA2[L2,,,]<-inter
inter<-curr.detA[L1]
curr.detA[L1]<-curr.detA[L2]
curr.detA[L2]<-inter

inter<-curr.S[L1,,]
curr.S[L1,,]<-curr.S[L2,,]
curr.S[L2,,]<-inter
inter<-curr.iS[L1,,]
curr.iS[L1,,]<-curr.iS[L2,,]
curr.iS[L2,,]<-inter

inter<-curr.m[L1,]
curr.m[L1,]<-curr.m[L2,]
curr.m[L2,]<-inter
inter<-curr.eta[L1,]
curr.eta[L1,]<-curr.eta[L2,]
curr.eta[L2,]<-inter
inter<-curr.dm[L1,,]
curr.dm[L1,,]<-curr.dm[L2,,]
curr.dm[L2,,]<-inter
inter<-curr.dm2[L1,,,]
curr.dm2[L1,,,]<-curr.dm2[L2,,,]
curr.dm2[L2,,,]<-inter}