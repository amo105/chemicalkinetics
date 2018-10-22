nu<-3

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

