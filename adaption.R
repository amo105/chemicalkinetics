

adaption<-function(rhos,u,a.star=0.5,iterations,zeta=0.6,c=1){

K<-length(rhos)+1

T<-1
for(jj in 2:K){
T[jj]<-prod(exp(exp(rhos[1:(jj-1)])))}

new.rhos<-c()
for(jj in 1:(K-1)){
new.rhos[jj]<-rhos[jj]+(c*((iterations+1)^(-zeta)))*(min(1,exp((u(jj+1)-u(jj))*((1/T[jj])-(1/T[jj+1]))))-a.star)
new.rhos[jj]<-ifelse(new.rhos[jj]<(-0.8),new.rhos[jj],-0.8)}

new.rhos}

rhos2T<-function(rhos){

K<-length(rhos)+1

T<-1
for(jj in 2:K){
T[jj]<-prod(exp(exp(rhos[1:(jj-1)])))}

T}