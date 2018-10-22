#library(tseries)


load("ode_data.RDATA")

resp<-read.matrix(file="resp_form.txt",header=FALSE)[,-1]
MARKER<-rep(marker,3)

N<-6
Ni<-table(marker)
kk<-3
n<-sum(Ni)*3
sNi<-sum(Ni)
mm<-93

y<-log(resp)

where<-rep(0,sNi)
where2<-rep(0,sNi)
for(i in 1:sNi){
if(any(y[i,]==(-Inf))){
where[i]<-1
where2[i]<-(1:kk)[y[i,]==(-Inf)]}}
where3<-rep(0,n)
where3[1:109]<-ifelse(where2==1,1,0)
where3[110:218]<-ifelse(where2==2,1,0)
where3[219:327]<-ifelse(where2==3,1,0)

c<-log(0.01)

G<-matrix(0,nrow=sNi,ncol=mm)
G[1:91,1:91]<-diag(91)
G[92:107,74:89]<-diag(16)
G[108:109,92:93]<-diag(2)
X<-kronecker(diag(kk),G)

ids<-c()
for(i in 1:5){
ids[i]<-min((1:109)[marker==i])}

DATASET<-as.matrix(dataset)

dataset2<-dataset[c(1:91,108,109),]

t.max<-max(dataset$time)
t.min<-min(dataset$time)
#t.ran<-t.max-t.min
t.ran<-1

T.array<-matrix(Inf,ncol=sNi,nrow=sNi)
for(i in 1:sNi){
for(j in 1:sNi){
if(marker[i]==marker[j]){
T.array[i,j]<-abs((log(dataset$time[i])-log(dataset$time[j]))/t.ran)^2}}}

dataset2<-dataset2[,-3]
minss<-apply(dataset2[,-6],2,min)
maxss<-apply(dataset2[,-6],2,max)
minss[6]<-0
maxss[6]<-1
dataset2[,6]<-log(dataset2[,6])
dataset3<-0*dataset2
for(i in 1:6){
dataset3[,i]<-(dataset2[,i]-minss[i])/(maxss[i]-minss[i])}
dataset4<-matrix(0,ncol=6,nrow=dim(dataset3)[1])
dataset4[,1]<-dataset3[,1]
dataset4[,2]<-dataset3[,2]
dataset4[,3]<-dataset3[,3]
dataset4[,4]<-dataset3[,4]
dataset4[,5]<-dataset3[,5]
dataset4[,6]<-dataset3[,6]
dataset3<-dataset4

A.array1<-matrix(Inf,ncol=mm,nrow=mm)
for(i in 1:mm){
for(j in 1:mm){
A.array1[i,j]<-sum(abs(dataset3[i,-6]-dataset3[j,-6])^2)}}

A.array2<-matrix(Inf,ncol=mm,nrow=mm)
for(i in 1:mm){
for(j in 1:mm){
A.array2[i,j]<-abs(dataset3[i,6]-dataset3[j,6])^2}}








