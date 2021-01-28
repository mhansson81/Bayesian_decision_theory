library("LearnBayes")
S=as.matrix(read.csv("stockvalues.txt"))
n_st=ncol(S) #number of stocks
n_days=nrow(S) #number of days

Z=matrix(nrow=n_days-1,ncol=n_st) #matrix for log(return)
Z_mean=numeric(n_st)
for (i in 1:{n_days-1}) {
  Z[i,]=log(S[i+1,]/S[i,])
}
for (i in 1:n_st) {
  Z_mean[i]=mean(Z[,i]) #gamma value for multinormal dist
}
Z_cov=cov(Z) #covariance of log(Z)

#1a ##########################################
S=1000
n=100
V <- rmnorm(S, n*Z_mean, n*Z_cov) #generating 1000 V-vectors

#1b ##########################################
#u=function computing utility
u=function(w,k,V) {
  s=matrix(nrow=nrow(V),ncol=1,0)
  for (i in 1:length(w)) {
    s=s+w[i]*exp(V[,i])
  }
  u=(1-s^-k)/k
  return(u)
}

#Case 1 - equal weight
w=rep(1/n_st,n_st) #equal weight
k=c(-.5,.5,1.5) #k-values
u1=u(w,k[1],V) #u for k1
E1=sum(u1)/S #Expected u when simulating S V-vectors
u2=u(w,k[2],V)
E2=sum(u2)/S
u3=u(w,k[3],V)
E3=sum(u3)/S

#Case 2 - "best" stock only
E_s=matrix(ncol=n_st,nrow=3)

for (i in 1:n_st) {
  w=numeric(n_st) #all w=0 except for one stock (at a time)
  w[i]=1
  u1i=u(w,k[1],V) #u for k1 and w_i
  E_s[1,i]=sum(u1i)/S
  u2i=u(w,k[2],V)
  E_s[2,i]=sum(u2i)/S
  u3i=u(w,k[3],V)
  E_s[3,i]=sum(u3i)/S
}

#1c ##########################################
Ef <- function(w3,k0,V0,S0) {
  -sum(u(c(w3,1-w3),k0,V0))/S0
}
knew=c(-.5,1.5) #k-values of interest
Vnew=V[,3:4] #V vectors with only S3 and S4

w3_1=nlm(Ef,.5,V0=Vnew,k0=knew[1],S0=S) #optimizing over w3 -->w4=1-w3
w3_2=nlm(Ef,.5,V0=Vnew,k0=knew[2],S0=S)











