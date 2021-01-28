library("LearnBayes")
library("Rlab")
MK=as.matrix(read.csv("matureknee.txt"))
IK=as.matrix(read.csv("immatureknee.txt"))
n=length(MK)+length(IK)

#Merge data to one matrix
K=matrix(nrow=n,ncol=2)
K[,1]=rbind(MK,IK)
K[1:{n/2},2]=rep(1,n/2)
K[{n/2+1}:n,2]=rep(0,n/2)

x=K[,1] #age
y=K[,2] #report, Mature=1, Immature=0

#2a
####################################
#Lf_log = -log (Likelihood)
Lf_log <- function(z) { 
  a <- z[1]
  b <- z[2]
  -sum(y*log(exp(a+b*(x-18))/(1+exp(a+b*(x-18))))+(1-y)*log((1-exp(a+b*(x-18))/(1+exp(a+b*(x-18))))))
  }
max_Lf=nlm(Lf_log,c(1,1)) #minimizing(-log(Lf))
a_max=max_Lf$estimate[1] 
b_max=max_Lf$estimate[2]
x0=seq(12,26,0.1)
p=exp(a_max+b_max*(x0-18))/(1+exp(a_max+b_max*(x0-18))) #calculating logistic regression curve
plot(x,y,xlab='age')
lines(x0,p,type="l",col = "blue",lty="dashed")

#2b plotting posterior
####################################
Lf_log2 <- function(a,b,x,y) {sum(y*log(exp(a+b*(x-18))/(1+exp(a+b*(x-18))))+(1-y)*log((1-exp(a+b*(x-18))/(1+exp(a+b*(x-18))))))}
s=0
a=seq(-.5,2,length.out = 21) #grid 21x21
b=seq(.5,3,length.out = 21)
Lf=matrix(nrow=length(a),ncol=length(b))
k=0
#for-loop for combinations of a and b
for (ai in a) {
  k=k+1
  l=0
  for (bi in b) {
    l=l+1
    Lf[k,l]=exp(Lf_log2(ai,bi,x,y)) #Lf(a,b)
  }
}
image(a,b,Lf)


#2c - cost diff when a=a_max, b=b_max
####################################
cost <- function(a,b,c,mu,alp) {
  B=10
  g <- function(x) {dgamma(x-14,alp,alp/(mu-14))} #age distribution
  f <- function(x) {1/(exp(-(a+b*(x-18)))+1)} #f(x)
  if (c==1) {
    c_low <- function(x) {B}
    c_up <- function(x) {1}
  }
  else {
    c_low <- function(x) {B*(18-x)}
    c_up <- function(x) {x-18}
  }
  f_g_c_low <- function(x) {f(x)*g(x)*c_low(x)} #cost function for x<=18
  f_g_c_up <- function(x) {f(x)*g(x)*c_up(x)} #cost function for x>18
  c_c=integrate(f_g_c_low,0,18)$value #cost of misclass to child
  c_a=integrate(f_g_c_up,18,200)$value #cost of misclass to adult
return(c_c-c_a)
}

q1=matrix(nrow=2,ncol=3) #cost diff for different combinations of age and c(x)
mu=c(18.5,19.5,20.58)
alp=c(3,6,3)
for (i in 1:2) {
  for (j in 1:3) {
    q1[i,j]=cost(a_max,b_max,i,mu[j],alp[j])
  }
}
q1

#2d - cost diff when averaging over a and b
####################################
q2=matrix(nrow=2,ncol=3) #cost diff for different combinations of age and c(x)
mu=c(18.5,19.5,20.58)
alp=c(3,6,3)

for (i in 1:2) {
  for (j in 1:3) {
    n=0
    q_ab=0
    for (a in seq(-.5,2,.01)) {
      for (b in seq(.5,3,.01)) {
      n=n+1
      q_ab=q_ab+cost(a,b,i,mu[j],alp[j]) #sum for cost diff for all comb of a and b
      }
    }
    q2[i,j]=q_ab/n  #resulting cost diff for a given age distr and c(x)
  }
}
q2







