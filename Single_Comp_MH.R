###Implementation of single component MH

###Proposals Binomial(n,0.5) for the
#first component and a uniform(0,1) for the second component


n=5; al=2;bt=1;nsim=20000
x=c();y=c();x[1]=1;y[1]=0.5
for(j in 2:nsim){
  #first component update
  z1=rbinom(1,n,0.5)
  num=(y[j-1]^z1)*((1-y[j-1])^(n-z1))
  den=(y[j-1]^x[j-1])*((1-y[j-1])^(n-x[j-1]))
  prob1=min(num/den,1)
  binom=rbinom(1,1,prob1)
  x[j]=ifelse(binom==1,z1,x[j-1])
  #second component update
  z2=runif(1)
  num= (z2^(x[j-1]+al-1))*((1-z2)^(n-x[j-1]+bt-1))
  den= (y[j-1]^(x[j-1]+al-1))*((1-y[j-1])^(n-x[j-1]+bt-1))
  prob2=min(num/den,1)
  binom=rbinom(1,1,prob2)
  y[j]=ifelse(binom==1,z2,y[j-1])}
z=cbind(x,y)
hist(x); hist(y)

##Change point problem
###simulate data
n=150; k=50; mu=4;lambda=15
x=c(rpois(k,mu),rpois((n-k),lambda))
hist(x)


####Priors: Gamma(a1,b1), Gamma(a2,b2),discrete uniform
#prior parameters
a=10;b=1;
##proposals: Gamma, Gamma, discrete unif
rdunif=function(nsim,a,b){floor(runif(nsim,a,b+1))}
#test=rdunif(500,0,10)
#table(test)
n=length(x)
u=c(); v=c();w=c()
u[1]=1;v[1]=2;w[1]=round(n/2)
for(j in 2:10000){
  #update for first component
  u[j]=rgamma(1,shape=a+sum(x[1:w[j-1]]), rate=(w[j-1]+b))
  #update for second component
  v[j]=rgamma(1,shape=a+sum(x[-(1:w[j-1])]), rate=((n-w[j-1]+b)))
  #update for third component
  #proposal 
  q=rdunif(1,1,n); cr=w[j-1]; mu=u[j-1];la=v[j-1]
  sterm1=(sum(x[1:q])-sum(x[1:cr])); 
  sterm2=sum(x[-(1:q)])-sum(x[-(1:cr)])
  sterm3=((q-cr)*mu)+((cr-q)*la)
  logp1=(sterm1*log(mu))+(sterm2*log(la))-sterm3
  rho=min(exp(logp1),1)
  binom=rbinom(1,1,rho)
  w[j]=ifelse(binom==1,q,w[j-1])}
chain=cbind(u,v,w)
chain.burn=chain[-(1:5000),]
estimates=colMeans(chain.burn)
estimates
###Alternative way: frequentist approach,
#any arbitrary starting point
sp=125
imu1=mean(x[1:sp]); imu2=mean(x[-(1:sp)])
loss=c()
for(j in 1:(n-1)){
  t=j
  loss[j]=sum((x[1:t]-imu1)^2)+sum((x[-(1:t)]-imu2)^2)}

plot(loss,type="l")
which.min(loss)

















































































