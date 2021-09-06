### cauchy distribution

x=c(-4.8,-2.8,-1.35,-0.02,0.70,0.98,2.92,5.50)
cauchy.l=function(alpha){
  y=-sum(log(0.01+(x-alpha)^2))
  return(y)
}
#cauchy.l(0.5)

al=seq(-6,6,length.out=500)
y=sapply(as.list(al),cauchy.l)
plot(al,y)

###implement Simulated Annealing to find the MLE of alpha

set.seed(123)
lik=cauchy.l
r=0.5
alpha.old=0
N=15000
Tt=function(t){y=1/((1+t)^{1/2});return(y)}
updates=c()

for(k in 1:N){
  #Step 1:
  alpha.can=runif(1,alpha.old-r,alpha.old+r)
  ###Evaluate the objective function at the old 
  #and the candidate value
  f.can=lik(alpha.can); f.old=lik(alpha.old)
  ###Step 2
  if(f.can>f.old){alpha.new=alpha.can}else{rho=exp((f.can-f.old)/Tt(k))
  b=rbinom(1,1,rho)
  if(b==1){alpha.new=alpha.can}
  if(b==0){alpha.new=alpha.old}}
  updates=c(updates,alpha.new)
  alpha.old=alpha.new}

plot(1:N, updates,type="l")


####Find the optimizer of the function
#h(x)=(cos(50x)+sin(50x))^2
##x>0
h=function(x){(cos(50*x)+sin(50*x))^2}

x=seq(0,1,length.out=100)
y=h(x)
plot(x,y,type="l")


set.seed(123)
lik=h
r=0.5
alpha.old=0
N=15000
Tt=function(t){y=1/((1+t)^2);return(y)}
updates=c()

for(k in 1:N){
  #Step 1:
  alpha.can=runif(1,max(0,alpha.old-r),alpha.old+r)
  ###Evaluate the objective function at the old 
  #and the candidate value
  f.can=lik(alpha.can); f.old=lik(alpha.old)
  ###Step 2
  if(f.can>f.old){alpha.new=alpha.can}else{rho=exp((f.can-f.old)/Tt(k))
  b=rbinom(1,1,rho)
  if(b==1){alpha.new=alpha.can}
  if(b==0){alpha.new=alpha.old}}
  updates=c(updates,alpha.new)
  alpha.old=alpha.new}

plot(1:N, updates,type="l")


