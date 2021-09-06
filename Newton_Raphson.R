## Newton raphson for logistic regression 
###Generate data
set.seed(123)
n=300; x=rnorm(n); b=c(0.25,0.75)
X=cbind(1,x)
invlogit=function(X,b){p=exp(X%*%b)/(1+exp(X%*%b))
return(p)}
p=invlogit(X,b)
y=c(); for(i in 1:n){y[i]=rbinom(1,1,p[i])}

###use (x,y) to estimate the b coefficients.

#for any matrix H: H^-1=solve(H)
p=dim(X)[2] ##measure the dimension of the model 
old_hb=rep(0,p) ##arbitrary starting point
n=length(y)

V=matrix(0,n,n)
tol=1e-5
maxits=5000
its=0; err=10 

while(err>tol & its<maxits){
  old_hp=invlogit(X,old_hb)
  diag(V)=old_hp*(1-old_hp)
  H=-t(X)%*%V%*%X; Hinv=solve(H)
  S=t(X)%*%(y-old_hp)
  new_hb=old_hb-Hinv%*%S
  err=max(abs(new_hb-old_hb))
  its=its+1
  old_hb=new_hb
}
new_hb


for(j in 1:100){
  old_hp=invlogit(X,old_hb)
  diag(V)=old_hp*(1-old_hp)
  H=-t(X)%*%V%*%X; Hinv=solve(H)
  S=t(X)%*%(y-old_hp)
  new_hb=old_hb-Hinv%*%S
  err=max(abs(new_hb-old_hb))
  its=its+1
  old_hb=new_hb
}




#new_hb
