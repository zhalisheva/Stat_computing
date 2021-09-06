###EM Algorithm for mixture of normals

###simulate data
n=500; pr= 0.4; mu1=-2; mu2=2; si1=0.25; si2=0.25
z1=rbinom(n,1,pr); z2=1-z1
x1=rnorm(n,mu1,si1);x2=rnorm(n,mu2,si2)
x=(z1*x1)+(z2*x2);
hist(x)


f=function(x,mu,si){
  fx=(1/(si*sqrt(2*pi)))*exp(-(x-mu)^2/(2*si^2))
  return(fx)}
f(3,0,1)

n=length(x)
oldpr=0.5; oldmu1=0; oldmu2=0.5; oldsi1sq=1;oldsi2sq=1

maxits=1000
its=0
tol=1e-5
err=100
while(err>tol & its<maxits){
  num=(oldpr*f(x,oldmu1,sqrt(oldsi1sq)))
  den=(oldpr*f(x,oldmu1,sqrt(oldsi1sq)))+((1-oldpr)*f(x,oldmu2,sqrt(oldsi2sq)))
  oldez1=num/den
  #oldez1  
  oldez2=1-oldez1
  newpr=sum(oldez1)/n
  newmu1=sum(oldez1*x)/sum(oldez1)
  newmu2=sum(oldez2*x)/sum(oldez2)
  newsi1sq=sum(oldez1*((x-newmu1)^2))/sum(oldez1)
  newsi2sq=sum(oldez2*((x-newmu2)^2))/sum(oldez2)
  err=max(abs(c(newpr-oldpr, newmu1-oldmu1,newmu2-oldmu2,
                newsi1sq-oldsi1sq,newsi2sq-oldsi2sq)))
  its=its+1
  oldpr=newpr; oldmu1=newmu1;oldmu2=newmu2;
  oldsi1sq=newsi1sq; oldsi2sq=newsi2sq}

newpr
newmu1
newmu2
newsi1sq
newsi2sq

esi1=sqrt(newsi1sq)
esi1

