## -----------------------------------------------------------------------------
h=rep(0,10000)
for (i in 1:10000)
{
  x=sample(1:100,100)
  for (j in 1:100)
  {
    if (x[j]==j)
      h[i]=h[i]+1
  }
}

## -----------------------------------------------------------------------------
(ex=mean(h))
(var=var(h))

## -----------------------------------------------------------------------------
my.sample<-function()
{
  sample(c(1,2,4,10),size=1e5,replace=TRUE,prob=1:4/10)
}
data<-my.sample()
d<-as.vector(table(data))
p<-1:4/10
d/1e5/p

## -----------------------------------------------------------------------------
set.seed(12345)
n<-1000
u<-runif(n)
x<-log(2*u)
for (i in 1000)
{
  if(u[i]>1/2)
  {
    x[i]=-log(2-2*u[i])
  }
}

## -----------------------------------------------------------------------------
hist(x,prob=TRUE, main=expression(f(x)==e^{-abs(x)}/2))
y<-seq(-10,10,0.01)
lines(y,exp(-abs(y))/2)

## -----------------------------------------------------------------------------
Beta<-function(a,b,n)
{
  x<-rep(0,n)
  k<-0
  while(k<n)
  {
    u<-runif(1)
    y<-runif(1)
    if(y^(a-1)*(1-y)^(b-1)>u)
    {
      k=k+1
      x[k]<-y
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
x<-Beta(3,2,1000)
hist(x,prob=TRUE, main=expression(beta(3,2)))
y<-seq(-10,10,0.01)
lines(y,gamma(5)/gamma(2)/gamma(3)*y^2*(1-y))

## -----------------------------------------------------------------------------
ek<-function(n)
{
  x<-rep(0,n)
  for (i in 1:n)
  {
    x[i]=u[3]
    u<-runif(3,-1,1)
    if(abs(u[3])>abs(u[2])&&abs(u[3])>abs(u[1]))
    {
      x[i]=u[2]
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
x<-ek(1000)
hist(x,prob=TRUE,ylim=c(0,0.8))
y<-seq(-1,1,0.01)
lines(y,3*(1-y^2)/4)

## -----------------------------------------------------------------------------
rho<-c(0.5,0.8,1)
variance<-function(r)
{
  p<-rep(0,100)
  for (i in 1:100)
  {
    m<-1e6
    x<-runif(m,0,1/r)
    Y<-runif(m,0,pi/2)
    pihat<-2*r/mean(sin(Y)>x)
    p[i]=pihat
  }
  return(var(p))
}

## -----------------------------------------------------------------------------
(c(variance(0.5),variance(0.8),variance(1)))

## -----------------------------------------------------------------------------
x1<-rep(0,100)
x2<-rep(0,100)
for (i in 1:100)
{
  m<-1e6
  U<-runif(m,0,1)
  U1<-U[1:m/2]
  x1[i]<-mean(exp(U))
  x2[i]<-(mean(exp(U1))+mean(exp(1-U1)))/2
}
((var(x1)-var(x2))/var(x1))

## -----------------------------------------------------------------------------
x<-seq(1,10,0.01)
y<-x^2*exp(-x^2/2)/sqrt(2*pi)
plot(x,y,type='l',ylim=c(0,1))
lines(x,2*dnorm(x,1),lty=2)
lines(x,dgamma(x-1,3/2,2),lty=3)
legend("topright",inset=0.02,legend=c("g(x)","f1","f2"),lty=1:3)

## -----------------------------------------------------------------------------
plot(x,y/(dgamma(x-1,1.5,2)),type="l",lty=3,ylab="")
lines(x,y/(2*dnorm(x,1)),lty=2)
legend("topright",inset=0.02,legend=c("f1","f2"),lty=2:3)

## -----------------------------------------------------------------------------
e1<-rep(0,1000)
for (i in 1:1000) 
{
x<-sqrt(rchisq(10000,1))+1
f<-2*dnorm(x,1)
g<-x^2*exp(-x^2/2)/sqrt(2*pi)
e1[i]<-mean(g/f)
}
e2<-rep(0,1000)
for (i in 1:1000) 
{
x<-rgamma(10000,1.5,2)+1
f<-dgamma(x-1,1.5,2)
g<-x^2*exp(-x^2/2)/sqrt(2*pi)
e2[i]<-mean(g/f)
}

## -----------------------------------------------------------------------------
c(mean(e1),mean(e2))
c(var(e1),var(e2))

## -----------------------------------------------------------------------------
M<-10000
k<-5
m<-M/k
e<-rep(0,k)
v<-rep(0,k)
g<-function(x)
   exp(-x)/(1+x^2)
f<-function(x) 
   k/(1-exp(-1))*exp(-x)

## -----------------------------------------------------------------------------
for (i in 1:k)
{
  u<-runif(m,(i-1)/k,i/k)
  x<--log(1-(1-exp(-1))*u)
  d<-g(x)/f(x)
  e[i]<-mean(d)
  v[i]<-var(d)
}
sum(e);
mean(v);
sqrt(mean(v))

## -----------------------------------------------------------------------------
M<-10000
k<-1
e<-0
v<-0
u<-runif(M)
x<--log(1-(1-exp(-1))*u)
d<-g(x)/f(x)
(e<-mean(d))
(v<-var(d))#variance of estimator without classification
(sqrt(v))

## -----------------------------------------------------------------------------
n<-20
t<-qt(0.975,df=n-1)
L<-rep(0,10000)
U<-rep(0,10000)
for (i in 1:10000)
{
  x<-rchisq(n,2)
  L[i]<-mean(x)-t*sd(x)/sqrt(n)
  U[i]<-mean(x)+t*sd(x)/sqrt(n)
}
sum(L<2&U>2)
mean(L<2&U>2)

## -----------------------------------------------------------------------------
n<-20
UV<-rep(0,10000)
for (i in 1:10000)
{
  x<-rchisq(n,2)
  UV[i]<-(n-1)*var(x)/qchisq(0.05,n-1)
}
sum(UV>4)
mean(UV>4)

## -----------------------------------------------------------------------------
n<-20
t<-qt(0.975,df=n-1)
L<-rep(0,10000)
U<-rep(0,10000)
for (i in 1:10000)
{
  x<-rchisq(n,1)
  L[i]<-mean(x)-t*sd(x)/sqrt(n)
  U[i]<-mean(x)+t*sd(x)/sqrt(n)
}
sum(L>1|U<1)
mean(L>1|U<1)

## -----------------------------------------------------------------------------
n<-20
t<-qt(0.975,df=n-1)
L<-rep(0,10000)
U<-rep(0,10000)
for (i in 1:10000)
{
  x<-runif(n,0,2)
  L[i]<-mean(x)-t*sd(x)/sqrt(n)
  U[i]<-mean(x)+t*sd(x)/sqrt(n)
}
sum(L>1|U<1)
mean(L>1|U<1)

## -----------------------------------------------------------------------------
n<-20
t<-qt(0.975,df=n-1)
L<-rep(0,10000)
U<-rep(0,10000)
for (i in 1:10000)
{
  x<-rexp(n,1)
  L[i]<-mean(x)-t*sd(x)/sqrt(n)
  U[i]<-mean(x)+t*sd(x)/sqrt(n)
}
sum(L>1|U<1)
mean(L>1|U<1)

## -----------------------------------------------------------------------------
set.seed(12345)
M<-1000; Bonf<-matrix(0,M,3);BH<-matrix(0,M,3);
for (i in 1:M)
{
  p_value<-rep(0,1000)
  p_value[1:950]<-runif(950)
  p_value[951:1000]<-rbeta(50,0.1,1)
  p_bonf<-p.adjust(p_value,method="bonferroni")
  p_BH<-p.adjust(p_value,method="BH")
  bonf<-rep(0,1000); ##拒绝第i个假设，记为1
  bh<-rep(0,1000); ##拒绝第i个假设，记为1
  for (j in 1:1000)
  {
    if (p_bonf[j]<0.1) 
    {
      bonf[j]<-1
    }
    if (p_BH[j]<0.1) 
    {
      bh[j]<-1
    }
  }
  FWER1<-(sum(bonf[1:950])>0)
  FWER2<-(sum(bh[1:950])>0)
  FDR1<-sum(bonf[1:950])/sum(bonf)
  FDR2<-sum(bh[1:950])/sum(bh)
  TPR1<-sum(bonf[951:1000])/50
  TPR2<-sum(bh[951:1000])/50
  Bonf[i, ]=c(FWER1,FDR1,TPR1)
  BH[i, ]=c(FWER2,FDR2,TPR2)       
}

## -----------------------------------------------------------------------------
output<-matrix(c(mean(Bonf[,1]),mean(Bonf[,2]),mean(Bonf[,3]),mean(BH[,1]),mean(BH[,2]),mean(BH[,3])),nrow=2,ncol=3,byrow=TRUE,dimnames = list(c("Bonferroni","BH"), c("FWER","FDR","TPR")))
output

## -----------------------------------------------------------------------------
B<-1000;m<-1000;n=5;
Bias<-numeric(m)
SE<-numeric(m)
for (i in 1:m)
{
  x<-rexp(n,2)
  lambda<-1/mean(x)
  lambdastar<-numeric(B)
  for (b in 1:B)
  {
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }
  Bias[i]<-mean(lambdastar)-lambda
  SE[i]<-sd(lambdastar)
}
d1<-c(mean(Bias),2/(n-1),mean(SE),2*n/(n-1)/sqrt(n-2))

## -----------------------------------------------------------------------------
B<-1000;m<-1000;n=10;
Bias<-numeric(m)
SE<-numeric(m)
for (i in 1:m)
{
  x<-rexp(n,2)
  lambda<-1/mean(x)
  lambdastar<-numeric(B)
  for (b in 1:B)
  {
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }
  Bias[i]<-mean(lambdastar)-lambda
  SE[i]<-sd(lambdastar)
}
d2<-c(mean(Bias),2/(n-1),mean(SE),2*n/(n-1)/sqrt(n-2))

## -----------------------------------------------------------------------------
B<-1000;m<-1000;n=20;
Bias<-numeric(m)
SE<-numeric(m)
for (i in 1:m)
{
  x<-rexp(n,2)
  lambda<-1/mean(x)
  lambdastar<-numeric(B)
  for (b in 1:B)
  {
    xstar<-sample(x,replace=TRUE)
    lambdastar[b]<-1/mean(xstar)
  }
  Bias[i]<-mean(lambdastar)-lambda
  SE[i]<-sd(lambdastar)
}
d3<-c(mean(Bias),2/(n-1),mean(SE),2*n/(n-1)/sqrt(n-2))

## -----------------------------------------------------------------------------
output<-matrix(c(d1,d2,d3),nrow=3,ncol=4,byrow=TRUE,dimnames = list(c("n=5","n=10","n=20"), c("Bias(Bootstrap)","Bias(Theoretical)","SE(Bootstrap)","SE(Theoretical)")))
output

## -----------------------------------------------------------------------------
set.seed(12345)
library(bootstrap)
x<-law$LSAT
y<-law$GPA
n<-length(x)
correlation<-cor(x,y)
R<-numeric(1000)
Tobs<-numeric(1000)
SED<-numeric(1000)
for (b in 1:1000)
{
  index<-sample(1:n,replace=TRUE)
  xstar<-x[index]
  ystar<-y[index]
  R[b]<-cor(xstar,ystar)
  RR<-numeric(1000)
  for (i in 1:1000)
  {
    index<-sample(1:n,replace=TRUE)
    xre<-xstar[index]
    yre<-ystar[index]
    RR[i]<-cor(xre,yre)
  }
  SED[b]<-sd(RR)
}
se<-sd(R)
Tobs=(R-correlation)/SED
t1<-quantile(Tobs,0.025)
t2<-quantile(Tobs,0.975)
L<-correlation-t2*se
U<-correlation-t1*se

## -----------------------------------------------------------------------------
L
U

## -----------------------------------------------------------------------------
library(boot)
set.seed(12345)
d<-c(3,5,7,18,43,85,91,98,100,130,230,487)
m<-1e3
boot.mean <- function(x,i) return(mean(x[i]))
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m)
{
    de <- boot(data=d,statistic=boot.mean, R = 999)
    ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
    ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
    ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5]
}


## -----------------------------------------------------------------------------
output<-matrix(c(mean(ci.norm[,1]),mean(ci.norm[,2]),mean(ci.basic[,1]),mean(ci.basic[,2]),mean(ci.perc[,1]),mean(ci.perc[,2]),mean(ci.bca[,1]),mean(ci.bca[,2])),nrow=4,ncol=2,byrow=TRUE,dimnames = list(c("norm","basic","perc","BCa"), c("Lower bound","Upper bound")))
output

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
SIGMA<-cov(scor)
lambda<-eigen(SIGMA)$values
theta<-max(lambda)/sum(lambda)

## -----------------------------------------------------------------------------
set.seed(12345)
n<-nrow(scor)
es<-as.numeric(n)
for (i in 1:n)
{
  scor_j<-scor[-i,]
  sigma_j<-cov(scor_j)
  lambda_j<-eigen(sigma_j)$values
  es[i]<-max(lambda_j)/sum(lambda_j)
}
je<-mean(es)
bias<-(n-1)*(je-theta)
sd<-sqrt((n-1)*var(es))
output<-matrix(c(bias,sd),nrow=2,ncol=1,byrow=TRUE,dimnames = list(c("bias","sd")))
output

## -----------------------------------------------------------------------------
library(DAAG)
set.seed(12345)
attach(ironslag)
n<-length(magnetic) #in DAAG ironslag
M<-n*(n-1)/2
e1 <- e2 <- e3 <- e4 <- numeric(M)
# for n-fold cross validation
# fit models on leave-one-out samples
for (i in 1:(n-1)) 
{
  for (j in (i+1):n)
  {
    index<-c(i,j)
    y <- magnetic[-index]
    x <- chemical[-index]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[index]
    e1[(2*n-i)*(i-1)/2+j-i]<- sum((magnetic[index] - yhat1)^2)
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[index] +J2$coef[3] * chemical[index]^2
    e2[(2*n-i)*(i-1)/2+j-i]<- sum((magnetic[index] - yhat2)^2)
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[index]
    yhat3 <- exp(logyhat3)
    e3[(2*n-i)*(i-1)/2+j-i]<- sum((magnetic[index] - yhat3)^2)
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[index])
    yhat4 <- exp(logyhat4)
    e4[(2*n-i)*(i-1)/2+j-i]<- sum((magnetic[index] - yhat4)^2)
  }
}

## -----------------------------------------------------------------------------
c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
attach(chickwts)
x1 <- sort(as.vector(weight[feed == "soybean"]))
x3 <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
f<-function(x,xs)
{
  num<-length(xs)
  numx<-numeric(num)
  for (i in 1:(num))
  {
     numx[i]<-sum(x<=xs[i])/length(x)
  }
  return(numx)
}
g<-function(y,ys)
{
  n0<-length(y)
  n1<-length(ys)
  numy<-numeric(n1)
  for (i in 1:(n1))
  {
    numy[i]<-sum(y<=ys[i])/n0
  }
  return(numy)
}
W<-function(x,y)
{
  z<-c(x,y)
  m<-length(x)
  n<-length(y)
  N<-m+n
  W=numeric(10000)
  for (i in 1:10000)
  {
  index<-sample(1:N)
  Z<-z[index]
  xp<-Z[1:m]
  yp<-Z[(m+1):N]
  W[i]<-m*n*(sum((f(x,xp)-g(y,xp))^2)+sum((f(x,yp)-g(y,yp))^2))/((m+n)^2)
  }
  wobs<-m*n*(sum((f(x,x)-g(y,x))^2)+sum((f(x,y)-g(y,y))^2))/((m+n)^2)
  p_value<-sum(W>wobs)/10000
  return(p_value)
}

## -----------------------------------------------------------------------------
W(x1,x3)

## -----------------------------------------------------------------------------
attach(chickwts)
x2 <- sort(as.vector(weight[feed == "sunflower"]))
W(x2,x3)

## -----------------------------------------------------------------------------
count5 <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(max(c(outx, outy)))
}
countp<-function(x,y)
{
  z<-c(x,y)
  n<-length(x)
  N<-length(z)
  stats<-replicate(199,expr={
    index<-sample(1:N)
    k1<-index[1:n]
    k2<-index[(n+1):N]
    count5(z[k1],z[k2])
  })
  stat<-count5(x,y)
  as<-c(stats,stat)
  return(list(estimate=stat,p=mean(as>=stat)))
}

## -----------------------------------------------------------------------------
set.seed(2023)
n1 <- 20
n2 <- 40
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x<-rnorm(n1,mu1,sigma1)
y<-rnorm(n2,mu2,sigma2)
countp(x,y)

## -----------------------------------------------------------------------------
set.seed(2023)
sigma1<-1
sigma2<-2
x<-rnorm(n1,mu1,sigma1)
y<-rnorm(n2,mu2,sigma2)
countp(x,y)

## -----------------------------------------------------------------------------
alpha<-function(N,b1,b2,b3,f0)
{
  x1<-rpois(N,1)
  x2<-rexp(N,1)
  x3<-rbinom(N,1,0.5)
  g<-function(alpha)
  {
      tmp<-exp(-alpha-b1*x1-b2*x2-b3*x3)
      p<-1/(1+tmp)
      mean(p)-f0
  }
  solution<-uniroot(g,c(-20,0))
  return(solution$root)
}

## -----------------------------------------------------------------------------
set.seed(12345)
alpha1<-alpha(1e6,0,1,-1,0.1)
alpha2<-alpha(1e6,0,1,-1,0.01)
alpha3<-alpha(1e6,0,1,-1,0.001)
alpha4<-alpha(1e6,0,1,-1,0.0001)
c(alpha1,alpha2,alpha3,alpha4)

## -----------------------------------------------------------------------------
f0<-seq(0.0001,0.5,0.01)
n<-length(f0)
y<-numeric(n)
for (i in 1:n)
{
    y[i]<-alpha(1e6,0,1,-1,f0[i])
}
plot(-log(f0),y)

## -----------------------------------------------------------------------------
M<-function(N,x0,sigma)
{
  x<-numeric(N)
  x[1]=x0
  u<-runif(N)
  r=0
  for (i in 2:N)
  {
    y<-rnorm(1,x[i-1],sigma)
    alpha<-exp(abs(x[i-1])-abs(y))
    if(u[i]<=alpha)
      x[i]<-y
    else
    {
      x[i]=x[i-1]
      r=r+1
    }
  }
  return(list(x=x,r=r))
}

## -----------------------------------------------------------------------------
N<-8000
sigma<-c(0.05,0.5,1,2,5,10)
x0<-rnorm(1)
M1<-M(N,x0,sigma[1])
M2<-M(N,x0,sigma[2])
M3<-M(N,x0,sigma[3])
M4<-M(N,x0,sigma[4])
M5<-M(N,x0,sigma[5])
M6<-M(N,x0,sigma[6])

## -----------------------------------------------------------------------------
print(c(M1$r,M2$r,M3$r,M4$r,M5$r,M6$r)/N)

## -----------------------------------------------------------------------------
B<-1000
MN1<-M1$x[(B+1):N]
MN2<-M2$x[(B+1):N]
MN3<-M3$x[(B+1):N]
MN4<-M4$x[(B+1):N]
MN5<-M5$x[(B+1):N]
MN6<-M6$x[(B+1):N]
par(mfrow=c(2,3))
plot(MN1,type='o')
plot(MN2,type='o')
plot(MN3,type='o')
plot(MN4,type='o')
plot(MN5,type='o')
plot(MN6,type='o')

## -----------------------------------------------------------------------------
N<-8000
B<-1000
X<-matrix(0,N,2)
x1<-x2<-0
X[1,]<-c(x1,x2)
for (i in 2:N)
{
  x2<-X[i-1,2]
  x1<-rnorm(1,0.9*x2,0.19)
  x2<-rnorm(1,0.9*x1,0.19)
  X[i,]<-c(x1,x2)
}

## -----------------------------------------------------------------------------
XN<-X[(B+1):N,]
x<-XN[,1]
y<-XN[,2]
L<-lm(y~x)

## -----------------------------------------------------------------------------
summary(L)

## -----------------------------------------------------------------------------
plot(x,type='l',col='blue',ylim=c(-2,2))
par(new=T)   
plot(y,type='l',col='red',ylab='random numbers',ylim=c(-2,2))
axis(side=4)
legend('topright',legend=c('x','y'),lty=c(1,2),bty='n')

## -----------------------------------------------------------------------------
library(coda)
f <- function(x, sigma)
{
  if (any(x < 0)) 
    return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

## -----------------------------------------------------------------------------
chain<-function(sigma,m,x0)
{
  x <- numeric(m)
  x[1] <- x0
  k <- 0
  u <- runif(m)
  for (i in 2:m) 
  {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) 
      x[i] <- y 
    else 
    {
      x[i] <- xt
      k <- k+1 #y is rejected
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) 
{
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
sigma<-4
x0<-c(1/sigma^2,1/sigma,sigma,sigma^2)
k<-4
m<-2000
X<-matrix(0,k,m)
for (i in 1:k)
{
  X[i,]<-chain(sigma,m,x0[i])
}
psi<-t(apply(X,1,cumsum))
for (i in 1:nrow(psi))
  psi[i,]<-psi[i,]/(1:ncol(psi))
rhat<-Gelman.Rubin(psi)
rhat

## -----------------------------------------------------------------------------
X1<-as.mcmc(X[1,])
X2<-as.mcmc(X[2,])
X3<-as.mcmc(X[3,])
X4<-as.mcmc(X[4,])
Y<-mcmc.list(X1,X2,X3,X4)
print(gelman.diag(Y))

## -----------------------------------------------------------------------------
gelman.plot(Y)

## -----------------------------------------------------------------------------
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)

## -----------------------------------------------------------------------------
mle<-function(lambda)
{
  n<-length(v)
  b=0
  for (i in 1:n)
  {
    b=b+(v[i]-u[i])/(exp(lambda*(v[i]-u[i]))-1)
  }
  b=b-sum(u)
  return(b)
}

## -----------------------------------------------------------------------------
uniroot(mle,c(0,1))

## -----------------------------------------------------------------------------
EM<-function(x)
{
  n<-length(v)
  t<-0
  for (i in 1:n)
  {
    t=t+((u[i]+1/x)*exp(-x*u[i])-(v[i]+1/x)*exp(-x*v[i]))/(exp(-x*u[i])-exp(-x*v[i]))
  }
  return(n/t)
}

## -----------------------------------------------------------------------------
x<-0.1
y<-EM(x)
t<-(y-x)
n=0
while(abs(t)>1e-5)
{
  x<-y
  y<-EM(x)
  t<-(y-x)
  n=n+1
}
c(y,n)

## -----------------------------------------------------------------------------
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

## -----------------------------------------------------------------------------
library(boot)
B<-A+2
sb<-solve.game(B)
sa<-solve.game(A)
sa$v
sb$v

## -----------------------------------------------------------------------------
round(cbind(sb$x,sb$y),7)

## -----------------------------------------------------------------------------
round(cbind(sa$x,sa$y),7)

## -----------------------------------------------------------------------------
round(sb$x*61,7)

## -----------------------------------------------------------------------------
round(sa$x*61,7)

## -----------------------------------------------------------------------------
d1<-data.frame(a=integer(),b=logical())
dim(d1)

## -----------------------------------------------------------------------------
d2<-data.frame(row.names = 1:3)
dim(d2)

## -----------------------------------------------------------------------------
d3<-data.frame()
dim(d3)

## -----------------------------------------------------------------------------
GR<-function(a,b,n,N)
{
  x<-y<-rep(0,N)
  x[1]<-rbinom(1,n,0.5)
  y[1]<-rbeta(1,x[1]+a,n-x[1]+b)
  for (i in 2:N)
  {
    x[i]<-rbinom(1,n,y[i-1])
    y[i]<-rbeta(1,x[i]+a,n-x[i]+b)
  }
  m<-matrix(c(x,y),nrow=2,byrow=TRUE)
  return(m)
}

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// Function to perform the Gibbs sampling
// [[Rcpp::export]]
NumericMatrix GC(double a, double b, int n, int N) {
  
  // Initialize vectors to store samples
  NumericVector x(N);
  NumericVector y(N);
  
  // Initialize random number generator seed
  srand(time(NULL));
  
  // Initialize x with initial value
  x[0] = R::rbinom(n,0.5);
  y[0] = R::rbeta(x[0]+a,n-x[0]+b);
  
  // Perform Gibbs sampling
  for (int i = 1; i < N; i++) {
    
    x[i] = R::rbinom(n,y[i-1]);
    y[i] = R::rbeta(x[i]+a,n-x[i]+b);}
    
  // Return stored samples as matrix
  NumericMatrix samples(N, 2);
  samples(_, 0) = x;
  samples(_, 1) = y;
  
  return samples;
}')

## -----------------------------------------------------------------------------
library(microbenchmark)
set.seed(12345)
a=1
b=2
n=10
N=10000
ts<-microbenchmark(gr=GR(a,b,n,N),gc=GC(a,b,n,N))

## -----------------------------------------------------------------------------
summary(ts)

