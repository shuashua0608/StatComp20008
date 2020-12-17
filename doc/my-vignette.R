## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ridge)
data(Hald)
cement<-data.frame(Hald)
#use Hald's cement data to carry out linear regression
lm.sol<-lm(y~.,data=cement)
#get fitted value
yfit<-predict(lm.sol)
yreal<-cement[,1]


## -----------------------------------------------------------------------------
plot(x=yfit, y=yreal, xlab="Model-predicted y value", ylab="Actual y value", xlim=c(60,120),ylim=c(60,120))
abline(a=0,b=1,col="grey")

## -----------------------------------------------------------------------------
e=residuals(lm.sol)
e.std=rstandard(lm.sol)
plot(e~yfit)
plot(e.std~yfit)


## -----------------------------------------------------------------------------
 xtable::xtable(head(cement))

## -----------------------------------------------------------------------------
n<-100 
u<-runif(n)
a<-2
b<-2
x<-b/(u^(1/a))## the inverse function of F(x)
hist(x,prob=T)
y<-seq(2,20,0.01)
lines(y,a*b^a/y^(a+1))

## -----------------------------------------------------------------------------
n<-100
y<-numeric(n)
k<-0 #counter for accepted
j<-0#iterations
while(k<n){
  x<-runif(1)#random variable from U(0,1)
  u1<-runif(1,-1,1)
  u2<-runif(1,-1,1)
  u3<-runif(1,-1,1)# generate U1,U2,U3 from U(-1,1)
  j<-j+1
  if(abs(u1)<=abs(u3)&&abs(u2)<=abs(u3))##deliver u2
    {
    if(abs(u2)<=x){
      k<-k+1
      y[k]<-u2
    }
  }
  else  ##deliver u3
    {   
    if(abs(u3)<=x){
      k<-k+1
      y[k]<-u3
    }
  }
}
hist(y,prob=T)## y:the simulated random sample with length 10000.

## -----------------------------------------------------------------------------
n<-100
u<-runif(n)
beta=2
r=4
x<-beta/(1-u)^(1/r)-beta
hist(x,prob=T,main=expression(f(x)==r*beta^r/(beta+x)^(r+1)))
y<-seq(0,20,0.01)
lines(y,r*beta^r/(beta+y)^(r+1))

## -----------------------------------------------------------------------------
set.seed(12345)
m <- 1e2; 
t <- runif(m, min=0, max=pi/3)##Generate random numbers~U(0,pi/3)
hat.f <- mean(sin(t)) * pi/3 ##estimate of f(t)
print(c(hat.f,0.5))

## -----------------------------------------------------------------------------
m<-1e2
set.seed(12345)
#####The antithetic variate approach ########

t1<-runif(m/2) ##Generate random numbers~U(0,1)
theta1=1/2*((exp(t1))+(exp(1-t1)))##the antithetic estimator

#####The simple Monte Carlo method ##########

t2<-runif(m)
theta2<-exp(t2)##simple MC estimator

##compare the results of antithetic variate approach, simple MC method and true value

print(c(mean(theta1),mean(theta2),exp(1)-1))

####variance reduction######
(var(theta2)-var(theta1))/var(theta2)*100  

## -----------------------------------------------------------------------------
antithetic<-function(m=1e2,anti=T){
        u<-runif(m/2)
        if(!anti) v<-runif(m/2) else ###simple MC method
                v<-1-u###antithetic variate approach
        t<-c(u,v)
        theta<-exp(t)
        mean(theta)##the estimator
}
n<-100
MC1<-MC2<-numeric(n)
for(i in 1:n){
        MC1[i]<-antithetic(m=1e2)
        MC2[i]<-antithetic(m=1e2,anti=F)
}
print(100*(var(MC2)-var(MC1))/var(MC2))###variance reduction


## -----------------------------------------------------------------------------
m<-1e2
g<-function(x){
  exp(-x^2/2)*x^2/sqrt(2*pi)*(x>1)
  }
  

#####using f1######
x<-rweibull(m,2,sqrt(2))
gf<-g(x)/dweibull(x,2,sqrt(2))
theta1<-mean(gf)
se1<-sd(gf)



####using f2######
x<-rgamma(m,3,2)
gf<-g(x)/dgamma(x,3,2)
theta2<-mean(gf)
se2<-sd(gf)


####PS:Using Normal ditribution#######
x<-rnorm(m)
gf<-g(x)/dnorm(x)
theta3<-mean(gf)
se3<-sd(gf)

####In one simulation
####theta.hat<-(0.3977415 0.4062168 0.383093)
####se<-(0.3589789 0.3008697 1.096234)
theta.hat<-c(theta1,theta2,theta3)

se<-c(se1,se2,se3)
rbind(theta.hat,se)


## -----------------------------------------------------------------------------
t<-seq(1,10,0.1)
g<-exp(-t^2/2)*t^2/sqrt(2*pi)
f1<-dweibull(t,2,sqrt(2))
f2<-dgamma(t,3,2)
 
####get the curve of g,f1 and f2####
plot(t,g,type="l",col="black",main="compare g(x), f1(x) and f2(x) ")   
lines(t,f1,col="red")  
lines(t,f2,col="green")  
legend("topright",legend =c('g(t)','f1(t)',"f2(t)") ,lty=1,col=c("black","red","green")) 


## -----------------------------------------------------------------------------
t<-seq(1,10,0.1)
 g<-exp(-t^2/2)*t^2/sqrt(2*pi)
 f1<-dweibull(t,2,sqrt(2))
 f2<-dgamma(t,3,2)
r1<-g/f1
r2<-g/f2
plot(t,r1,col="red")
points(t,r2,col="green")
title(main="ratio function")

## -----------------------------------------------------------------------------

###example 5.10##
m<-1e2
g<-function(x)
  exp(-x-log(1+x^2))*(x > 0)*(x < 1)
f3<-function(x)
  exp(-x)/(1-exp(-1))

##### using f3 ####
t <- runif(m) #inverse transform method for f3
x <- -log(1-t*(1 - exp(-1)))
gf <- g(x) / ((exp(-x) / (1-exp(-1))))
theta1<- mean(gf)
se1<- sd(gf)


###Stratified importance sampling #####
k<-5     ##number of stratum
r<-m/k   ##replicates per stratum
T2<-numeric(k)
var<-numeric(k)
for (i in 1:k) {
  u<-  runif(r,(i-1)/k,i/k)
  x<- -log(1-u*(1-exp(-1)))
#### estimated mean in each subinterval
  T2[i]<-mean(g(x)/f3(x)) 
  
####estimated variance in each subinterval
  var[i]<-sd(g(x)/f3(x))^2 
}
theta2<-mean(T2) 
se2<-sqrt(sum(var)/r) 

###compare###
print(c(theta1,theta2))
print(c(se1,se2))

## -----------------------------------------------------------------------------
mu<-1;
sigma<-1
n<-20
alpha<-0.05
UCL<-replicate(10,expr = {
  x<-rlnorm(n,mu,sigma)
  y<-log(x)
  abs(sqrt(n/var(y))*(mean(y)-mu))
})
sum(UCL<qt(1-alpha/2,n-1))
mean(UCL<qt(1-alpha/2,n-1)) 

## -----------------------------------------------------------------------------
####Example 6.4 ： N(0,sigma^2=4)#####
n<-20
alpha<-0.05

####example 6.4
UCL.norm<-replicate(10,expr={
  x<-rnorm(n,mean=0,sd=2)
  ((mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n)>0)&&(mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n)<0))####the mean of x=0
})
sum(UCL.norm==TRUE)###check if the CI covers 0

####Chi分布
###Using t-interval
UCL.chisq<-replicate(10,expr={
  y<-rchisq(n,df=2)
  ((mean(y)+qt(1-alpha/2,n-1)*sd(y)/sqrt(n)>2)&&(mean(y)-qt(1-alpha/2,n-1)*sd(y)/sqrt(n)<2))#### the mean of x=2
})
sum(UCL.chisq==TRUE) ###check if the CI covers 2


###Using example 6.4
UCL.chi<-replicate(10,expr={
  y<-rchisq(n,df=2)
  (n-1)*var(y)/qchisq(alpha,df=n-1)
})
sum(UCL.chi>4)  ###check if CI covers 4(sigma^2)

## -----------------------------------------------------------------------------

n<-c(10,20,30,50,100,500)##sample size
cv<-qnorm(0.975,0,sqrt(6*(n-2)/((n+1)*(n+3))))

skewness<-function(x){
  barx<-mean(x)
  s2<-mean((x-barx)^2)
  s3<-mean((x-barx)^3)
  return(s3/s2^1.5)
}

r1<-numeric(length(n))
##store results for beta distribution
r2<-numeric(length(n))
##store results for t distribution

B<-1e2 ##repeat each simulation
a<-function(a){
for(i in 1:length(n)){
  test<-numeric(B)
  compare<-numeric(B)
  for(j in 1:B){
    t<-seq(0,1,length=n[i])
    ## generate samples from beta distribution##
    x<-rbeta(t,a,a)
    ##obtain samples from t distribution
    y<-rt(n[i],n[i]-1)
    test[j]<-as.integer(abs(skewness(x))>=cv[i])
    compare[j]<-as.integer(abs(skewness(y))>=cv[i])
  }
  r1[i]<-mean(test)
  r2[i]<-mean(compare)
}
  r<-cbind(r1,r2)
  colnames(r)<-c("beta","t")
  rownames(r)<-c("10","20","30","50","100","500")
  return(r)
}
a(0.05)
a(0.5)
a(1)




## -----------------------------------------------------------------------------
### Example 6.16 ###

sigma1<-1
sigma2<-1.5
u1<-u2<-0

##Count Five test##
Countfive<-function(x,y){
  sx<-x-mean(x)
  sy<-y-mean(y)
  out1<-sum(sx>max(sy))+sum(sx<min(sy))
  out2<-sum(sy>max(sx))+sum(sy<min(sx))
  ##return 1##
  return(as.integer(max(c(out1,out2))>5))
}

m<-1e2
power1<-function(n1,n2){
  ### Countfive test ##
  r<-replicate(m,expr={
  x<-rnorm(n1,u1,sd=sigma1)
  y<-rnorm(n2,u2,sd=sigma2)
  x<-x-mean(x)
  y<-y-mean(y)
  Countfive(x,y)
  })
  return(mean(r))
}

power2<-function(n1,n2){
  ###F test ##
  r<-replicate(m,expr={
    x<-rnorm(n1,u1,sd=sigma1)
    y<-rnorm(n2,u2,sd=sigma2)
    ftest<-var.test(x,y)$p.value
    return(ftest)
    })
  sum(r>0.055)/m
}
data1<-c(power1(20,20),power1(20,30),power1(50,50),power1(100,100),power1(200,200),power1(500,500))
data2<-c(power2(20,20),power2(20,30),power2(50,50),power2(100,100),power2(200,200),power2(500,500))
d<-matrix(cbind(data1,data2),nrow=2,byrow=T)
rownames(d)<-c("Countfive test","F test")
colnames(d)<-c("(20,20)","(20,30)","(50,50)","(100,100)","(200,200)","(500,500)")
d


## -----------------------------------------------------------------------------
###Example 6.8###
library(MASS)
n <- c(10,20,50,100) #sample size

mardia <- function(x) {
  mx<- mean(x)
  n <- nrow(x)
  m <- 1/n^2*sum(((x-mx)%*%solve(cov(x))%*%t(x-mx))^3)
  return(m)
}

reject <- numeric(length(n))
m <- 1e2
s <- matrix(c(1,0,0,1),2,2)

r<-function(d){
for (i in 1:length(n)) {
  test <- numeric(m)  #test decisions
  for (j in 1:m) {
    x <- mvrnorm(n[i],c(0,0),s)##generate x 
    cv<-qchisq(0.975,d*(d+1)*(d+2)/6)#crit.value
    test[j] <- as.integer(abs(mardia(x)) >= 6*cv/n[i])
  }
  reject[i] <- mean(test)
}
  return(reject)
}
##set d=2 and d=3##
##The results##
r(2)
r(3)


## ----warning=F----------------------------------------------------------------
###In Example 7.2, dataset: law ###
data(law,package="bootstrap")
n<-nrow(law)
l<-law$LSAT
gpa<-law$GPA

cor.jack<-numeric(n)
for( i in 1:n){
        cor.jack[i]<-cor(l[-i],gpa[-i])
}
###estimate bias###
bias<-(n-1)*(mean(cor.jack)-cor(l,gpa))
###estimate standard error###
se<-sqrt((n-1)*mean((cor.jack-mean(cor.jack))^2))
print(cbind(rawdata=cor(l,gpa),bias=bias,se=se))



## ----warning=F, eval=FALSE----------------------------------------------------
#  library(boot)
#  data(aircondit,package="boot")
#  
#  hour<-aircondit$hours
#  
#  ###use meantime function to obtain the estimate for mean and sd
#  meantime<-function(d,index){
#          m<-mean(d[index])
#          n<-length(index)
#          s<-(n-1)*var(d[index])/n^2
#          c(m,s)
#  }
#  
#  boot.air<-boot(hour,meantime,R=5)
#  print(boot.air)
#  ###using boot.ci function to obtain CI by different method###
#  boot.ci(boot.air,conf=0.95,type=c("norm","basic","perc","bca"))
#  

## -----------------------------------------------------------------------------
data(scor,package="bootstrap")
n<-nrow(scor)
lambda<-eigen(cov(scor))$values
mec<-scor[,1]
vec<-scor[,2]
alg<-scor[,3]
ana<-scor[,4]
sta<-scor[,5]
###the estimate###
theta<-lambda[1]/sum(lambda)

lambda.jack<-numeric(n)
for( i in 1:n){
        scor.hat<-cbind(mec[-i],vec[-i],alg[-i],ana[-i],sta[-i])
        lambda.hat<-eigen(cov(scor.hat))$values
        lambda.jack[i]<-lambda.hat[1]/sum(lambda.hat)
}
###estimate bias###
bias<-(n-1)*(mean(lambda.jack)-theta)
###estimate standard error###
se<-sqrt((n-1)*mean((lambda.jack-mean(lambda.jack))^2))
print(cbind(rawdata=theta,bias=bias,se=se))

## -----------------------------------------------------------------------------
###the 4 models are from example 7.17###
library(DAAG)
library(boot)
attach(ironslag)
n<-length(magnetic)

t<- seq(10, 40, .1)
###model 1: linear
mod1 <- lm(magnetic ~ chemical) 

###model 2:quadratic 
mod2 <- lm(magnetic ~ chemical + I(chemical^2))


###model 3: exponential
mod3 <- lm(log(magnetic) ~ chemical) 


###model 4:log-log
mod4 <- lm(log(magnetic) ~ log(chemical)) 

####leave-two-out###
###assign n*n matrix for each time 's error###
###each time randomly pick i and j### 
###take the average of the two prediction errors###
e1<-e2<-e3<-e4<-matrix(nrow=n,ncol=n)
I<-rep(1,n^2)
for(i in 1:n){
        y<-magnetic[-i]
        x<-chemical[-i]
        for(j in 1:n){
                y<-y[-j]
                x<-x[-j]
               ###linear model### 
                J1<-lm(y~x)
                yhat11<-J1$coef[1]+J1$coef[2]*chemical[i]
                yhat12<-J1$coef[1]+J1$coef[2]*chemical[j]
                e1[i,j]<-0.5*((magnetic[i]-yhat11)^2+(magnetic[j]-yhat12)^2)
               ###quadratic model###
                J2 <- lm(y ~ x + I(x^2))
                yhat21<-J2$coef[1]+J2$coef[2]*chemical[i]+J2$coef[3]*chemical[i]^2
                 yhat22<-J2$coef[1]+J2$coef[2]*chemical[j]+J2$coef[3]*chemical[j]^2
                 e2[i,j]<-0.5*((magnetic[i]-yhat21)^2+(magnetic[j]-yhat22)^2)
              ###exponential model###
                 J3 <- lm(log(y) ~ x)
                 yhat31 <- exp(J3$coef[1] + J3$coef[2] * chemical[i])
                 yhat32 <- exp(J3$coef[1] + J3$coef[2] * chemical[j])
                 e3[i,j]<-0.5*((magnetic[i]-yhat31)^2+(magnetic[j]-yhat32)^2)
              ###log-log model###
                 J4 <- lm(log(y) ~ log(x))
                 yhat41 <- exp(J4$coef[1] + J4$coef[2] * log(chemical[i]))
                 yhat42 <- exp(J4$coef[1] + J4$coef[2] * log(chemical[j]))
                 e4[i,j]<-0.5*((magnetic[i]-yhat41)^2+(magnetic[j]-yhat42)^2)
 
        }
}
###estimate for leave-two-out ###
E<-function(matr){
        sum(colSums(matr))/(n*n)
}

leaveoneout<-c(19.55644, 17.85248 ,18.44188 ,20.45424)
leavetwoout<-c(E(e1),E(e2),E(e3),E(e4))
cbind(leaveoneout,leavetwoout)


## -----------------------------------------------------------------------------
set.seed(12345)

# Count Five test
countfive<-function(x, y) {
 ##centralize x and y##
sx<-x - mean(x)
sy<-y - mean(y)
out1<-sum(sx > max(sy)) + sum(sx< min(sy))
out2 <- sum(sy> max(sx)) + sum(sy < min(sx))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(out1, out2)) > 5))
}
# Count Five test permutation
count_permutation<-function(z) {
n <-length(z)
k<-1:(n/2)
x <-z[k]
y <-z[-k]
sx <-x - mean(x)
sy<-y - mean(y)
out1<-sum(sx > max(sy)) + sum(sx< min(sy))
out2 <- sum(sy> max(sx)) + sum(sy < min(sx))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(out1, out2)) > 5))
}

per = function(z,m) {
  n <- length(z)
  r<- numeric(m)
  for (i in 1: m){
      p = sample(1:n ,n ,replace = FALSE)
      r[i] = count_permutation(z[p])
  }
  return(mean(r))
}              


u1<-u2<-0
s1<-s2<-1
m<-1e2

compare<-function(n1,n2){
p1 = mean(replicate(m, expr={
x = rnorm(n1, u1, s1)
y = rnorm(n2, u2, s2)
x = x - mean(x) #centered by sample mean
y = y - mean(y)
countfive(x, y)
}))
p2 = mean(replicate(m, expr={
x = rnorm(n1, u1, s1)
y = rnorm(n2, u2, s2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
per(z,m) 
})<0.05)
p<-c(p1,p2)
return(p)
}
p<-cbind(compare(20,20),compare(20,30),compare(20,50))
rownames(p)<-c("Countfive","permutation")
colnames(p)<-c("(20,20)","(20,30)","(20,50)")
p



## ----eval=FALSE---------------------------------------------------------------
#  ###First load packages we need###
#  ####For NN test###
#  library(RANN)
#  ####For energy test###
#  library(energy)
#  library(boot)
#  ###For ball test###
#  library(Ball)

## ----eval=FALSE---------------------------------------------------------------
#  ###From the slides: Tn function###
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]
#    n2 <- sizes[2]
#    n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ]
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5)
#    i2 <- sum(block2 > n1+.5)
#    return((i1 + i2) / (k * n))
#  }
#  
#  set.seed(12345)
#  m<-100
#  k<-3
#  n1<-n2<-50
#  n <- n1+n2
#  N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#     boot.obj <- boot(data=z,statistic=Tn,R=999,
#     sim = "permutation", sizes = sizes,k=k)
#     ts <- c(boot.obj$t0,boot.obj$t)
#     p.value <- mean(ts>=ts[1])
#     list(statistic=ts[1],p.value=p.value)
#  }
#   p.values <- matrix(NA,m,3)
#   power.compare<-function(mu1,mu2,s1,s2,a){
#     x <- rnorm(n1,mu1,s1)
#     y <- rnorm(n2,mu2,s2)
#     z <- c(x,y)
#     for(i in 1:m){
#     ### t distribution
#     p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#     p.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
#  }
#   pow<- colMeans(p.values<a)
#   names(pow)<-c("NN","energy","Ball")
#   return(pow)
#   }
#  
#  ###Unequal variances and equal expectations###
#  power.compare(0,0,1,1.5,0.055)
#  ###Unequal variances and unequal expectations###
#  power.compare(0.5,0,1,1.5,0.02)
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  m<-100
#  k<-3
#  n1<-n2<-50
#  n <- n1+n2
#  N = c(n1,n2)
#    eqdist.nn <- function(z,sizes,k){
#     boot.obj <- boot(data=z,statistic=Tn,R=999,
#     sim = "permutation", sizes = sizes,k=k)
#     ts <- c(boot.obj$t0,boot.obj$t)
#     p.value <- mean(ts>=ts[1])
#     list(statistic=ts[1],p.value=p.value)
#  }
#   pt.values <- matrix(NA,m,3)
#  ###t distribution###
#   for(i in 1:m){
#     ### t distribution
#     x <- rt(n1,df=1)
#     y <- rt(n2,df=1)
#     z <- c(x,y)
#     pt.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     pt.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#     pt.values[i,3]<-bd.test(x,y,R=999,seed=i*12345)$p.value
#  }
#   alpha.t <- 0.2
#   pow.t <- colMeans(pt.values<alpha.t)
#   names(pow.t)<-c("NN","energy","Ball")
#   pow.t

## ----eval=FALSE---------------------------------------------------------------
#  ### Bimodel distribution###
#   pb.values <- matrix(NA,m,3)
#   for(i in 1:m){
#     ### t distribution
#     x <- 0.5*rnorm(n1,0,1)+0.5*rnorm(n1,1,1)
#     y <- 0.5*rnorm(n1,0,1.5)+0.5*rnorm(n1,1,1.5)
#     z <- c(x,y)
#     pb.values[i,1] <- eqdist.nn(z,N,k)$p.value
#     pb.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#     pb.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
#  }
#   alpha.b <- 0.2
#   pow.b <- colMeans(pb.values<alpha.b)
#   names(pow.b)<-c("NN","energy","Ball")
#   pow.b

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  n1<-10
#  n2<-100
#  n <- n1+n2
#  N = c(n1,n2)
#  m<-50
#  power.compare<-function(mu1,mu2,s1,s2,a){
#     x <- rnorm(n1,mu1,s1)
#     y <- rnorm(n2,mu2,s2)
#     z <- c(x,y)
#     for(i in 1:m){
#     ### t distribution
#     pt.values[i,1] <- eqdist.nn(z,N,3)$p.value
#     pt.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
#     pt.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
#  }
#   alpha <- a
#   pow<- colMeans(pt.values<a)
#   names(pow)<-c("NN","energy","Ball")
#   return(pow)
#  }
#  power.compare(0,0,1,1.5,0.1)

## ----eval=FALSE---------------------------------------------------------------
#  rv.metro<-function(s,x0,N){
#    ##s: different sigma for different variances##
#    ##x0:initial value##
#    ##N:length of the chain##
#    x<-numeric(N)
#    u<-runif(N)
#    x[1]<-x0
#    k<-0
#    for(i in 2:N){
#      y<-rnorm(1,x[i-1],s)
#      if(u[i]<=exp(abs(x[i-1])-abs(y)))
#        x[i]<-y else{
#          x[i]<-x[i-1]
#          k<-k+1
#        }
#    }
#    return(list(x=x,k=k))
#  }
#  s<-c(0.1,0.5,1,5)
#  x0<-10
#  N<-4000
#  
#  r1<-rv.metro(s[1],x0,N)
#  r2<-rv.metro(s[2],x0,N)
#  r3<-rv.metro(s[3],x0,N)
#  r4<-rv.metro(s[4],x0,N)
#  ##display the numbers of rejected points for different variances##
#  dat<-c(r1$k,r2$k,r3$k,r4$k)
#  ##display the rejection rate##
#  dat/N
#  ##the acceptance rate##
#  (N-dat)/N

## ----eval=FALSE---------------------------------------------------------------
#  par(mfrow = c(2,2))
#  plot(r1$x,ylab='value of x',type='l',main='s=0.1')
#  plot(r2$x,ylab='value of x',type='l',main='s=0.5')
#  plot(r3$x,ylab='value of x',type='l',main='s=1')
#  plot(r4$x,ylab='value of x',type='l',main='s=5')

## ----eval=F-------------------------------------------------------------------
#  ###Gelman-Rubin function from example 9.8###
#  
#  GR<-function(p){
#    ###p[i,j] represents the statistic psi(X[i,1:j])###
#    p<-as.matrix(p)
#    n<-ncol(p)
#    k<-nrow(p)
#  
#    p.B<-rowMeans(p)
#    ###between-sequence variance###
#    B<-n*var(p.B)
#    p.w<-apply(p,1,"var")
#    ###within sample variance###
#    W<-mean(p.w)
#    ###upper bound for psi###
#    upper<-(n-1)*W/n+B/n
#    ###G-R statistic###
#    rhat<-upper/W
#    rhat
#  }
#  
#  ###M-H sampler###
#  chain<-function(s,N,X0){
#    x<-numeric(N)
#    x[1]<-X0
#    u<-runif(N)
#  
#    for(i in 2:N){
#      ##candidate points##
#      y<-rnorm(1,x[i-1],s)
#      ##target density##
#      r1<-dnorm(x[i-1],y,s)*exp(-abs(y))/2
#      r2<-dnorm(y,x[i-1],s)*exp(-abs(x[i-1]))/2
#      r<-r1/r2
#      if(u[i]<=r) x[i]<-y else
#        x[i]<-x[i-1]
#    }
#    return(x)
#  }
#  
#  s<-1##sigma chosen for proposal distribution
#  k<-4 ##4 chains
#  n<-150 ##length of chain
#  burn<-10 ##burn-in length
#  
#  a<-c(-10,-5,5,10)
#  X<-matrix(0,nrow=k,ncol=n)
#  for(i in 1:k){
#    X[i,]<-chain(s,n,a[i])
#  }
#  ##diagnostic statistics##
#  psi<-t(apply(X,1,cumsum))
#  for(i in 1:nrow(psi)){
#    psi[i, ]<-psi[i,]/(1:ncol(psi))
#  }
#  print(GR(psi))
#  
#  ##psi for four chains##
#  par(mfrow=c(2,2))
#  for( i in 1:k)
#    plot(psi[i,(burn+1):n],type='l',xlab=i,ylab=bquote(psi))
#  
#  par(mfrow=c(1,1))
#  ##R-hat##
#  rhat<-numeric(n)
#  for(i in (burn+1):n){
#    rhat[i]<-GR(psi[,1:i])
#  }
#  plot(rhat[(burn+1):n],type='l',xlab="",ylab="R")
#  abline(h=1.2,lty=2)
#  
#  

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
for( i in 1:length(k)){
  a<-seq(0,sqrt(k[i]),0.001)
  r1<-pt(a*sqrt((k[i]-1)/(k[i]-a^2)),df=k[i]-1)
  r2<-pt(a*sqrt(k[i]/(k[i]+1-a^2)),df=k[i])
  plot(a,r1-r2,type="l",col="black",xlab='x',ylab='y')
  lines(a,r2,col="red")
  abline(0,0)
}

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
##location of the root ##
root<-numeric(length(k))
## the value of the function evaluated at that point##
value<-numeric(length(k))
for( i in 1:length(k)){
  out<-uniroot(function(a){
    pt(a*sqrt((k[i]-1)/(k[i]-a^2)),df=k[i]-1)-pt(a*sqrt(k[i]/(k[i]+1-a^2)),df=k[i])},lower=1,upper=2)
  root[i]<-out$root
  value[i]<-out$f.root
}
r<-cbind(root,value)
rownames(r)<-k
r


## -----------------------------------------------------------------------------
## nA=nAA+nAO ##
nA<-444
##nB=nBB+nBO##
nB<-132
nOO<-361
nAB<-63
n<-nA+nB+nOO+nAB
##rhat##
r0<-sqrt(nOO/n)
n1<-nA/n
n2<-nB/n
##Solve equations to get p0 and q0##
##p^2+2pr=n1##
##q^2+2qr=n2##
p0<- -r0+sqrt(r0^2+n1)
q0<- -r0+sqrt(r0^2+n2)

##set max iterations to be 30##
iter<-30
p<-p0
q<-q0
dat<-matrix(data=NA,nrow=iter,ncol=3)

for(i in 1:iter){
  a<-nA*(1-p-q)/(2-p-2*q)
  b<-nB*(1-p-q)/(2-q-2*p)
  ##using matrix to solve p1 and p2##
  n<-matrix(c(2*nA+2*nOO+nAB+2*b,2*nA+nAB-2*a,2*nB+nAB-2*b,2*nB+2*nOO+nAB+2*a),nrow=2,byrow=T)
  rn<-matrix(c(2*nA+nAB-2*a,2*nB+nAB-2*b),nrow=2)

  result<-solve(n,rn)
  p<-result[1]
  q<-result[2]
  r<-1-q-p
  likelihood<-nA*log(p^2+2*p*r)+nB*log(q^2+2*q*r)+2*nOO*log(r)+nAB*log(2*p*q)
  
  dat[i,1]<-p
  dat[i,2]<-q
  dat[i,3]<-likelihood
}
dat



## -----------------------------------------------------------------------------
data(mtcars)
###lapply###
formulas <- list(
  mpg ~ disp, 
  mpg ~ I(1 / disp), 
  mpg ~ disp + wt, 
  mpg ~ I(1 / disp) + wt
)
model1<-lapply(formulas,lm,data=mtcars)
model1
###loops###
loops<-function(x,lm){
  out<-vector("list",length(x))
  for(i in seq_along(x)){
    out[[i]]<-lm(x[[i]],data=mtcars)
  }
  out
}
loops(formulas,lm)

## -----------------------------------------------------------------------------
 trials <- replicate( 
   100, 
   t.test(rpois(10, 10), rpois(7, 10)),
   simplify = FALSE
)
###sapply###
###anonymous function###
p1<-sapply(trials,function(x)x$p.value)
p1
###extra challenge###
###using "[[ "###
p2<-sapply(trials,"[[","p.value")
p2

## -----------------------------------------------------------------------------
##x:list/matrix##
##fun:function we need###
##output:type of output##
lapply1<-function(x,fun,output){
  my<-Map(fun,x)
  vapply(my,function(x)x,output)
}
lapply1(mtcars,mean,numeric(1))
lapply1(mtcars,sd,double(1))

## -----------------------------------------------------------------------------
library(Rcpp)
set.seed(12345)
rv<-function(s,x0,N){
  ##s: different sigma for different variances##
  ##x0:initial value##
  ##N:length of the chain##
  x<-numeric(N)
  u<-runif(N)
  x[1]<-x0

  for(i in 2:N){
    y<-rnorm(1,x[i-1],s)
    if(u[i]<=exp(abs(x[i-1])-abs(y)))
      x[i]<-y else{
        x[i]<-x[i-1]
      }
  }
  return(x)
}

## -----------------------------------------------------------------------------
sourceCpp(
  code= '
    #include<cmath>
    #include<Rcpp.h>
    using namespace Rcpp;
    
    double abs(double x){
      double y=0.5*exp(-fabs(x));
      return y;
    }
    // [[Rcpp::export]]
    NumericVector rvC(double s,double x0,int N) {
      NumericVector x(N);
      NumericVector u=runif(N);
      x[0]=x0;
      for(int i=1;i<N; i++){
        NumericVector y=rnorm(1,x[(i-1)],s);
        if(u[i]<=abs(y[0])/abs(x[(i-1)])){
          x[i]=y[0];
        }
        else{
          x[i]=x[(i-1)];
        }
      }
      return x;
    }
  '
)

## -----------------------------------------------------------------------------
s<-c(0.05,0.1)
x0<-20
N<-200
r<-matrix(data=NA,ncol=2,nrow=N)
rC<-matrix(data=NA,ncol=2,nrow=N)
for(i in 1:length(s)){
  r[ ,i]<-rv(s[i],x0,N)
  rC[ ,i]<-rvC(s[i],x0,N)
}
      

## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
for(i in 1:2){
  plot(r[,i],ylab='value of x',type='l',main = bquote(sigma == .(s[i])))
}


## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
for(i in 1:2){
  plot(rC[,i],ylab='value of x',type='l',main = bquote(sigma == .(s[i])))
}

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
for(i in 1:2){
  qqplot(r[,i],rC[,i],xlab="rv.metro",ylab="C.metro",
         main=bquote(sigma == .(s[i])))
  abline(0,1,col='red')
}

## -----------------------------------------------------------------------------
library(microbenchmark)
microbenchmark(t1=rv(s[1],x0,N),t2=rvC(s[1],x0,N))
microbenchmark(t1=rv(s[2],x0,N),t2=rvC(s[2],x0,N))

