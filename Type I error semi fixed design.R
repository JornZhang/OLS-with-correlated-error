rm(list = ls())

library("sandwich")
library("lmtest")

alpha=0.05

z_alpha=qnorm(alpha, mean = 0, sd=1, lower.tail = F)

n=60

d=2


beta=c(0,1)


run_boxplot=1000 # number of people
run_cp=10000 # each people run regression run_cp times to calculate one coverage probability










############ V is AR1  #########################
rho=0.9
V=matrix(0, nrow=n, ncol=n)
for (g in 1:n) {
  for (l in 1:n) {
    if(g==l){
      V[g,l]=1
    } else{
      V[g,l]=(rho)^(abs(g-l))
    }
  }
}
Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, V, Lambda.f0.5)
#######################################################

















############ V has noc equal-sized clusters  #########################
noc.true=10
#rho=rep(0.9, times=noc.true)
rho=seq(from=0, to=0.9, length.out=noc.true)
V=matrix(0, nrow=n, ncol=n)
for (g in 1:noc.true) {
  Vcluster=matrix(rho[g], (n/noc.true), (n/noc.true)) + (1-rho[g])*diag((n/noc.true))
  V[((g-1)*(n/noc.true)+1): (g*(n/noc.true)), ((g-1)*(n/noc.true)+1): (g*(n/noc.true))] = Vcluster
}

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, V, Lambda.f0.5, Vcluster)
#######################################################




















################# V(1-rho) #######################
rho=0.9
V=matrix(rho, nrow = n, ncol = n) + (1-rho)*diag(n)
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Q, V, Lambda0.5, Lambda.f0.5)
###################################################################











################### V is factor model from Cholesky #####################################
K=2 # the number of factors in the correlated error term
alpha.fac=0.9
# generate A
set.seed(420)
Atmp = matrix(rnorm(n*K,0,1),n)
Lt = chol(t(Atmp)%*%Atmp)
A=sqrt((alpha.fac*n)/K)*(Atmp %*% solve(Lt))

# check the orthogonality of A columns
t(A)%*%A

V=A%*%t(A)+(1-alpha.fac)*diag(n)
sum(diag(V)) # check tr(V)=n

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, Lambda.f0.5, Atmp, Lt, A, K, V, alpha.fac)
###############################################################################





















############# V is permutated AR1 #############################################
rho=0.9
V=matrix(0, nrow=n, ncol=n)
for (g in 1:n) {
  for (l in 1:n) {
    if(g==l){
      V[g,l]=1
    } else{
      V[g,l]=(rho)^(abs(g-l))
    }
  }
}
set.seed(419)
index.V=sample(1:n)
V=V[index.V, index.V]

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, V, Lambda.f0.5, index.V)
##############################################################













###################### V is factor model from eigen-value decomposition ##########
K=2 # the number of factors in the correlated error term
alpha.fac=0.9
# generate A
set.seed(420)
Atmp = matrix(rnorm(n*K,0,1),n)
Qtmp=eigen(t(Atmp)%*%Atmp)$vectors
A=sqrt(alpha.fac * n / sum(diag((t(Atmp)%*%Atmp))))*(Atmp%*%Qtmp)

# check the orthogonality of A columns
t(A)%*%A

V=A%*%t(A)+(1-alpha.fac)*diag(n)
sum(diag(V)) # check tr(V)=n

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, Lambda.f0.5, Atmp, Qtmp, A, K, V, alpha.fac)
######################################################################################













################## V has three different blocks ################################
rho=0.9
V1=matrix(rho, nrow = (n/4), ncol = (n/4)) + (1-rho)*diag((n/4)) # first diagonal blk is 1rho
V2=matrix(0, nrow=(n/4), ncol=(n/4))
for (g in 1:(n/4)) {
  for (l in 1:(n/4)) {
    if(g==l){
      V2[g,l]=1
    } else{
      V2[g,l]=(rho)^(abs(g-l))
    }
  }
}
set.seed(419)
index.V2=sample(1:(n/4))
V2=V2[index.V2, index.V2] # second diagonal is permuted AR1(0.9)
rho=-0.9 
V3=matrix(0, nrow=(n/2), ncol=(n/2))
for (g in 1:(n/2)) {
  for (l in 1:(n/2)) {
    if(g==l){
      V3[g,l]=1
    } else{
      V3[g,l]=(rho)^(abs(g-l))
    }
  }
} # third diagonal blk is AR1(-0.9)
V=matrix(0,n,n)
V[1:(n/4),1:(n/4)]=V1
V[( (n/4) + 1):(n/2),( (n/4) + 1):(n/2)]=V2
V[( (n/2) + 1):(n),( (n/2) + 1):(n)]=V3

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, V, Lambda.f0.5, index.V2, V1, V2, V3)
##############################################################


















##################### V has two layers #############################
rho0=0.4
rho1=0.6
rho2=0.7
rho3=0.8
rho4=0.9
V=matrix(rho0, n, n)
V[1:(n/4),1:(n/4)]=matrix(rho1, nrow = (n/4), ncol = (n/4)) + (1-rho1)*diag((n/4))
V[((n/4)+1):(n/2), ((n/4)+1):(n/2)]=matrix(rho2, nrow = (n/4), ncol = (n/4)) + (1-rho2)*diag((n/4))
V[((n/2)+1):(3*n/4), ((n/2)+1):(3*n/4)]=matrix(rho3, nrow = (n/4), ncol = (n/4)) + (1-rho3)*diag((n/4))
V[((3*n/4)+1):n, ((3*n/4)+1):n]=matrix(rho4, nrow = (n/4), ncol = (n/4)) + (1-rho4)*diag((n/4))

Lambda=diag(c(eigen(V)$values))
Q=eigen(V)$vectors
Lambda0.5=diag(c(sqrt(eigen(V)$values)))
V0.5=Q%*%Lambda0.5%*%t(Q)
# check whether it is correct
print(max(abs(V0.5%*%V0.5-V)))

Lambda.f0.5=diag(c((eigen(V)$values)^(-0.5)))
V.f0.5=Q%*%Lambda.f0.5%*%t(Q)
print(max(abs(V.f0.5%*%V%*%V.f0.5 - diag(n))))
rm(Lambda, Lambda0.5, Q, V, Lambda.f0.5, rho0, rho1, rho2, rho3, rho4)
##############################################################################################




























reject.ols=matrix(0, nrow = run_cp, ncol = run_boxplot)
reject.gls=matrix(0, nrow = run_cp, ncol = run_boxplot)
reject.hac=matrix(0, nrow = run_cp, ncol = run_boxplot)
reject.lz=matrix(0, nrow = run_cp, ncol = run_boxplot)
reject.ehw=matrix(0, nrow = run_cp, ncol = run_boxplot)


set.seed(424)
Xboxplot=matrix(rnorm(n*d*run_cp, mean=0, sd=1), nrow=n, ncol=(d*run_cp))




for (j in 1:(run_boxplot)) {
  
  for (i in 1:run_cp) {
    
  
  
  X=Xboxplot[,((i-1)*d+1):(i*d)]
  #X=matrix(runif(n*d, min = -sqrt(3), max = sqrt(3)), n,d) # uniform distribution
  #X=matrix( (2*rbinom(n*d,1,0.5)-1) , n,d) # symmetric bernoulli
  
  
  #w=rnorm(n)
  #w=runif(n, min = -sqrt(3), max = sqrt(3))
  #w=2*rbinom(n,1,0.5)-1
  w=c(rnorm((n/4)), runif((n/2), min = -sqrt(3), max = sqrt(3)), (2*rbinom((n/4),1,0.5)-1))
  
  
  eps=V0.5%*%w
  y=X%*%beta+eps
  
  
  ######### ols ###############################
  fit.ols=lm(y~X-1)
  z.ols=summary(fit.ols)$coefficients["X1", "t value"]
  if (z.ols>z_alpha){
    reject.ols[i,j]=1
  }
  #################################################
  
  
  
  
  
  ############# gls ######################################################
  fit.gls=lm((V.f0.5%*%y)~(V.f0.5%*%X)-1)
  z.gls=summary(fit.gls)$coefficients["V.f0.5 %*% X1", "t value"]
  if (z.gls>z_alpha){
    reject.gls[i,j]=1
  }
  #################################################
  
  
  
  
  
  
  
  
  
  
  ######### Apply HAC s.e. adjustment #######################################
  # Inf is calculating p value based on N(0,1), if df is a positive integer, p is calculated by t(df).
  # vcov is the covariance matrix adjustment for fit.ols
  # by default of kernHAC, it computes a quadratic spectral kernel HAC estimator with VAR(1) prewhitening and automatic bandwidth selection based on an AR(1) approximation.
  #z.hac=coeftest(fit.ols, df=Inf, vcov = kernHAC(fit.ols))["X1","z value"] 
  # by default of NeweyWest, it calculate automatic bandwidth with prewhitening.
  z.hac=coeftest(fit.ols, df=Inf, vcov = NeweyWest(fit.ols, lag=floor(n^(1/4)), prewhite = F, adjust = T))["X1","z value"]
  if (z.hac>z_alpha){
    reject.hac[i,j]=1
  }
  #########################################################################
  
  
  
  
  
  
  ############## Liang-Zeger with noc=10 clusters ###############################
  noc.fit=10 # note that if noc.fit is small, its type I error is much larger than 0.05.
  eps.hat=as.matrix(fit.ols$residuals)
  meat.lz=matrix(0, d, d)
  for (g in 1:noc.fit) {
    meat.lz= meat.lz + t(t(eps.hat[((g-1)*(n/noc.fit)+1): (g*(n/noc.fit)),]) %*% X[((g-1)*(n/noc.fit)+1): (g*(n/noc.fit)), ])%*%(t(eps.hat[((g-1)*(n/noc.fit)+1): (g*(n/noc.fit)),]) %*% X[((g-1)*(n/noc.fit)+1): (g*(n/noc.fit)), ])
  }
  vcov.lz=solve(t(X)%*%X) %*% meat.lz %*% solve(t(X)%*%X)
  z.lz=fit.ols$coefficients["X1"]/sqrt(vcov.lz[1,1])
  if (z.lz>z_alpha){
    reject.lz[i,j]=1
  }
  ##########################################################################
  
  
  
  
  
  ########### EHW ############################
  vcov.ehw=solve(t(X)%*%X) %*% (t(X) %*% diag(fit.ols$residuals)%*%diag(fit.ols$residuals) %*% X) %*% solve(t(X)%*%X)
  z.ehw=fit.ols$coefficients["X1"]/sqrt(vcov.ehw[1,1])
  if (z.ehw>z_alpha){
    reject.ehw[i,j]=1
  }
  #############################################
  
  
  
  print(i)
 }

  print(j)
}


















cp_ols=apply(reject.ols, 2, mean)
cp_gls=apply(reject.gls, 2, mean)
cp_hac=apply(reject.hac, 2, mean)
cp_lz=apply(reject.lz, 2, mean)
cp_ehw=apply(reject.ehw, 2, mean)

boxplot(cbind(cp_ols, cp_gls, cp_hac, cp_lz, cp_ehw),
        names=c("OLS", "GLS", "HAC", "LZ 10", "EHW"),
        main=paste("Type I error n=",n," d=",d))
abline(h=alpha, col=2)
abline(h=alpha+0.005, col=2)

mean(reject.ols)
mean(reject.gls)
mean(reject.hac)
mean(reject.lz)
mean(reject.ehw)


save(list = c("reject.ehw", "reject.gls", "reject.hac", "reject.lz", "reject.ols",
     "V0.5", "Xboxplot", "alpha", "beta", "d", "n", "run_boxplot", "run_cp", "z_alpha"), 
     file = "Type I error permuted AR1 V five normal based methods n=60.Rdata")








