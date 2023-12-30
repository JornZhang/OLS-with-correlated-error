rm(list = ls())

library("sandwich")
library("lmtest")
#devtools::install_url("https://cran.r-project.org/src/contrib/Archive/gaoptim/gaoptim_1.1.tar.gz")
library("gaoptim")

alpha=0.05

z_alpha=qnorm(alpha, mean = 0, sd=1, lower.tail = F)

n=100

d=2


h=c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5)

beta1=h/sqrt(n-d)




run=200








###### lihua rank test function ###############################################################
# note that m=(1/alpha)-1, and n should satisfy n \geq d*m; n/(1+m) is an integer.
# The output is delta.star, and (eta.star, Pi_1 eta.star, ..., Pi_m eta.star) ready to be multiplied by y
lihua=function(n, d, alpha, X){
  m=(1/alpha)-1
  tperm=n/(m+1)
  BX=matrix(0, nrow = n, ncol = d*m)
  Pi_mtX = rbind(X[(n-m*tperm+1):n,], X[1:(n-m*tperm),])
  BX[,1:d]=X-Pi_mtX
  for (g in 1:(m-1)) {
    BX[,(g*d+1):((g+1)*d)]=rbind(X[(n-g*tperm+1):n, ], X[1:(n-g*tperm), ])-Pi_mtX
  }
  
  eta.tilde=as.matrix(lm(BX[,1]~BX[,2:(d*m)]-1)$residuals)
  delta.star=as.numeric(sqrt(t(eta.tilde)%*%eta.tilde))
  eta.star=eta.tilde/delta.star
  
  Eta.total=matrix(0, nrow=n, ncol=(m+1))
  Eta.total[,1]=eta.star
  for (g in 1:m) {
    Eta.total[,(g+1)] = c(eta.star[(g*tperm+1):n,1], eta.star[1:(g*tperm),1]) 
  }
  return(list(delta.star=delta.star, Eta.total=Eta.total))
}
#######################################################################################








############# use gaoptim to find the best ordering #####################################
# input X, n=nrow(X), d=ncol(X), rounds and popSize are for GA algorithm.
# output is the best ordering to have largest delta.star and the corresponding
lihua_GA <- function(X, n, d, alpha, ga_obj = NULL, rounds = 100, popSize = 10, ...){
  if (is.null(ga_obj)){
    objective_fun <- function(ordering){
      lihua(n, d, alpha, X[ordering, ])$delta.star
    }
    ga_obj <- gaoptim::GAPerm(objective_fun, n, popSize = popSize)
  }
  ga_obj$evolve(rounds)
  best_ordering <- ga_obj$bestIndividual()
  return(list(ordering = best_ordering,
              ga_obj = ga_obj))
}
########################################################################################


















############### bootstrap #########################################
bspair <- function(y, X){
  index.bs <- sample(n, replace=T)
  return(lm(y[index.bs,]~X[index.bs,]-1)$coef)
}
#######################################################################










########### FL ###########################################
# y=X[,1] beta1 + X[,-1] beta[-1] + epsilon and assume epsilon is exchangable
# Input is fit.nuisance=lm(y~X[,-1]-1). Permute fit.nuisance$residuals and add them back to fitZ$fitted.values
# Another input is the X=(X[,1],X[,-1])
# Output is t statistics from permuted model
FLt=function(X, fit.nuisance){
  y_FLperm=sample(fit.nuisance$residuals) + fit.nuisance$fitted.values
  fit.permfull=lm(y_FLperm~X-1)
  return(summary(fit.permfull)$coefficients[, "t value"])
}
####################################################################





















































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
rho=0.9
noc.true=10
V=matrix(0, nrow=n, ncol=n)
Vcluster=matrix(rho, (n/noc.true), (n/noc.true)) + (1-rho)*diag((n/noc.true))
for (g in 1:noc.true) {
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
Atmp = matrix(rnorm(n*K,0,10),n)
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










############ V has noc equal-sized clusters with different rho's  #########################
noc.true=5
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



























power_ols=matrix(0,nrow = run,ncol=length(h))
power_gls=matrix(0,nrow = run,ncol=length(h))
power_hac=matrix(0,nrow = run,ncol=length(h))
power_lz=matrix(0,nrow = run,ncol=length(h))
power_ehw=matrix(0,nrow = run,ncol=length(h))
power_lihua=matrix(0,nrow = run,ncol=length(h))
power_bspair=matrix(0,nrow = run,ncol=length(h))
power_FL=matrix(0,nrow = run,ncol=length(h))
#power_lz1=matrix(0,nrow = run,ncol=length(h))


for (j in 1:length(h)) {
  
  beta=c(beta1[j], 1)
  
  for (i in 1:run) {
    
    #X=matrix(rnorm(n*d, mean=0, sd=1), n,d)
    X=matrix(runif(n*d, min = -sqrt(3), max = sqrt(3)), n,d) # uniform distribution
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
      power_ols[i,j]=1
    }
    #################################################
    
    
    
    
    
    ############# gls ######################################################
    fit.gls=lm((V.f0.5%*%y)~(V.f0.5%*%X)-1)
    z.gls=summary(fit.gls)$coefficients["V.f0.5 %*% X1", "t value"]
    if (z.gls>z_alpha){
      power_gls[i,j]=1
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
      power_hac[i,j]=1
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
      power_lz[i,j]=1
    }
    ##########################################################################
    
    
    
    
    
    
    
    
    ############## EHW SE ###############################
    vcov.ehw=solve(t(X)%*%X) %*% (t(X) %*% diag(fit.ols$residuals)%*%diag(fit.ols$residuals) %*% X) %*% solve(t(X)%*%X)
    z.ehw=fit.ols$coefficients["X1"]/sqrt(vcov.ehw[1,1])
    if (z.ehw>z_alpha){
      power_ehw[i,j]=1
    }
    ##########################################################################
    
   
    
    
    
    
    
    
    
    
    
    ############################ Apply genetic optimization ##########################
    #select.lihua_delta.star=replicate(1000, lihua.permX(n, d, alpha, X))
    #lihua_optm.index=as.vector(select.lihua_delta.star[2:(n+1),which.max(select.lihua_delta.star[1,])])
    lihua_optm.index=lihua_GA(X, n, d, alpha, ga_obj = NULL, rounds = 400, popSize = 10)$ordering
    Xoptim_lihua=X[lihua_optm.index,]
    lihua(n, d, alpha, Xoptim_lihua)$delta.star
    Eta.total.lihua_optim=lihua(n, d, alpha, Xoptim_lihua)$Eta.total
    yoptim_lihua=y[lihua_optm.index,]
    S_lihua=t(yoptim_lihua)%*%Eta.total.lihua_optim # we test one side alternative, so no median procedure is needed.
    R0_lihua=rank(-S_lihua)[1]
    if ( (R0_lihua/(1/alpha)) <= alpha ){
      power_lihua[i,j] = 1
    }
    #######################################################
    
    
    
    
    
    
    
    
    ############# paired bootstrap ######################################
    bspair_coef=replicate(1000, bspair(y,X))
    # use bootstrap estimated betahats to adjust standard error
    #z.bspair=fit.ols$coefficients["X1"]/sd(bspair_coef[1,])
    #if (z.bspair>z_alpha){
    #  power_bspair[i,j]=1
    #}
    #use bootstrap estimated betahats to construct CI's
    ci_bspair.oneside=quantile(bspair_coef[1,], 0.05)
    if (0<ci_bspair.oneside){
        power_bspair[i,j]=1
      }
    #####################################################################
    
    
    
    
    
    
    
    ################ FL ##############################
    fit.nuisance=lm(y~X[,-1]-1)
    FLrun=1000
    FLtperm=replicate(FLrun, FLt(X,fit.nuisance = fit.nuisance))
    if (sum(FLtperm[1,]>=z.ols)/FLrun<alpha){
      power_FL[i,j]=1
    }
    #################################################
    
    
    
    ############## Liang-Zeger with 1 cluster ###############################
    # vcov.lz1=solve(t(X)%*%X) %*% t(X)%*%eps.hat %*%t(eps.hat)%*%X %*% solve(t(X)%*%X)
    # z.lz1=fit.ols$coefficients["X1"]/sqrt(vcov.lz1[1,1])
    # if (z.lz1>z_alpha){
    #  power_lz1[i,j]=1
    # }
    ##########################################################################
    
    
    
    
    print(i)
  }
  print(j)
}




plot(x=h, y=apply(power_ols, 2, mean), type = "b", pch=19, col=1, 
     ylab = " ", main = paste("Comparison of Power n=",n," d=",d, " rep=", run),
     yaxt="n")
lines(x=h,y=apply(power_gls, 2, mean), type="b",pch=1,col=2)
lines(x=h,y=apply(power_hac, 2, mean), type="b",pch=19,col=3)
lines(x=h,y=apply(power_lz, 2, mean), type="b",pch=19,col=4)
lines(x=h,y=apply(power_ehw, 2, mean), type="b",pch=19,col=5)
lines(x=h,y=apply(power_lihua, 2, mean), type="b",pch=10,col=6)
lines(x=h,y=apply(power_bspair, 2, mean), type="b",pch=10,col=7)
lines(x=h,y=apply(power_FL, 2, mean), type="b",pch=11,col=rgb(0.9, 0.5,0.1))
abline(h=alpha, col=rgb(0.5, 0.5, 0.5), lty=3)
legend("bottomright", legend = c("OLS", "GLS", "HAC", "LZ 10", "EHW", "Lihua", "bootstrap", "FL"),
       pch = c(19,1,19,19,19,10,10,11), lty = c(1,1,1,1,1,1,1,1), col=c(1,2,3,4,5,6,7,rgb(0.9, 0.5,0.1)) )
ytick<-c(0.05, 0.2, 0.5, 1)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1], ytick,  
     labels = ytick, srt = 0, pos = 2, xpd = TRUE)




apply(power_ols, 2, mean)
apply(power_gls, 2, mean) - apply(power_ols, 2, mean)
apply(power_hac, 2, mean) - apply(power_ols, 2, mean)
apply(power_lz, 2, mean) - apply(power_ols, 2, mean)
apply(power_ehw, 2, mean) - apply(power_ols, 2, mean)
apply(power_lihua, 2, mean) - apply(power_ols, 2, mean)
apply(power_bspair, 2, mean) - apply(power_ols, 2, mean)
apply(power_FL, 2, mean) - apply(power_ols, 2, mean)




plot(x=h, y=apply(power_gls, 2, mean) - apply(power_ols, 2, mean), 
     type = "b", pch=1, col=2, ylim = c(-0.1, 0.6),
     ylab = " ", main = paste("Power difference compared with OLS n=",n," d=",d, " rep=", run),
     yaxt="n")
abline(h=0, col=rgb(0.5, 0.5, 0.5), lty=3)
lines(x=h,y=apply(power_hac, 2, mean) - apply(power_ols, 2, mean), type="b",pch=19,col=3)
lines(x=h,y=apply(power_lz, 2, mean) - apply(power_ols, 2, mean), type="b",pch=19,col=4)
lines(x=h,y=apply(power_ehw, 2, mean) - apply(power_ols, 2, mean), type="b",pch=19,col=5)
lines(x=h,y=apply(power_lihua, 2, mean) - apply(power_ols, 2, mean), type="b",pch=10,col=6)
lines(x=h,y=apply(power_bspair, 2, mean) - apply(power_ols, 2, mean), type="b",pch=10,col=7)
lines(x=h,y=apply(power_FL, 2, mean) - apply(power_ols, 2, mean), type="b",pch=10,col=rgb(0.9, 0.5,0.1))
legend("topright", legend = c("GLS", "HAC", "LZ 10", "EHW", "Lihua", "bootstrap", "FL"),
       pch = c(1,19,19,19,10,10,11), lty = c(1,1,1,1,1,1,1), col=c(2,3,4,5,6,7,rgb(0.9, 0.5,0.1)) )
ytick<-c(-0.1, 0, 0.2, 0.5)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1], ytick,  
     labels = ytick, srt = 0, pos = 2, xpd = TRUE)




save(list = c("power_ols", "power_gls", "power_hac", "power_lz", "power_ehw", "power_lihua", "power_bspair", "power_FL",
              "V0.5", "V.f0.5",
              "alpha", "d", "h", "n", "noc.fit", "run", "z_alpha"),
     file = "Power Comparison four methods V is blk diag with different rhos and non identical w.Rdata")

load("Power Comparison four methods V is factor model K=2 n=100.Rdata")
