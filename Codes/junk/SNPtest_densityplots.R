## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
# setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
library(lmm)
library(fda.usc)
library(parallel)

source('misc_functions.R')

set.seed(11092015)
p = 10
sig = 1
nboot = 1e2
n.iter = 1e2
beta = c(0,rep(.2,2),rep(0,p-3)) # change this so that length 9
beta.names = "x2x3"
X.names = c('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10')

D = matrix(c(9,   4.8, 0.6, 0,
             4.8, 4,   1,   0,
             0.6, 1,   1,   0,
             0,   0,   0,   0),
           ncol=4, byrow=T)

## set up data
ni = 4
m = 250
n = ni*m
# sdm = seq(1,10, by=.1)
sdm = 1:24
p = length(beta)
pZ = 4
X = cbind(1,matrix(runif(n*(p-1), -2, 2), ncol=p-1))
Z = X[,1:pZ]
subj = rep(1:m, rep(ni,m))

# construct error covariance matrix
W.true = matrix(0, ncol=n, nrow=n)
for(im in 1:m){
  iind = ((im-1)*ni+1):(im*ni)
  Zi = Z[iind,]
  W.true[iind,iind] = sig * diag(ni) + Zi %*% D %*% t(Zi)
}
W.true = solve(W.true)

## Generate samples
set.seed(04182017)
y = rep(0,n)
for(im in 1:m){
  iind = ((im-1)*ni+1):(im*ni)
  y[iind] = my.mvrnorm(1, X[iind,] %*% beta, W.true[iind,iind])
}

## full model constants
lmmfull = ecmeml1(y=y, subj=subj, pred=X, xcol=1:p, zcol=1:pZ)

W = matrix(0, nrow=n, ncol=n)
for(im in 1:m){
  iind = ((im-1)*ni+1):(im*ni)
  Zi = Z[iind,]
  W[iind,iind] = with(lmmfull, solve(sigma2 * diag(ni) + Zi %*% psi %*% t(Zi)))
}
beta.hat = lmmfull$beta
H = lmmfull$cov.beta %*% t(X) %*% W
r = y - X %*% beta.hat; r = r - mean(r)

plotfun = function(sd){
  SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
  
  # depth model selection
  beta.mat = matrix(0 ,nrow=nboot, ncol=p)
  for(i in 1:nboot){
    iresid = as.matrix(sd * rep(rnorm(m), rep(ni,m)) * r, ncol=1)
    beta.mat[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
  }
  
  beta.mat1 = matrix(0 ,nrow=nboot, ncol=p)
  for(i in 1:nboot){
    iresid = as.matrix(sd * rep(rnorm(m), rep(ni,m)) * r, ncol=1)
    beta.mat1[i,] = as.numeric(beta.hat) + as.numeric(H %*% iresid)
  }
  SSPmat.d[,p+1] = mdepth.RP(beta.mat1, beta.mat)$dep
  
  ## now marginal models
  for(j in 1:p){
    
    jbeta.mat = matrix(0 , nrow=nboot, ncol=p)
    for(i in 1:nboot){
      iresid = sd * rep(rnorm(m), rep(ni,m))* (y - X[,-j] %*% beta.hat[-j])
      iresid = iresid - mean(iresid)
      jbeta.mat[i,-j] = beta.hat[-j] + H[-j,] %*% iresid
    }
    
    SSPmat.d[,j] = mdepth.RP(jbeta.mat, beta.mat)$dep
  }
  
  ## plot to check
  plot(density(SSPmat.d[,p+1]), xlim=c(0,.5), ylim=c(0,15), lwd=2)
  for(i in 1:p){ lines(density(SSPmat.d[,i]), col="red", lwd=2)}
  for(i in 2:3){ lines(density(SSPmat.d[,i]), col="blue", lwd=2)}
  lines(density(SSPmat.d[,p+1]), lwd=2)
  legend("topright", c("Non-zero indices","Zero indices","Full model"), col=c("blue","red","black"), lty=1, lwd=2)
}

plotfun(.01)
plotfun(1)
plotfun(10)
