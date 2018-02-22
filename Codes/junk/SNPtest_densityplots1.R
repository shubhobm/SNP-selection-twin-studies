## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(fda.usc)
library(parallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 250 # number of qqfamilies
rho = 0.7
nboot = 1e3

## preapre data structure
FAMID = rep(1:nfamily, rep(4,nfamily))*10
INDIV = rep(1:4, nfamily)
ID = FAMID + INDIV
FTYPE = 1

# pedigree
v = c(0,0,1,1)
pedigree = data.frame(FAMID=FAMID, ID=ID,
                      PID=(FAMID+2)*v,
                      MID=(FAMID+1)*v)

## generate X data family-wise
MAF = c(0.2,0.4,0.25,0.4,0.25)
n.block = c(6,4,30,6,4)
p = sum(n.block)
n = nfamily*4
p.causal = 20

## generate two error matrices
Kmat = matrix(c(1,1,0.5,0.5,
                1,1,0.5,0.5,
                0.5,0.5,1,0,
                0.5,0.5,0,1),nrow=4,byrow=T)
Env = matrix(1, nrow=4, ncol=4)

set.seed(05302017)
signal = 0.05
X = simgen(LD=rep(rho, 5), MAF=MAF,
           n.block=n.block, n.person=nfamily)
X = data.frame(X)
y = simphe.r(gendat=X, h2=c(signal,signal,0,signal,signal),
             MAF=MAF, n.block=n.block, nfam=nfamily,
             sigma2=c(6,5,4), sigmae2=c(4,3,2), type="MZ", r=.5)
beta = y$beta
y = y$Y[,3]; y = y - mean(y)

datamat = data.frame(cbind(y,X))
form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))
active.ind = c(1,7,41,47)
Kron1 = kronecker(diag(nfamily), Kmat)
Kron2 = kronecker(diag(nfamily), Env)
mod = regress(form, ~Kron1+Kron2, data=datamat)

# save relevant elements from the model in a list
# and delete full momodel object to save space
model.list = list(
  beta.hat = mod$beta,
  beta.cov = mod$beta.cov,
  W = mod$W,
  X1 = mod$X,
  H = with(mod, beta.cov %*% t(X) %*% W),
  r = with(mod, model$y - fitted))
mod = model.list

plotfun = function(sdn, Emap, ...){
  beta.hat = mod$beta.hat
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X1
  H = mod$H
  r = mod$r
  n = nrow(X1)
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = depth.dist(score.mat, score.mat, Emap)

  wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  loopfun = function(j){
    set.seed(j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = depth.dist(jbeta.mat, beta.mat, Emap)
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  ## plot to check
  pdf(paste0("plot_h",signal,"_tau",10*sdn,".pdf"), width=7, height=4)
  plot(density(depth.full), lwd=2, main=paste("6h =",signal,", s =",sdn), xlab=Emap, ...)
  for(i in 1:p){ lines(density(depth.mat[,i]), col="red", lwd=2)}
  for(i in active.ind){ lines(density(depth.mat[,i]), col="blue", lwd=2)}
  lines(density(depth.full), lwd=2)
  legend("topright", c("Non-zero indices","Zero indices","Full model"), col=c("blue","red","black"), lty=1, lwd=2)
  dev.off()
}

plotfun(.2, "E1", xlim=c(0,0.05), ylim=c(0,4e2))
plotfun(.3, "E1", xlim=c(0,0.05), ylim=c(0,4e2))
plotfun(.6, "E1", xlim=c(0,0.05), ylim=c(0,4e2))
plotfun(1, "E1", xlim=c(0,0.05), ylim=c(0,4e2))
plotfun(2, "E1", xlim=c(0,0.05), ylim=c(0,4e2))
