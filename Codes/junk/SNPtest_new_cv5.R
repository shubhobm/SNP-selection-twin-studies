## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
# setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)
library(caret)

source('misc_functions.R')
source('simgen.r')

nfamily = 1e2 # number of qqfamilies

## generate X data family-wise
MAF = c(0.2,0.4,0.25)
n.block = c(6,4,40)
X = simgen(LD=c(0.9,0.8,0), MAF=MAF, n.block=n.block, n.person=nfamily)
p = ncol(X)
X = data.frame(X)

## generate two error matrices
Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1), nrow=4,ncol=4)

## looping for multiple instances of the analysis
nrep = 1e2
TPTN.mat = matrix(0, nrow=nrep, ncol=2)
set.seed(11162016)

for(rep in 1:nrep){
  ## generate y
  y = simphe(gendat=X, h2=c(.2,.2,0), MAF=MAF, n.block=n.block,
             sigma2=c(6,5,4), sigmae2=c(4,3,2))
  beta = y$beta
  y = y$Y[,3]
  datamat = data.frame(cbind(y,X))
  form = as.formula(paste0("y~", paste(names(X), collapse="+")))

  p = ncol(X)
  p.causal = 10
  n = nrow(X)
  sdn.vec = n^seq(.2, .5, by=.01)
  nboot = 1e3
  RMSE.vec = rep(0, length(sdn.vec))
  
  ## first get full model to calculate index sets for each sdn
  Kron = kronecker(diag(nfamily), Kmat)
  mod = regress(form, ~Kron, data=datamat)
  beta.hat = mod$beta
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X
  n1 = nrow(X1)
  H = beta.cov %*% t(X1) %*% W
  r = with(mod, y - fitted)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
  resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat[,-1], score.mat[,-1])$dep
  Cn.full = mean(depth.full)
  ## we are going to recycle this score matrix for different values of bootstrap standard deviation
  
  ## loop over the parameters to save on memory:
  ## recycle memory by not storing H, covariance matrix etc for truncated models
  loopfun = function(j){
    
    set.seed(rep*j)
    ## calculate quantities for truncated model
    Hj = H[-(j+1),] - outer(beta.cov[j+1,-(j+1)], as.numeric(W %*% X1[,j+1]))
    rj = r + X1[,j+1] * beta.hat[j+1]
    jresid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    jdepth.mat = matrix(0, ncol=length(sdn.vec), nrow=nboot)
    for(i in 1:length(sdn.vec)){
      sdn = sdn.vec[i]
      jbeta.mat = matrix(0, ncol=p+1, nrow=nboot)
      jbeta.mat[,-(j+1)] = matrix(beta.hat[-(j+1)], nrow=nboot, ncol=p, byrow=T) + sdn * jscore.mat
      jdepth.mat[,i] = mdepth.RP(matrix(beta.hat[-1], nrow=nboot, ncol=p, byrow=T) + sdn*score.mat[,-1],
                                 jbeta.mat[,-1])$dep
    }
    
    # return mean depth of truncated model, for all values of bootstrap sd
    Cn.j = apply(jdepth.mat, 2, mean)
    Cn.j
  }
  
  Cn.mat1 = mclapply(1:p, loopfun, mc.cores=10)
  # Cn.mat1 = lapply(1:p, loopfun)
  Cn.mat.full = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
  
  ## save best index list
  index.list = list()
  for(i in 1:length(sdn.vec)){
    index.list[[i]] = which(Cn.mat.full[,i] < Cn.full)
  }
  
  ## now do cross validation
  # 5-fold CV
  folds = createFolds(1:nfamily, k=5, list=F)
  
  for(foldnum in 1:5){
    # train model on fold data
    train = which(folds != foldnum)
    ftrain = length(train)
    train = rep(4*train, rep(4,ftrain)) - rep(c(3,2,1,0), ftrain)
    Kron = kronecker(diag(ftrain), Kmat)
    ntrain = length(train)
    mod = regress(form, ~Kron, data=datamat[train,])
    
    # full model quantities
    beta.hat = mod$beta
    beta.cov = mod$beta.cov
    W = mod$W
    X1 = mod$X
    n1 = nrow(X1)
    H = beta.cov %*% t(X1) %*% W
    r = with(mod, y[train] - fitted)
    
    ## matrix of full model bootstrap betas
    beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
    resid.mat = matrix(rnorm(ntrain*nboot), ncol=nboot) * matrix(r, nrow=ntrain, ncol=nboot, byrow=F)
    score.mat = t(H %*% resid.mat)
    depth.full = mdepth.RP(score.mat[,-1], score.mat[,-1])$dep
    Cn.full = mean(depth.full)
    ## we are going to recycle this score matrix for different values of bootstrap standard deviation
    
    ## loop over the parameters to save on memory:
    ## recycle memory by not storing H, covariance matrix etc for truncated models
    loopfun = function(j){
      
      ## calculate quantities for truncated model
      Hj = H[-(j+1),] - outer(beta.cov[j+1,-(j+1)], as.numeric(W %*% X1[,j+1]))
      rj = r + X1[,j+1] * beta.hat[j+1]
      jresid.mat = matrix(rnorm(ntrain*nboot), ncol=nboot) * matrix(rj, nrow=ntrain, ncol=nboot, byrow=F)
      jscore.mat = t(Hj %*% jresid.mat)
      
      ## calculate Cn for truncated model, for a range of bootstrap variances
      jdepth.mat = matrix(0, ncol=length(sdn.vec), nrow=nboot)
      for(i in 1:length(sdn.vec)){
        sdn = sdn.vec[i]
        jbeta.mat = matrix(0, ncol=p+1, nrow=nboot)
        jbeta.mat[,-(j+1)] = matrix(beta.hat[-(j+1)], nrow=nboot, ncol=p, byrow=T) + sdn * jscore.mat
        jdepth.mat[,i] = mdepth.RP(matrix(beta.hat[-1], nrow=nboot, ncol=p, byrow=T) + sdn*score.mat[,-1],
                                   jbeta.mat[,-1])$dep
      }
      
      # return mean depth of truncated model, for all values of bootstrap sd
      Cn.j = apply(jdepth.mat, 2, mean)
      Cn.j
    }
    
    # Cn.mat1 = mclapply(1:p, loopfun)
    # system.time(
    Cn.mat1 <- mclapply(1:p, loopfun, mc.cores=12)
    # )
    Cn.mat = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
    
    ## find RMSE
    for(i in 1:length(sdn.vec)){
      which.final = which(Cn.mat[,i] < Cn.full)
      ibeta = beta.hat[c(1,which.final+1)]
      RMSE.vec[i] = RMSE.vec[i] + sum((y[-train] - cbind(1, as.matrix(X[-train,which.final])) %*% ibeta)^2)
    }
  }
  
  # find index set for best model
  RMSE.vec = RMSE.vec/n
  best.index = index.list[[which.min(RMSE.vec)]]
  TP.num = (sum(best.index %in% 1:6)>0) + (sum(best.index %in% 7:10)>0)
  TN.num = p - p.causal - sum(best.index %in% 11:50)
  TPTN.mat[rep,] = c(TP.num, TN.num)
  cat(paste("Replication",rep,":",
            round(TP.num/2,2),
            round(TN.num/(p-p.causal),2),"done\n"))
}

apply(TPTN.mat,2,mean)/c(2,p-p.causal)
apply(TPTN.mat,2,sd)/c(2,p-p.causal)
# save(TPTN.mat, file="TPTNmat.Rda")
rm(list=ls())
