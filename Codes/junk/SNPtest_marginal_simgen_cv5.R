## Finding causal SNPs from twin data
## large simulation study with 500 families, 1000 SNPs
rm(list=ls())
setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)
library(caret)

source('misc_functions.R')
source('simgen.r')

nfamily = 1e2 # number of qqfamilies
p = 1e2 # number of SNPs

## generate X data family-wise
X = simgen()
p = ncol(X)

## generate two error matrices
Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1), nrow=4,ncol=4)

# other constants
Kron = kronecker(diag(nfamily), Kmat)
X = data.frame(X)

## looping for multiple instances of the analysis
nrep = 1e2
TPTN.mat = matrix(0, nrow=nrep, ncol=2)
set.seed(11162016)

for(rep in 1:nrep){
  ## generate y
  y = simphe(gendat=X, h2=rep(.1,6))
  b0.index = y[[2]]; p.causal = length(b0.index)
  y = (y[[1]])[,1]
  datamat = data.frame(cbind(y,X))
  form = as.formula(paste0("y~", paste(names(X), collapse="+")))
  
  p = ncol(X)
  n = nrow(X)
  sdn.vec = n^seq(.2, .5, by=.01)
  nboot = 1e3
  RMSE.vec = rep(0, length(sdn.vec))
  
  ## 5-fold CV
  folds = createFolds(1:nfamily, k=5, list=F)
  
  for(foldnum in 1:5){
    # train model on fold data
    train = which(folds != foldnum)
    ftrain = length(train)
    train = rep(4*train, rep(4,ftrain)) - rep(c(3,2,1,0), ftrain)
    Kron = kronecker(diag(ftrain), Kmat)
    ntrain = length(train)
    form = as.formula(paste0("y~", paste(names(X), collapse="+")))
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
        jdepth.mat[,i] = mdepth.RP(jbeta.mat[,-1],
                                   matrix(beta.hat[-1], nrow=nboot, ncol=p, byrow=T) + sdn*score.mat[,-1])$dep
      }
      
      # return mean depth of truncated model, for all values of bootstrap sd
      Cn.j = apply(jdepth.mat, 2, mean)
      Cn.j
    }
    
    Cn.mat1 = lapply(1:p, loopfun)
    # system.time(
    # Cn.mat1 <- mclapply(1:p, loopfun, mc.cores=12)
    # )
    Cn.mat = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
    
    ## find RMSE
    index.list = list()
    for(i in 1:length(sdn.vec)){
      which.final = which(Cn.mat[,i] < Cn.full)
      ibeta = beta.hat[c(1,which.final+1)]
      RMSE.vec[i] = RMSE.vec[i] + mean((y[-train] - cbind(1, as.matrix(X[-train,which.final])) %*% ibeta)^2)
      index.list[[i]] = which.final
    }
  }
  
  # find index set for best model
  best.index = index.list[[which.min(RMSE.vec)]]
  TP.num = sum(best.index %in% b0.index)
  TN.num = p - p.causal - sum(best.index %in% (1:p)[-b0.index])
  TPTN.mat[rep,] = c(TP.num, TN.num)
  cat(paste("Replication",rep,
            round(TP.num/p.causal,2),
            round(TN.num/(p-p.causal),2),"done\n"))
}

# save(TPTN.mat, file="TPTNmat.Rda")
rm(list=ls())
