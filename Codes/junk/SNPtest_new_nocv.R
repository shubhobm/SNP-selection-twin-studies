## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
# setwd('D:/Study/My projects/SNP-selection-twin-studies/Codes')

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
MAF = c(0.2,0.4,0.25,0.4,0.25)
n.block = c(6,4,6,4,30)
X = simgen(LD=c(0.7,0.7,.7,.7,0), MAF=MAF, n.block=n.block, n.person=nfamily)
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
  y = simphe(gendat=X, h2=c(.01,.01,.01,.01,0), MAF=MAF, n.block=n.block,
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
  
  # 5-fold CV
  folds = createFolds(1:nfamily, k=5, list=F)
  
  p.causal = 20
  n = nrow(X)
  sdn.vec = n^seq(.01, .5, by=.01)
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
  Cn.mat1 = mclapply(1:p, loopfun, mc.cores=16)
  # Cn.mat1 = lapply(1:p, loopfun)
  Cn.mat.full = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
  
  ## save best index list... top 10 (at most) predictors that have less depth than Cn.full
  index.list = list()
  for(i in 1:length(sdn.vec)){
    Cn.full.i = Cn.mat.full[,i]
    order.i = order(Cn.full.i[which(Cn.full.i<Cn.full)], decreasing=F)
    if(length(order.i)>=10){
      index.list[[i]] = order.i[1:20]
    } else{
      index.list[[i]] = order.i
    }
  }
  
  ## now variable selection
  cnt.list = list()
  for(j in 1:p.causal){
    cnt.list[[j]] = as.numeric(lapply(index.list, function(x) j %in% x))
  }
  pos.mat = matrix(unlist(cnt.list),nrow=p.causal, byrow=TRUE)
  
  cnt.list = list()
  for(j in (p.causal+1):p){
    cnt.list[[j]] = as.numeric(lapply(index.list, function(x) j %in% x))
  }
  neg.mat = matrix(unlist(cnt.list),nrow=p-p.causal, byrow=TRUE)
  total.cnt = c(rowSums(pos.mat),rowSums(neg.mat))
  best.index = which(total.cnt>25)
  
  TP.num = (sum(best.index %in% 1:6)>0) + (sum(best.index %in% 7:10)>0)
  TN.num = p - p.causal - sum(best.index %in% 11:50)
  TPTN.mat[rep,] = c(TP.num, TN.num)
  cat(paste("Replication",rep,":",
            round(TP.num/2,2),
            round(TN.num/(p-p.causal),2),"done\n"))
}

apply(TPTN.mat,2,mean)/c(2,p-p.causal)
  # best.index = which(total.cnt>.9*length(sdn.vec))
  best.index = order(total.cnt,decreasing=T)[1:10]
  
  TP.num = (sum(best.index %in% 1:6)>0) + (sum(best.index %in% 7:10)>0) +
    (sum(best.index %in% 11:16)>0) + (sum(best.index %in% 17:20)>0)
  TN.num = p - p.causal - sum(best.index %in% 21:p)
  TPTN.mat[rep,] = c(TP.num, TN.num)
  cat(paste("Replication",rep,":",
            round(TP.num/4,2),
            round(TN.num/(p-p.causal),2),"done\n"))
}

apply(TPTN.mat,2,mean)/c(4,p-p.causal)
apply(TPTN.mat,2,sd)/c(2,p-p.causal)
# save(TPTN.mat, file="TPTNmat.Rda")
rm(list=ls())
