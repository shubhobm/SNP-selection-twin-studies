## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
setwd('D:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)
library(caret)

source('misc_functions.R')
source('simgen.r')

nfamily = 1e2 # number of qqfamilies

## generate X data family-wise
MAF = c(0.2,0.4,0.25,0.4,0.25)
n.block = c(6,4,30,6,4)
p = sum(n.block)
n = nfamily*4
p.causal = 20

## generate two error matrices
Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1), nrow=4,ncol=4)

## looping for multiple instances of the analysis
get.outputs = function(signal){
  nrep = 1e2
  TPTN.mat = matrix(0, nrow=nrep, ncol=6)
  TPTN.mat.BIC = TPTN.mat
  set.seed(11162016)
  
  pb = txtProgressBar(0,nrep)
  for(rep in 1:nrep){
    ## generate y
    set.seed(rep)
    X = simgen(LD=c(0.7,0.7,.7,.7,0.7), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    X = data.frame(X)
    y = simphe(gendat=X, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block,
               sigma2=c(6,5,4), sigmae2=c(4,3,2))
    beta = y$beta
    y = y$Y[,3]
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))
    
    # best index for linear model BIC
    lmod = lm(form, data=datamat)
    lmod.step = step(lmod, trace=F, k = log(n))
    best.index.BIC = which(names(datamat)[-1] %in% names(lmod.step$coef))
    # TP.num.BIC = c((sum(best.index.BIC %in% 1:6)>0), (sum(best.index.BIC %in% 7:10)>0),
    #                (sum(best.index.BIC %in% 41:46)>0), (sum(best.index.BIC %in% 47:50)>0))
    TP.num.BIC = c((1 %in% best.index.BIC), (7 %in% best.index.BIC),
                   (41 %in% best.index.BIC), (47 %in% best.index.BIC))
    TN.num.BIC = p - p.causal - sum(best.index.BIC %in% 21:p)
    TPTN.mat.BIC[rep,] = c(TP.num.BIC, sum(TP.num.BIC), TN.num.BIC)
    
    sdn.vec = n^seq(.02, .5, by=.02)
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
    beta.mat = matrix(0, nrow=nboot, ncol=p)
    resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
    score.mat = t(H %*% resid.mat)
    depth.full = mdepth.RP(score.mat, score.mat)$dep
    Cn.full = mean(depth.full)
    ## we are going to recycle this score matrix for different values of bootstrap standard deviation
    
    ## loop over the parameters to save on memory:
    ## recycle memory by not storing H, covariance matrix etc for truncated models
    loopfun = function(j){
      
      set.seed(rep*j)
      ## calculate quantities for truncated model
      Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
      rj = r + X1[,j] * beta.hat[j]
      jresid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(rj, nrow=n, ncol=nboot, byrow=F)
      jscore.mat = t(Hj %*% jresid.mat)
      
      ## calculate Cn for truncated model, for a range of bootstrap variances
      jdepth.mat = matrix(0, ncol=length(sdn.vec), nrow=nboot)
      for(i in 1:length(sdn.vec)){
        sdn = sdn.vec[i]
        jbeta.mat = matrix(0, ncol=p, nrow=nboot)
        jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
          sdn * jscore.mat
        jdepth.mat[,i] = mdepth.RP(jbeta.mat,
                                   matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat)$dep
      }
      
      # return mean depth of truncated model, for all values of bootstrap sd
      Cn.j = apply(jdepth.mat, 2, mean)
      Cn.j
    }
    
    Cn.mat1 = mclapply(1:p, loopfun, mc.cores=16)
    # Cn.mat1 = lapply(1:p, loopfun)
    Cn.mat.full = matrix(unlist(Cn.mat1), ncol=length(sdn.vec), byrow=T)
    
    ## save best index list... top 10 (at most) predictors that have less depth than Cn.full
    index.list = list()
    for(i in 1:length(sdn.vec)){
      Cn.full.i = Cn.mat.full[,i]
      which.less.i = which(Cn.full.i<Cn.full)
      order.i = order(Cn.full.i[which.less.i], decreasing=F)
      if(length(order.i)>=10){
        index.list[[i]] = which.less.i[order.i]
      } else{
        index.list[[i]] = which.less.i[order.i]
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
    best.index = which(total.cnt>.9*length(sdn.vec))
    # best.index = order(total.cnt,decreasing=T)[1:10]
    
    # TP.num = c((sum(best.index %in% 1:6)>0), (sum(best.index %in% 7:10)>0),
    #            (sum(best.index %in% 41:46)>0), (sum(best.index %in% 47:50)>0))
    TP.num = c((1 %in% best.index), (7 %in% best.index),
    (41 %in% best.index), (47 %in% best.index))
    TN.num = p - p.causal - sum(best.index %in% 21:p)
    TPTN.mat[rep,] = c(TP.num, sum(TP.num), TN.num)
    # cat(paste("Replication",rep,":",
    #           round(TP.num/4,2),
    #           round(TN.num/(p-p.causal),2),"done\n"))
    setTxtProgressBar(pb,rep)
  }
  close(pb)
  c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,p-p.causal), apply(TPTN.mat,2,sd)/c(rep(1,4),4,p-p.causal),
    apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,p-p.causal), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,p-p.causal))
}

get.outputs(signal=.05)
get.outputs(signal=.02)
get.outputs(signal=.01)
# get.outputs(signal=.005)
# get.outputs(signal=.002)
# get.outputs(signal=.001)
get.outputs(signal=0)
rm(list=ls())
