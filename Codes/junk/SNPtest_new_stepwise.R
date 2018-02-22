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

## stepwise depth-model selection function
step.depth = function(mod, sdn, nboot=1e3){
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X
  n1 = nrow(X1)
  H = beta.cov %*% t(X1) %*% W
  r = mod$model$y - mod$fitted
  n = nrow(X1)
  p1 = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p1)
  resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  active.set = 1:p1
  p0 = 1
  while(p1>p0){
    
    p1 = length(active.set)
    wildmat = matrix(rnorm(n*nboot), ncol=nboot)
    loopfun = function(j){
      set.seed(j)
      ## calculate quantities for truncated model
      Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
      rj = r + X1[,j] * beta.hat[j]
      jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
      jscore.mat = t(Hj %*% jresid.mat)
      
      ## calculate Cn for truncated model, for a range of bootstrap variances
      jbeta.mat = matrix(0, ncol=p1, nrow=nboot)
      jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p1-1, byrow=T) +
        sdn * jscore.mat
      jdepth.vec = mdepth.RP(jbeta.mat,
                             matrix(beta.hat, nrow=nboot, ncol=p1, byrow=T)+ sdn*score.mat)$dep
      # return mean depth of truncated model, for all values of bootstrap sd
      mean(jdepth.vec)
    }
    
    Cn.vec = as.numeric(lapply(1:p1, loopfun))
    
    ## update all concerned variables
    active.set = which(Cn.vec < Cn.full)
    p0 = length(active.set)
    r = r + X1[,-active.set] %*% as.matrix(beta.hat[-active.set],ncol=1)
    beta.hat = beta.hat[active.set]
    if(p0<2){
      break
    } # if too less variables then terminate
    
    
    beta.cov = beta.cov[active.set,active.set]
    X1 = X1[,active.set]
    if(is.numeric(X1)){
      X1 = as.matrix(X1, ncol=1)
    }
    resid.mat = matrix(rnorm(n*nboot), ncol=nboot) *
      matrix(r, nrow=n, ncol=nboot, byrow=F)
    H = beta.cov %*% t(X1) %*% W
    score.mat = t(H %*% resid.mat)
    depth.full = mdepth.RP(score.mat, score.mat)$dep
    Cn.full = mean(depth.full)
  }
  
  z = which(mod$beta %in% beta.hat)
  if(p0==0){
    z = NULL
  }
  z
}

## looping for multiple instances of the analysis
get.outputs = function(signal,sdn){
  nrep = 1e2
  set.seed(11162016)
  
  loopfun = function(rep){
    library(regress)
    library(gap)
    library(fda.usc)
    library(parallel)
    library(caret)
    source('misc_functions.R')
    source('simgen.r')
    
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
    TPTN.vec = rep(0,14)
    
    # best index for linear model BIC
    lmod = lm(form, data=datamat)
    lmod.step = step(lmod, trace=F, k = log(n))
    best.index.BIC = which(names(datamat)[-1] %in% names(lmod.step$coef))
    # TP.num.BIC = c((sum(best.index.BIC %in% 1:6)>0), (sum(best.index.BIC %in% 7:10)>0),
    #                (sum(best.index.BIC %in% 41:46)>0), (sum(best.index.BIC %in% 47:50)>0))
    TP.num.BIC = c((1 %in% best.index.BIC), (7 %in% best.index.BIC),
                   (41 %in% best.index.BIC), (47 %in% best.index.BIC))
    TN.num.BIC = p - p.causal - sum(best.index.BIC %in% 11:40)
    TPTN.vec[1:7] = c(TP.num.BIC, sum(TP.num.BIC), TN.num.BIC, length(best.index.BIC))
    
    Kron = kronecker(diag(nfamily), Kmat)
    mod = regress(form, ~Kron, data=datamat)
    best.index = step.depth(mod, sdn=sdn)
    TP.num = c((1 %in% best.index), (7 %in% best.index),
               (41 %in% best.index), (47 %in% best.index))
    # TP.num = c((sum(best.index %in% 1:6)>0), (sum(best.index %in% 7:10)>0),
    #            (sum(best.index %in% 41:46)>0), (sum(best.index %in% 47:50)>0))
    TN.num = p - p.causal - sum(best.index %in% 11:40)
    TPTN.vec[8:14] = c(TP.num, sum(TP.num), TN.num, length(best.index))
    
    ## recycle memory by not storing H, covariance matrix etc for truncated models
    # cat(paste("Replication",rep,":",
    #           round(TP.num/4,2),
    #           round(TN.num/(p-p.causal),2),"done\n"))
    setTxtProgressBar(pb,rep)
    TPTN.vec
  }
  
  pb = txtProgressBar(0,nrep)
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(16,nrep))
  close(pb)
  
  TPTN.mat1 = matrix(unlist(TPTN.mat1), ncol=14, byrow=T)
  TPTN.mat.BIC = TPTN.mat1[,1:7]
  TPTN.mat = TPTN.mat1[,8:14]
  z = c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,p-p.causal,1), apply(TPTN.mat,2,sd)/c(rep(1,4),4,p-p.causal,1),
        apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,p-p.causal,1), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,p-p.causal,1))
  matrix(z, nrow=4, byrow=T)
}

outs = matrix(0,nrow=20,ncol=7)
for(i in 1:20){
  z = get.outputs(signal=.02,sdn=.9+i/10)
  outs[i,] = z[1,]
}
save(outs,file="outs.rda")
pdf('plot.pdf',5,5)
plot(outs[,5], type='l', ylim=c(0,1), lwd=2, col="blue",xlab="sd",main="TP")
abline(h = z[3,5], col="red")
dev.off()

pdf('plot1.pdf',5,5)
plot(outs[,6], type='l', ylim=c(0,1), lwd=2, col="blue",xlab="sd",main="TN")
abline(h = z[3,6], col="red")
dev.off()

rm(list=ls())
