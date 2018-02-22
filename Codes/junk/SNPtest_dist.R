## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
# setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 250 # number of qqfamilies

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
      jdepth.vec
    }
    
    depth.mat = lapply(1:p1, loopfun)
    depth.mat = matrix(unlist(depth.mat), ncol=p1, byrow=F)
    
    ## update all concerned variables
    tail.probs = apply(depth.mat,2,function(x) (ecdf(x))(Cn.full))
    active.set = which(tail.probs > 0.6)
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

step1.depth = function(mod, sdn, adj, nboot=1e3){
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
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    jdepth.vec = mdepth.RP(jbeta.mat,
                           matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat)$dep
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  ## update all concerned variables
  # tail.probs = apply(depth.mat,2,function(x) (ecdf(x))(Cn.full))
  full.ecdf = ecdf(depth.full)
  tail.probs = apply(depth.mat,2,function(x) 1-full.ecdf(median(x)))
  tail.probs.adj = p.adjust(tail.probs, method=adj)
  which(tail.probs.adj > median(tail.probs.adj))
}

## looping for multiple instances of the analysis
get.outputs = function(signal, sdn, rho, adj, nrep){
  set.seed(11162016)
  
  loopfun = function(rep){
    library(regress)
    library(gap)
    library(fda.usc)
    library(parallel)
    source('misc_functions.R')
    source('simgen.r')
    
    ## generate y
    set.seed(rep)
    X = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    X = data.frame(X)
    y = simphe(gendat=X, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=c(6,5,4), sigmae2=c(4,3,2))
    beta = y$beta
    y = y$Y[,3]
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))
    TPTN.vec = rep(0,21)
    
    # best index for linear model BIC
    active.ind = c(1,7,41,47)
    # active.ind = c(1:10,41:50)
    parent.ind = c(4*(1:nfamily)-3, 4*(1:nfamily)-2)
    lmod = lm(form, data=datamat[parent.ind,])
    lmod.step = step(lmod, trace=F, k = log(n))
    best.BIC = which(names(datamat)[-1] %in% names(lmod.step$coef))
    # TP.BIC = c((sum(best.BIC %in% 1:6)>0), (sum(best.BIC %in% 7:10)>0),
    #                (sum(best.BIC %in% 41:46)>0), (sum(best.BIC %in% 47:50)>0))
    TP.BIC = c((1 %in% best.BIC), (7 %in% best.BIC),
               (41 %in% best.BIC), (47 %in% best.BIC))
    TN.BIC = sum((1:p)[-best.BIC] %in% (1:p)[-active.ind])
    TPTN.vec[1:7] = c(TP.BIC, sum(TP.BIC), TN.BIC, length(best.BIC))
    
    # # adjusted p-values from univariate models
    Kron = kronecker(diag(nfamily), Kmat)
    # p.vec = rep(0,p)
    # jfun = function(j){
    #   form.j = as.formula(paste0("y~", paste(names(X)[j], collapse="+")))
    #   mod.j = regress(form.j, ~Kron, data=datamat)
    #   p.vec[j] = 1 - pnorm(abs(with(mod.j, beta[2]/beta.se[2])))
    # }
    # system.time(p.vec <- as.numeric(lapply(1:p, jfun)))
    # best.mc = which(p.adjust(p.vec, 'hochberg') < 0.05)
    # # TP.BIC = c((sum(best.BIC %in% 1:6)>0), (sum(best.BIC %in% 7:10)>0),
    # #                (sum(best.BIC %in% 41:46)>0), (sum(best.BIC %in% 47:50)>0))
    # TP.mc = c((1 %in% best.mc), (7 %in% best.mc),
    #                (41 %in% best.mc), (47 %in% best.mc))
    # TN.mc = sum((1:p)[-best.mc] %in% (1:p)[-active.ind])
    # TPTN.vec[8:14] = c(TP.mc, sum(TP.mc), TN.mc, length(best.mc))
    
    mod = regress(form, ~Kron, data=datamat)
    best.index = step1.depth(mod, sdn=sdn, adj=adj)
    TP.num = c((1 %in% best.index), (7 %in% best.index),
               (41 %in% best.index), (47 %in% best.index))
    # TP.num = c((sum(best.index %in% 1:6)>0), (sum(best.index %in% 7:10)>0),
    #            (sum(best.index %in% 41:46)>0), (sum(best.index %in% 47:50)>0))
    TN.num = sum((1:p)[-best.index] %in% (1:p)[-active.ind])
    TPTN.vec[15:21] = c(TP.num, sum(TP.num), TN.num, length(best.index))
    TPTN.vec
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(16,nrep))
  TPTN.mat1 = matrix(unlist(TPTN.mat1), ncol=21, byrow=T)
  TPTN.mat.BIC = TPTN.mat1[,1:7]
  TPTN.mat.mc = TPTN.mat1[,8:14]
  TPTN.mat = TPTN.mat1[,15:21]
  z = c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat,2,sd)/c(rep(1,4),4,46,1),
        apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,46,1),
        apply(TPTN.mat.mc,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.mc,2,sd)/c(rep(1,4),4,46,1)
  )
  matrix(z, nrow=6, byrow=T)
}

##### Output function
outfun = function(signal, rho, adj, nrep, filename){
  outs = matrix(0,nrow=20,ncol=7)
  for(i in 1:20){
    z = get.outputs(signal, sdn=.1+i/10, rho, adj, nrep)
    cat(paste(i,"Done\n"))
    outs[i,] = z[1,]
  }
  
  pdf(filename,10,5)
  par(mfrow=c(1,2))
  outs[,5]
  plot(outs[,5], type='l', ylim=c(0,1), lwd=2, xlab="sd",main="TP")
  abline(h = z[5,5], col="blue")
  abline(h = z[3,5], col="red")
  legend('bottomright',c("MC","BIC","BS"),col=c("blue","red","black"), lwd=2)
  
  plot(outs[,6], type='l', ylim=c(0,1), lwd=2, xlab="sd",main="TN")
  abline(h = z[5,6], col="blue")
  abline(h = z[3,6], col="red")
  legend('bottomright',c("MC","BIC","BS"),col=c("blue","red","black"), lwd=2)
  par(mfrow=c(1,1))
  dev.off()
}

nrep = 1e1
# outfun(signal=.01, rho=.1, adj='fdr', nrep, filename='plot_pt01_rho1_fdr.pdf')
outfun(signal=.01, rho=.7, adj='fdr', nrep, filename='plot_pt01_rho7_fdr_big.pdf')
# outfun(signal=0, rho=.1, adj='fdr', nrep, filename='plot_0_rho1_fdr.pdf')
outfun(signal=0, rho=.7, adj='fdr', nrep, filename='plot_0_rho7_fdr_big.pdf')
