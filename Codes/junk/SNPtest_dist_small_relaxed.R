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
               MAF=MAF, n.block=n.block,
               sigma2=c(6,5,4), sigmae2=c(4,3,2))
    beta = y$beta
    y = y$Y[,3]
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))
    TPTN.vec = rep(0,21)
    
    # best index for linear model BIC
    # active.ind = c(1,7,41,47)
    active.ind = c(1:10,41:50)
    parent.ind = c(4*(1:nfamily)-3, 4*(1:nfamily)-2)
    lmod = lm(form, data=datamat[parent.ind,])
    lmod.step = step(lmod, trace=F, k = log(n))
    best.BIC = which(names(datamat)[-1] %in% names(lmod.step$coef))
    TP.BIC = c((sum(best.BIC %in% 1:6)>0), (sum(best.BIC %in% 7:10)>0),
                   (sum(best.BIC %in% 41:46)>0), (sum(best.BIC %in% 47:50)>0))
    TN.BIC = sum((1:p)[-best.BIC] %in% (1:p)[-active.ind])
    TPTN.vec[1:7] = c(TP.BIC, sum(TP.BIC), TN.BIC, length(best.BIC))
    
    # # adjusted p-values from univariate models
    Kron = kronecker(diag(nfamily), Kmat)
    mod = regress(form, ~Kron, data=datamat)
    best.index = step1.depth(mod, sdn=sdn, adj=adj)
    TP.num = c((sum(best.index %in% 1:6)>0), (sum(best.index %in% 7:10)>0),
               (sum(best.index %in% 41:46)>0), (sum(best.index %in% 47:50)>0))
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
  z = c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,30,1), apply(TPTN.mat,2,sd)/c(rep(1,4),4,30,1),
        apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,30,1), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,30,1),
        apply(TPTN.mat.mc,2,mean)/c(rep(1,4),4,30,1), apply(TPTN.mat.mc,2,sd)/c(rep(1,4),4,30,1)
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
outfun(signal=.01, rho=.1, adj='fdr', nrep, filename='plot_pt01_rho1_fdr_relaxed.pdf')
outfun(signal=.01, rho=.7, adj='fdr', nrep, filename='plot_pt01_rho7_fdr_relaxed.pdf')
outfun(signal=0, rho=.1, adj='fdr', nrep, filename='plot_0_rho1_fdr_relaxed.pdf')
outfun(signal=0, rho=.7, adj='fdr', nrep, filename='plot_0_rho7_fdr_relaxed.pdf')

## plotting function
plot(density(depth.full), xlim=c(0,.4), ylim=c(0,40),
     main=paste("sd = ",sdn))
for(i in 1:p){
  lines(density(depth.mat[,i]),col=adjustcolor(grey(.5),alpha.f=.2),lwd=2)
}
for(i in active.ind){
  lines(density(depth.mat[,i]),col="blue",lwd=2)
  # abline(v=median(depth.mat[,i]),lty=2)
}
lines(density(depth.full), lwd=2,col="red")
