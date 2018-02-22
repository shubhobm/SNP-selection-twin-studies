## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
# setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(RFGLS)
library(fda.usc)
library(parallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 500 # number of qqfamilies

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
get.outputs.depth = function(signal, sdn, rho, adj, nrep){
  set.seed(11162016)
  
  loopfun = function(rep){
    
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
    
    #################################
    # best index for all-variable mixed model
    #################################
    active.ind = c(1,7,41,47)
    Kron = kronecker(diag(nfamily), Kmat)
    mod = regress(form, ~Kron, data=datamat)
    best.index = step1.depth(mod, sdn=sdn, adj=adj)
    TP.num = c((1 %in% best.index), (7 %in% best.index),
               (41 %in% best.index), (47 %in% best.index))
    TN.num = 46 - sum((1:p)[best.index] %in% (1:p)[-active.ind])
    c(TP.num, sum(TP.num), TN.num, length(best.index))
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(16,nrep))
  TPTN.mat = matrix(unlist(TPTN.mat1), ncol=7, byrow=T)
  z = c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,46,1),
        apply(TPTN.mat,2,sd)/c(rep(1,4),4,46,1))
  matrix(z, nrow=2, byrow=T)
}

get.outputs.others = function(signal, rho, adj, nrep){
  set.seed(11162016)
  
  loopfun = function(rep){
    
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
    TPTN.vec = rep(0,14)
    
    #################################
    # best index for linear model BIC
    #################################
    active.ind = c(1,7,41,47)
    parent.ind = c(4*(1:nfamily)-3, 4*(1:nfamily)-2)
    lmod = lm(form, data=datamat[parent.ind,])
    lmod.step = step(lmod, trace=F, k = log(n))
    best.BIC = which(names(datamat)[-1] %in% names(lmod.step$coef))
    TP.BIC = c((1 %in% best.BIC), (7 %in% best.BIC),
               (41 %in% best.BIC), (47 %in% best.BIC))
    TN.BIC = 46 - sum((1:p)[best.BIC] %in% (1:p)[-active.ind])
    TPTN.vec[1:7] = c(TP.BIC, sum(TP.BIC), TN.BIC, length(best.BIC))
    
    #################################
    # best index for RFGLS+multiple correction
    #################################
    # make data
    pheno = data.frame(cbind(FAMID,ID,FTYPE,INDIV,y))
    geno = data.frame(t(X))
    colnames(geno) = ID
    
    # model
    mod.gls = gls.batch(
      phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
      theta=NULL,snp.names=row.names(geno),
      input.mode=c(1,2,3),pediheader=FALSE, 
      pedicolname=c("FAMID","ID","PID","MID"),
      phen="y",covars=NULL,med=c("UN","VC"),
      outfile=NULL,col.names=TRUE,return.value=TRUE,
      covmtxfile.out=NULL,covmtxparams.out=NULL,
      sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
    
    # adjust p-values
    pval.adj = p.adjust(mod.gls$pval, method=adj)
    
    # find best index
    best.RF = which(pval.adj<.05)
    TP.num.RF = c((1 %in% best.RF), (7 %in% best.RF),
                  (41 %in% best.RF), (47 %in% best.RF))
    TN.num.RF = 46 - sum((1:p)[best.RF] %in% (1:p)[-active.ind])
    TPTN.vec[8:14] = c(TP.num.RF, sum(TP.num.RF), TN.num.RF, length(best.RF))
    
    # return output
    TPTN.vec
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(16,nrep))
  TPTN.mat = matrix(unlist(TPTN.mat1), ncol=14, byrow=T)
  TPTN.mat.BIC = TPTN.mat[,1:7]
  TPTN.mat.RF = TPTN.mat[,8:14]
  z = c(apply(TPTN.mat.RF,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.RF,2,sd)/c(rep(1,4),4,46,1),
        apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,46,1))
  matrix(z, nrow=4, byrow=T)
}

##### Output function
outfun = function(signal, rho, adj, nrep, filename){
  out.list = list()
  out.list[[1]] = get.outputs.others(signal, rho, adj, nrep)
  
  for(i in 1:20){
    out.list[[i+1]] = get.outputs.depth(signal, sdn=.1+i/10, rho, adj, nrep)
  }
  
  save(out.list, file=filename)
}

nrep = 1e2
outfun(signal=.01, rho=.1, adj='fdr', nrep, filename='rho1_pt01.Rda')
outfun(signal=.01, rho=.7, adj='fdr', nrep, filename='rho7_pt01.Rda')
outfun(signal=0, rho=.1, adj='fdr', nrep, filename='rho1_zero.Rda')
outfun(signal=0, rho=.7, adj='fdr', nrep, filename='rho7_zero.Rda')
