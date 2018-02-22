## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
# setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(RFGLS)
library(fda.usc)
library(parallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 250 # number of qqfamilies

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

## looping for multiple instances of the analysis
get.outputs.depth = function(signal, sdn.vec, rho, q, mult, nrep){
  set.seed(11162016)
  
  loopfun = function(rep){
    
    ## generate y
    set.seed(rep)
    X = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    X = data.frame(X)
    y = simphe.r(gendat=X, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=c(6,5,4), sigmae2=c(4,3,2), type="MZ", r=.5)
    beta = y$beta
    y = y$Y[,3]; y = y - mean(y)
    
    # generate test data
    Xtest = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    ytest = simphe.r(gendat=Xtest, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=c(6,5,4), sigmae2=c(4,3,2), type="MZ", r=.5)
    ytest = ytest$Y[,3]; ytest = ytest - mean(ytest)
    
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))

    #################################
    # best index for all-variable mixed model
    #################################
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
    
    ## now cycle through range of sdn values
    nsdn = length(sdn.vec)
    TPTN.sdn = matrix(0, nrow=nsdn, ncol=8)
    
    for(i in 1:nsdn){
      best.index = step11.depth(mod, sdn=sdn.vec[i], q=q, mult=mult)[[2]]
      TP.num = c((1 %in% best.index), (7 %in% best.index),
                 (41 %in% best.index), (47 %in% best.index))
      TN.num = 46 - sum((1:p)[best.index] %in% (1:p)[-active.ind])
      test.error = sum(diag(cov(matrix(ytest, nrow=nfamily, byrow=F))))
      if(length(best.index)>0){
        err = as.numeric(ytest - as.matrix(Xtest[,best.index]) %*% mod$beta.hat[best.index])
        test.error = sum(diag(cov(matrix(err, nrow=nfamily, byrow=F))))
        # test.error = sum(diag(matrix(err, nrow=nfamily, byrow=F)^2))
      }
      TPTN.sdn[i,] = c(TP.num, sum(TP.num), TN.num, length(best.index), test.error)
    }
    
    best.sdn = which.min(TPTN.sdn[,8])
    TPTN.sdn[best.sdn,-8]
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
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
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
  TPTN.mat = matrix(unlist(TPTN.mat1), ncol=14, byrow=T)
  TPTN.mat.BIC = TPTN.mat[,1:7]
  TPTN.mat.RF = TPTN.mat[,8:14]
  z = c(apply(TPTN.mat.RF,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.RF,2,sd)/c(rep(1,4),4,46,1),
        apply(TPTN.mat.BIC,2,mean)/c(rep(1,4),4,46,1), apply(TPTN.mat.BIC,2,sd)/c(rep(1,4),4,46,1))
  matrix(z, nrow=4, byrow=T)
}

##### Output function
outfun = function(rho, adj, nrep, q, mult, filename){
  out.list = list()
  for(i in 1:length(h.vec)){
    z = rbind(get.outputs.others(signal=h.vec[i], rho, adj, nrep),
              get.outputs.depth(signal=h.vec[i], sdn.vec, rho, q, mult, nrep))
    out.list[[i]] = z
  }
  save(out.list, file=filename)
}

nrep = 100
sdn.vec = seq(.1,1,by=.05)
h.vec = c(.1, .05, .02, .01, 0)
adj = 'fdr'
outfun(rho=.7, adj=adj, nrep, q=.9, mult=.8, filename='all_q9_mult8.Rda')
outfun(rho=.7, adj=adj, nrep, q=.5, mult=.8, filename='all_q5_mult8.Rda')
outfun(rho=.7, adj=adj, nrep, q=.2, mult=.8, filename='all_q2_mult8.Rda')
outfun(rho=.7, adj=adj, nrep, q=.1, mult=.8, filename='all_q1_mult8.Rda')
outfun(rho=.7, adj=adj, nrep, q=.05, mult=.8, filename='all_q05_mult8.Rda')
