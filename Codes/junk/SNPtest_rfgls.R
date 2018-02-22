## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
# setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(parallel)
library(RFGLS)

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

## looping for multiple instances of the analysis
get.outputs = function(signal, rho, adj, nrep){
  set.seed(11162016)
  
  loopfun = function(rep){
    
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
    
    # make data
    pheno = data.frame(cbind(FAMID,ID,FTYPE,INDIV,y))
    geno = data.frame(t(X))
    colnames(geno) = ID
    
    mod.gls = gls.batch(
      phenfile=pheno,genfile=data.frame(t(geno)),pedifile=pedigree,
      theta=NULL,snp.names=row.names(geno),
      input.mode=c(1,2,3),pediheader=FALSE, 
      pedicolname=c("FAMID","ID","PID","MID"),
      sep.phe=" ",sep.gen=" ",sep.ped=" ",
      phen="y",covars=NULL,med=c("UN","VC"),
      outfile=NULL,col.names=TRUE,return.value=TRUE,
      covmtxfile.out=NULL,covmtxparams.out=NULL,
      sizeLab=NULL,Mz=NULL,Bo=NULL,Ad=NULL,Mix=NULL,indobs=NULL)
    
    # adjust p-values
    pval.adj = p.adjust(mod.gls$pval, method=adj)
    
    # best index for linear model BIC
    active.ind = c(1,7,41,47)
    best.index = which(pval.adj<.05)
    TP.num = c((1 %in% best.index), (7 %in% best.index),
               (41 %in% best.index), (47 %in% best.index))
    TN.num = sum((1:p)[-best.index] %in% (1:p)[-active.ind])
    c(TP.num, sum(TP.num), TN.num, length(best.index))
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  TPTN.mat1 = mclapply(1:nrep, loopfun, mc.cores=min(16,nrep))
  TPTN.mat = matrix(unlist(TPTN.mat1), ncol=7, byrow=T)
  z = c(apply(TPTN.mat,2,mean)/c(rep(1,4),4,46,1),
        apply(TPTN.mat,2,sd)/c(rep(1,4),4,46,1))
  matrix(z, nrow=2, byrow=T)
}

##### Output function
get.outputs(signal=.5, rho=.1, adj="fdr", nrep=1e2)
get.outputs(signal=.05, rho=.1, adj="fdr", nrep=1e2)
get.outputs(signal=.01, rho=.1, adj="fdr", nrep=1e2)

get.outputs(signal=5., rho=.7, adj="fdr", nrep=1e2)
get.outputs(signal=.05, rho=.7, adj="fdr", nrep=1e2)
get.outputs(signal=.01, rho=.7, adj="fdr", nrep=1e2)