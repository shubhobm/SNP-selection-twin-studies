## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
# library(regress)
# library(gap)
library(RFGLS)
library(parallel)
library(locfdr)
library(doParallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 250 # number of families

## preapre data structure
FAMID = rep(1:nfamily, rep(4,nfamily))*10
INDIV = rep(1:4, nfamily)
ID = FAMID + INDIV
FTYPE = 1

# pedigree
v = c(1,1,0,0)
pedigree = data.frame(FAMID=FAMID, ID=ID,
                      PID=(FAMID+4)*v,
                      MID=(FAMID+3)*v)


## generate X data family-wise
MAF = c(0.2,0.4,0.25,0.4,0.25)
n.block = c(6,4,30,6,4)
p = sum(n.block)
n = nfamily*4
p.causal = 4

## generate two error matrices
Kmat = matrix(c(1,1,0.5,0.5,
                1,1,0.5,0.5,
                0.5,0.5,1,0,
                0.5,0.5,0,1),nrow=4,byrow=T)
Env = matrix(1, nrow=4, ncol=4)

get.outputs.others = function(signal, rho, adj, nrep, filename=NULL){
  set.seed(11162016)
  
  loopfun = function(rep){
    library(RFGLS)
    library(locfdr)
    source('misc_functions.R')
    source('simgen.r')
    
    nfamily = 250 # number of families
    
    ## preapre data structure
    FAMID = rep(1:nfamily, rep(4,nfamily))*10
    INDIV = rep(1:4, nfamily)
    ID = FAMID + INDIV
    FTYPE = 1
    
    # pedigree
    v = c(1,1,0,0)
    pedigree = data.frame(FAMID=FAMID, ID=ID,
                          PID=(FAMID+4)*v,
                          MID=(FAMID+3)*v)
    
    
    ## generate X data family-wise
    MAF = c(0.2,0.4,0.25,0.4,0.25)
    n.block = c(6,4,30,6,4)
    p = sum(n.block)
    n = nfamily*4
    p.causal = 4
    
    ## generate two error matrices
    Kmat = matrix(c(1,1,0.5,0.5,
                    1,1,0.5,0.5,
                    0.5,0.5,1,0,
                    0.5,0.5,0,1),nrow=4,byrow=T)
    Env = matrix(1, nrow=4, ncol=4)
    
    
    out.mat = matrix(0, nrow=4, ncol=2)
    
    ## generate y
    set.seed(rep)
    X = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    X = data.frame(X)
    y = simphe(gendat=X, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=4, sigmae2=2, type="MZ", r=.5)
    beta = y$beta
    y = y$Y[,1]; y = y - mean(y)
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))

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
    
    # do FDR correction and find best index
    fdr = locfdr(mod.gls$t.stat, plot=0)$fdr
    best.RF = which(fdr<1)
    
    active0 = c(1,7,41,47)
    active1 = c(1:10, 41:50)
    TP0 = sum(active0 %in% best.RF)
    TN0 = 46 - sum((1:p)[best.RF] %in% (1:p)[-active0])
    TP1 = sum(c(sum(1:6 %in% best.RF)>0,
                sum(7:10 %in% best.RF)>0,
                sum(41:46 %in% best.RF)>0,
                sum(47:50 %in% best.RF)>0))
    TN1 = 30 - sum((1:p)[best.RF] %in% (1:p)[-active1])
    out.list = c(TP0/4, TN0/46, TP1/4, TN1/30)
    
    # return output
    out.list
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  all.list = foreach(i=1:nrep) %dopar% loopfun(i)
  # all.list = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
  if(is.null(filename)){
    matrix(unlist(all.list), nrow=nrep, ncol=4, byrow=T)
  } else{
    save(all.list, file=filename)
  }
}


##### Output function

nrep = 1e2
active0 = c(1,7,41,47)
active1 = c(1:10, 41:50)

all = get.outputs.others(signal=.1, rho=.7, adj, nrep=20)
apply(all,2,mean)
# all = get.outputs.others(signal=.1, rho=.7, adj, nrep=20, filename="test.Rda")
# get.outputs.others(signal=.05, rho=.7, adj, nrep, filename="others1_h05_rho07.Rda")
# get.outputs.others(signal=.02, rho=.7, adj, nrep, filename="others1_h02_rho07.Rda")
# get.outputs.others(signal=.01, rho=.7, adj, nrep, filename="others1_h01_rho07.Rda")
# get.outputs.others(signal=0, rho=.7, adj, nrep, filename="others1_h00_rho07.Rda")

# get.outputs.others(signal=.1, rho=.1, adj, nrep, filename="others1_h10_rho01.Rda")
# get.outputs.others(signal=.05, rho=.1, adj, nrep, filename="others1_h05_rho01.Rda")
# get.outputs.others(signal=.02, rho=.1, adj, nrep, filename="others1_h02_rho01.Rda")
# get.outputs.others(signal=.01, rho=.1, adj, nrep, filename="others1_h01_rho01.Rda")
# get.outputs.others(signal=0, rho=.1, adj, nrep, filename="others1_h00_rho01.Rda")


