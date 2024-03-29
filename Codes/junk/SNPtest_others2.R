## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
# setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(RFGLS)
library(fda.usc)
library(parallel)
library(bigstep)

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
n.block = c(6,4,80,6,4)
p = sum(n.block)
n = nfamily*4
p.causal = 4

## generate two error matrices
Kmat = matrix(c(1,1,0.5,0.5,
                1,1,0.5,0.5,
                0.5,0.5,1,0,
                0.5,0.5,0,1),nrow=4,byrow=T)
Env = matrix(1, nrow=4, ncol=4)

get.outputs.others = function(signal, rho, adj, nrep, filename){
  set.seed(11162016)
  
  loopfun = function(rep){
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
    # best index for linear model mBIC2
    #################################
    parent.ind = c(4*(1:nfamily)-3, 4*(1:nfamily)-2)
    lmod.step = selectModel(X, y, p=p, verbose=F, crit="mbic2")
    best.BIC = which(names(datamat)[-1] %in% lmod.step)
    TP0 = sum(active0 %in% best.BIC)
    TN0 = 96 - sum((1:p)[best.BIC] %in% (1:p)[-active0])
    TP1 = sum(c(sum(1:6 %in% best.BIC)>0,
                         sum(7:10 %in% best.BIC)>0,
                         sum(91:96 %in% best.BIC)>0,
                         sum(97:100 %in% best.BIC)>0))
    TN1 = 80 - sum((1:p)[best.BIC] %in% (1:p)[-active1])
    out.mat[,1] = c(TP0/4, TN0/96, TP1/4, TN1/80)

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
    
    # do for FDR correction and find best index
    alpha.vec = .05*(1:p)/p
    pval = mod.gls$pval
    pval.order = order(pval, decreasing = F)
    which.less = which(pval[pval.order] < alpha.vec)
    best.RF = integer(0)
    if(length(which.less)>0){
      best.RF = pval.order[1:max(which.less)]
    }
    # best.RF = which(pval < 0.05)
    
    # pval.adj = p.adjust(mod.gls$pval, method=adj)
    # best.RF = which(pval.adj<.05)
    
    TP0 = sum(active0 %in% best.RF)
    TN0 = 96 - sum((1:p)[best.RF] %in% (1:p)[-active0])
    TP1 = sum(c(sum(1:6 %in% best.RF)>0,
                sum(7:10 %in% best.RF)>0,
                sum(91:96 %in% best.RF)>0,
                sum(97:100 %in% best.RF)>0))
    TN1 = 80 - sum((1:p)[best.RF] %in% (1:p)[-active1])
    out.mat[,2] = c(TP0/4, TN0/96, TP1/4, TN1/80)
    
    # return output
    out.mat
  }
  
  # TPTN.mat1 = lapply(1:nrep, loopfun)
  all.list = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
  save(all.list, file=filename)
}


##### Output function

nrep = 1e2
active0 = c(1,7,91,97)
active1 = c(1:10, 91:100)
adj = "fdr"

# get.outputs.others(signal=.1, rho=.7, adj, nrep, filename="others1_h10_rho07.Rda")
# get.outputs.others(signal=.05, rho=.7, adj, nrep, filename="others1_h05_rho07.Rda")
# get.outputs.others(signal=.02, rho=.7, adj, nrep, filename="others1_h02_rho07.Rda")
# get.outputs.others(signal=.01, rho=.7, adj, nrep, filename="others1_h01_rho07.Rda")
# get.outputs.others(signal=0, rho=.7, adj, nrep, filename="others1_h00_rho07.Rda")

# get.outputs.others(signal=.1, rho=.1, adj, nrep, filename="others1_h10_rho01.Rda")
# get.outputs.others(signal=.05, rho=.1, adj, nrep, filename="others1_h05_rho01.Rda")
# get.outputs.others(signal=.02, rho=.1, adj, nrep, filename="others1_h02_rho01.Rda")
# get.outputs.others(signal=.01, rho=.1, adj, nrep, filename="others1_h01_rho01.Rda")
# get.outputs.others(signal=0, rho=.1, adj, nrep, filename="others1_h00_rho01.Rda")

nrep = 1e3
# get.outputs.others(signal=.1, rho=.7, adj, nrep, filename="others2_h10_rho07_big.Rda")
get.outputs.others(signal=.07, rho=.7, adj, nrep, filename="others2_h07_rho07_big.Rda")
# get.outputs.others(signal=.05, rho=.7, adj, nrep, filename="others2_h05_rho07_big.Rda")
get.outputs.others(signal=.03, rho=.7, adj, nrep, filename="others2_h03_rho07_big.Rda")
# get.outputs.others(signal=.02, rho=.1, adj, nrep, filename="others2_h02_rho01_big.Rda")
# get.outputs.others(signal=.01, rho=.7, adj, nrep, filename="others2_h01_rho07_big.Rda")
# get.outputs.others(signal=0, rho=.7, adj, nrep, filename="others2_h00_rho07_big.Rda")

# get.outputs.others(signal=.1, rho=.1, adj, nrep, filename="others2_h10_rho01_big.Rda")
get.outputs.others(signal=.07, rho=.1, adj, nrep, filename="others2_h07_rho01_big.Rda")
# get.outputs.others(signal=.05, rho=.1, adj, nrep, filename="others2_h05_rho01_big.Rda")
get.outputs.others(signal=.03, rho=.1, adj, nrep, filename="others2_h03_rho01_big.Rda")
# get.outputs.others(signal=.02, rho=.1, adj, nrep, filename="others2_h02_rho01_big.Rda")
# get.outputs.others(signal=.01, rho=.1, adj, nrep, filename="others2_h01_rho01_big.Rda")
# get.outputs.others(signal=0, rho=.1, adj, nrep, filename="others2_h00_rho01_big.Rda")



