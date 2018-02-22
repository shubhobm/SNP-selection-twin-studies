## select optimum sd using train-test split, and do it over values of t as well
## use covariate information
rm(list=ls())
setwd("/share/bulk/gedi/majum010")
library(regress)
library(fda.usc)
library(bigstep)
library(parallel)

source('functions.R')
source('misc_functions.R')

## save data dictionary
datadir = "/share/bulk/gedi/data/GEDI/gedi5_EA_autosomes_none-missing/R/"
load("SNP_Dictionary.Rda")
gene.info.df = read.csv("geneinfo.csv")
names(gene.info.df)[1] = "genename"

## function to extract only essential elements from model
compress.model = function(mod){
  list(
    beta.hat = mod$beta,
    beta.cov = mod$beta.cov,
    W = mod$W,
    X1 = mod$X,
    H = with(mod, beta.cov %*% t(X) %*% W),
    r = with(mod, model$y - fitted))
}

# initialize quantities
sdn.vec = seq(.5,2,by=.1)
q.vec = rev(5:9)/10
t.vec = seq(.8, .1, by=-.05)
Emap = "E2"

# Emap = "E1"
# t.vec = exp(-(1:5))
set.seed(11162016)
nq = length(q.vec)
nt = length(t.vec)

## will take families with MZ and DZ twins and non-NA CON_FAC for all individuals
# which families have size 4?
pheno = read.csv("/share/bulk/gedi/mbmiller/phenotypes/20110508_McGue/GEDI_Pheno_5_8_2011.csv")
pheno = pheno[which(pheno$CON_FAC != -99),]
pheno$FAMILYID = floor(pheno$ID / 100)
pheno.MZ = pheno[which(pheno$FTYPE %in% 1:2),]
Fam.freq = as.data.frame(table(pheno.MZ$FAMILYID))
Fam.size4 = as.numeric(paste(Fam.freq[ which(Fam.freq$Freq == 4), 1]))

# function to get model elements for a single gene
do.mbic2 = function(gene){
  # isolate SNPs in the named gene
  gene.info = gene.info.df[which(gene.info.df$genename==gene),]
  gene.df =  snp.df[with(snp.df, which( chr == as.numeric(gene.info[2]) &
                                          pos_B37 >= as.numeric(gene.info[3]) &
                                          pos_B37 <= as.numeric(gene.info[4]))),]
  
  # find out 4-size families in z
  filename = unique(gene.df$file)
  z = readRDS(paste0(datadir,filename))
  gene.z = z[,which(colnames(z) %in% gene.df$snp)]
  ID.z = as.numeric(paste(rownames(gene.z)))
  FAMID.z = floor(ID.z / 100)
  Fam.freq.z = as.data.frame(table(FAMID.z))
  Fam.size4.z = as.numeric(paste(Fam.freq.z[ which(Fam.freq.z$Freq == 4), 1]))
  
  # get common family IDs in z and pheno.MZ
  common.FAMID = Fam.size4[ which( Fam.size4 %in% Fam.size4.z )]
  X.MZ4 = gene.z[which(FAMID.z %in% common.FAMID),]
  X.MZ4 = X.MZ4[order(as.numeric(rownames(X.MZ4))),]
  pheno.MZ4 = pheno.MZ[which( pheno.MZ$FAMILYID %in% common.FAMID),]
  pheno.MZ4 = pheno.MZ4[order(as.numeric(pheno.MZ4$ID)),]
  
  # put data in shape
  X = cbind(X.MZ4, pheno.MZ4[,c(4,5,7,8)]) # add covariate data to X
  y = pheno.MZ4$CON_FAC
  p = ncol(X)

  ## break into train-test
  set.seed(04052017)
  nfamily = nrow(X)/4
  train = as.numeric(sample(nfamily, ceiling(.75*nfamily), replace=F))
  ftrain = length(train)
  train = rep(4*train, rep(4,ftrain)) - rep(c(3,2,1,0), ftrain)
  ytrain = scale(y[train], scale=F)
  Xtrain = X[train,]

  # select variables
  lmod.step = selectModel(Xtrain, ytrain, p=p, verbose=F, crit="mbic2")
  indices = union(which(names(Xtrain) %in% lmod.step), p-4 + 1:4) # keep covariates always

  ## now build a model on selected indices and get prediction error
  # specify error structure
  K1 = matrix(c(1,1,.5,.5,
                1,1,.5,.5,
                .5,.5,1,0,
                .5,.5,0,1), nrow=4,ncol=4)
  K2 = matrix(c(1,.5,.5,.5,
                1,.5,.5,.5,
                .5,.5,1,0,
                .5,.5,0,1), nrow=4,ncol=4)
  Env = matrix(1, nrow=4, ncol=4)

  # construct pedigree matrix
  V1 = kronecker(diag(ftrain), K1)
  for(f in 1:ftrain){
    if((pheno.MZ4[train,])$FTYPE[4*f-3]==2){
      inds = 4*(f-1) + 1:4
      V1[inds,inds] = K2
    }
  }
  V2 = kronecker(diag(ftrain), Env)

  df = data.frame(cbind(ytrain, Xtrain[,indices]))
  names(df)[1] = "y"
  mod.select = regress(y~.-1, ~V1+V2, data=df)
  model.list = compress.model(mod.select)
  err = as.numeric(y[-train] - as.matrix(X[-train,indices]) %*% model.list$beta.hat)
  test.error = sum(diag(cov(matrix(err, nrow=nfamily-ftrain, byrow=F))))
  
  ## return
  list(names(Xtrain)[indices], test.error)
}

loopfun = function(i){
  igene = paste(gene.info.df$genename[i])
  z = do.mbic2(igene)
  cat(paste(i,igene,"Done\n"))
  z
}
model.list = mclapply(1:nrow(gene.info.df), loopfun, mc.cores=5)
save(model.list, file=paste0("outputMZDZ_mbic2.Rda"))

# function to get testing error for a value of sdn
# loopfun = function(sdn){
# 	best.ind = step11.depth(model.list, sdn=sdn, adj=adj)
# 	err = sum((scale(y[-train]) - as.matrix(X[-train, best.ind]) %*% mod$beta[best.ind])^2)
# 	list(best.ind=best.ind, err=err)
# }
# adj = "fdr"
# err.vec = mclapply(seq(.1,1,by=.05), loopfun, mc.cores=8)
# save(err.vec, "errPKNOX2.Rda")

