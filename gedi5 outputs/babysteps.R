rm(list=ls())
setwd("/share/bulk/gedi/majum010")
library(regress)
library(fda.usc)
#library(caret)
library(parallel)
source('functions.R')

snp.df = readRDS(paste0(datadir,"map_G1.rds"))
snp.df$file = "G1.rds"
for(i in 2:100){
	ifile = paste0("map_G",i,".rds")
	ipath = paste0(datadir, ifile)
	i.df = readRDS(ipath)
	i.df$file = paste0("G",i,".rds")
	snp.df = rbind(snp.df, i.df)
}
save(snp.df, file="SNP_Dictionary.Rda")

pheno = read.csv("/share/bulk/gedi/mbmiller/phenotypes/20110508_McGue/GEDI_Pheno_5_8_2011.csv")
z = readRDS("/share/bulk/gedi/data/GEDI/gedi5_EA_autosomes_none-missing/R/G1.rds")

## will take families with MZ twins and non-NA CON_FAC for all individuals
# which families have size 4?
pheno = pheno[which(pheno$CON_FAC != -99),]
pheno$FAMILYID = floor(pheno$ID / 100)
pheno.MZ = pheno[which(pheno$FTYPE == 1),]
Fam.freq = as.data.frame(table(pheno.MZ$FAMILYID))
Fam.size4 = as.numeric(paste(Fam.freq[ which(Fam.freq$Freq == 4), 1]))

# find out 4-size families in z
ID.z = as.numeric(paste(rownames(z)))
FAMID.z = floor(ID.z / 100)
Fam.freq.z = as.data.frame(table(FAMID.z))
Fam.size4.z = as.numeric(paste(Fam.freq.z[ which(Fam.freq.z$Freq == 4), 1]))

# get common family IDs in z and pheno
common.FAMID = Fam.size4[ which( Fam.size4 %in% Fam.size4.z )]
X.MZ4 = z[which(FAMID.z %in% common.FAMID),]
pheno.MZ4 = pheno.MZ[which( pheno.MZ$FAMILYID %in% common.FAMID),]

# now model
y = pheno.MZ4$CON_FAC; y = scale(y)
X = X.MZ4[,1:10]
df = data.frame(cbind(y,X))
nfamily = nrow(df)/4
Kmat = matrix(c(1,.5,.5,.5,
                .5,1,.5,.5,
                .5,.5,1,0,
                .5,.5,0,1), nrow=4,ncol=4)
Kron = kronecker(diag(nfamily), Kmat)
#mod = regress(y~.-1, ~Kron, data=df)

## 5-fold cv
k = 8
folds = createFolds(1:nfamily, k=k, list=F)
loopfun = function(foldnum){
	train = which(folds != foldnum)
	ftrain = length(train)
    train = rep(4*train, rep(4,ftrain)) - rep(c(3,2,1,0), ftrain)
    Kron = kronecker(diag(ftrain), Kmat)
    ntrain = length(train)
    mod = regress(y~.-1, ~Kron, data=df[train,])
	mod
}
system.time(model.list <- mclapply(1:k, loopfun, mc.cores=min(k,8)))

loopfun1 = function(i){
	# train model for the i-th fold
	train = which(folds != i)
	ftrain = length(train)
    train = rep(4*train, rep(4,ftrain)) - rep(c(3,2,1,0), ftrain)
    Kron = kronecker(diag(ftrain), Kmat)
    ntrain = length(train)
    mod.i = regress(y~.-1, ~Kron, data=df[train,])

# now loop over all sd values
	err.fold = rep(0, 10)
	for(j in 1:10){
		sdn.j = .1 + .5*j
		select.ind = step1.depth(mod.i, sdn.j, adj)
		err.fold[j] = sum(y[-train]^2)
		if(length(select.ind) > 0){
			err.fold[j] = sum((y[-train] - as.matrix(X[-train, select.ind]) %*% mod.i$beta[select.ind])^2)
		}
	}
	err.fold
}

cv.err = matrix(unlist(mclapply(1:k, loopfun1, mc.cores=min(k,8))), nrow=k, byrow=F)

# final model
best.sd = .1 + .5*which.min(colSums(cv.err))
mod = regress(y~.-1, ~Kron, data=df)
select.ind = step1.depth(mod, best.sd, adj)
select.ind