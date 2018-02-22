

setwd('C:/Study/UMN files/fellowship/IDF 2015')
library(regress)

X = read.csv('Genotype_500fam.txt')
ymat = read.csv('5phenotype_500fam.txt')
y = ymat[,1]

Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1),
              nrow=4,ncol=4)
Kron = kronecker(diag(500), Kmat)

datamat = data.frame(cbind(y,X))
form = as.formula(paste0("y~", paste(names(X), collapse="+")))

system.time(lmod <- regress(form, data=datamat))
with(lmod, View(cbind(y,fitted,predicted)))

system.time(mixedmod <- regress(form, ~Kron, data=datamat))

library(gap)
library(fda.usc)
k2 <- kin.morgan(l51)$kin.matrix*2
k2[1:10,1:10]

p = ncol(X)
n = nrow(X)
sdn = n^.1
nboot = 1e3
  
# depth model selection

## full models
SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)

form = as.formula(paste0("y~", paste(names(X), collapse="+")))
mod = regress(form, ~Kron, data=datamat)
beta = with(mod, beta)
H = with(mod, beta.cov %*% t(X) %*% W)
g = BLUP(mod)$Mean
r = with(mod, y - fitted - g)

beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
for(i in 1:nboot){
  iresid = as.matrix(sdn * (rep(rnorm(n/4), rep(4,n/4))*g + rnorm(n)*r), ncol=1)
  beta.mat[i,] = as.numeric(beta) + as.numeric(H %*% iresid)
}
SSPmat.d[,p+1] = mdepth.TD(beta.mat[,-1], beta.mat[,-1])$dep

print("Calculate Hessians and score matrices")
Const.H = list()
Const.beta = list()
Const.g = list()
Const.r = list()
for(j in 1:p){
  jform = as.formula(paste0("y~", paste(names(X[,-j]), collapse="+")))
  jmod = regress(jform, ~Kron, data=datamat)
  Const.beta[[j]] = with(jmod, beta)
  Const.H[[j]] = with(jmod, beta.cov %*% t(X) %*% W)
  Const.g[[j]] = BLUP(jmod)$Mean
  Const.r[[j]] = with(jmod, y - fitted - Const.g[[j]])
  print(paste("Variable", j, "done."))
}
  
#   # matrix of full model estimates: Fn is approx by Fn^b1
#   beta.mat = matrix(0, nrow=nboot, ncol=p)
#   for(i in 1:nboot){
#     beta.mat[i,] = beta.hat + P %*% (sdn[nsd]*rnorm(n)*r)
#   }
  
# get truncated estimates: beta_alpha approx by beta_alpha^b
print("Do Monte-Carlo simulation")
for(j in 1:p){
  
  jbeta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
  for(i in 1:nboot){
    iresid = sdn * (rep(rnorm(n/4), rep(4,n/4))* Const.g[[j]] + rnorm(n)* Const.r[[j]])
    jbeta.mat[i,-(j+1)] = Const.beta[[j]] + Const.H[[j]] %*% iresid
  }
  
  SSPmat.d[,j] = mdepth.TD(jbeta.mat[,-1], beta.mat[,-1])$dep
  print(paste("Variable", j, "done."))
}

pVal = rep(1, p+1)
for(i in 1:p){
  pVal[i] = t.test(SSPmat.d[,i], SSPmat.d[,p+1], paired=TRUE)$p.value
}

Cn.frame = data.frame(DroppedVar = c(paste("-", names(X)), "<none>"),
                      Cn = apply(SSPmat.d, 2, mean),
                      pValue = pVal)
Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
row.names(Cn.frame) = NULL
Cn.frame
