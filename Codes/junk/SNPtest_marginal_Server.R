rm(list=ls())
# setwd('C:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)

X = read.csv('Genotype_500fam.txt')
ymat = read.csv('5phenotype_500fam.txt')
y = ymat[,3]

Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1),
              nrow=4,ncol=4)
Kron = kronecker(diag(500), Kmat)

datamat = data.frame(cbind(y,X))
form = as.formula(paste0("y~", paste(names(X), collapse="+")))

# system.time(lmod <- regress(form, data=datamat))
# with(lmod, View(cbind(y,fitted,predicted)))
# 
# system.time(mixedmod <- regress(form, ~Kron, data=datamat))

k2 <- kin.morgan(l51)$kin.matrix*2
k2[1:10,1:10]

p = ncol(X)
n = nrow(X)
sdn.vec = n^seq(.05, .25, by=.01)
nboot = 1e3
  
# depth model selection

## full models
SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)

# train full model
set.seed(10252016)
train = sample(1:500, 400, replace=F)
train = rep(4*train, rep(4,400)) - rep(c(3,2,1,0), 400)
ntrain = length(train)
Kron = kronecker(diag(400), Kmat)
form = as.formula(paste0("y~", paste(names(X), collapse="+")))
mod = regress(form, ~Kron, data=datamat[train,])

# full model quantities
beta = with(mod, beta)
H = with(mod, beta.cov %*% t(X) %*% W)
# g = BLUP(mod)$Mean
# r = with(mod, y - fitted - g)
r = with(mod, y[train] - fitted)

# truncated model quantities
Const.H = list()
Const.beta = list()
Const.g = list()
Const.r = list()
for(j in 1:p){
  # jform = as.formula(paste0("y~", paste(names(X[,-j]), collapse="+")))
  # jmod = regress(jform, ~Kron, data=datamat)
  Const.beta[[j]] = beta[-j]
  Const.H[[j]] = with(mod, beta.cov[-j,-j] %*% t(X[,-j]) %*% W)
  # Const.g[[j]] = BLUP(jmod)$Mean
  # Const.r[[j]] = with(mod, y - X[,-j] %*% beta[-j] - Const.g[[j]])
  Const.r[[j]] = with(mod, y[train] - X[,-j] %*% beta[-j])
  # print(paste("Variable", j, "done."))
}
  
#   # matrix of full model estimates: Fn is approx by Fn^b1
# get truncated estimates: beta_alpha approx by beta_alpha^b
RMSE = rep(0, length(sdn.vec))
Cn.frame.list = list()
loopfun = function(sdn){
  ## full model bootstrap
  beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
  for(i in 1:nboot){
    # iresid = as.matrix(sdn * (rep(rnorm(n/4), rep(4,n/4))*g + rnorm(n)*r), ncol=1)
    iresid = sdn * as.matrix(rnorm(ntrain)*r, ncol=1)
    beta.mat[i,] = as.numeric(beta) + as.numeric(H %*% iresid)
  }
  SSPmat.d[,p+1] = mdepth.RP(beta.mat[,-1], beta.mat[,-1])$dep
  
  ## Truncated models bootstrap
  for(j in 1:p){
    
    jbeta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
    for(i in 1:nboot){
      # iresid = sdn * (rep(rnorm(n/4), rep(4,n/4))* Const.g[[j]] + rnorm(n)* Const.r[[j]])
      iresid = sdn * rnorm(ntrain)* Const.r[[j]]
      jbeta.mat[i,-(j+1)] = Const.beta[[j]] + Const.H[[j]] %*% iresid
    }
    
    SSPmat.d[,j] = mdepth.RP(jbeta.mat[,-1], beta.mat[,-1])$dep
    print(paste("Variable", j, "done."))
  }
  
  Cn.frame = data.frame(DroppedVar = c(paste("-", names(X)), "<none>"),
                        Cn = apply(SSPmat.d, 2, mean))
  Cn.frame = Cn.frame[with(Cn.frame, order(Cn)),]
  row.names(Cn.frame) = NULL

  # build final model
  noneCn = Cn.frame$Cn[which(Cn.frame$DroppedVar == "<none>")]
  which.final = which(apply(SSPmat.d, 2, mean) < noneCn)
  names.final = "1"
  if(length(which.final)>0){
    names.final = names(X)[which.final]
  }
  
  form.final = as.formula(paste("y ~", paste(names.final, collapse="+")))
  mod.final = regress(form.final, ~Kron, data=datamat[train,])
  
  # make predictions
  sum((y[-train] - as.matrix(cbind(1,X[-train,which.final])) %*% mod.final$beta)^2)
}


RMSE = as.numeric(mclapply(sdn, loopfun, mc.cores=12))
(RMSE <- sqrt(RMSE/400))


# plot errors
lmod.full = lm(y~., data=datamat[train,])
lmod.best = step(lmod.full)
RMSE.lin = sqrt(mean((y[-train] - predict(lmod.best, newdata=X[-train,]))^2))
RMSE.lin.full = sqrt(mean((y[-train] - predict(lmod.full, newdata=X[-train,]))^2))
RMSE.full = sqrt(mean((y[-train] - as.matrix(cbind(1,X[-train,])) %*% mod$beta)^2))

pdf('errorplot.pdf', width=6, height=6)
plot(RMSE~log(sdn.vec,n), type='b', lwd=2, xlab="Tuning parameter for D1",
     ylim=c(min(RMSE)-.01, max(RMSE)+.01))
abline(h=RMSE.full, lty=2, lwd=2)
abline(h=RMSE.lin, lwd=2, col="red")
abline(h=RMSE.lin.full, lty=2, lwd=2, col="red")
# legend('bottomleft',
#        c("Full LM","Reduced LM","Full LMM","Reduced LMM"),
#        col=c("red","red","black","black"),
#        lty=c(2,1,2,1),
#        lwd=2)
# save(Cn.frame.list,file="SNP5.Rda")
