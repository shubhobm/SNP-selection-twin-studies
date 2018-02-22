setwd('C:/Study/UMN files/fellowship/IDF 2015')
library(regress)

X = read.csv('Genotype_500fam.txt')
ymat = read.csv('5phenotype_500fam.txt')
y = ymat[,1]
set.seed(10202015)

Kmat = matrix(c(1,0,.5,.5,
                0,1,.5,.5,
                .5,.5,1,.5,
                .5,.5,.5,1),
              nrow=4,ncol=4)
Kron = kronecker(diag(400), Kmat)

datamat = data.frame(cbind(y,X))
train = sample(1:500, 400, replace=F)
train = rep(4*train, rep(4,400)) - rep(c(3,2,1,0), 400)

form = as.formula(paste0("y~", paste(names(X), collapse="+")))

system.time(lmod <- regress(form, data=datamat))
#with(lmod, View(cbind(y,fitted,predicted)))

library(gap)
library(fda.usc)
k2 <- kin.morgan(l51)$kin.matrix*2
k2[1:10,1:10]

p = ncol(X)
n = nrow(X)
sdn.vec = n^seq(.05, .25, by=.01)
nboot = 1e3
  
## Calculate fixed quantities
print("Calculate Hessians and score matrices")
Const.H = list()
Const.beta = list()
Const.g = list()
Const.r = list()
for(j in 1:p){
  jform = as.formula(paste0("y~", paste(names(X[,-j]), collapse="+")))
  jmod = regress(jform, ~Kron, data=datamat[train,])
  Const.beta[[j]] = with(jmod, beta)
  Const.H[[j]] = with(jmod, beta.cov %*% t(X) %*% W)
  Const.g[[j]] = BLUP(jmod)$Mean
  Const.r[[j]] = with(jmod, y[train] - fitted - Const.g[[j]])
  print(paste("Variable", j, "done."))
}

## full models
form = as.formula(paste0("y~", paste(names(X), collapse="+")))
mod = regress(form, ~Kron, data=datamat[train,])
beta = with(mod, beta)
H = with(mod, beta.cov %*% t(X) %*% W)
g = BLUP(mod)$Mean
r = with(mod, y[train] - fitted - g)

RMSE = rep(0, length(sdn.vec))
Cn.frame.list = list()
isd = 1
for(sdn in sdn.vec){
  
  SSPmat.d = matrix(0, nrow=nboot, ncol=p+1)
  
  # depth model selection
  beta.mat = matrix(0 ,nrow=nboot, ncol=p+1)
  for(i in 1:nboot){
    iresid = as.matrix(sdn * (rep(rnorm(400), rep(4,400))*g + rnorm(1600)*r), ncol=1)
    beta.mat[i,] = as.numeric(beta) + as.numeric(H %*% iresid)
  }
  SSPmat.d[,p+1] = mdepth.TD(beta.mat[,-1], beta.mat[,-1])$dep
  
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
      iresid = sdn * (rep(rnorm(400), rep(4,400))* Const.g[[j]] + rnorm(1600)* Const.r[[j]])
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
  Cn.frame.list[[isd]] = Cn.frame
  
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
  RMSE[isd] = sum((y[-train] - as.matrix(cbind(1,X[-train,which.final])) %*% mod.final$beta)^2)
  print(paste("  Iteration", isd, "done."))
  isd = isd+1
}

RMSE = sqrt(RMSE/400)
Cn.frame.list[[which.min(RMSE)]]

# best linear model prediction
lmod.full = lm(y~., data=datamat[train,])
lmod.best = step(lmod.full)
RMSE.lin = sqrt(mean((y[-train] - predict(lmod.best, newdata=X[-train,]))^2))
RMSE.lin.full = sqrt(mean((y[-train] - predict(lmod.full, newdata=X[-train,]))^2))
RMSE.full = sqrt(mean((y[-train] - as.matrix(cbind(1,X[-train,])) %*% mod$beta)^2))

pdf('')
plot(RMSE~log(sdn.vec,n), type='b', lwd=2, xlab="Tuning parameter for D1",
     ylim=c(min(RMSE)-.01, max(RMSE)+.01))
abline(h=RMSE.full, lty=2, lwd=2)
abline(h=RMSE.lin, lwd=2, col="red")
abline(h=RMSE.lin.full, lty=2, lwd=2, col="red")
legend('bottomleft',
       c("Full LM","Reduced LM","Full LMM","Reduced LMM"),
       col=c("red","red","black","black"),
       lty=c(2,1,2,1),
       lwd=2)
save(Cn.frame.list,file="SNP5.Rda")
