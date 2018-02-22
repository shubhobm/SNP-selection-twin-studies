## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
rm(list=ls())
# setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')
library(regress)
library(gap)
library(fda.usc)
library(parallel)

source('misc_functions.R')
source('simgen.r')

nfamily = 250 # number of families

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

## looping for multiple instances of the analysis
get.outputs.depth = function(signal, rho,
                             sdn.vec, q.vec, t.vec,
                             nrep, filename){
  set.seed(11162016)
  nq = length(q.vec)
  nt = length(t.vec)
  
  loopfun = function(rep){
    
    ## generate y
    set.seed(1e3*rep)
    X = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    X = data.frame(X)
    y = simphe.r(gendat=X, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=4, sigmae2=2, type="MZ", r=.5)
    beta = y$beta
    y = y$Y[,1]; y = y - mean(y)
    
    # generate test data
    Xtest = simgen(LD=rep(rho, 5), MAF=MAF,
               n.block=n.block, n.person=nfamily)
    ytest = simphe.r(gendat=Xtest, h2=c(signal,signal,0,signal,signal),
               MAF=MAF, n.block=n.block, nfam=nfamily,
               sigma2=4, sigmae2=2, type="MZ", r=.5)
    ytest = ytest$Y[,1]; ytest = ytest - mean(ytest)
    
    datamat = data.frame(cbind(y,X))
    form = as.formula(paste0("y~", paste(names(X), collapse="+"), "-1"))

    #################################
    # best index for all-variable mixed model
    #################################
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
    TP0.list = list() # detect actual causal SNPs
    TN0.list = list()
    TP1.list = list() # detect any of correlated SNPs
    TN1.list = list()
    err.list = list()
    
    for(isdn in 1:nsdn){
      qt.list = step.depth.qt(mod, sdn=sdn.vec[isdn], q.vec, t.vec)
      TP0.mat = matrix(0, nrow=nq, ncol=nt)
      TN0.mat = TP0.mat
      TP1.mat = TP0.mat
      TN1.mat = TP0.mat
      err.mat = TP0.mat
      
      # deal with index sets one by one
      for(i in 1:nq){
        for(j in 1:nt){
          best.index = qt.list[[i]][[j]]
          
          # get relevant quantities
          TP0.mat[i,j] = sum(active0 %in% best.index)
          TN0.mat[i,j] = 46 - sum((1:p)[best.index] %in% (1:p)[-active0])
          TP1.mat[i,j] = sum(c(sum(1:6 %in% best.index)>0,
                               sum(7:10 %in% best.index)>0,
                               sum(41:46 %in% best.index)>0,
                               sum(47:50 %in% best.index)>0))
          TN1.mat[i,j] = 30 - sum((1:p)[best.index] %in% (1:p)[-active1])
          
          test.error = sum(diag(cov(matrix(ytest, nrow=nfamily, byrow=F))))
          if(length(best.index)>0){
            err = as.numeric(ytest - as.matrix(Xtest[,best.index]) %*% mod$beta.hat[best.index])
            test.error = sum(diag(cov(matrix(err, nrow=nfamily, byrow=F))))
            # test.error = sum(diag(matrix(err, nrow=nfamily, byrow=F)^2))
          }
          err.mat[i,j] = test.error
        }
      }
      
      # now store
      TP0.list[[isdn]] = TP0.mat
      TN0.list[[isdn]] = TN0.mat
      TP1.list[[isdn]] = TP1.mat
      TN1.list[[isdn]] = TN1.mat
      err.list[[isdn]] = err.mat
    }
    
    # now select best sdn for every (q,t)
    best.TP0 = matrix(0, nrow=nq, ncol=nt)
    best.TN0 = best.TP0
    best.TP1 = best.TP0
    best.TN1 = best.TP0
    
    
    for(i in 1:nq){
      for(j in 1:nt){
        err.vec = lapply(err.list, function(x) x[i,j])
        best.sdn = which.min(err.vec)
        best.TP0[i,j] = (TP0.list[[best.sdn]])[i,j]/4
        best.TN0[i,j] = (TN0.list[[best.sdn]])[i,j]/46
        best.TP1[i,j] = (TP1.list[[best.sdn]])[i,j]/4
        best.TN1[i,j] = (TN1.list[[best.sdn]])[i,j]/30
      }
    }
    
    # return
    rbind(best.TP0, best.TN0, best.TP1, best.TN1)
  }
  
  # all.list = lapply(1:nrep, loopfun)
  all.list = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
  save(all.list, file=filename)
}

##### Output function

nrep = 1e3
active0 = c(1,7,41,47)
active1 = c(1:10, 41:50)
sdn.vec = seq(.3,2,by=.05)
h.vec = c(.1, .05, .02, .01, 0)
q.vec = c(.95, .9, .5, .2, .1)
t.vec = rev(1:9)/10
get.outputs.depth(signal=.1, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="qt_h10_rho7_big.Rda")
get.outputs.depth(signal=.05, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="qt_h05_rho7_big.Rda")
get.outputs.depth(signal=.02, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="qt_h02_rho7_big.Rda")
get.outputs.depth(signal=.01, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="qt_h01_rho7_big.Rda")
get.outputs.depth(signal=0, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="qt_h00_rho7_big.Rda")

