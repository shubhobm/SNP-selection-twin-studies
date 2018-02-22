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

nfamily = 250 # number of qqfamilies

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
                             Emap, sdn.vec, q.vec, t.vec,
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
    result.list0 = list()
    result.list1 = list() # list to store results
    # lists to store results
    best.result = matrix(0, nrow=nt, ncol=8)

    for(isdn in 1:nsdn){
      # qt.list = step.nodepth.qtgamma.fdr(mod, sdn=sdn.vec[isdn], Emap=Emap,
      #                                q.vec, t.vec, nboot=1e3)
      qt.list.both = step.nodepth.qtgamma.both(mod, sdn=sdn.vec[isdn], Emap=Emap,
                                         q.vec, t.vec, nboot=1e3)
      
      ###############################################
      ### first element of output: BH corrected list
      qt.list0 = qt.list.both[[1]]
      qt.list1 = qt.list.both[[2]]
      result.mat0 = matrix(0, nrow=5, ncol=nt)
      result.mat1 = result.mat0
      
      process = function(qt.list, result.mat){
        # deal with index sets one by one
        for(j in 1:nt){
          
          ## best index set is the indices that get selected in 4 of 5 quantiles checked
          out.mat = matrix(0, nrow=nq, ncol=p)
          for(i in 1:nq){
            out.mat[i, qt.list[[i]][[j]] ] = 1
          }
          best.index = which(colSums(out.mat) == nq)
          
          # get relevant quantities
          result.mat[1,j] = sum(active0 %in% best.index)
          result.mat[2,j] = 46 - sum((1:p)[best.index] %in% (1:p)[-active0])
          result.mat[3,j] = sum(c(sum(1:6 %in% best.index)>0,
                                  sum(7:10 %in% best.index)>0,
                                  sum(41:46 %in% best.index)>0,
                                  sum(47:50 %in% best.index)>0))
          result.mat[4,j] = 30 - sum((1:p)[best.index] %in% (1:p)[-active1])
          
          test.error = sum(diag(cov(matrix(ytest, nrow=nfamily, byrow=F))))
          if(length(best.index)>0){
            err = as.numeric(ytest - as.matrix(Xtest[,best.index]) %*% mod$beta.hat[best.index])
            test.error = sum(diag(cov(matrix(err, nrow=nfamily, byrow=F))))
            # test.error = sum(diag(matrix(err, nrow=nfamily, byrow=F)^2))
          }
          result.mat[5,j] = test.error
        }
        result.mat
      }

      # now store
      result.list0[[isdn]] = process(qt.list0, result.mat0)
      result.list1[[isdn]] = process(qt.list1, result.mat1)
    }
    
    # now select best sdn for every t
    for(j in 1:nt){
      err.vec = lapply(result.list0, function(x) x[5,j])
      best.sdn = which.min(err.vec)
      best.result[j,1:4] = (result.list0[[best.sdn]])[-5,j]/c(4,46,4,30)

      err.vec = lapply(result.list1, function(x) x[5,j])
      best.sdn = which.min(err.vec)
      best.result[j,5:8] = (result.list1[[best.sdn]])[-5,j]/c(4,46,4,30)
    }
    
    best.result
  }
  
  # all.list = lapply(1:nrep, loopfun)
  all.list = mclapply(1:nrep, loopfun, mc.cores=min(8,nrep))
  save(all.list, file=filename)
}

##### Output function

active0 = c(1,7,41,47)
active1 = c(1:10, 41:50)
sdn.vec = seq(.3,1,by=.05)
h.vec = c(.1, .05, .02, .01, 0)
q.vec = rev(5:9)/10
t.vec = exp(-(1:5))

Emap = "E2"
nrep = 1e2
t.vec = seq(.9, .5, by=-.025)
get.outputs.depth(signal=.07, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h07_rho7_both.Rda")
# get.outputs.depth(signal=.05, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h05_rho7_big.Rda")
# get.outputs.depth(signal=.03, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h03_rho7_big.Rda")
# get.outputs.depth(signal=.02, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h02_rho7_big.Rda")
# get.outputs.depth(signal=.01, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h01_rho7_big.Rda")
# get.outputs.depth(signal=0, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h00_rho7_big.Rda")
# get.outputs.depth(signal=.1, rho=.7, Emap, sdn.vec, q.vec, t.vec, nrep, filename="E2_h10_rho7_big.Rda")

# nrep = 1e3
# get.outputs.depth(signal=.07, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h07_rho7_big.Rda")
# get.outputs.depth(signal=.05, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h05_rho7_big.Rda")
# get.outputs.depth(signal=.03, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h03_rho7_big.Rda")
# get.outputs.depth(signal=.02, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h02_rho7_big.Rda")
# get.outputs.depth(signal=.01, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h01_rho7_big.Rda")
# get.outputs.depth(signal=0, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h00_rho7_big.Rda")
# get.outputs.depth(signal=.1, rho=.7, sdn.vec, q.vec, t.vec, nrep, filename="E2_h10_rho7_big.Rda")
