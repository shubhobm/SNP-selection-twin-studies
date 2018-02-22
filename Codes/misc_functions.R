## misc_functions: miscellaneous functions for ease of use
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

ecmeml1 = function (y, subj, pred, xcol, zcol, vmax, occ, start, maxits = 1000, 
                    eps = 1e-04) {
  tmp <- table(subj)
  m <- length(tmp)
  nmax <- max(tmp)
  ntot <- length(y)
  pcol <- ncol(pred)
  q <- length(zcol)
  p <- length(xcol)
  {
    if (missing(vmax)) {
      vmax <- diag(rep(1, nmax))
      occ <- integer(ntot)
      iflag <- as.integer(1)
    }
    else iflag <- as.integer(0)
  }
  storage.mode(vmax) <- "double"
  {
    if (!missing(start)) {
      beta <- start$beta
      sigma2 <- start$sigma2
      xi <- start$psi/start$sigma2
      storage.mode(beta) <- "double"
      storage.mode(xi) <- "double"
      storage.mode(sigma2) <- "double"
      sflag <- as.integer(1)
    }
    else {
      beta <- numeric(p)
      xi <- matrix(0, q, q)
      sigma2 <- 0
      sflag <- as.integer(0)
    }
  }
  # if(trace==F){
  #   cat("Performing ECME...")
  # }
  now <- proc.time()
  err <- 0
  tmp <- .Fortran("ecmeml", ntot, as.integer(subj), m, ist = integer(m), 
                  ifin = integer(m), as.integer(occ), nmax, vmax, w = array(0, 
                                                                            c(nmax, nmax, m)), vinv = array(0, c(nmax, nmax, 
                                                                                                                 m)), pcol, as.double(pred), q, as.integer(zcol), 
                  ztvinv = array(0, c(q, nmax, m)), ztvinvz = array(0, 
                                                                    c(q, q, m)), iflag = iflag, err = as.integer(err), 
                  msg = integer(1), u = array(0, c(q, q, m)), iter = integer(1), 
                  sflag, sigma2 = sigma2, p, as.integer(xcol), beta = beta, 
                  as.double(y), delta = rep(0, ntot), xtw = matrix(0, p, 
                                                                   nmax), xtwx = matrix(0, p, p), xtwy = numeric(p), 
                  xtwxinv = matrix(0, p, p), wkqq1 = matrix(0, q, q), wkqq2 = matrix(0, 
                                                                                     q, q), xi = xi, wkqnm = array(0, c(q, nmax, m)), 
                  b = matrix(0, q, m), cvgd = integer(1), obeta = rep(0, 
                                                                      p), oxi = matrix(0, q, q), maxits = as.integer(maxits), 
                  llvec = numeric(as.integer(maxits)), eps = as.double(eps), 
                  PACKAGE = "lmm")
  clock <- proc.time() - now
  # if(trace==F){
  #   cat("\n")
  # }
  iter <- tmp$iter
  msg <- tmp$msg
  {
    if (msg == 1) 
      warning("Supplied V <- i matrix is not positive definite")
    else if (msg == 2) 
      warning("GLS failed for start vals, t(X)%*%inv(V)%*%X not full rank")
    else if (msg == 3) 
      warning("Inadequate information to obtain starting value of psi")
    else if (msg == 4) 
      warning("Value of psi became non-pos.def. during iterations")
    else if (msg == 5) 
      warning("t(X)%*%W%*%X became non-pos.def. during iterations")
  }
  llvec <- tmp$llvec[1:iter]
  converged <- tmp$cvgd == as.integer(1)
  cov.beta <- tmp$xtwxinv * tmp$sigma2
  b.hat <- tmp$b
  cov.b <- tmp$u * tmp$sigma2
  psi <- tmp$xi * tmp$sigma2
  beta <- tmp$beta
  if (!is.null(dimnames(pred)[[2]])) {
    colnames <- dimnames(pred)[[2]]
    names(beta) <- colnames[xcol]
    dimnames(psi) <- list(colnames[zcol], colnames[zcol])
  }
  list(beta = beta, sigma2 = tmp$sigma2, psi = psi, converged = converged, 
       iter = iter, loglik = llvec, cov.beta = cov.beta, b.hat = b.hat, 
       cov.b = cov.b)
}

depth.dist = function(x, xx, Emap){
  
  # get center and scale factors
  mu = as.numeric(apply(xx,2,mean))
  sig = as.numeric(apply(xx,2,sd))
  
  xs = scale(x, mu, sig)
  
  # calculate evaluation map
  Eval = 1/(1+rowSums(xs^2))
  if(Emap=="E2"){
    Eval = exp(-rowSums(abs(xs)))
  }
  Eval
}

step11.depth = function(mod, sdn, q=0.05, mult=.5, nboot=1e3){
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta.hat
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X1
  H = mod$H
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  loopfun = function(j){
    set.seed(j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = mdepth.RP(jbeta.mat, beta.mat)$dep
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  ## plot to check
  # plot(density(depth.full), xlim=c(0,.5), ylim=c(0,20))
  # for(i in 1:ncol(depth.mat)){ lines(density(depth.mat[,i]), col="red")}
  # lines(density(depth.full), lwd=2)
  # 
  ## update all concerned variables
  # tail.probs = apply(depth.mat,2,function(x) (ecdf(x))(Cn.full))
  full.ecdf = ecdf(depth.full)
  tail.probs = apply(depth.mat,2,function(x) full.ecdf(quantile(x,q)))
  # tail.probs.adj = p.adjust(tail.probs, method=adj)
  # which(tail.probs.adj < median(tail.probs.adj))
  list(evalue=tail.probs,
       select=which(tail.probs < q*mult))
}

step.depth.qt = function(mod, sdn, q.vec, t.vec, nboot=1e3){
  
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta.hat
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X1
  H = mod$H
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  loopfun = function(j){
    set.seed(j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = mdepth.RP(jbeta.mat, beta.mat)$dep
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    # evalues.order = order(evalues, decreasing = F)
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      jlist[[j]] = which(evalues < qi*tj)
    }
    qt.list[[i]] = jlist
  }
  
  qt.list
}

step.depth.fdr = function(mod, sdn, q.vec, t.vec, nboot=1e3){
  
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta.hat
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X1
  H = mod$H
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rnorm(nfam*nboot), rep(4, nfam*nboot)), ncol=nboot, byrow=F)
  loopfun = function(j){
    set.seed(j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = mdepth.RP(jbeta.mat, beta.mat)$dep
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    eval.adj = 1 - p.adjust(1-evalues, "fdr")
    # evalues.order = order(evalues, decreasing = F)
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      jlist[[j]] = which(eval.adj < qi*tj)
    }
    qt.list[[i]] = jlist
  }
  
  qt.list
}

step.nodepth.qtgamma0 = function(mod, sdn, Emap, q.vec, t.vec, nboot=1e3){
  
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta.hat
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X1
  H = mod$H
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  # depth.full = mdepth.RP(score.mat, score.mat)$dep
  depth.full = depth.dist(score.mat, score.mat, Emap)
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  loopfun = function(j){
    set.seed(1e3*j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    # jdepth.vec = mdepth.RP(jbeta.mat, beta.mat)$dep
    jdepth.vec = depth.dist(jbeta.mat, beta.mat, Emap)
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))

    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      jlist[[j]] = which(evalues < qi*tj)
    }
    qt.list[[i]] = jlist
  }
  
  qt.list
}

step.nodepth.qtgamma = function(mod, sdn, Emap, q.vec, t.vec, nboot=1e3, indices=NULL){
  
  ## first get full model to calculate index sets for each sdn
  if(is.null(indices)){
    indices = 1:ncol(mod$X1) # select over all predictors if indices is not specified
  }
  p = length(indices)
  
  beta.hat = mod$beta.hat[indices,]
  beta.cov = mod$beta.cov[indices,indices]
  W = mod$W
  X1 = mod$X1[,indices]
  H = mod$H[indices,]
  r = mod$r
  n = nrow(X1)
  nfam = n/4

  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H[indices,] %*% resid.mat)
  depth.full = depth.dist(score.mat, score.mat, Emap)
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  loopfun = function(j){
    set.seed(1e3*j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = depth.dist(jbeta.mat, beta.mat, Emap)
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(indices, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      jlist[[j]] = which(evalues < qi*tj)
    }
    jlist[[nt+1]] = evalues
    qt.list[[i]] = jlist
  }
  
  qt.list
}

step.nodepth.qtgamma.fdr = function(mod, sdn, Emap, q.vec, t.vec, nboot=1e3, indices=NULL){
  
  ## first get full model to calculate index sets for each sdn
  if(is.null(indices)){
    indices = 1:ncol(mod$X1) # select over all predictors if indices is not specified
  }
  p = length(indices)
  
  beta.hat = mod$beta.hat[indices,]
  beta.cov = mod$beta.cov[indices,indices]
  W = mod$W
  X1 = mod$X1[,indices]
  H = mod$H[indices,]
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H[indices,] %*% resid.mat)
  depth.full = depth.dist(score.mat, score.mat, Emap)
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  loopfun = function(j){
    set.seed(1e3*j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = depth.dist(jbeta.mat, beta.mat, Emap)
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(indices, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    evalues.order = order(evalues, decreasing = F)
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      thres = qi*tj*(1:p)/p
      which.less = which(evalues[evalues.order] < thres)
      jlist[[j]] = integer(0)
      if(length(which.less)>0){
        jlist[[j]] = sort(evalues.order[ 1: max(which.less)])
      }
    }
    jlist[[nt+1]] = evalues
    qt.list[[i]] = jlist
  }
  
  qt.list
}

step.nodepth.qtgamma.both = function(mod, sdn, Emap, q.vec, t.vec, nboot=1e3, indices=NULL){
  
  ## first get full model to calculate index sets for each sdn
  if(is.null(indices)){
    indices = 1:ncol(mod$X1) # select over all predictors if indices is not specified
  }
  p = length(indices)
  
  beta.hat = mod$beta.hat[indices,]
  beta.cov = mod$beta.cov[indices,indices]
  W = mod$W
  X1 = mod$X1[,indices]
  H = mod$H[indices,]
  r = mod$r
  n = nrow(X1)
  nfam = n/4
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  # wildmat1 = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat1 = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  resid.mat = wildmat1 * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H[indices,] %*% resid.mat)
  depth.full = depth.dist(score.mat, score.mat, Emap)
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rgamma(nfam*nboot, shape=1), rep(4, nfam*nboot)), ncol=nboot, byrow=F) - 1
  loopfun = function(j){
    set.seed(1e3*j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    beta.mat = matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    # jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T)
    # jbeta.mat = jbeta.mat + sdn*score.mat
    jdepth.vec = depth.dist(jbeta.mat, beta.mat, Emap)
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(indices, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  full.ecdf = ecdf(depth.full)
  
  nq = length(q.vec)
  nt = length(t.vec)
  qt.list0 = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    evalues.order = order(evalues, decreasing = F)
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      thres = qi*tj*(1:p)/p
      which.less = which(evalues[evalues.order] < thres)
      jlist[[j]] = integer(0)
      if(length(which.less)>0){
        jlist[[j]] = sort(evalues.order[ 1: max(which.less)])
      }
    }
    jlist[[nt+1]] = evalues
    qt.list0[[i]] = jlist
  }
  
  qt.list1 = list()
  
  # get selected index set corresponding to (q,t)
  for(i in 1:nq){
    qi = q.vec[i]
    evalues = apply(depth.mat,2,function(x) full.ecdf(quantile(x,qi)))
    
    jlist = list()
    for(j in 1:nt){
      tj = t.vec[j]
      jlist[[j]] = which(evalues < qi*tj)
    }
    jlist[[nt+1]] = evalues
    qt.list1[[i]] = jlist
  }
  
  list(qt.list0, qt.list1)
}




