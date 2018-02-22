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
  resid.mat = matrix(rep(rnorm(nfam*nboot), rep(4,nfam*nboot)),
                     ncol=nboot, byrow=F) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  # wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  wildmat = matrix(rep(rnorm(nfam*nboot), rep(4,nfam*nboot)),
                   ncol=nboot, byrow=F)
  
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

step1.depth = function(mod, sdn, adj, nboot=1e3){
  ## first get full model to calculate index sets for each sdn
  beta.hat = mod$beta
  beta.cov = mod$beta.cov
  W = mod$W
  X1 = mod$X
  n1 = nrow(X1)
  H = beta.cov %*% t(X1) %*% W
  r = mod$model$y - mod$fitted
  n = nrow(X1)
  p = ncol(X1)
  
  ## matrix of full model bootstrap betas
  beta.mat = matrix(0, nrow=nboot, ncol=p)
  resid.mat = matrix(rnorm(n*nboot), ncol=nboot) * matrix(r, nrow=n, ncol=nboot, byrow=F)
  score.mat = t(H %*% resid.mat)
  depth.full = mdepth.RP(score.mat, score.mat)$dep
  Cn.full = mean(depth.full)
  
  wildmat = matrix(rnorm(n*nboot), ncol=nboot)
  loopfun = function(j){
    set.seed(j)
    ## calculate quantities for truncated model
    Hj = H[-j,] - outer(beta.cov[j,-j], as.numeric(W %*% X1[,j]))
    rj = r + X1[,j] * beta.hat[j]
    jresid.mat = wildmat * matrix(rj, nrow=n, ncol=nboot, byrow=F)
    jscore.mat = t(Hj %*% jresid.mat)
    
    ## calculate Cn for truncated model, for a range of bootstrap variances
    jbeta.mat = matrix(0, ncol=p, nrow=nboot)
    jbeta.mat[,-j] = matrix(beta.hat[-j], nrow=nboot, ncol=p-1, byrow=T) +
      sdn * jscore.mat
    jdepth.vec = mdepth.RP(jbeta.mat,
                           matrix(beta.hat, nrow=nboot, ncol=p, byrow=T)+ sdn*score.mat)$dep
    # return mean depth of truncated model, for all values of bootstrap sd
    jdepth.vec
  }
  
  depth.mat = lapply(1:p, loopfun)
  depth.mat = matrix(unlist(depth.mat), ncol=p, byrow=F)
  
  ## update all concerned variables
  # tail.probs = apply(depth.mat,2,function(x) (ecdf(x))(Cn.full))
  full.ecdf = ecdf(depth.full)
  tail.probs = apply(depth.mat,2,function(x) 1-full.ecdf(median(x)))
  tail.probs.adj = p.adjust(tail.probs, method=adj)
  list(full.dist=depth.full,
		drop.dists=depth.mat,
		select.ind=which(tail.probs.adj > median(tail.probs.adj)))
}

## creating folds for cross-validation. from caret package
createFolds = function (y, k = 10, list = TRUE, returnTrain = FALSE) 
{
    if (class(y)[1] == "Surv") 
        y <- y[, "time"]
    if (is.numeric(y)) {
        cuts <- floor(length(y)/k)
        if (cuts < 2) 
            cuts <- 2
        if (cuts > 5) 
            cuts <- 5
        breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
        y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
        y <- factor(as.character(y))
        numInClass <- table(y)
        foldVector <- vector(mode = "integer", length(y))
        for (i in 1:length(numInClass)) {
            min_reps <- numInClass[i]%/%k
            if (min_reps > 0) {
                spares <- numInClass[i]%%k
                seqVector <- rep(1:k, min_reps)
                if (spares > 0) 
                  seqVector <- c(seqVector, sample(1:k, spares))
                foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
            }
            else {
                foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                  size = numInClass[i])
            }
        }
    }
    else foldVector <- seq(along = y)
    if (list) {
        out <- split(seq(along = y), foldVector)
        names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
            sep = "")
        if (returnTrain) 
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
}