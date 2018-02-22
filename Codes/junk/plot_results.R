## Finding causal SNPs from twin data
## simulation study with relaxed definition of false positives
## see Frommelet et al CSDA paper
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')

get.plots = function(filename){
  load(filename)
  sd.vec = .1 + (1:20)/10
  
  # pdf(paste0('plot_',filename,'.pdf'),6,3)
  par(mfrow=c(1,2))
  
  TP.vec = unlist(lapply(1:20, function(i){ out.list[[i+1]][1,5]}))
  plot(sd.vec, TP.vec, lwd=2, type='l', col="red", ylim=c(0,1), main="TP", xlab="sd", ylab="TP")
  abline(h=out.list[[1]][1,5], col="blue", lwd=2)
  abline(h=out.list[[1]][3,5], col="black", lwd=2)
  
  TN.vec = unlist(lapply(1:20, function(i){ out.list[[i+1]][1,6]}))
  plot(sd.vec, TN.vec, lwd=2, type='l', col="red", ylim=c(0,1), main="TN", xlab="sd", ylab="TN")
  abline(h=out.list[[1]][1,6], col="blue", lwd=2)
  abline(h=out.list[[1]][3,6], col="black", lwd=2)
  par(mfrow=c(1,1))
  
  # dev.off()
}

get.plots('rho7_pt01.Rda')
get.plots('rho7_zero.Rda')
