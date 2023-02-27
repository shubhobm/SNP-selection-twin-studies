rm(list=ls())
library(locfdr)

setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')

# load data
load("../gedi5 outputs/Current/outputMZDZ_rfgls.Rda")
genes = paste(read.csv("../gedi5 outputs/geneinfo.csv")[,1])
names(model.list) = genes

for(i in 1:length(genes)){
  gf = locfdr(head(model.list[[i]], -4)$t.stat)
  print(genes[i])
  print(sum(gf$fdr < 0.05)) # LFDR correction
  print(sum(p.adjust(head(model.list[[i]], -4)$pval, method='hochberg')< 0.05)) # BH correction
}
