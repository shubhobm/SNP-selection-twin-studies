setwd("d:/Study/My projects/SNP-selection-twin-studies/gedi5 outputs/")

library(fda.usc)
library(GenomicFeatures)
source('functions.R')

## takes a while
tx = makeTxDbFromUCSC(
  genome="hg19",
  tablename="knownGene",
  transcript_ids=NULL,
  circ_seqs=DEFAULT_CIRC_SEQS,
  url="http://genome.ucsc.edu/cgi-bin/",
  goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
  taxonomyId=NA,
  miRBaseBuild=NA)

## function to get outputs from a gene model
do.gene = function(gene){
  load(paste0("model",gene,".Rda"))
  best.list = step11.depth(model.list, sdn=sdn, q=q, mult=mult)
  load("SNP_Dictionary.Rda")
  gene.info.df = read.csv("geneinfo.csv")
  gene.info = gene.info.df[which(gene.info.df[,1]==gene),]
  gene.df =  snp.df[with(snp.df, which( chr == as.numeric(gene.info[2]) &
                                          pos_B37 >= as.numeric(gene.info[3]) &
                                          pos_B37 <= as.numeric(gene.info[4]))),]
  
  gene.df$evalue = NA
  gene.df$evalue[which(gene.df$snp %in% rownames(model.list$beta.hat))] = best.list$evalue
  
  gene.df$beta = NA
  gene.df$beta[which(gene.df$snp %in% rownames(model.list$beta.hat))] = model.list$beta.hat
  gene.df = gene.df[complete.cases(gene.df),]
  
  # make object for gene
  SNPdata = with(gene.df, data.frame(
    ASSOC = sign(beta),
    SNP.NAME = snp,
    LOC = pos_B37,
    SS.PVAL = evalue
  ))
  
  tx1 = tx
  seqlevels(tx1) = paste0("chr", gene.df$chr[1])
  # plotfun(SNPdata,tx1,gene)
  ## now plot
  attach(SNPdata)
  pdf(paste0("plot",gene,".pdf"),8,4.5)
  
  ##### Plotting starts
  defaultPar = par()
  par(mai = c(1.32,0.82,0.52,0.42))
  pos1 = gene.info$end
  pos0 = gene.info$start
  plot((1-SS.PVAL)~LOC, col="red", pch=19, xaxt='n',ylim=c(0,1), xlim=c(pos0,pos1),bty="n",
       main=gene, ylab="1 - e-value", xlab="")
  abline(h=1-q*mult, lty=2)
  abline(h=0)
  colvec = ifelse(SS.PVAL < q*mult, "red", "black")
  for(i in 1:nrow(SNPdata)){
    segments(LOC[i],0,LOC[i],(1-SS.PVAL)[i],col="blue")
  }
  
  buffer = .00*(pos1-pos0)
  baseseq = seq(pos0-buffer, pos1+buffer, by=10)
  text("|", col=grey(.8), x=baseseq, y=-.05, pos=1, xpd=TRUE, cex=1.5)
  text("|",col=colvec, x=LOC, y=-.05, pos=1, xpd=TRUE, cex=1.5)
  
  # plot exon regions
  col.bases = rep("lightgreen", length(baseseq))
  ex = exons(tx1, columns=c("EXONID", "TXNAME"))
  which.exons = which(ex@ranges@start >= pos0 & ex@ranges@start <= pos1)
  start.exons = ex@ranges@start[which.exons]
  end.exons = ex@ranges@start[which.exons] + ex@ranges@width[which.exons] - 1
  for(i in 1:length(which.exons)){
    col.bases[which(baseseq >= start.exons[i] & baseseq <= end.exons[i])] = "blue"
  }
  # which.blue = which(col.bases == "blue")
  text("|", col=col.bases, x=baseseq, y=-.2, pos=1, xpd=TRUE, cex=1.5)
  # text("|", col=col.bases[which.blue], x=baseseq[which.blue], y=-.2, pos=1, xpd=TRUE, cex=1)
  
  #text("SNP",x=min(LOC)-2.5e4, y=-.1, xpd=TRUE)
  text(c(pos0, paste("Position, in","chr",gene.df$chr[1]), pos1),
       x=c(pos0, (pos0+pos1)/2, pos1),
       y=-.35, pos=1, xpd=TRUE)
  legend("topright", c("Imp. SNP", "Unimp. SNP","Exon","Intron"),
         col=c("red","black","blue","lightgreen"), lwd=2, lty=1)
  par(defaultPar)
  ##### Plotting done
  
  dev.off()
  detach(SNPdata)
  
  # return output data frame
  SNPdata
}

## analyze all genes and get table of non-zero SNPs
set.seed(04192017)
q = .9
mult = .5
sdn = 30
geneinfo = read.csv("geneinfo.csv")
SNPdata.list = lapply(paste(geneinfo[,1]), do.gene)
SNP.df = data.frame(gene=paste(geneinfo[,1]),
                    total.SNP=rep(NA,9),
                    detect.SNP=rep(NA,9),
                    detect.SNP.name=rep(NA,9))
for(i in 1:9){
  SNPdata = SNPdata.list[[i]]
  SNP.df[i,2] = sum(!is.na(SNPdata$SS.PVAL))
  best.index = which(SNPdata$SS.PVAL < q*mult)
  SNP.df[i,3] = length(best.index)
  z = with(SNPdata, paste0(SNP.NAME,"(",ifelse(ASSOC==1,"+","-"),")"))
  SNP.df[i,4] = paste(z[best.index],collapse=", ")
}

mytable = xtable::xtable(SNP.df)
print(mytable, include.rownames = F)
