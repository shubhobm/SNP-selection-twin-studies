##################################################
#Simulation code to generate data on correlated SNPs in nuclear #families and simulate multivariate traits conditional on few 
#  causal variants
#  Reference: Basu S, Zhang Y, Ray D, Miller MB, Iacono WG, #McGue M (2013), A Rapid Genome-wide Gene-based Association #Tests with Multivariate Traits. Human Heredity 76(2):53-63. 
#
#############################################################
## LD is a vector of correlation among snps; each entry is for 1 block
## MAF is a vector of MAFs, each entry is for 1 block
## n.block is a vector the number of snps for each block
## n.person is the total number of each family member
## this code is adaptive to any number of blocks

simgen<-function(LD=c(0.9,0.8,0.9,0.8,0.9,0.8),
                 MAF=c(0.2,0.4,0.25,0.3,0.1,0.15),
                 n.block=c(6,4,6,4,6,4),
                 n.person=100
                 ){
  library(mvtnorm)
    xmat0<-NULL
for (block in 1:length(LD)){
    k<-n.block[block]
    Vg<-diag(NA,n.block[block],n.block[block])
    for (row in 1:k) for (col in 1:k){
      a<-ifelse(row==col,0,1)
      Vg[row,col]=LD[block]^a
    } 

## Generate genotype of Mother's and Father's
    mo<-fa<-ofs1<-ofs2<-matrix(NA,nrow=n.person,ncol=k) 
    z<-rmvnorm(n.person,rep(0,k),Vg)
    p<-MAF[block]
    mo<-apply(z,2,function(y)ifelse(y<=qnorm(p^2),2,ifelse(y>qnorm(1-(1-p)^2),0,1)))
    z<-rmvnorm(n.person,rep(0,k),Vg)
    fa<-apply(z,2,function(y)ifelse(y<=qnorm(p^2),2,ifelse(y>qnorm(1-(1-p)^2),0,1)))
## offspring 1    
    for (snp in 1:k){      
      ofs1[which(mo[,snp]==0 & fa[,snp]==0),snp]=0
      ofs1[which((fa[,snp]==0 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==0)) ,snp]=1
      ofs1[which(fa[,snp]==2 & mo[,snp]==2),snp]=2
      ofs1[which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0)),snp] = sample(c(0,1),size=length(which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0))),prob=c(0.5,0.5),replace=T)
      ofs1[which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1)),snp] = sample(c(1,2),size=length(which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1))),prob=c(0.5,0.5),replace=T)
      ofs1[which(fa[,snp]==1 & mo[,snp]==1),snp] = sample(c(0,1,2),size=length(which(fa[,snp]==1 & mo[,snp]==1)),prob=c(0.25,0.5,0.25),replace=T)
    }
## offspring 2
    for (snp in 1:k){      
      ofs2[which(mo[,snp]==0 & fa[,snp]==0),snp]=0
      ofs2[which((fa[,snp]==0 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==0)) ,snp]=1
      ofs2[which(fa[,snp]==2 & mo[,snp]==2),snp]=2
      ofs2[which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0)),snp] = sample(c(0,1),size=length(which((fa[,snp]==0 & mo[,snp]==1) | (fa[,snp]==1 & mo[,snp]==0))),prob=c(0.5,0.5),replace=T)
      ofs2[which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1)),snp] = sample(c(1,2),size=length(which((fa[,snp]==1 & mo[,snp]==2) | (fa[,snp]==2 & mo[,snp]==1))),prob=c(0.5,0.5),replace=T)
      ofs2[which(fa[,snp]==1 & mo[,snp]==1),snp] = sample(c(0,1,2),size=length(which(fa[,snp]==1 & mo[,snp]==1)),prob=c(0.25,0.5,0.25),replace=T)
    }
## putting snps together
    xmat<-matrix(NA, nrow=4*n.person,ncol=k)
    xmat[seq(1,4*n.person,4),]=mo
    xmat[seq(2,4*n.person,4),]=fa
    xmat[seq(3,4*n.person,4),]=ofs1
    xmat[seq(4,4*n.person,4),]=ofs2
    xmat0<-cbind(xmat0,xmat)  
  } 
    return(xmat0)
}

### to generate 3 correlated traits 
# sigma2 is a vector of polygenic variance, each entry for 1 trait
# sigmae2 is a vector of enviromental variance
# rho =( rho12,rho13,rho23,rho14,rho24,rho34): correlation between environmental effect of phenotype i with phenotype j of the same individual
# r is within family correlation for each trait; I assumed it's the same for all the traits
# causal SNPs is NOT randomly selected; it is the first snp in the block
# if partial=0,causal SNPs affecting all phenotypes; if partial=a vector, causal SNP is only affecting some of the phenotypes
# type can be "DZ","MZ","ADOPT"
############################################################
# we also need the parameters for genotypes
# MAF: MAF for each block
# n.block: block size
# h2 is the vector of genetic variance for the causal SNP in each block; it must have length as number of blocks (same length as MAF)


simphe<-function(gendat=gen,
                 b0=0,
                 sigma2=c(60,50,40),
                 sigmae2=c(40,30,20),
                 nfam=100,
                 r=0.5,
                 rho=c(0,0,0),
                 h2=c(0,0,0,0,0,0),
                 partial=0,
                 MAF=c(0.2,0.4,0.25,0.3,0.1,0.15),
                 n.block=c(6,4,6,4,6,4),
                 rhoAB=c(0.33,0.33,0.33),
                 type=c("DZ")){
   library(mvtnorm)
   library(Matrix)
   ntrait=length(sigma2)
   if(type=="DZ") phi<-matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0,0.5,0.5,0,1),nrow=4,byrow=T) else
   if(type=="MZ") phi<-matrix(c(1,1,0.5,0.5,1,1,0.5,0.5,0.5,0.5,1,0,0.5,0.5,0,1),nrow=4,byrow=T) else
   if (type=="ADOPT") phi<-diag(4)
    I=diag(4)
    #W=r*matrix(1,4,4)+(1-r)*I
    W=I
### first generate Random effect term, which is polygenic effect+ environment effect
### construction of covariance matrix
	Vmat<-matrix(0,4*ntrait,4*ntrait)
## first put in diagonals
   	 lists.diag<-list()
   	 for (b in 1:ntrait){
	     lists.diag[[b]]<-sigma2[b]*phi+sigmae2[b]*W
}
	Vmat<-bdiag(lists.diag)

## then put in off-diagonal entries
	if (ntrait>1){
	rho.index<-1
   	for (b in 2:ntrait){
	g=1
	while (g<b){
  #           Vmat[(4*(g-1)+1):(4*g),(4*(b-1)+1):(4*b)]<-sqrt(sigma2[b]*sigma2[g])*phi+rho[rho.index]*sqrt(sigmae2[b]*sigmae2[g])*I
         rhop=sqrt(sigma2[b]*sigma2[g])*rhoAB[rho.index]/sqrt((sigma2[b]+sigmae2[b])*(sigma2[g]+sigmae2[g]))
        Vmat[(4*(g-1)+1):(4*g),(4*(b-1)+1):(4*b)]<-rhop*sqrt((sigma2[b]+sigmae2[b])*(sigma2[g]+sigmae2[g]))*phi+rho[rho.index]*sqrt(sigmae2[b]*sigmae2[g])*I
	    rho.index=rho.index+1   
	     g=g+1
       }
	}
}
    Vmat<-as.matrix(forceSymmetric(Vmat))
   #print(Vmat)
    e<-rmvnorm(nfam,rep(0,ntrait*4),Vmat)
   RE<-matrix(0,ncol=ntrait, nrow=4*nfam)
   for (t in 1:ntrait) RE[,t]<-array(t(e[,(4*(t-1)+1):(4*t)]))
     
## additive model
## the first SNP of each block is the causal SNP
   h2<-as.matrix(h2)
   MAF<-as.matrix(MAF)
   b<-c(sqrt(h2/(2*MAF*(1-MAF))))
    nsnp<-dim(gendat)[2]
   beta<-rep(0,nsnp)
   b.index<-b0.index<-1
   if (length(n.block)>1){
   for (bb in 1:(length(n.block)-1)){
     b.index<-b.index+n.block[bb]
     b0.index<-c(b0.index,b.index)
   }}else b0.index=1
   # print(b0.index)
   beta[b0.index]<-b
   beta<-as.matrix(beta)
   mu=as.matrix(gendat)%*%beta+b0
   if (prod(partial)==0) Y=apply(RE,2,function(y)y+mu) 
else if(prod(partial)>0){
   Y1<-apply(as.matrix(RE[,partial]),2,function(y)y+mu)  
   Y<-cbind(Y1,RE[,-partial])
 }
   return(list(Y=Y,beta=beta))     
   }
   
  
simphe.r <-function(gendat=gen,
                 b0=0,
                 sigma2=c(60,50,40),
                 sigmae2=c(40,30,20),
                 nfam=100,
                 r=0.5,
                 rho=c(0,0,0),
                 h2=c(0,0,0,0,0,0),
                 partial=0,
                 MAF=c(0.2,0.4,0.25,0.3,0.1,0.15),
                 n.block=c(6,4,6,4,6,4),
                 rhoAB=c(0.33,0.33,0.33),
                 type=c("DZ")){
  library(mvtnorm)
  library(Matrix)
  ntrait=length(sigma2)
  if(type=="DZ") phi<-matrix(c(1,0.5,0.5,0.5,0.5,1,0.5,0.5,0.5,0.5,1,0,0.5,0.5,0,1),nrow=4,byrow=T) else
    if(type=="MZ") phi<-matrix(c(1,1,0.5,0.5,1,1,0.5,0.5,0.5,0.5,1,0,0.5,0.5,0,1),nrow=4,byrow=T) else
      if (type=="ADOPT") phi<-diag(4)
  I=diag(4)
  W=r*matrix(1,4,4)+(1-r)*I

  ### first generate Random effect term, which is polygenic effect+ environment effect
  ### construction of covariance matrix
  Vmat<-matrix(0,4*ntrait,4*ntrait)
  ## first put in diagonals
  lists.diag<-list()
  for (b in 1:ntrait){
    lists.diag[[b]]<-sigma2[b]*phi+sigmae2[b]*W
  }
  Vmat<-bdiag(lists.diag)
  
  ## then put in off-diagonal entries
  if (ntrait>1){
    rho.index<-1
    for (b in 2:ntrait){
      g=1
      while (g<b){
        #           Vmat[(4*(g-1)+1):(4*g),(4*(b-1)+1):(4*b)]<-sqrt(sigma2[b]*sigma2[g])*phi+rho[rho.index]*sqrt(sigmae2[b]*sigmae2[g])*I
        rhop=sqrt(sigma2[b]*sigma2[g])*rhoAB[rho.index]/sqrt((sigma2[b]+sigmae2[b])*(sigma2[g]+sigmae2[g]))
        Vmat[(4*(g-1)+1):(4*g),(4*(b-1)+1):(4*b)]<-rhop*sqrt((sigma2[b]+sigmae2[b])*(sigma2[g]+sigmae2[g]))*phi+rho[rho.index]*sqrt(sigmae2[b]*sigmae2[g])*I
        rho.index=rho.index+1   
        g=g+1
      }
    }
  }
  Vmat<-as.matrix(forceSymmetric(Vmat))
  #print(Vmat)
  e<-rmvnorm(nfam,rep(0,ntrait*4),Vmat)
  RE<-matrix(0,ncol=ntrait, nrow=4*nfam)
  for (t in 1:ntrait) RE[,t]<-array(t(e[,(4*(t-1)+1):(4*t)]))
  
  ## additive model
  ## the first SNP of each block is the causal SNP
  h2<-as.matrix(h2)
  MAF<-as.matrix(MAF)
  b<-c(sqrt(h2/(2*MAF*(1-MAF))))
  nsnp<-dim(gendat)[2]
  beta<-rep(0,nsnp)
  b.index<-b0.index<-1
  if (length(n.block)>1){
    for (bb in 1:(length(n.block)-1)){
      b.index<-b.index+n.block[bb]
      b0.index<-c(b0.index,b.index)
    }}else b0.index=1
  # print(b0.index)
  beta[b0.index]<-b
  beta<-as.matrix(beta)
  mu=as.matrix(gendat)%*%beta+b0
  if (prod(partial)==0) Y=apply(RE,2,function(y)y+mu) 
  else if(prod(partial)>0){
    Y1<-apply(as.matrix(RE[,partial]),2,function(y)y+mu)  
    Y<-cbind(Y1,RE[,-partial])
  }
  return(list(Y=Y,beta=beta))     
}

	



