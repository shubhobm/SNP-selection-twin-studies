## generate_table.R
## generate latex table from simulation results
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes/E1E2_nosignal_outputs')

# make q and t vectors
h.vec = c("10", "07", "05", "03", "02", "01", "00")
nh = length(h.vec)

TP1.df = matrix(NA, nrow=nh, ncol=13)
TN1.df = matrix(NA, nrow=nh, ncol=13)

for(i in 1:nh){
  h = h.vec[i]

  # add BIC and RFGLS values
  load(paste0("others1_h",h,"_rho07_nosignal.Rda"))
  mean.mat = round(apply(simplify2array(all.list), 1:2, mean),2)
  
  TP1.df[i,1:2] = mean.mat[1,]
  TN1.df[i,1:2] = mean.mat[2,]
  
  # add evalue values
  load(paste0("E1_h",h,"_rho7_nosignal.Rda"))
  mean.mat = round(apply(simplify2array(all.list), 1:2, mean),2)
  
  TP1.df[i,3:7] = mean.mat[,1]
  TN1.df[i,3:7] = mean.mat[,2]
  
  load(paste0("E2_h",h,"_rho7_nosignal.Rda"))
  mean.mat = round(apply(simplify2array(all.list), 1:2, mean),2)
  
  TP1.df[i,8:13] = mean.mat[2*(1:6)-1,1]
  TN1.df[i,8:13] = mean.mat[2*(1:6)-1,2]
  
}

## relaxed definition
table.df = t(matrix(paste0(TP1.df, "/", TN1.df), ncol=13, byrow=F))
c1 = c("mBIC2", "RFGLS+BH",NA,NA,"E1",rep(NA,4),"E2",rep(NA,3))
c2 = c("","",paste0("t = exp(",1:5,")"), paste0("t = ",rev(seq(.5,.8,by=.06))))
table.df = cbind(c1,c2,table.df)

# insert h column
# table.df = cbind(paste0('h=',c(10,7,5,3,2,1,0)), table.df)

# make table with heading
mytable = xtable::xtable(table.df)
names(mytable) = c("Method","",paste("h =",c(10,7,5,3,2,1,0)))
print(mytable, include.rownames = F)