## generate_table.R
## generate latex table from simulation results
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')

# make q and t vectors
h.vec = c("10", "05", "02", "01", "00")
h.starts = 5*(0:4)

TP0.df = matrix(NA, nrow=25, ncol=6)
TN0.df = matrix(NA, nrow=25, ncol=6)
TP1.df = matrix(NA, nrow=25, ncol=6)
TN1.df = matrix(NA, nrow=25, ncol=6)

for(i in 1:5){
  h = h.vec[i]
  index.i = h.starts[i] + 1:5
  
  # add BIC and RFGLS values
  load(paste0("others1_h",h,"_rho07_big.Rda"))
  mean.mat = round(apply(simplify2array(all.list), 1:2, mean),2)
  
  TP0.df[index.i[3],1:2] = mean.mat[1,]
  TN0.df[index.i[3],1:2] = mean.mat[2,]
  TP1.df[index.i[3],1:2] = mean.mat[3,]
  TN1.df[index.i[3],1:2] = mean.mat[4,]
  
  # add evalue values
  load(paste0("qtg_h",h,"_rho7_big.Rda"))
  mean.mat = round(apply(simplify2array(all.list), 1:2, mean),2)
  
  TP0.df[index.i,3:6] = mean.mat[1:5, 2:5]
  TN0.df[index.i,3:6] = mean.mat[6:10, 2:5]
  TP1.df[index.i,3:6] = mean.mat[11:15, 2:5]
  TN1.df[index.i,3:6] = mean.mat[16:20, 2:5]
  
}

table.df = matrix(paste0(TP0.df, "/", TN0.df), ncol=6, byrow=F)
table.df[which(table.df=="NA/NA", arr.ind=T)] = " "

# insert h column
table.df = cbind(NA, table.df)
table.df[h.starts+3,1] = paste0('h=',c(10,5,2,1,0))

# insert q column
q.column = rep(c(0.95, 0.9, 0.5, 0.2, 0.1), 5)
table.df = cbind(table.df[,1:3], q.column, table.df[,4:7])

# make table with heading
mytable = xtable::xtable(table.df)
names(mytable) = c(" ","mBIC2","RFGLS+BH", "q", paste0('t=',as.numeric(t.vec)/10))
print(mytable, include.rownames = F)

table.df = matrix(paste0(TP1.df, "/", TN1.df), ncol=6, byrow=F)
table.df[which(table.df=="NA/NA", arr.ind=T)] = " "

# insert h column
table.df = cbind(NA, table.df)
table.df[h.starts+3,1] = paste0('h=',c(10,5,2,1,0))

# insert q column
q.column = rep(c(0.9, 0.5, 0.2, 0.1, 0.05), 5)
table.df = cbind(table.df[,1:3], q.column, table.df[,4:7])

# make table with heading
mytable = xtable::xtable(table.df)
names(mytable) = c(" ","BIC","RFGLS+FDR", "q", paste0('t=',as.numeric(t.vec)/10))
print(mytable, include.rownames = F)
