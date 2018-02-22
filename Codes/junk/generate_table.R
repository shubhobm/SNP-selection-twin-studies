## generate_table.R
## generate latex table from simulation results
rm(list=ls())
setwd('d:/Study/My projects/SNP-selection-twin-studies/Codes')

# make q and t vectors
q.vec = c("9", "5", "2", "1", "05")
t.vec = c("8", "7", "6", "5")
h.starts = 5*(0:4)

TP.df = matrix(NA, nrow=25, ncol=6)
for(i in 1:5){
  for(j in 1:4){
    q = q.vec[i]
    t = t.vec[j]
    load(paste0("Varcriterion/all_q",q,"_mult",t,".Rda"))
    if(j == 1){
      TP.df[h.starts+3,1:2] = matrix(unlist(mapply(function(x){c(x[1,5],x[3,5])}
                                                      ,out.list)),
                                        ncol=2, byrow=T)
    }
    TP.df[h.starts+i,2+j] = unlist(mapply(function(x){x[5,5]} ,out.list))
  }
}

TN.df = matrix(NA, nrow=25, ncol=6)
for(i in 1:5){
  for(j in 1:4){
    h.starts = 5*(0:4)
    q = q.vec[i]
    t = t.vec[j]
    load(paste0("all_q",q,"_mult",t,".Rda"))
    if(j == 1){
      TN.df[h.starts+3,1:2] = matrix(unlist(mapply(function(x){c(x[1,6],x[3,6])}
                                                      ,out.list)),
                                        ncol=2, byrow=T)
    }
    TN.df[h.starts+i,2+j] = unlist(mapply(function(x){x[5,6]} ,out.list))
  }
}

table.df = matrix(paste0(round(TP.df,2), "/", round(TN.df,2)), ncol=6, byrow=F)
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
