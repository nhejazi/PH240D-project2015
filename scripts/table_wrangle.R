
library(xtable)

proj_dir <- "/Users/nimahejazi/github-repos/bmcsa-project240D/data/"
setwd(proj_dir)

data.1 = read.table("diff_exp/DE.benzene.txt", header=TRUE, 
                    col.names=c('Gene','Ratio','Perm P-Value','Perm Q-Value'))

data.2 = read.table("diff_exp/DE.nickel.txt", header=TRUE, 
                    col.names=c('logFC','AveExpr','t','P.Value','adj.P.Val','B'))
data.2 <- data.2[-c(6,7),]

data.3 = read.table("diff_exp/DE.smoke.txt", header=TRUE, 
                    col.names=c('logFC','AveExpr','t','P.Value','adj.P.Val','B'))
data.3 <- data.3[-c(2,7,10),]
data.3 <- data.3[-6,]

data.auc <- read.csv("AUC_table.csv")
data.auc <- data.auc[-c(1,8,9),-c(1,2)]
colnames(data.auc) <- c("Benzene","Nickel","PAD","Smoking",
                        "Arsenic","Stress","Arthritis")
rownames(data.auc) <- c('ACSL1 & AQP9','NFKB1 & IFNB1','PRG2 & ACSL1',
                        'PRG2 & CLEC5A','NFKB1 & CLEC5A','ACSL1 & CLEC5A')
xtable(data.1)
xtable(data.2)
xtable(data.3)
xtable(data.auc)

