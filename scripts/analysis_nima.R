### ================== ###
### Nima Hejazi        ###
### Public Health 240D ###
### Group Project      ###
### Dec. 03, 2015      ###
### ================== ###

# analysis for genes on Rheumatoid Arthritis and Stress
rm(list=ls())
set.seed(0) # for reproducibility across group analyses
data_dir = '/Users/nimahejazi/github-repos/bmcsa-project240D/data/' # change to run
stress_data = "final_stress_genes.txt" # change to run
rheumatoid_data = "final_ra_genes.txt" # change to run
if (getwd() != data_dir) { setwd(data_dir) }

stress_data = read.table(stress_data, sep="", fill=FALSE, 
                         strip.white=TRUE)

stress_pairs1.v1 = as.data.frame(t(rbind(stress_data["ACSL1", ], 
                                         stress_data["AQP9", ])))
stress_pairs1.v2 = as.data.frame(t(rbind(stress_data["ACSL1", ], 
                                         stress_data["AQP9.2", ])))
stress_pairs2 = as.data.frame(t(rbind(stress_data["NFKB1", ], 
                                      stress_data["IFNB1", ])))
stress_pairs3.v1 = as.data.frame(t(rbind(stress_data["PRG2", ], 
                                         stress_data["ACSL1", ])))
stress_pairs3.v2 = as.data.frame(t(rbind(stress_data["PRG2.2", ], 
                                         stress_data["ACSL1", ])))
stress_pairs4.v1 = as.data.frame(t(rbind(stress_data["PRG2", ], 
                                         stress_data["CLEC5A", ])))
stress_pairs4.v2 = as.data.frame(t(rbind(stress_data["PRG2.2", ], 
                                         stress_data["CLEC5A", ])))
stress_pairs5 = as.data.frame(t(rbind(stress_data["NFKB1", ], 
                                      stress_data["CLEC5A", ])))
stress_pairs6 = as.data.frame(t(rbind(stress_data["ACSL1", ], 
                                      stress_data["CLEC5A", ])))

stress_outcome <- replace(as.vector(rep(1,length(colnames(stress_data)))), 
                          grep('Non',colnames(stress_data)), 0)
if (length(stress_outcome) != length(colnames(stress_data))) { 
  print('PROBLEM WITH OUTCOMES - STOP')}


rheumatoid_data = read.table(rheumatoid_data, sep="", fill=FALSE, 
                             strip.white=TRUE)

rheumatoid_pairs1.v1 = as.data.frame(t(rbind(rheumatoid_data["ACSL1", ], 
                                             rheumatoid_data["AQP9", ])))
rheumatoid_pairs1.v2 = as.data.frame(t(rbind(rheumatoid_data["ACSL1", ], 
                                             rheumatoid_data["AQP9.2", ])))
rheumatoid_pairs1.v3 = as.data.frame(t(rbind(rheumatoid_data["ACSL1", ], 
                                             rheumatoid_data["AQP9.3", ])))
rheumatoid_pairs2 = as.data.frame(t(rbind(rheumatoid_data["NFKB1", ], 
                                          rheumatoid_data["IFNB1", ])))
rheumatoid_pairs3.v1 = as.data.frame(t(rbind(rheumatoid_data["PRG2", ], 
                                             rheumatoid_data["ACSL1", ])))
rheumatoid_pairs3.v2 = as.data.frame(t(rbind(rheumatoid_data["PRG2.2", ], 
                                             rheumatoid_data["ACSL1", ])))
rheumatoid_pairs4.v1 = as.data.frame(t(rbind(rheumatoid_data["PRG2", ], 
                                             rheumatoid_data["CLEC5A", ])))
rheumatoid_pairs4.v2 = as.data.frame(t(rbind(rheumatoid_data["PRG2.2", ], 
                                             rheumatoid_data["CLEC5A", ])))
rheumatoid_pairs5 = as.data.frame(t(rbind(rheumatoid_data["NFKB1", ], 
                                          rheumatoid_data["CLEC5A", ])))
rheumatoid_pairs6 = as.data.frame(t(rbind(rheumatoid_data["ACSL1", ], 
                                          rheumatoid_data["CLEC5A", ])))

rheumatoid_outcome <- replace(as.vector(rep(1,length(colnames(rheumatoid_data)))), 
                              grep('.Control',colnames(rheumatoid_data)), 0)
if (length(rheumatoid_outcome) != length(colnames(rheumatoid_data))) { 
  print('PROBLEM WITH OUTCOMES - STOP')}

rm(data_dir, stress_data, rheumatoid_data)



# Setting up Super Learner
library(ROCR)
library(cvAUC)
library(SuperLearner)
SL.lib <- c("SL.glm","SL.stepAIC","SL.bayesglm","SL.randomForest","SL.gam","SL.mean")
CV_folds = 10


# Applying Super Learner to "Stress" data
Y = stress_outcome
levels(Y) = list("0"="CON", "1"="DIS")
Y = as.vector(Y)
Y = as.numeric(Y) # numeric vector of our zero-one outcomes


# Xfinal=data.frame(X=normal[,c(18,23)]) # this should be an nx2 data frame where the two columns
# # are the expressions/counts of the 2 genes
# 
# fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
# fld=fit.test.final$fold
# predsY.final=fit.test.final$SL.predict
# fold=rep(NA,length(predsY.final))
# for(k in 1:V) {
#   ii=unlist(fld[k])
#   fold[ii]=k
# }
# 
# ciout.final=ci.cvAUC(predsY.final, Y, folds = fold)
# ciout.final
# 
# # graph the ROC curve
# txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
# pred <- prediction(predsY.final,Y)
# perf1 <- performance(pred, "sens", "spec")
# plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
#      main="ROC Curve, Normalized Counts")
# text(0.6,0.4,txt)
# abline(0,1)


# TWO GENES IN THE SL: ACSL1, AQP9 
set.seed(0) # for reproducibility across group analyses
X = stress_pairs1.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p1.ciout.final.1 <- stress.ciout.final.a
} else {
  stress.p1.ciout.final.1 <- stress.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = stress_pairs1.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p1.ciout.final.2 <- stress.ciout.final.a
} else {
  stress.p1.ciout.final.2 <- stress.ciout.final.b
}


# TWO GENES IN THE SL: NFKB1, IFNB1 
set.seed(0) # for reproducibility across group analyses
X = stress_pairs2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p2.ciout.final <- stress.ciout.final.a
} else {
  stress.p2.ciout.final <- stress.ciout.final.b
}


# TWO GENES IN THE SL: PRG2,ACSL1 
set.seed(0) # for reproducibility across group analyses
X = stress_pairs3.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p3.ciout.final.1 <- stress.ciout.final.a
} else {
  stress.p3.ciout.final.1 <- stress.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = stress_pairs3.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p3.ciout.final.2 <- stress.ciout.final.a
} else {
  stress.p3.ciout.final.2 <- stress.ciout.final.b
}


# TWO GENES IN THE SL: PRG2, CLEC5A 
set.seed(0) # for reproducibility across group analyses
X = stress_pairs4.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p4.ciout.final.1 <- stress.ciout.final.a
} else {
  stress.p4.ciout.final.1 <- stress.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = stress_pairs4.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p4.ciout.final.2 <- stress.ciout.final.a
} else {
  stress.p4.ciout.final.2 <- stress.ciout.final.b
}


# TWO GENES IN THE SL: NFKB1, CLEC5A 
set.seed(0) # for reproducibility across group analyses
X = stress_pairs5
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p5.ciout.final <- stress.ciout.final.a
} else {
  stress.p5.ciout.final <- stress.ciout.final.b
}


# TWO GENES IN THE SL: ACSL1, CLEC5A
set.seed(0) # for reproducibility across group analyses
X = stress_pairs6
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
stress.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
stress.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("stress.ciout.final.a") == TRUE) {
  stress.p6.ciout.final <- stress.ciout.final.a
} else {
  stress.p6.ciout.final <- stress.ciout.final.b
}


# Applying Super Learner to "Rheumatoid Arthritis" data
Y = rheumatoid_outcome
levels(Y) = list("0"="CON", "1"="DIS")
Y = as.vector(Y)
Y = as.numeric(Y) # numeric vector of our zero-one outcomes


# TWO GENES IN THE SL: ACSL1, AQP9 
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs1.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p1.ciout.final.1 <- rheum.ciout.final.a
} else {
  rheum.p1.ciout.final.1 <- rheum.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs1.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p1.ciout.final.2 <- rheum.ciout.final.a
} else {
  rheum.p1.ciout.final.2 <- rheum.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs1.v3
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p1.ciout.final.3 <- rheum.ciout.final.a
} else {
  rheum.p1.ciout.final.3 <- rheum.ciout.final.b
}


# TWO GENES IN THE SL: NFKB1, IFNB1 
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p2.ciout.final <- rheum.ciout.final.a
} else {
  rheum.p2.ciout.final <- rheum.ciout.final.b
}


# TWO GENES IN THE SL: PRG2,ACSL1 
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs3.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p3.ciout.final.1 <- rheum.ciout.final.a
} else {
  rheum.p3.ciout.final.1 <- rheum.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs3.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p3.ciout.final.2 <- rheum.ciout.final.a
} else {
  rheum.p3.ciout.final.2 <- rheum.ciout.final.b
}


# TWO GENES IN THE SL: PRG2, CLEC5A 
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs4.v1
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p4.ciout.final.1 <- rheum.ciout.final.a
} else {
  rheum.p4.ciout.final.1 <- rheum.ciout.final.b
}

set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs4.v2
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p4.ciout.final.2 <- rheum.ciout.final.a
} else {
  rheum.p4.ciout.final.2 <- rheum.ciout.final.b
}


# TWO GENES IN THE SL: NFKB1, CLEC5A 
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs5
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p5.ciout.final <- rheum.ciout.final.a
} else {
  rheum.p5.ciout.final <- rheum.ciout.final.b
}


# TWO GENES IN THE SL: ACSL1, CLEC5A
set.seed(0) # for reproducibility across group analyses
X = rheumatoid_pairs6
fit.pairs.SL <- CV.SuperLearner(Y, X, family = binomial(), 
                                SL.library = SL.lib, V = CV_folds)
fld = fit.pairs.SL$fold
predsY.final = fit.pairs.SL$SL.predict

fold = rep(NA, nrow(X))
for(k in 1:CV_folds) {
  ii = unlist(fld[k])
  fold[ii] = k
}
rheum.ciout.final.a = ci.cvAUC(predsY.final, Y, folds = fold)
rheum.ciout.final.b = cvAUC(predsY.final, Y, folds = fold)
if (exists("rheum.ciout.final.a") == TRUE) {
  rheum.p6.ciout.final <- rheum.ciout.final.a
} else {
  rheum.p6.ciout.final <- rheum.ciout.final.b
}


stress.AUC <- rbind(max(stress.p1.ciout.final.1$cvAUC,stress.p1.ciout.final.2$cvAUC),
                    max(stress.p3.ciout.final.1$cvAUC,stress.p3.ciout.final.2$cvAUC),
                    max(stress.p4.ciout.final.1$cvAUC,stress.p4.ciout.final.2$cvAUC),
                    stress.p2.ciout.final$cvAUC,stress.p5.ciout.final$cvAUC,
                    stress.p6.ciout.final$cvAUC)

rheum.AUC <- rbind(max(rheum.p1.ciout.final.1$cvAUC,rheum.p1.ciout.final.2$cvAUC,rheum.p1.ciout.final.3$cvAUC),
                   max(rheum.p3.ciout.final.1$cvAUC,rheum.p3.ciout.final.2$cvAUC),
                   max(rheum.p4.ciout.final.1$cvAUC,rheum.p4.ciout.final.2$cvAUC),
                   rheum.p2.ciout.final$cvAUC,rheum.p5.ciout.final$cvAUC,
                   rheum.p6.ciout.final$cvAUC)

AUC.table <- cbind(stress.AUC,rheum.AUC)

rm(list= ls()[!(ls() %in% c('AUC.table'))])

colnames(AUC.table) <- c('stress_cvAUC','RA_cvAUC')
rownames(AUC.table) <- c('ACSL1 & AQP9','NFKB1 & IFNB1','PRG2 & ACSL1',
                         'PRG2 & CLEC5A','NFKB1 & CLEC5A','ACSL1 & CLEC5A')
(AUC.table <- as.data.frame(AUC.table))
# xtable(AUC.table)


#EndScript
