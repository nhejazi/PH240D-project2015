# Libraries
library(ROCR)
library(cvAUC)
library(SuperLearner)

SL.library <- c("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.randomForest", "SL.gam", "SL.mean")
V=10

#### NICKEL DATA #### 
nickel <- read.csv("~/Documents/Fall 2015/240D/newnickel.genes.txt", sep="")
nickelT <- t(nickel)

nickel.binary <- c(rep(0,10),rep(1,8))
Y=matrix(nickel.binary,nrow=18,ncol=1)



##### ACSL1, AQP9  
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(6,2)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve ## AUC 1
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


###### NFKB1, IFNB1
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(5,8)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### PRG2.2, ACSL1  
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(3,6)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### PRG2.3, ACSL1  
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(4,6)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### PRG2.2, CLEC5A
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(3,7)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### PRG2.3, CLEC5A  
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(4,7)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### NFKB1, CLEC5A 
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(5,7)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)



##### ACSL1, CLEC5A  
set.seed(77)
Xfinal=data.frame(X=nickelT[,c(6,7)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)



######## PAD DATA ########

pad <- read.csv("~/Downloads/final.pad.genes.txt", sep="")
padT <- t(pad)

pad.binary <- c(rep(1,19),rep(0,18))
Y=matrix(pad.binary,nrow=37,ncol=1)



##### ACSL1, AQP9 
set.seed(77)
Xfinal=data.frame(X=padT[,c(2,1)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve ## AUC .121
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### NFKB1, IFNB1
set.seed(77)
Xfinal=data.frame(X=padT[,c(5,10)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### PRG2, ACSL1  
set.seed(77)
Xfinal=data.frame(X=padT[,c(6,2)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### PRG2.2, ACSL1  
set.seed(77)
Xfinal=data.frame(X=padT[,c(7,2)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### PRG2.3, ACSL1  
set.seed(77)
Xfinal=data.frame(X=padT[,c(8,2)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### PRG2, CLEC5A 
set.seed(77)
Xfinal=data.frame(X=padT[,c(6,4)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### PRG2.2, CLEC5A 
set.seed(77)
Xfinal=data.frame(X=padT[,c(7,4)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)


##### PRG2.3, CLEC5A 
set.seed(77)
Xfinal=data.frame(X=padT[,c(8,4)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)

##### NFKB1, CLEC5A    
set.seed(1)
Xfinal=data.frame(X=padT[,c(5,5)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)



##### ACSL1, CLEC5A  
set.seed(77)
Xfinal=data.frame(X=padT[,c(2,4)]) # this should be an nx2 data frame where the two columns

fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout=ci.cvAUC(predsY.final, Y, folds = fold)

out.final=cvAUC(predsY.final, Y, folds = fold)
out.final$cvAUC

# graph the ROC curve
txt=paste("AUC = ",round(out.final$cvAUC,2))
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)
