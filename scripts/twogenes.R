
SL.library <- c("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.randomForest", 
                "SL.gam", "SL.mean")

library(ROCR)
library(cvAUC)
library(SuperLearner)

V=10
Xfinal=data.frame(X=normal[,c(18,23)]) # this should be an nx2 data frame where the two columns
# are the expressions/counts of the 2 genes
Y=Benzene3
levels(Y)=list("0"="Ctr", "1"="<1ppm")
Y=as.vector(Y)
Y=as.numeric(Y)# numeric vector of our zero-one outcomes

SL.library <- c("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.randomForest", "SL.gam", "SL.mean")
fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout.final=ci.cvAUC(predsY.final, Y, folds = fold)
ciout.final

# graph the ROC curve
txt=paste("AUC = ",round(ciout$cvAUC,2),",  95% CI = ",round(ciout$ci[1],2),"-",round(ciout$ci[2],2),sep="")
pred <- prediction(predsY.final,Y)
perf1 <- performance(pred, "sens", "spec")
plot(1-slot(perf1,"x.values")[[1]],slot(perf1,"y.values")[[1]],type="s",xlab="1-Specificity",ylab="Sensitivity",
     main="ROC Curve, Normalized Counts")
text(0.6,0.4,txt)
abline(0,1)