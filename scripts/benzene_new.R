
####################

# LASSO for Benzene

benzene.new <- read.table("~/Downloads/final.benzene.genes.txt", quote="\"", comment.char="")
benzene2 <- as.matrix(benzene.new)

Y <- as.matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 
                 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0,0, 0, 0, 0, 
                 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0), ncol = 1)


library(glmnet)
LASSO = glmnet(benzene2, Y, family = c("binomial"), alpha = 1, nlambda = 101, lambda = (0:100)/200, standardize = TRUE, intercept = TRUE )

plot(NULL, NULL, xlim = c(0,2), ylim = c(min(LASSO$beta), max(LASSO$beta)), xlab = "lambda", ylab = "coefficients", "main" = "Benzene LASSO")
for(i in 1:10){
  lines(rev(0:100)/10, LASSO$beta[i,],col = rainbow(10)[i])
}


# Last pair is INFB1, CLEC5A

library(SamplerCompare)
LASSO.mse.learning<- matrix(,nrow=101,ncol=2) 
colnames(LASSO.mse.learning) = c("lambda", "mse")
for (j in 1:101){
  LASSO.mse.learning[j,2] <- twonorm((Y - benzene2 %*% LASSO$beta[,(j)]))
  LASSO.mse.learning[j,1] <- 101-j #lambdastar
}
min(LASSO.mse.learning)
plot(LASSO.mse.learning, main = "Benzene LASSO MSE", col = "blue", ylab = "MSE")

# minimum lasso
LASSO.mse.learning
LASSO$beta[,48] 


lasso.result <- matrix(,nrow=length(Y),ncol = 1)
coef <- coef(LASSO)
for(i in 1:length(Y)){
  lasso.result[i] <- coef[1,48] + coef[2,48]*benzene2[i,1]
}

lasso.pair <- matrix(,nrow=length(Y),ncol = 1)
coef <- coef(LASSO)
for(i in 1:length(Y)){
  lasso.pair[i] <- coef[1,90] + coef[2,90]*benzene2[i,1] + coef[3,90]*benzene2[i,2]
} 



##### to get folds
set.seed(0)
V = 10
fit.test.final=CV.SuperLearner(Y,benzene2,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,length(predsY.final))
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}


out.INFB1 <- cvAUC(lasso.result, Y, label.ordering = NULL, folds = fold)
out.INFB1$cvAUC
# 0.8533058

out.pair <- cvAUC(lasso.pair, Y, label.ordering = NULL, folds = fold)
out.pair$cvAUC
# 0.8842975

# SuperLearner with determined pairs

V=10
set.seed(0)
Xfinal=data.frame(X=benzene2[,1])
SL.library <- c("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.randomForest", "SL.gam", "SL.mean")
fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout.final.last=ci.cvAUC(predsY.final, Y, folds = fold)
ciout.final.last$cvAUC
# 0.9138095


V=10
set.seed(0)
Xfinal=data.frame(X=benzene2[,c(1,2)])
SL.library <- c("SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.randomForest", "SL.gam", "SL.mean")
fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}



ciout.final.pair=ci.cvAUC(predsY.final, Y, folds = fold)
ciout.final.pair$cvAUC
# 0.9113095

# INFB1 had really high differential expression

