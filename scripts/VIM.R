# TWO GENES IN THE SL: ACSL1, AQP9 

V=10
set.seed(0)
Xfinal=data.frame(X=normal[,c(1,2)])
fit.test.final=CV.SuperLearner(Y,Xfinal,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final$fold
predsY.final=fit.test.final$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}

ciout.final=ci.cvAUC(predsY.final, Y, folds = fold)
ciout.final # 0.92

# TWO GENES IN THE SL: NFKB1, IFNB1 

set.seed(0)
Xfinal1=data.frame(X=normal[,c(12,18)])
fit.test.final1=CV.SuperLearner(Y,Xfinal1,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final1$fold
predsY.final1=fit.test.final1$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout.final1=ci.cvAUC(predsY.final1, Y, folds = fold)
ciout.final1 # 0.939881


# TWO GENES IN THE SL: PRG2,ACSL1 

set.seed(0)
Xfinal3=data.frame(X=normal[,c(21,1)])
fit.test.final3=CV.SuperLearner(Y,Xfinal3,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final3$fold
predsY.final3=fit.test.final3$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout.final3=ci.cvAUC(predsY.final3, Y, folds = fold)
ciout.final3 # 0.920476


# TWO GENES IN THE SL: PRG2, CLEC5A 
set.seed(0)
Xfinal4=data.frame(X=normal[,c(6,21)])
fit.test.final4=CV.SuperLearner(Y,Xfinal4,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final4$fold
predsY.final4=fit.test.final4$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout.final4=ci.cvAUC(predsY.final4, Y, folds = fold)
ciout.final4 # 0.9382143

# TWO GENES IN THE SL: NFKB1, CLEC5A 
set.seed(0)
Xfinal5=data.frame(X=normal[,c(18,6)])
fit.test.final5=CV.SuperLearner(Y,Xfinal5,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final5$fold
predsY.final5=fit.test.final5$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout.final5=ci.cvAUC(predsY.final5, Y, folds = fold)
ciout.final5 #0.9057143


# TWO GENES IN THE SL: ACSL1, CLEC5A

set.seed(0)
Xfinal6=data.frame(X=normal[,c(1,6)])
fit.test.final6=CV.SuperLearner(Y,Xfinal6,family=binomial(),SL.library=SL.library,V=V)
fld=fit.test.final6$fold
predsY.final6=fit.test.final6$SL.predict
fold=rep(NA,77)
for(k in 1:V) {
  ii=unlist(fld[k])
  fold[ii]=k
}
ciout.final6=ci.cvAUC(predsY.final6, Y, folds = fold)
ciout.final6 #0.9150595



