library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.csv", header=T)
td<-topdose[,-1]
attach(td)
ids<-names(td)
ids = ids[-1]

leaps.all<-regsubsets(nextDayCFU ~ ., data=td, nbest=3, nvmax=10, really.big=T)
#shows the best models and which variables to include
summary(leaps.all)
#shows the value for each parameter in the model along with intercept
vcov(leaps.all, 50) ##the second parameter corresponds to model number
#coef
coef(leaps.all, 1:2)
#scales=bic, Cp, and adjr2, r2
#plot(leaps.all,scale="r2")
plot(leaps.all,scale="adjr2")
plot(leaps.all,scale="bic")
