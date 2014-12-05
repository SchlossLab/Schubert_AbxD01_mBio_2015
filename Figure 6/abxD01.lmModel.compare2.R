library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.19otus.rfnegpos.csv", header=T)
td<-topdose[,-1]
attach(td)
ids<-names(td)
ids = ids[-1]

#this will make 50 models--the 5 best for each model with 1 to 10 variables
leaps<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10)
#shows the best models and which variables to include
summary(leaps)
#shows the value for each parameter in the model along with intercept
vcov(leaps, 50) ##the second parameter corresponds to model number, so 50 is the max/last model
#coef
coef(leaps, 1:2)
#scales=bic, Cp, and adjr2, r2
plot(leaps,scale="r2")
plot(leapsf,scale="adjr2")
plot(leaps,scale="bic")

leapsf<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10, method="forward")
leapsb<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10, method="backward")
leapss<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10, method="seqrep")

summary.leaps<-summary(leaps)
which.max(summary.leaps$adjr2) #outputted 28, which is model number
summary.leaps$which[28,]

lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)

summary.leapss<-summary(leapss)
which.max(summary.leapss$adjr2) 
summary.leapss$which[26,]

lmadjr2<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)
summary(lmadjr2)


####################
##testing the model with 3 OTUs against the new titration

lm3<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020, data=td)
summary(lm3)
newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
predictlm3 <- as.data.frame(predict.lm(lm3, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)

SSres = sum((actual-predictlm3)^2)
rsq = 1-(SSres/SStot)
rsq


####################
##testing the model with 5 OTUs against the new titration
lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)
newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
predictlm5 <- as.data.frame(predict.lm(lm5, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)

SSres = sum((actual-predictlm5)^2)
rsq = 1-(SSres/SStot)
rsq




####################
##testing the best model against the new titration
lmadjr2<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)

newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
predictlmadjr2 <- as.data.frame(predict.lm(lmadjr2, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)

SSres = sum((actual-predictlmadjr2)^2)
rsq = 1-(SSres/SStot)
rsq

####################
##try leaving out OTU7 or OTU6 because they have the highest correlation with each other
leapsno7<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10)
summary(leapsno7)
plot(leapsno7,scale="bic")
leapssno7<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00006 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10, method="seqrep")
plot(leapssno7,scale="bic")

leapsno6<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10)
summary(leapsno6)
plot(leapsno6,scale="bic")
leapssno6<-regsubsets(nextDayCFU ~ Otu00002 + Otu00003 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00283 + Otu00431, data=td, nbest=3, nvmax=10, method="seqrep")
plot(leapssno6,scale="bic")


