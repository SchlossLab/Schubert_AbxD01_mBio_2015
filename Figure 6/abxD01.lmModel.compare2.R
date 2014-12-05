library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.19otus.rfnegpos.csv", header=T)
td<-topdose[,-1]
attach(td)
ids<-names(td)
ids = ids[-1]

#this will make 30 models--the 3 best for each model with 1 to 10 variables
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


#which is the best model? using method=exhaustive, based on adjusted R^2
summary.leaps<-summary(leaps)
which.max(summary.leaps$adjr2) #outputted 28, which is model number
summary.leaps$which[28,]

lmB<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)
summary(lmB)


####################
##testing the model with 3 OTUs against the new titration
lm3<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020, data=td)
summary(lm3)
p <- 3  #number of parameters

newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
n <- dim(newtit)[1] #number of samples

predictlm3 <- as.data.frame(predict.lm(lm3, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlm3)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlm3)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)


####################
##testing the model with 5 OTUs against the new titration
lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)
p <- 5  #number of parameters

newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
n <- dim(newtit)[1] #number of samples

predictlm5 <- as.data.frame(predict.lm(lm5, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlm5)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlm5)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)


######################
##Compare the lm3 to lm5, can use the anova() function because these two models are nested
anova(lm3, lm5)


####################
##testing the model with 5 OTUs against the delayed data
lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)
p <- 5  #number of parameters

delay2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.delay.day0.logtrans.filter16mintotal.19otus.csv", header=T)
delay2<-delay2[,-1] #remove sample names
actual.2<-as.data.frame(delay2[,1])
delay<-delay2[,-1]
n <- dim(delay)[1] #number of samples

predictlm5.2 <- as.data.frame(predict.lm(lm5, newdata=delay))

ybar = colMeans(actual.2)[1]
SStot = sum((actual.2-ybar)^2)
SSres = sum((actual.2-predictlm5.2)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual.2, predictlm5.2)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)


####################
##testing the model with 3 OTUs against the delayed data
lm3<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020, data=td)
summary(lm3)
p <- 3  #number of parameters

delay2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.delay.day0.logtrans.filter16mintotal.19otus.csv", header=T)
delay2<-delay2[,-1] #remove sample names
actual<-as.data.frame(delay2[,1])
delay<-delay2[,-1]
n <- dim(delay)[1] #number of samples

predictlm3delay <- as.data.frame(predict.lm(lm3, newdata=delay))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlm3delay)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlm3delay)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)


####################
##testing the best model against the delay data
lmB<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)
p <- 10  #number of parameters

delay2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.delay.logtrans.filter16mintotal.19otus.csv", header=T)
delay2<-delay2[,-1] #remove sample names
actual<-as.data.frame(delay2[,1])
delay<-delay2[,-1]
n <- dim(delay)[1] #number of samples

predictlmBdelay <- as.data.frame(predict.lm(lmB, newdata=delay))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlmBdelay)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlmBdelay)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)


####################
##testing the best model against the titration data
lmB<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)
p <- 10  #number of parameters

newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
n <- dim(newtit)[1] #number of samples

predictlmBtit <- as.data.frame(predict.lm(lmB, newdata=newtit))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlmBtit)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlmBtit)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)

act_cfu <- 10^(res$actual)
pred_cfu <-  10^(res$predict)
res_cfu<-as.data.frame(cbind(act_cfu, pred_cfu))
names(res_cfu) <- c( "actual", "predict")
plot(res_cfu$actual, res_cfu$predict, log="xy")


####################
##testing the best model against the METRO delay
lmB<-lm(nextDayCFU ~ Otu00002 + Otu00006 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00120 + Otu00283, data=td)
p <- 10  #number of parameters

mdelay2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.metrodelay.logtrans.filter16mintotal.19otus.csv", header=T)
mdelay2<-mdelay2[,-1] #remove sample names
actual<-as.data.frame(mdelay2[,1])
mdelay<-mdelay2[,-1]
n <- dim(mdelay)[1] #number of samples

predictlmBmdel <- as.data.frame(predict.lm(lmB, newdata=mdelay))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlmBmdel)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlmBmdel)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)

####################
##testing the 5 OTU model against the METRO delay
lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)
p <- 5  #number of parameters

mdelay2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.metrodelay.day0.logtrans.filter16mintotal.19otus.csv", header=T)
mdelay2<-mdelay2[,-1] #remove sample names
actual<-as.data.frame(mdelay2[,1])
mdelay<-mdelay2[,-1]
n <- dim(mdelay)[1] #number of samples

predictlm5mdel <- as.data.frame(predict.lm(lm5, newdata=mdelay))

ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)
SSres = sum((actual-predictlm5mdel)^2)
rsq = 1-(SSres/SStot)
rsq

#adjusted r^2
numer <- ((1-rsq)*(n-1))
denom <- (n-p-1)
adjr2 <- (1 - numer/denom)
adjr2

res<-cbind(actual, predictlm5mdel)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)




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


