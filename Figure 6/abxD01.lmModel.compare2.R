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
