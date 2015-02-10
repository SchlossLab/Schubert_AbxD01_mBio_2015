
############################
# Uncurated Approach 2/9/15, using new shared file with cutoff of 0.005 for at least one average by experimental group.
############################
library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.avg0.005.logtrans.csv", header=T)

actual <-NULL
actual<-as.data.frame(topdose[,2]) #save the actual results in new df
row.names(actual) <- topdose[,1] #save the group names
td<-topdose[,-1] 
attach(td)
#detach(td)
ids<-names(td)
ids = ids[-1] #ids of OTUs in topdose

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(which(ids == "Otu00006"))
leaps.build<-regsubsets(nextDayCFU ~ ., data=td, nbest=5, nvmax=10, force.in=inGroup, force.out=outGroup, really.big=T)

leaps.topdose.uncurated2.20150209 <- leaps.build

leaps.plots(leaps.topdose.uncurated.20150208, 4, 10)
models8<- getModels(leaps.build, 8)
models9<- getModels(leaps.build, 9)

#Test the best model based on the leaps analysis in a linear model trained on the topdose data.
# This model is close to the best by BIC
topdose_results <- NULL
lm8.3_7_17_20_39_42_134_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00017 + Otu00020 + Otu00039 + Otu00042 + Otu00134 + Otu00181, data=td)
lm8.3_7_17_20_39_42_134_181.results <- lm_Analysis_Tests(lm8.3_7_17_20_39_42_134_181, actual)
lm8.3_7_17_20_39_42_134_181.rsqs <- RSQcomparisons(lm8.3_7_17_20_39_42_134_181.results, "lm8.3_7_17_20_39_42_134_181")
topdose_results <- rbind(topdose_results, lm8.3_7_17_20_39_42_134_181.rsqs)

lm8.3_7_17_20_39_42_147_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00017 + Otu00020 + Otu00039 + Otu00042 + Otu00147 + Otu00181, data=td)
lm8.3_7_17_20_39_42_147_181.results <- lm_Analysis_Tests(lm8.3_7_17_20_39_42_147_181, actual)
lm8.3_7_17_20_39_42_147_181.rsqs <- RSQcomparisons(lm8.3_7_17_20_39_42_147_181.results, "lm8.3_7_17_20_39_42_147_181")
topdose_results <- rbind(topdose_results, lm8.3_7_17_20_39_42_147_181.rsqs)

lm8.3_7_13_20_39_42_134_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00042 + Otu00134 + Otu00181, data=td)
lm8.3_7_13_20_39_42_134_181.results <- lm_Analysis_Tests(lm8.3_7_13_20_39_42_134_181, actual)
lm8.3_7_13_20_39_42_134_181.rsqs <- RSQcomparisons(lm8.3_7_13_20_39_42_134_181.results, "lm8.3_7_13_20_39_42_134_181")
topdose_results <- rbind(topdose_results, lm8.3_7_13_20_39_42_134_181.rsqs)

topdose_results_uncurated_20150208 <- topdose_results

detach(td)

####################### now try to use toptit data
####################### 

toptit<-read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.avg0.001.csv", header=T)

actual <-NULL
actual<-as.data.frame(toptit[,2]) #save the actual results in new df
row.names(actual) <- toptit[,1] #save the group names
toptit<-toptit[,-1] 
attach(toptit)
#detach(toptit)
ids<-names(toptit)
ids = ids[-1] #ids of OTUs in toptit

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(which(ids == "Otu00006"))
leaps.build<-regsubsets(nextDayCFU ~ ., data=toptit, nbest=5, nvmax=15, force.in=inGroup, force.out=outGroup, really.big=T)

leaps.toptit.uncurated.20150208 <- leaps.build

leaps.plots(leaps.build, 6, 15)
models10<- getModels(leaps.build, 10)
models9<- getModels(leaps.build, 9)


#Test the best model based on the leaps analysis in a linear model trained on the toptit data.
# This model is close to the best by BIC
toptit_results <- NULL
# lm10.3_4_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptit)
# lm10.3_4_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.3_4_7_13_15_20_23_39_42_181, actual)
# lm10.3_4_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.3_4_7_13_15_20_23_39_42_181.results, "lm10.3_4_7_13_15_20_23_39_42_181")
# toptit_results <- rbind(toptit_results, lm10.3_4_7_13_15_20_23_39_42_181.rsqs)
# 
# lm10.3_7_13_15_20_23_39_42_108_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00108 + Otu00181, data=toptit)
# lm10.3_7_13_15_20_23_39_42_108_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_20_23_39_42_108_181, actual)
# lm10.3_7_13_15_20_23_39_42_108_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_20_23_39_42_108_181.results, "lm10.3_7_13_15_20_23_39_42_108_181")
# toptit_results <- rbind(toptit_results, lm10.3_7_13_15_20_23_39_42_108_181.rsqs)
# 
# lm10.3_7_13_15_18_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptit)
# lm10.3_7_13_15_18_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_18_20_23_39_42_181, actual)
# lm10.3_7_13_15_18_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_18_20_23_39_42_181.results, "lm10.3_7_13_15_18_20_23_39_42_181")
# toptit_results <- rbind(toptit_results, lm10.3_7_13_15_18_20_23_39_42_181.rsqs)
# 
# # lm10.3_7_13_15_20_23_39_42_76_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00076 + Otu00181, data=toptit)
# # lm10.3_7_13_15_20_23_39_42_76_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_20_23_39_42_76_181, actual)
# # lm10.3_7_13_15_20_23_39_42_76_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_20_23_39_42_76_181.results, "lm10.3_7_13_15_20_23_39_42_76_181")
# # toptit_results <- rbind(toptit_results, lm10.3_7_13_15_20_23_39_42_76_181.rsqs)
# #gave some crazy error: t 0x10e728b50. This is a serious error. This application, or a library it uses, is using an invalid context  and is thereby contributing to an overall degradation of system stability and reliability. This notice is a courtesy: please fix this problem. It will become a fatal error in an upcoming update.
# #s: invalid context 0x10e728b50. This is a serious error. This application, or a library it uses, is using an invalid context  and is thereby contributing to an overall degradation of system stability and reliability. This notice is a courtesy: please fix this problem. It will become a fatal error in an upcoming update.
# 
# lm10.2_3_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00002 + Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptit)
# lm10.2_3_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.2_3_7_13_15_20_23_39_42_181, actual)
# lm10.2_3_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.2_3_7_13_15_20_23_39_42_181.results, "lm10.2_3_7_13_15_20_23_39_42_181")
# toptit_results <- rbind(toptit_results, lm10.2_3_7_13_15_20_23_39_42_181.rsqs)


#All of the models with 10 variables included these core 9:
lm9.3_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptit)
lm9.3_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm9.3_7_13_15_20_23_39_42_181, actual)
lm9.3_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm9.3_7_13_15_20_23_39_42_181.results, "lm9.3_7_13_15_20_23_39_42_181")
toptit_results <- rbind(toptit_results, lm9.3_7_13_15_20_23_39_42_181.rsqs)

toptit_results_uncurated_2015028 <- toptit_results

detach(toptit)

####################### now try to use toptitdel data
####################### 

toptitdel<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.toptitdel.noNewUntr.regression.logtrans.filter16mintot.csv", header=T)

actual <-NULL
actual<-as.data.frame(toptitdel[,2]) #save the actual results in new df
row.names(actual) <- toptitdel[,1] #save the group names
toptitdel<-toptitdel[,-1] 
attach(toptitdel)
#detach(toptitdel)
ids<-names(toptitdel)
ids = ids[-1] #ids of OTUs in toptitdel

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(which(ids == "Otu00006"))
leaps.toptitdel.uncurated <- regsubsets(nextDayCFU ~ ., data=toptitdel, nbest=5, nvmax=18, force.in=inGroup, force.out=outGroup, really.big=T)

leaps.toptitdel.uncurated.20150208 <- leaps.toptitdel.uncurated

leaps.plots(leaps.build, 6, 15)
models10<- getModels(leaps.build, 10)
models9<- getModels(leaps.build, 9)


#Test the best model based on the leaps analysis in a linear model trained on the toptitdel data.
# This model is close to the best by BIC
toptitdel_results <- NULL
# lm10.3_4_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptitdel)
# lm10.3_4_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.3_4_7_13_15_20_23_39_42_181, actual)
# lm10.3_4_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.3_4_7_13_15_20_23_39_42_181.results, "lm10.3_4_7_13_15_20_23_39_42_181")
# toptitdel_results <- rbind(toptitdel_results, lm10.3_4_7_13_15_20_23_39_42_181.rsqs)
# 
# lm10.3_7_13_15_20_23_39_42_108_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00108 + Otu00181, data=toptitdel)
# lm10.3_7_13_15_20_23_39_42_108_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_20_23_39_42_108_181, actual)
# lm10.3_7_13_15_20_23_39_42_108_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_20_23_39_42_108_181.results, "lm10.3_7_13_15_20_23_39_42_108_181")
# toptitdel_results <- rbind(toptitdel_results, lm10.3_7_13_15_20_23_39_42_108_181.rsqs)
# 
# lm10.3_7_13_15_18_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptitdel)
# lm10.3_7_13_15_18_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_18_20_23_39_42_181, actual)
# lm10.3_7_13_15_18_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_18_20_23_39_42_181.results, "lm10.3_7_13_15_18_20_23_39_42_181")
# toptitdel_results <- rbind(toptitdel_results, lm10.3_7_13_15_18_20_23_39_42_181.rsqs)
# 
# # lm10.3_7_13_15_20_23_39_42_76_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00076 + Otu00181, data=toptitdel)
# # lm10.3_7_13_15_20_23_39_42_76_181.results <- lm_Analysis_Tests(lm10.3_7_13_15_20_23_39_42_76_181, actual)
# # lm10.3_7_13_15_20_23_39_42_76_181.rsqs <- RSQcomparisons(lm10.3_7_13_15_20_23_39_42_76_181.results, "lm10.3_7_13_15_20_23_39_42_76_181")
# # toptitdel_results <- rbind(toptitdel_results, lm10.3_7_13_15_20_23_39_42_76_181.rsqs)
# #gave some crazy error: t 0x10e728b50. This is a serious error. This application, or a library it uses, is using an invalid context  and is thereby contributing to an overall degradation of system stability and reliability. This notice is a courtesy: please fix this problem. It will become a fatal error in an upcoming update.
# #s: invalid context 0x10e728b50. This is a serious error. This application, or a library it uses, is using an invalid context  and is thereby contributing to an overall degradation of system stability and reliability. This notice is a courtesy: please fix this problem. It will become a fatal error in an upcoming update.
# 
# lm10.2_3_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00002 + Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptitdel)
# lm10.2_3_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm10.2_3_7_13_15_20_23_39_42_181, actual)
# lm10.2_3_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm10.2_3_7_13_15_20_23_39_42_181.results, "lm10.2_3_7_13_15_20_23_39_42_181")
# toptitdel_results <- rbind(toptitdel_results, lm10.2_3_7_13_15_20_23_39_42_181.rsqs)


#All of the models with 10 variables included these core 9:
lm9.3_7_13_15_20_23_39_42_181<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00181, data=toptitdel)
lm9.3_7_13_15_20_23_39_42_181.results <- lm_Analysis_Tests(lm9.3_7_13_15_20_23_39_42_181, actual)
lm9.3_7_13_15_20_23_39_42_181.rsqs <- RSQcomparisons(lm9.3_7_13_15_20_23_39_42_181.results, "lm9.3_7_13_15_20_23_39_42_181")
toptitdel_results <- rbind(toptitdel_results, lm9.3_7_13_15_20_23_39_42_181.rsqs)

toptitdel_results_uncurated_2015028 <- toptitdel_results

detach(toptitdel)





############################
#Curated Approach 2/5/15
############################

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.avg0.001.csv", header=T)
td<-topdose[,-1] 

#fit the randomforest model
rf.model <- randomForest(nextDayCFU~., 
                         data = td,  outscale=TRUE,
                         importance=TRUE, proximity=TRUE,
                         keep.forest=TRUE, ntree=5000
)
plot(rf.model)
print(rf.model) # % Var explained: 87.73


#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
par(mfrow=c(1, 1)) 
varImpPlot(rf.model, type=1)
imp<-importance(rf.model)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.001.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.model, td, Otu00013)

############ Then used these values in combination with their correlation values with subsequent C. difficile CFU. 
############ Put those in the following candidatePool.csv:

library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.candidatePool.csv", header=T)

actual <-NULL
actual<-as.data.frame(topdose[,2]) #save the actual results in new df
row.names(actual) <- topdose[,1] #save the group names
td<-topdose[,-1] 
attach(td)
#detach(td)
ids<-names(td)
ids = ids[-1] #ids of OTUs in topdose

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(NULL)
leaps.build<-regsubsets(nextDayCFU ~ ., data=td, nbest=5, nvmax=10, force.in=inGroup, force.out=outGroup)
#leaps.build
#leaps.build.no6

leaps.topdose.curated.20150208 <- leaps.build

leaps.plots(leaps.topdose.curated.20150208, 4, 10)
models7<- getModels(leaps.build, 7)

#Test the best model based on the leaps analysis in a linear model trained on the topdose data.
# This model is close to the best by BIC, but 13 is more commonly found in top models
topdose_results <- NULL
# lm_3_7_13_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00042, data=td)
# lm_3_7_13_20_39_42.results <- lm_Analysis_Tests(lm_3_7_13_20_39_42, actual)
# lm_3_7_13_20_39_42.rsqs <- RSQcomparisons(lm_3_7_13_20_39_42.results, "lm_3_7_13_20_39_42")
# topdose_results <- rbind(topdose_results, lm_3_7_13_20_39_42.rsqs)
# 
# #This one is the best by BIC, limiting parameter size
# lm_3_7_9_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00009 + Otu00020 + Otu00039 + Otu00042, data=td)
# lm_3_7_9_20_39_42.results <- lm_Analysis_Tests(lm_3_7_9_20_39_42, actual)
# lm_3_7_9_20_39_42.rsqs <- RSQcomparisons(lm_3_7_9_20_39_42.results, "lm_3_7_9_20_39_42")
# topdose_results <- rbind(topdose_results, lm_3_7_9_20_39_42.rsqs)

# This model is also good by BIC and Cp
lm7.3_7_13_15_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00039 + Otu00042, data=td)
lm7.3_7_13_15_20_39_42.results <- lm_Analysis_Tests(lm7.3_7_13_15_20_39_42, actual)
lm7.3_7_13_15_20_39_42.rsqs <- RSQcomparisons(lm7.3_7_13_15_20_39_42.results, "lm7.3_7_13_15_20_39_42")
topdose_results <- rbind(topdose_results, lm7.3_7_13_15_20_39_42.rsqs)

# lm7.3_7_9_20_39_42_86<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00009 + Otu00020 + Otu00039 + Otu00042 + Otu00086, data=td)
# lm7.3_7_9_20_39_42_86.results <- lm_Analysis_Tests(lm7.3_7_9_20_39_42_86, actual)
# lm7.3_7_9_20_39_42_86.rsqs <- RSQcomparisons(lm7.3_7_9_20_39_42_86.results, "lm7.3_7_9_20_39_42_86")
# topdose_results <- rbind(topdose_results, lm7.3_7_9_20_39_42_86.rsqs)
# 
# lm7.3_7_9_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00009 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=td)
# lm7.3_7_9_20_23_39_42.results <- lm_Analysis_Tests(lm7.3_7_9_20_23_39_42, actual)
# lm7.3_7_9_20_23_39_42.rsqs <- RSQcomparisons(lm7.3_7_9_20_23_39_42.results, "lm7.3_7_9_20_23_39_42")
# topdose_results <- rbind(topdose_results, lm7.3_7_9_20_23_39_42.rsqs)
# 
# lm7.3_7_9_13_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00009 + Otu00020 + Otu00013 + Otu00039 + Otu00042, data=td)
# lm7.3_7_9_13_20_39_42.results <- lm_Analysis_Tests(lm7.3_7_9_13_20_39_42, actual)
# lm7.3_7_9_13_20_39_42.rsqs <- RSQcomparisons(lm7.3_7_9_13_20_39_42.results, "lm7.3_7_9_13_20_39_42")
# topdose_results <- rbind(topdose_results, lm7.3_7_9_13_20_39_42.rsqs)
# 
# lm7.3_7_9_12_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00009 + Otu00020 + Otu00012 + Otu00039 + Otu00042, data=td)
# lm7.3_7_9_12_20_39_42.results <- lm_Analysis_Tests(lm7.3_7_9_12_20_39_42, actual)
# lm7.3_7_9_12_20_39_42.rsqs <- RSQcomparisons(lm7.3_7_9_12_20_39_42.results, "lm7.3_7_9_12_20_39_42")
# topdose_results <- rbind(topdose_results, lm7.3_7_9_12_20_39_42.rsqs)

lm5.3_7_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00020 + Otu00039 + Otu00042, data=td)
lm5.3_7_20_39_42.results <- lm_Analysis_Tests(lm5.3_7_20_39_42, actual)
lm5.3_7_20_39_42.rsqs <- RSQcomparisons(lm5.3_7_20_39_42.results, "lm5.3_7_20_39_42")
topdose_results <- rbind(topdose_results, lm5.3_7_20_39_42.rsqs)

topdose_results_curated_20150208 <- topdose_results

anova(lm5.3_7_20_39_42, lm_3_7_13_20_39_42) # 0.0006315 ***
anova(lm5.3_7_20_39_42, lm_3_7_9_20_39_42) # 0.0004592 ***

anova(lm5.3_7_20_39_42, lm7.3_7_13_15_20_39_42) # 0.0003349 ***
anova(lm5.3_7_20_39_42, lm7.3_7_9_20_39_42_86) # 0.0003624 ***
anova(lm5.3_7_20_39_42, lm7.3_7_9_20_23_39_42) # 0.0004339 ***
anova(lm5.3_7_20_39_42, lm7.3_7_9_13_20_39_42) # 0.0004527 ***
anova(lm5.3_7_20_39_42, lm7.3_7_9_12_20_39_42) # 0.0006398 ***

anova(lm_3_7_13_20_39_42, lm7.3_7_13_15_20_39_42) # 0.03673 *


detach(td)

############### Now want to revise the model to incorporate information from the titration data. 
############### Combined the topdose and titration samples. 

toptit <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.avg0.001.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained: 87.85

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
par(mfrow=c(1, 1)) 
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.001.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

############ Then used these values in combination with their correlation values with subsequent C. difficile CFU. 
############ Put those in the following candidatePool.csv:

library(leaps)

toptit<-read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.avg0.001.candidatePool.csv", header=T)

actual <-NULL
actual<-as.data.frame(toptit[,2]) #save the actual results in new df
row.names(actual) <- toptit[,1] #save the group names
toptit<-toptit[,-1] 
attach(toptit)
#detach(toptit)
ids<-names(toptit)
ids = ids[-1] #ids of OTUs in topdose

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(NULL)
leaps.toptit<-regsubsets(nextDayCFU ~ ., data=toptit, nbest=5, nvmax=15, force.in=inGroup, force.out=outGroup)

leaps.toptit.curated.20150208 <- leaps.toptit

leaps.plots(leaps.toptit, 6, 12)
models10 <- getModels(leaps.toptit, 10)
models9 <- getModels(leaps.toptit, 9)
models8 <- getModels(leaps.toptit, 8)


#Test the best model based on the leaps analysis in a linear model trained on the topdose data.
toptit_results <- NULL
lm10.3_7_13_15_18_20_23_39_42_76<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00076, data=toptit)
lm10.3_7_13_15_18_20_23_39_42_76.results <- lm_Analysis_Tests(lm10.3_7_13_15_18_20_23_39_42_76, actual)
lm10.3_7_13_15_18_20_23_39_42_76.rsqs <- RSQcomparisons(lm10.3_7_13_15_18_20_23_39_42_76.results, "lm10.3_7_13_15_18_20_23_39_42_76")
toptit_results <- rbind(toptit_results, lm10.3_7_13_15_18_20_23_39_42_76.rsqs)

lm10.3_7_13_15_18_19_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00019 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm10.3_7_13_15_18_19_20_23_39_42.results <- lm_Analysis_Tests(lm10.3_7_13_15_18_19_20_23_39_42, actual)
lm10.3_7_13_15_18_19_20_23_39_42.rsqs <- RSQcomparisons(lm10.3_7_13_15_18_19_20_23_39_42.results, "lm10.3_7_13_15_18_19_20_23_39_42")
toptit_results <- rbind(toptit_results, lm10.3_7_13_15_18_19_20_23_39_42.rsqs)

lm10.3_7_11_13_15_18_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm10.3_7_11_13_15_18_20_23_39_42.results <- lm_Analysis_Tests(lm10.3_7_11_13_15_18_20_23_39_42, actual)
lm10.3_7_11_13_15_18_20_23_39_42.rsqs <- RSQcomparisons(lm10.3_7_11_13_15_18_20_23_39_42.results, "lm10.3_7_11_13_15_18_20_23_39_42")
toptit_results <- rbind(toptit_results, lm10.3_7_11_13_15_18_20_23_39_42.rsqs)

# lm10.3_4_7_13_15_20_23_39_42_76<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00076, data=toptit)
# lm10.3_4_7_13_15_20_23_39_42_76.results <- lm_Analysis_Tests(lm10.3_4_7_13_15_20_23_39_42_76, actual)
# lm10.3_4_7_13_15_20_23_39_42_76.rsqs <- RSQcomparisons(lm10.3_4_7_13_15_20_23_39_42_76.results, "lm10.3_4_7_13_15_20_23_39_42_76")
# toptit_results <- rbind(toptit_results, lm10.3_4_7_13_15_20_23_39_42_76.rsqs)

lm10.3_7_13_15_18_20_21_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm10.3_7_13_15_18_20_21_23_39_42.results <- lm_Analysis_Tests(lm10.3_7_13_15_18_20_21_23_39_42, actual)
lm10.3_7_13_15_18_20_21_23_39_42.rsqs <- RSQcomparisons(lm10.3_7_13_15_18_20_21_23_39_42.results, "lm10.3_7_13_15_18_20_21_23_39_42")
toptit_results <- rbind(toptit_results, lm10.3_7_13_15_18_20_21_23_39_42.rsqs)

lm9.3_7_13_15_18_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm9.3_7_13_15_18_20_23_39_42.results <- lm_Analysis_Tests(lm9.3_7_13_15_18_20_23_39_42, actual)
lm9.3_7_13_15_18_20_23_39_42.rsqs <- RSQcomparisons(lm9.3_7_13_15_18_20_23_39_42.results, "lm9.3_7_13_15_18_20_23_39_42")
toptit_results <- rbind(toptit_results, lm9.3_7_13_15_18_20_23_39_42.rsqs)

lm9.3_4_7_13_15_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm9.3_4_7_13_15_20_23_39_42.results <- lm_Analysis_Tests(lm9.3_4_7_13_15_20_23_39_42, actual)
lm9.3_4_7_13_15_20_23_39_42.rsqs <- RSQcomparisons(lm9.3_4_7_13_15_20_23_39_42.results, "lm9.3_4_7_13_15_20_23_39_42")
toptit_results <- rbind(toptit_results, lm9.3_4_7_13_15_20_23_39_42.rsqs)

lm9.3_7_13_15_20_23_39_42_76<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00076, data=toptit)
lm9.3_7_13_15_20_23_39_42_76.results <- lm_Analysis_Tests(lm9.3_7_13_15_20_23_39_42_76, actual)
lm9.3_7_13_15_20_23_39_42_76.rsqs <- RSQcomparisons(lm9.3_7_13_15_20_23_39_42_76.results, "lm9.3_7_13_15_20_23_39_42_76")
toptit_results <- rbind(toptit_results, lm9.3_7_13_15_20_23_39_42_76.rsqs)

lm8.3_7_13_15_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptit)
lm8.3_7_13_15_20_23_39_42.results <- lm_Analysis_Tests(lm8.3_7_13_15_20_23_39_42, actual)
lm8.3_7_13_15_20_23_39_42.rsqs <- RSQcomparisons(lm8.3_7_13_15_20_23_39_42.results, "lm8.3_7_13_15_20_23_39_42")
toptit_results <- rbind(toptit_results, lm8.3_7_13_15_20_23_39_42.rsqs)


anova(lm8.3_7_13_15_20_23_39_42, lm9.3_7_13_15_18_20_23_39_42) # 0.03106 *
anova(lm8.3_7_13_15_20_23_39_42, lm9.3_4_7_13_15_20_23_39_42) # 0.05798 .
anova(lm8.3_7_13_15_20_23_39_42, lm9.3_7_13_15_20_23_39_42_76) # 0.06302 .


anova(lm8.3_7_13_15_20_23_39_42, lm10.3_7_13_15_18_20_23_39_42_76) # 0.007971 **
anova(lm8.3_7_13_15_20_23_39_42, lm10.3_7_13_15_18_19_20_23_39_42) # 0.02669 *
anova(lm8.3_7_13_15_20_23_39_42, lm10.3_7_11_13_15_18_20_23_39_42) # 0.02709 *
anova(lm8.3_7_13_15_20_23_39_42, lm10.3_4_7_13_15_20_23_39_42_76) # 0.03208 *
anova(lm8.3_7_13_15_20_23_39_42, lm10.3_7_13_15_18_20_21_23_39_42) # 0.03775 *

toptit_results_curated_20150208 <- toptit_results

detach(toptit)

############### Now want to revise the model to incorporate information from the delay data. 
############### Combined the topdose, titration, and delay samples. 

toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.toptitdel.noNewUntr.regression.logtrans.filterAvg0.001.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained: 85.91

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
par(mfrow=c(1, 1)) 
varImpPlot(rf.toptitdel, type=1)
imp<-importance(rf.toptitdel)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.avg0.001.txt", sep="\t", row.names=T, col.names=T)

############ Then used these values in combination with their correlation values with subsequent C. difficile CFU. 
############ Put those in the following candidatePool.csv:

library(leaps)

toptitdel<-read.csv("~/Desktop/mothur/abxD01/model/shared.toptitdel.avg0.001.candidatePool.csv", header=T)

actual <-NULL
actual<-as.data.frame(toptitdel[,2]) #save the actual results in new df
row.names(actual) <- toptitdel[,1] #save the group names
toptitdel<-toptitdel[,-1] 
attach(toptitdel)
#detach(toptitdel)
ids<-names(toptitdel)
ids = ids[-1] #ids of OTUs in topdose

#Perform an exhaustive search for the best linear models with different numbers of parameters.
#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(NULL)
leaps.toptitdel<-regsubsets(nextDayCFU ~ ., data=toptitdel, nbest=5, nvmax=13, force.in=inGroup, force.out=outGroup)

leaps.toptitdel.curated.20150208 <- leaps.toptitdel

leaps.plots(leaps.toptitdel, 6, 13)
models10 <- getModels(leaps.toptitdel, 10)
models9 <- getModels(leaps.toptitdel, 9)


toptitdel_results <- NULL
lm10.3_4_7_8_13_17_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00007 + Otu00008 + Otu00013 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm10.3_4_7_8_13_17_20_23_39_42.results <- lm_Analysis_Tests(lm10.3_4_7_8_13_17_20_23_39_42, actual)
lm10.3_4_7_8_13_17_20_23_39_42.rsqs <- RSQcomparisons(lm10.3_4_7_8_13_17_20_23_39_42.results, "lm10.3_4_7_8_13_17_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm10.3_4_7_8_13_17_20_23_39_42.rsqs)

lm10.3_4_8_13_15_17_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00008 + Otu00013 + Otu00015 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm10.3_4_8_13_15_17_20_23_39_42.results <- lm_Analysis_Tests(lm10.3_4_8_13_15_17_20_23_39_42, actual)
lm10.3_4_8_13_15_17_20_23_39_42.rsqs <- RSQcomparisons(lm10.3_4_8_13_15_17_20_23_39_42.results, "lm10.3_4_8_13_15_17_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm10.3_4_8_13_15_17_20_23_39_42.rsqs)

lm10.3_4_13_15_17_20_23_39_42_65<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00013 + Otu00015 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00065, data=toptitdel)
lm10.3_4_13_15_17_20_23_39_42_65.results <- lm_Analysis_Tests(lm10.3_4_13_15_17_20_23_39_42_65, actual)
lm10.3_4_13_15_17_20_23_39_42_65.rsqs <- RSQcomparisons(lm10.3_4_13_15_17_20_23_39_42_65.results, "lm10.3_4_13_15_17_20_23_39_42_65")
toptitdel_results <- rbind(toptitdel_results, lm10.3_4_13_15_17_20_23_39_42_65.rsqs)

lm10.3_4_13_15_17_20_23_39_42_78<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00013 + Otu00015 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042 + Otu00078, data=toptitdel)
lm10.3_4_13_15_17_20_23_39_42_78.results <- lm_Analysis_Tests(lm10.3_4_13_15_17_20_23_39_42_78, actual)
lm10.3_4_13_15_17_20_23_39_42_78.rsqs <- RSQcomparisons(lm10.3_4_13_15_17_20_23_39_42_78.results, "lm10.3_4_13_15_17_20_23_39_42_78")
toptitdel_results <- rbind(toptitdel_results, lm10.3_4_13_15_17_20_23_39_42_78.rsqs)

lm10.3_4_8_13_17_18_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00008 + Otu00013 + Otu00017 + Otu00018 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm10.3_4_8_13_17_18_20_23_39_42.results <- lm_Analysis_Tests(lm10.3_4_8_13_17_18_20_23_39_42, actual)
lm10.3_4_8_13_17_18_20_23_39_42.rsqs <- RSQcomparisons(lm10.3_4_8_13_17_18_20_23_39_42.results, "lm10.3_4_8_13_17_18_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm10.3_4_8_13_17_18_20_23_39_42.rsqs)

lm9.3_4_8_13_17_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00008 + Otu00013 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm9.3_4_8_13_17_20_23_39_42.results <- lm_Analysis_Tests(lm9.3_4_8_13_17_20_23_39_42, actual)
lm9.3_4_8_13_17_20_23_39_42.rsqs <- RSQcomparisons(lm9.3_4_8_13_17_20_23_39_42.results, "lm9.3_4_8_13_17_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm9.3_4_8_13_17_20_23_39_42.rsqs)

lm9.3_4_13_15_17_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00013 + Otu00015 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm9.3_4_13_15_17_20_23_39_42.results <- lm_Analysis_Tests(lm9.3_4_13_15_17_20_23_39_42, actual)
lm9.3_4_13_15_17_20_23_39_42.rsqs <- RSQcomparisons(lm9.3_4_13_15_17_20_23_39_42.results, "lm9.3_4_13_15_17_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm9.3_4_13_15_17_20_23_39_42.rsqs)

lm8.3_4_13_17_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00004 + Otu00013 + Otu00017 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=toptitdel)
lm8.3_4_13_17_20_23_39_42.results <- lm_Analysis_Tests(lm8.3_4_13_17_20_23_39_42, actual)
lm8.3_4_13_17_20_23_39_42.rsqs <- RSQcomparisons(lm8.3_4_13_17_20_23_39_42.results, "lm8.3_4_13_17_20_23_39_42")
toptitdel_results <- rbind(toptitdel_results, lm8.3_4_13_17_20_23_39_42.rsqs)


anova(lm8.3_4_13_17_20_23_39_42, lm9.3_4_13_15_17_20_23_39_42) # 0.004704 **
anova(lm8.3_4_13_17_20_23_39_42, lm9.3_4_8_13_17_20_23_39_42) # 0.004491 **

anova(lm8.3_4_13_17_20_23_39_42, lm10.3_4_7_8_13_17_20_23_39_42) # 0.001618 **
anova(lm8.3_4_13_17_20_23_39_42, lm10.3_4_8_13_15_17_20_23_39_42) # 0.0002017 ***
anova(lm8.3_4_13_17_20_23_39_42, lm10.3_4_13_15_17_20_23_39_42_65) # 0.004751 **
anova(lm8.3_4_13_17_20_23_39_42, lm10.3_4_13_15_17_20_23_39_42_78) # 0.004979 **
anova(lm8.3_4_13_17_20_23_39_42, lm10.3_4_8_13_17_18_20_23_39_42) # 0.005226 **

anova(lm10.3_4_8_13_15_17_20_23_39_42, lm9.3_4_13_15_17_20_23_39_42) # 0.002663 **
anova(lm10.3_4_13_15_17_20_23_39_42_65, lm9.3_4_13_15_17_20_23_39_42) # 0.09881 .
anova(lm10.3_4_13_15_17_20_23_39_42_78, lm9.3_4_13_15_17_20_23_39_42) # 0.1048

anova(lm10.3_4_7_8_13_17_20_23_39_42, lm9.3_4_8_13_17_20_23_39_42) # 0.02864 *
anova(lm10.3_4_8_13_15_17_20_23_39_42, lm9.3_4_8_13_17_20_23_39_42) # 0.002787 **
anova(lm10.3_4_8_13_17_18_20_23_39_42, lm9.3_4_8_13_17_20_23_39_42) # 0.1174

toptitdel_results_curated_20150208 <- toptitdel_results

detach(toptitdel)
