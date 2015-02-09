############################
# Rf models 2/8/15
############################
library(randomForest)

topdose.avg0.001<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.avg0.001.csv", header=T)
topdose.avg0.001<-topdose.avg0.001[,-1]

#fit the randomforest model
topdose.avg0.001.rf <- randomForest(nextDayCFU~., 
                               data = topdose.avg0.001,  outscale=TRUE,
                               importance=TRUE, proximity=TRUE,
                               keep.forest=TRUE, ntree=5000
)
plot(topdose.avg0.001.rf)
print(topdose.avg0.001.rf)


#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(topdose.avg0.001.rf, type=1)
imp<-importance(topdose.avg0.001.rf)
most.imp<-rownames(imp[order(imp[, "%IncMSE"], decreasing=T)[1:20],])

#write.table(imp, file="~/Desktop/mothur/abxD01/rf/RF.regression.topdose2.filter.filter.txt", sep="\t", row.names=TRUE)

topdose.rf.titr <- rf.modelNewData(RFmodel = topdose.avg0.001.rf, descr = "Titration Data", newDataFile = "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv")

topdose.rf.delay <- rf.modelNewData(RFmodel = topdose.avg0.001.rf, descr = "Delay Data", newDataFile = "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv")



############################
# Uncurated Approach 2/6/15
############################
library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.avg0.001.csv", header=T)

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

leaps.topdose.uncurated.20150208 <- leaps.build

leaps.plots(leaps.build, 4, 10)
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

leaps.plots(leaps.build, 4, 10)
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


################################################################################################
# An attempt to build - revise - revise model using topdose - titration - delay data sets. 
# 2/5/15

library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.18otus.rfnegpos.csv", header=T)
topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.delay.logtrans.filter16mintot.select.csv", header=T)
#topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal2X.topdose2.shared.no6.csv", header=T)

actual <-NULL
actual<-as.data.frame(topdose[,2]) #save the actual results in new df
row.names(actual) <- topdose[,1] #save the group names
td<-topdose[,-1] 
attach(td)
#detach(td)
ids<-names(td)
ids = ids[-1] #ids of OTUs in topdose
#names(td)[2:(length(td))] <- paste0("OTU", substr(ids, 6, 8))

#inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
inGroup <- c(NULL)
outGroup <- c(NULL)
leaps.build<-regsubsets(nextDayCFU ~ ., data=td, nbest=5, nvmax=10, force.in=inGroup, force.out=outGroup)
plot(leaps.build, scale="adjr2", main="leaps.build")
plot(leaps.build, scale="bic", main="leaps.build")
plot(leaps.build, scale="Cp", main="leaps.build")

library(car)

abbrNames <- substr(leaps.build$xnames, 6, 8)
abbrNames <- abbrNames[-1]

minSubsetSize <- 4
subsets(leaps.build, names=abbrNames, statistic="adjr2", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12))
subsets(leaps.build, names=abbrNames, statistic="bic", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12))
subsets(leaps.build, names=abbrNames, statistic="cp", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12))
abline(h=c(1:11), lwd=1)
abline(v=c(5:11), lwd=1)
abline(0, 1, col="red") #for mallow's Cp "good" is Cp <= p, parameters

lm_3_7_13_15_20_39<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00039, data=td)
lm_3_7_13_15_20_39.results <- lm_Analysis_Tests(lm_3_7_13_15_20_39, actual)
lm_3_7_13_15_20_39.rsqs <- RSQcomparisons(lm_3_7_13_15_20_39.results, "lm_3_7_13_15_20_39")
compiled_results <- rbind(compiled_results,lm_3_7_13_15_20_39.rsqs)

lm_3_7_13_15_20_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00015 + Otu00020 + Otu00039 + Otu00042, data=td)
lm_3_7_13_15_20_39_42.results <- lm_Analysis_Tests(lm_3_7_13_15_20_39_42, actual)
lm_3_7_13_15_20_39_42.rsqs <- RSQcomparisons(lm_3_7_13_15_20_39_42.results, "lm_3_7_13_15_20_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_13_15_20_39_42.rsqs)

lm_3_7_13_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=td)
lm_3_7_13_20_23_39_42.results <- lm_Analysis_Tests(lm_3_7_13_20_23_39_42, actual)
lm_3_7_13_20_23_39_42.rsqs <- RSQcomparisons(lm_3_7_13_20_23_39_42.results, "lm_3_7_13_20_23_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_13_20_23_39_42.rsqs)
#adding 23 helped metro


lm_3_7_13_20_21_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042, data=td)
lm_3_7_13_20_21_23_39_42.results <- lm_Analysis_Tests(lm_3_7_13_20_21_23_39_42, actual)
lm_3_7_13_20_21_23_39_42.rsqs <- RSQcomparisons(lm_3_7_13_20_21_23_39_42.results, "lm_3_7_13_20_21_23_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_13_20_21_23_39_42.rsqs)
#adding 21, barely helped metro

lm_3_7_8_13_20_21_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042, data=td)
lm_3_7_8_13_20_21_23_39_42.results <- lm_Analysis_Tests(lm_3_7_8_13_20_21_23_39_42, actual)
lm_3_7_8_13_20_21_23_39_42.rsqs <- RSQcomparisons(lm_3_7_8_13_20_21_23_39_42.results, "lm_3_7_8_13_20_21_23_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_8_13_20_21_23_39_42.rsqs)
#adding 8, barely increased


lm_3_7_8_14_13_20_21_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008  + Otu00014 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042, data=td)
lm_3_7_8_14_13_20_21_23_39_42.results <- lm_Analysis_Tests(lm_3_7_8_14_13_20_21_23_39_42, actual)
lm_3_7_8_14_13_20_21_23_39_42.rsqs <- RSQcomparisons(lm_3_7_8_14_13_20_21_23_39_42.results, "lm_3_7_8_14_13_20_21_23_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_8_14_13_20_21_23_39_42.rsqs)
#add 14, helped

lm_3_7_8_14_13_20_21_23_39_42_47<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008  + Otu00014 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042  + Otu00047, data=td)
lm_3_7_8_14_13_20_21_23_39_42_47.results <- lm_Analysis_Tests(lm_3_7_8_14_13_20_21_23_39_42_47, actual)
lm_3_7_8_14_13_20_21_23_39_42_47.rsqs <- RSQcomparisons(lm_3_7_8_14_13_20_21_23_39_42_47.results, "lm_3_7_8_14_13_20_21_23_39_42_47")
compiled_results <- rbind(compiled_results,lm_3_7_8_14_13_20_21_23_39_42_47.rsqs)
#add 47


lm_3_8_14_13_20_21_23_39_47<-lm(nextDayCFU ~ Otu00003 + Otu00008  + Otu00014 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00047, data=td)
lm_3_8_14_13_20_21_23_39_47.results <- lm_Analysis_Tests(lm_3_8_14_13_20_21_23_39_47, actual)
lm_3_8_14_13_20_21_23_39_47.rsqs <- RSQcomparisons(lm_3_8_14_13_20_21_23_39_47.results, "lm_3_8_14_13_20_21_23_39_47")
compiled_results <- rbind(compiled_results,lm_3_8_14_13_20_21_23_39_47.rsqs)
#remove 42 

lm_3_7_8_14_13_20_21_23_39_47<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008  + Otu00014 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00047, data=td)
lm_3_7_8_14_13_20_21_23_39_47.results <- lm_Analysis_Tests(lm_3_7_8_14_13_20_21_23_39_47, actual)
lm_3_7_8_14_13_20_21_23_39_47.rsqs <- RSQcomparisons(lm_3_7_8_14_13_20_21_23_39_47.results, "lm_3_7_8_14_13_20_21_23_39_47")
compiled_results <- rbind(compiled_results,lm_3_7_8_14_13_20_21_23_39_47.rsqs)
#add back 7

lm_3_7_8_14_13_15_20_21_23_39_47<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008  + Otu00014 + Otu00013 + Otu00015 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00047, data=td)
lm_3_7_8_14_13_15_20_21_23_39_47.results <- lm_Analysis_Tests(lm_3_7_8_14_13_15_20_21_23_39_47, actual)
lm_3_7_8_14_13_15_20_21_23_39_47.rsqs <- RSQcomparisons(lm_3_7_8_14_13_15_20_21_23_39_47.results, "lm_3_7_8_14_13_15_20_21_23_39_47")
compiled_results <- rbind(compiled_results,lm_3_7_8_14_13_15_20_21_23_39_47.rsqs)
#add back 15
#STOPPED AT THIS POINT


lm_3_7_8_11_14_13_20_21_23_39_42_47<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00008 + Otu00011 + Otu00014 + Otu00013 + Otu00020 + Otu00021 + Otu00023 + Otu00039 + Otu00042 + Otu00047, data=td)
lm_3_7_8_11_14_13_20_21_23_39_42_47.results <- lm_Analysis_Tests(lm_3_7_8_11_14_13_20_21_23_39_42_47, actual)
lm_3_7_8_11_14_13_20_21_23_39_42_47.rsqs <- RSQcomparisons(lm_3_7_8_11_14_13_20_21_23_39_42_47.results, "lm_3_7_8_11_14_13_20_21_23_39_42_47")
compiled_results <- rbind(compiled_results,lm_3_7_8_11_14_13_20_21_23_39_42_47.rsqs)
#add 11


# delay_results <- NULL
# lm_8_11_14_36_39 <- lm(nextDayCFU ~ Otu00008 + Otu00011 + Otu00014 + Otu00036 + Otu00039, data=td)
# lm_8_11_14_36_39.results <- lm_Analysis_Tests(lm_8_11_14_36_39, actual)
# lm_8_11_14_36_39.rsqs <- RSQcomparisons(lm_8_11_14_36_39.results, "lm_8_11_14_36_39")
# delay_results <- rbind(delay_results,lm_8_11_14_36_39.rsqs)
# 
# lm_8_11_14_36_39_47 <- lm(nextDayCFU ~ Otu00008 + Otu00011 + Otu00014 + Otu00036 + Otu00039 + Otu00047, data=td)
# lm_8_11_14_36_39_47.results <- lm_Analysis_Tests(lm_8_11_14_36_39_47, actual)
# lm_8_11_14_36_39_47.rsqs <- RSQcomparisons(lm_8_11_14_36_39_47.results, "lm_8_11_14_36_39_47")
# delay_results <- rbind(delay_results,lm_8_11_14_36_39_47.rsqs)
# 
# lm_8_11_14_47 <- lm(nextDayCFU ~ Otu00008 + Otu00011 + Otu00014 + Otu00047, data=td)
# lm_8_11_14_47.results <- lm_Analysis_Tests(lm_8_11_14_47, actual)
# lm_8_11_14_47.rsqs <- RSQcomparisons(lm_8_11_14_47.results, "lm_8_11_14_47")
# delay_results <- rbind(delay_results,lm_8_11_14_47.rsqs)


lm_3_7_13_20_23_36_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00013 + Otu00020 + Otu00023 + Otu00036 + Otu00039 + Otu00042, data=td)
lm_3_7_13_20_23_36_39_42.results <- lm_Analysis_Tests(lm_3_7_13_20_23_36_39_42, actual)
lm_3_7_13_20_23_36_39_42.rsqs <- RSQcomparisons(lm_3_7_13_20_23_36_39_42.results, "lm_3_7_13_20_23_36_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_13_20_23_36_39_42.rsqs)
#added 36, dropped metro

lm_3_7_13_11_20_23_39_42<-lm(nextDayCFU ~ Otu00003 + Otu00007 + Otu00011 + Otu00013 + Otu00020 + Otu00023 + Otu00039 + Otu00042, data=td)
lm_3_7_13_11_20_23_39_42.results <- lm_Analysis_Tests(lm_3_7_13_11_20_23_39_42, actual)
lm_3_7_13_11_20_23_39_42.rsqs <- RSQcomparisons(lm_3_7_13_11_20_23_39_42.results, "lm_3_7_13_11_20_23_39_42")
compiled_results <- rbind(compiled_results,lm_3_7_13_11_20_23_39_42.rsqs)
#adding 11 didn't help

lm_3_7_13_15_20_39_120<-lm(nextDayCFU ~ Otu00007 + Otu00020 + Otu00039 + Otu00015 + Otu00003 + Otu00013 + Otu00120, data=td)
lm_3_7_13_15_20_39_120.results <- lm_Analysis_Tests(lm_3_7_13_15_20_39_120, actual)
lm_3_7_13_15_20_39_120.rsqs <- RSQcomparisons(lm_3_7_13_15_20_39_120.results, "lm_3_7_13_15_20_39_120")
compiled_results <- rbind(compiled_results,lm_3_7_13_15_20_39_120.rsqs)

del_8_11_14_36_47<-lm(nextDayCFU ~ Otu00008 + Otu00011 + Otu00014 + Otu00036 + Otu00047, data=td)
del_8_11_14_36_47.results <- lm_Analysis_Tests(del_8_11_14_36_47, actual)
del_8_11_14_36_47.rsqs <- RSQcomparisons(del_8_11_14_36_47.results, "del_8_11_14_36_47")
compiled_results <- rbind(compiled_results,del_8_11_14_36_47.rsqs)


compiled_results
detach(td)






###################################################################
#Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00053 + Otu00431
#which(ids == "Otu00003")
inGroup <- c(which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"), which(ids == "Otu00003"), which(ids == "Otu00015"), which(ids == "Otu00053"))
outGroup <- c(which(ids == "Otu00006"))
leaps.3_7_20_39_13_15_53.no6<-regsubsets(nextDayCFU ~ ., data=td, nbest=3, nvmax=10, force.in=inGroup, force.out=outGroup)
plot(leaps.3_7_20_39_13_15_53.no6, scale="adjr2", main="leaps.3_7_20_39_13_15_53.no6")
plot(leaps.3_7_20_39_13_15_53.no6, scale="bic", main="leaps.3_7_20_39_13_15_53.no6")
#leaps.13, force in otu13
#leaps.3 force in otu3
#leaps.3.no6
#leaps.3.no7
#leaps.no6
#leaps.no7
#leaps.3_7_15.no6
#leaps.3_7_15.no6_13
#leaps.3_7_20_39.no6
#leaps.3_7_20_39_15.no6_13
#leaps.3_7_20_39_13.no6_15
#leaps.7_20_39_13_15.no6
#leaps.3_7_20_39_13_15.no6
#leaps.7_20_39_13.no6_15
#leaps.7_20_39_13.no6
#leaps.3_7_20_39_13_15_53.no6
#shows the best models and which variables to include
summary(leaps)
#shows the value for each parameter in the model along with intercept
vcov(leaps, 50) ##the second parameter corresponds to model number, so 50 is the max/last model
#coef
coef(leaps, 1:2)

#Linear Models
lm_7_3_20_39_15<-lm(nextDayCFU ~ Otu00007 + Otu00003 + Otu00020 + Otu00039 + Otu00015, data=td)
lm_7_3_20_39_15.results <- lm_Analysis_Tests(lm_7_3_20_39_15, actual)
lm_7_3_20_39_15.rsqs <- RSQcomparisons(lm_7_3_20_39_15.results, "lm_7_3_20_39_15")
compiled_results <- lm_7_3_20_39_15.rsqs

lm_7_3_20_39_13<-lm(nextDayCFU ~ Otu00007 + Otu00003 + Otu00020 + Otu00039 + Otu00013, data=td)
lm_7_3_20_39_13.results <- lm_Analysis_Tests(lm_7_3_20_39_13, actual)
lm_7_3_20_39_13.rsqs <- RSQcomparisons(lm_7_3_20_39_13.results, "lm_7_3_20_39_13")
compiled_results <- rbind(compiled_results,lm_7_3_20_39_13.rsqs)

lm_7_3_20_39<-lm(nextDayCFU ~ Otu00007 + Otu00003 + Otu00020 + Otu00039, data=td)
lm_7_3_20_39.results <- lm_Analysis_Tests(lm_7_3_20_39, actual)
lm_7_3_20_39.rsqs <- RSQcomparisons(lm_7_3_20_39.results, "lm_7_3_20_39")
compiled_results <- rbind(compiled_results,lm_7_3_20_39.rsqs)

leaps.7_20_39_13.no6_15<-lm(nextDayCFU ~ Otu00007 + Otu00013 + Otu00020 + Otu00039, data=td)
leaps.7_20_39_13.no6_15.results <- lm_Analysis_Tests(leaps.7_20_39_13.no6_15, actual)
leaps.3_7_20_39_13_15_53.no6_15.rsqs <- RSQcomparisons(leaps.7_20_39_13.no6_15.results, "leaps.7_20_39_13.no6_15")
compiled_results <- rbind(compiled_results,leaps.7_20_39_13.no6_15.rsqs)

leaps.7_20_39_13_15.no6<-lm(nextDayCFU ~ Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00015, data=td)
leaps.7_20_39_13_15.no6.results <- lm_Analysis_Tests(leaps.7_20_39_13_15.no6, actual)
leaps.7_20_39_13_15.no6.rsqs <- RSQcomparisons(leaps.7_20_39_13_15.no6.results, "leaps.7_20_39_13_15.no6")
compiled_results <- rbind(compiled_results,leaps.7_20_39_13_15.no6.rsqs)

leaps.3_7_20_39_13_15.no6<-lm(nextDayCFU ~ Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00015 + Otu00003, data=td)
leaps.3_7_20_39_13_15.no6.results <- lm_Analysis_Tests(leaps.3_7_20_39_13_15.no6, actual)
leaps.3_7_20_39_13_15.no6.rsqs <- RSQcomparisons(leaps.3_7_20_39_13_15.no6.results, "leaps.3_7_20_39_13_15.no6")
compiled_results <- rbind(compiled_results,leaps.3_7_20_39_13_15.no6.rsqs)

leaps.3_7_20_39_13_15_53.no6<-lm(nextDayCFU ~ Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00015 + Otu00003 + Otu00053, data=td)
#leaps.3_7_20_39_13_15_53.no6.results <- lm_Analysis_Tests(leaps.3_7_20_39_13_15_53.no6, actual)
#leaps.3_7_20_39_13_15_53.no6.rsqs <- RSQcomparisons(leaps.3_7_20_39_13_15_53.no6.results, "leaps.3_7_20_39_13_15_53.no6")
#compiled_results <- rbind(compiled_results,leaps.3_7_20_39_13_15_53.no6.rsqs)

lm_3_7_13_15_20_39_11<-lm(nextDayCFU ~ Otu00007 + Otu00013 + Otu00020 + Otu00039 + Otu00015 + Otu00003 + Otu00011, data=td)
lm_3_7_13_15_20_39_11.results <- lm_Analysis_Tests(lm_3_7_13_15_20_39_11, actual)
lm_3_7_13_15_20_39_11.rsqs <- RSQcomparisons(lm_3_7_13_15_20_39_11.results, "lm_3_7_13_15_20_39_11")
compiled_results <- rbind(compiled_results,lm_3_7_13_15_20_39_11.rsqs)

lm_7_3_20_39_15_11<-lm(nextDayCFU ~ Otu00007 + Otu00020 + Otu00039 + Otu00015 + Otu00003 + Otu00011, data=td)
lm_7_3_20_39_15_11.results <- lm_Analysis_Tests(lm_7_3_20_39_15_11, actual)
lm_7_3_20_39_15_11.rsqs <- RSQcomparisons(lm_7_3_20_39_15_11.results, "lm_7_3_20_39_15_11")
compiled_results <- rbind(compiled_results,lm_7_3_20_39_15_11.rsqs)

compiled_results
detach(td)










############################################################
#***********************************************************
##testing the model with 5 OTUs against the new titration
lm5<-lm(nextDayCFU ~ Otu00007 + Otu00003 + Otu00020 + Otu00039 + Otu00013, data=td)
summary(lm5)

titr<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.noUntr.csv", header=T)
actual<-as.data.frame(titr[,2])
row.names(actual) <- titr[,1]
titr<-titr[,-1]
titr<-titr[,-1]
n <- dim(titr)[1] #number of samples

#predict C. difficile colonization using the model on the new titration data
predictlm5 <- as.data.frame(predict.lm(lm5, newdata=titr, se.fit=TRUE))

#res<-residuals(lm5)
#plot(predict(lm5), res, xlab="fitted values", ylab="residuals", ylim=max(abs(res)) * c(-1, 1))

#pred.w.plim <-predict(lm5, titr, interval = "prediction")
#pred.w.clim <-predict(lm5, titr, interval = "confidence")
#matplot(actual, cbind(pred.w.clim, pred.w.plim[, -1]), lty = c(1, 2, 2, 3, 3), type="l", ylab = "predicted y")

results<-cbind(actual, predictlm5$fit)
names(results) <- c( "actual", "predict")
plot(results$actual, results$predict)

ybar <- apply(results, 2, mean)
num<-sum((results$actual-ybar["actual"])*(results$predict-ybar["predict"]))
denA <- sum((results$actual-ybar["actual"])^2)
denB <- sum((results$predict-ybar["predict"])^2)

rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared
print(rsq)

results[(results$actual == 0 & results$predict < 4 & results$predict > 1.7),]
resid <- (results$actual - results$predict)
results <- cbind(results, resid)
absResid <- abs(results$actual - results$predict)
results <- cbind(results, absResid)
ord_residuals <- order(-results$absResid)
ord_residuals <- results[ord_residuals,]

plot(results$predict, results$resid)
write.table(ord_residuals, file="~/Desktop/mothur/abxD01/model/topdose5Model_on_newTitration.txt", row.names=TRUE, sep="\t")

######################
##Compare the lm3 to lm5, can use the anova() function because these two models are nested
anova(lm3, lm5)


####################
##testing the model with 5 OTUs against the delayed data
lm5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00283, data=td)
summary(lm5)
p <- 5  #number of parameters

delay<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.delay.day0.logtrans.filter16mintotal.19otus.csv", header=T)
actual.delay<-as.data.frame(delay[,2])
row.names(actual.delay) <- delay[,1]
delay<-delay[,-1] #remove sample names
delay<-delay[,-1]
n <- dim(delay)[1] #number of samples

predictlm5.delay <- as.data.frame(predict.lm(lm5, newdata=delay, se.fit=TRUE))

results.delay<-cbind(actual.delay, predictlm5.delay$fit)
names(results.delay) <- c( "actual", "predict")
plot(results.delay$actual, results.delay$predict)


ybar <- apply(results.delay, 2, mean)
num<-sum((results.delay$actual-ybar["actual"])*(results.delay$predict-ybar["predict"]))
denA <- sum((results.delay$actual-ybar["actual"])^2)
denB <- sum((results.delay$predict-ybar["predict"])^2)

rsq.delay <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared
print(rsq.delay)

#results.delay[(results.delay$actual == 0 & results.delay$predict < 4 & results.delay$predict > 1.7),]
resid <- (results.delay$actual - results.delay$predict)
results.delay <- cbind(results.delay, resid)
absResid <- abs(results.delay$actual - results.delay$predict)
results.delay <- cbind(results.delay, absResid)
ord_residuals.delay <- order(-results.delay$absResid)
ord_residuals.delay <- results.delay[ord_residuals.delay,]

plot(results.delay$predict, results.delay$resid)
write.table(ord_residuals.delay, file="~/Desktop/mothur/abxD01/model/topdose5Model_on_delay.txt", row.names=TRUE, sep="\t")

pred.w.plim <-predict(lm5, delay, interval = "prediction")
pred.w.clim <-predict(lm5, delay, interval = "confidence")
matplot(actual, cbind(pred.w.clim, pred.w.plim[, -1]), lty = c(1, 2, 2, 3, 3), type="l", ylab = "predicted y")


#use this! ALYXXXXXX
res<-cbind(actual, predictlm5)
names(res) <- c( "actual", "predict")
plot(res$actual, res$predict)

ybar <- apply(res, 2, mean)
num<-sum((res$actual-ybar["actual"])*(res$predict-ybar["predict"]))
denA <- sum((res$actual-ybar["actual"])^2)
denB <- sum((res$predict-ybar["predict"])^2)

rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared

# 
# 
# #predict the outcome of the testing data
# predicted_logtrans <- predict(model_logtrans, newdata=testing_logtrans[ ,-1])
# predicted2_logtrans <- predict(model_logtrans, newdata=testing2_logtrans[ ,-1])
# 
# # what is the proportion variation explained in the outcome of the testing data?
# # i.e., what is 1-(SSerror/SStotal)
# actual_logtrans <- testing_logtrans$nextDayCFU
# actual2_logtrans <- testing2_logtrans$nextDayCFU
# 
# ybar <- apply(res, 2, mean)
# num<-sum((res$actual-ybar["actual"])*(res$predict-ybar["predict"]))
# denA <- sum((res$actual-ybar["actual"])^2)
# denB <- sum((res$predict-ybar["predict"])^2)
# 
# rsq <- 1-(num/(denA*denB))
# 
# print(rsq)
# 


# ybar = colMeans(actual)[1]
# SStot = sum((actual-ybar)^2)
# SSres = sum((actual-predictlm3delay)^2)
# rsq = 1-(SSres/SStot)
# rsq

# 
# #adjusted r^2
# numer <- ((1-rsq)*(n-1))
# denom <- (n-p-1)
# adjr2 <- (1 - numer/denom)
# adjr2
# 
# 




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

newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.noUntr.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
n <- dim(newtit)[1] #number of samples

predictlm3 <- as.data.frame(predict(lm3, newdata=newtit, se.fit=TRUE))
predictlm3 <- predict(lm3, newdata=newtit, se.fit=TRUE)
res<-residuals(lm3)
plot(predict(lm3), res, xlab="fitted values", ylab="residuals", ylim=max(abs(res)) * c(-1, 1))

pred.w.plim <-predict(lm3, newtit, interval = "prediction")
pred.w.clim <-predict(lm3, newtit, interval = "confidence")
matplot(actual, cbind(pred.w.clim, pred.w.plim[, -1]), lty = c(1, 2, 2, 3, 3), type="l", ylab = "predicted y")
RsquareAdj(predictlm3)

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

