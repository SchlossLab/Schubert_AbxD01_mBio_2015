# ############################
# # Random forest models 2/18/15, using top 20 OTUs from toptitdel otus=82 (0.5% cutoff) from rf feature selection to put in a rf model
# ############################
#

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.noNewUntr.logtrans.20OTUs.csv", header=T)
topdose<-topdose[,-1] 

#fit the randomforest model
td.rf <- randomForest(nextDayCFU~., 
                      data = topdose,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(td.rf)
print(td.rf) # % Var explained: 88.5

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(td.rf, type=1)
#imp<-importance(td.rf)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)


# function to give me the prediction results on each data set
td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"

rf.predict(file=td, descr="td.rf on Topdose Data", rf.model=td.rf, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=titr, descr="td.rf on Titration Data", rf.model=td.rf, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=delay, descr="td.rf on Delay Data", rf.model=td.rf, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")



######## Now with the toptit data

toptit <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.noNewUntr.logtrans.20OTUs.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained:87.94

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")



######## Now with the toptitdel data


toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptitdel.noNewUntr.logtrans.20OTUs.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained:85.64

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptitdel, type=1)
#imp<-importance(rf.toptitdel)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.0.5p.txt", sep="\t", row.names=T, col.names=T)

rf.predict(file=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")
rf.predict(file=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel, title="RF model, toptitdel @ 0.5% cutoff: top 20 OTUs by RF")




############################
# # Random forest models 2/18/15, using OTUs at 0.5% cutoff
# ############################
#

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.noNewUntr.logtrans.0.5p.csv", header=T)
topdose<-topdose[,-1] 

#fit the randomforest model
td.rf <- randomForest(nextDayCFU~., 
                      data = topdose,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(td.rf)
print(td.rf) # % Var explained: 88.4

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(td.rf, type=1)
imp<-importance(td.rf)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.model, td, Otu00003) 


# function to give me the prediction results on each data set
td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"

rf.predict(file=td, descr="td.rf on Topdose Data", rf.model=td.rf, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=titr, descr="td.rf on Titration Data", rf.model=td.rf, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=delay, descr="td.rf on Delay Data", rf.model=td.rf, title="RF, 0.5% relabund cutoff, OTUs=82")



######## Now with the toptit data

toptit <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.noNewUntr.logtrans.0.5p.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained: 87.69

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit, title="RF, 0.5% relabund cutoff, OTUs=82")



######## Now with the toptitdel data


toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptitdel.noNewUntr.logtrans.0.5p.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained: 87.84

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptitdel, type=1)
imp<-importance(rf.toptitdel)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.0.5p.txt", sep="\t", row.names=T, col.names=T)

rf.predict(file=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel, title="RF, 0.5% relabund cutoff, OTUs=82")
rf.predict(file=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel, title="RF, 0.5% relabund cutoff, OTUs=82")





# ############################
# # Random forest models 2/18/15, using OTUs with filter16mintot
# ############################
#

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv", header=T)
topdose<-topdose[,-1] 

#fit the randomforest model
td.rf <- randomForest(nextDayCFU~., 
                      data = topdose,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(td.rf)
print(td.rf) # % Var explained: 88.4

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(td.rf, type=1)
imp<-importance(td.rf)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.model, td, Otu00003) 


# function to give me the prediction results on each data set
td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"

rf.predict(file=td, descr="td.rf on Topdose Data", rf.model=td.rf)
rf.predict(file=titr, descr="td.rf on Titration Data", rf.model=td.rf)
rf.predict(file=delay, descr="td.rf on Delay Data", rf.model=td.rf)



######## Now with the toptit data

toptit <- read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.toptit.noNewUntr.regression.logtrans.filter16mintot.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained: 87.

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit)
rf.predict(file=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit)
rf.predict(file=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit)



######## Now with the toptitdel data


toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.toptitdel.noNewUntr.regression.logtrans.filter16mintot.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained: 87.84

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptitdel, type=1)
#imp<-importance(rf.toptitdel)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.avg0.01.txt", sep="\t", row.names=T, col.names=T)

rf.predict(file=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel)
rf.predict(file=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel)
rf.predict(file=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel)


# ############################
# # Random forest models 2/16/15, using OTUs above the 1% cutoff 
# ############################
#

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.avg0.01.logtrans.csv", header=T)
topdose<-topdose[,-1] 

#fit the randomforest model
td.rf <- randomForest(nextDayCFU~., 
                      data = topdose,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(td.rf)
print(td.rf) # % Var explained: 87.68

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(td.rf, type=1)
imp<-importance(td.rf)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.model, td, Otu00003) 


# function to give me the prediction results on each data set
td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"

rf.predict(file=td, descr="td.rf on Topdose Data", rf.model=td.rf)
rf.predict(file=titr, descr="td.rf on Titration Data", rf.model=td.rf)
rf.predict(file=delay, descr="td.rf on Delay Data", rf.model=td.rf)



######## Now with the toptit data

toptit <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.avg0.01.logtrans.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained: 87.69

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
par(mfrow=c(1, 1)) 
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit)
rf.predict(file=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit)
rf.predict(file=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit)



######## Now with the toptitdel data


toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.avg0.01.logtrans.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained: 87.84

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
par(mfrow=c(1, 1)) 
varImpPlot(rf.toptitdel, type=1)
imp<-importance(rf.toptitdel)
write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel)
rf.predict(file=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel)
rf.predict(file=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel, title="RF, 1% relabund cutoff, OTUs=45")


