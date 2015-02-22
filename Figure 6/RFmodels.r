# ############################
# # Random forest models 2/20/15, using otus=50 (0.875% cutoff) 
# ############################
#

data<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.toptitdel.noNewUntr.logtrans.filter16mintot.grouped.csv", header=T)

means <- as.data.frame(matrix(nrow=nlevels(data$expgroup), ncol=0))
row.names(means) <- levels(data$expgroup)
for(i in 3:length(data)){
  
  df <- aggregate( data[,i] ~ data$expgroup, data=data, mean )
  names(df)[2] <- names(data)[i]
  
  means <- as.data.frame(cbind( means, df[,2]))
  names(means)[dim(means)[2]] <- names(data)[i]
  
}

means <- rbind(means, max=-1)
for(j in 1:length(means)){
  means["max", j] <- max(means[,j])
}

means <- rbind(means, relabund=-1)
for(k in 1:length(means)){
  means["relabund", k] <- (means[17,k])/1625
}

means <- rbind(means, percent=-1)
for(k in 1:length(means)){
  means["percent", k] <- (means[18,k])*100
}

meanOTU <- means[,-1] #remove the nextDayCFU column

set <- NULL
set <- seq(0,3, by=0.01)
numOTUxThreshold <- NULL
numOTUxThreshold <- as.data.frame(cbind(numOTUxThreshold, set))
names(numOTUxThreshold)[1] <- "threshold"
numOTUxThreshold <- cbind(numOTUxThreshold, numOTU=NA)
for( i in 1:(length(set)) ) {
  
  test <- meanOTU[ , meanOTU[19,] >= set[i] ]
  #test2 <- sum( meanOTU[19,] >= set[i] )
  
  numOTUxThreshold[i, 2] <- length(test)
  
}

plot(numOTUxThreshold$numOTU ~ numOTUxThreshold$threshold, xlab="Avg % RelAbund Threshold", ylab="# OTUs included")
abline(v=seq(0,3, by=.1), col="light gray")
abline(h=seq(1,300, by=10), col="light gray")

otusAtPercAbund <- function(percThresh, meanOTU){
  perc <- meanOTU[ , meanOTU[19,] >= percThresh ]
  otus <- names(perc)[-(length(perc))]
  return(otus)
}

otus <- otusAtPercAbund(1, meanOTU)

rf.results<-RF.analysis(dataFile = "~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.toptitdel.noNewUntr.logtrans.filter16mintot.grouped.csv", 
            otus = otus, plot.title = paste("RF model, 1% cutoff relabund otus=", length(otus)))

rf.results.perc5 <- rf.results
# rf.results.0.875perc
# rf.results.0.9perc
# rf.results.perc1
# rf.results.perc1.5
# rf.results.perc3
# rf.results.perc.5

# ############################
# # Random forest models 2/20/15, using otus=57 (0.75% cutoff) 
# ############################
#

library(randomForest)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/shared.topdose.noNewUntr.logtrans.0.75p.csv", header=T)
topdose<-topdose[,-1] 

#fit the randomforest model
td.rf <- randomForest(nextDayCFU~., 
                      data = topdose,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(td.rf)
print(td.rf) # % Var explained: 88.06

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(td.rf, type=1)
#imp<-importance(td.rf)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)


# function to give me the prediction results on each data set
td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"

rf.predict(file=td, descr="td.rf on Topdose Data", rf.model=td.rf, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=titr, descr="td.rf on Titration Data", rf.model=td.rf, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=delay, descr="td.rf on Delay Data", rf.model=td.rf, title="RF model, 0.75% cutoff, otus=57")



######## Now with the toptit data

toptit <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptit.noNewUntr.logtrans.0.75p.csv", header=T)
toptit <- toptit[,-1]

#fit the randomforest model
rf.toptit <- randomForest(nextDayCFU~., 
                          data = toptit,  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
)
plot(rf.toptit)
print(rf.toptit)
# % Var explained: 87.81

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptit, type=1)
imp<-importance(rf.toptit)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
#partialPlot(rf.toptit, toptit, Otu00013)

rf.predict(file=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit, title="RF model, 0.75% cutoff, otus=57")



######## Now with the toptitdel data


toptitdel <- read.csv("~/Desktop/mothur/abxD01/model/shared.toptitdel.noNewUntr.logtrans.0.75p.csv", header=T)
toptitdel <- toptitdel[,-1]

#fit the randomforest model
rf.toptitdel <- randomForest(nextDayCFU~., 
                             data = toptitdel,  outscale=TRUE,
                             importance=TRUE, proximity=TRUE,
                             keep.forest=TRUE, ntree=5000
)
plot(rf.toptitdel)
print(rf.toptitdel)
# % Var explained:85.98

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(rf.toptitdel, type=1)
#imp<-importance(rf.toptitdel)
#write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.0.5p.txt", sep="\t", row.names=T, col.names=T)

rf.predict(file=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel, title="RF model, 0.75% cutoff, otus=57")
rf.predict(file=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel, title="RF model, 0.75% cutoff, otus=57")
#0.951259




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


