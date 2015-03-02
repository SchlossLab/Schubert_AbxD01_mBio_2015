



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

perc<-1
otus <- otusAtPercAbund(perc, meanOTU)
dataFile <- "~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.toptitdel.noNewUntr.logtrans.filter16mintot.grouped.csv"
rf.results<-RF.analysis(dataFile=dataFile, 
            otus = otus, plot.title = paste0("RF model, ", perc, "% cutoff relabund otus=", length(otus)))
rf.results.perc1  <- rf.results
write.table(rf.results.perc1, file="~/Desktop/mothur/abxD01/rf/rfmodels.results.1percent.txt", sep="\t", row.names=T, col.names=T)
# rf.results.perc1
# rf.results.perc1.5
# rf.results.perc3, has otu3, but not otu 39
# rf.results.perc.5
# rf.results.perc5, doesn't contain otu3 or otu39
# rf.results.perc2, has both otu3 and 39

all <- read.csv(file=dataFile, header=T)

cols<- NULL
for(i in 1:(length(otus))){
  cols <-c(cols, which(names(all)==otus[i]))
}
cols <- c(which(names(all)=="Group"), which(names(all)=="nextDayCFU"), cols)

topdose <- all[1:99,cols]
titr <- all[100:153,cols]
del <- all[154:179, cols]
toptit <- all[1:153,cols]
titdel <- all[c(100:153, 154:179),cols]
topdel <- all[c(1:99, 154:179),cols]
toptitdel <- all[,cols]

tdOTUs <- c("Otu00007", "Otu00003", "Otu00015")
titrOTUs <- c("Otu00003", "Otu00013", "Otu00020")
delOTUs <- c("Otu00003", "Otu00039")

allcombos <- function(otuPool, data, dataName){
  resultsAll <- data.frame(matrix(ncol=5, nrow=0))
  colnames(resultsAll) <- c("model", "PercVarExpl", "rsq_td", "rsq_titr", "rsq_delay") 
  for(i in 1:(length(otuPool))){
    combos <- combn(otuPool, i)
    results <- data.frame(matrix(ncol=5, nrow=dim(combos)[2]))
    colnames(results) <- c("model", "PercVarExpl", "rsq_td", "rsq_titr", "rsq_delay")
    for(j in 1:dim(combos)[2]){
      formula <- makeFormula(combos[,j])
      subsets.rf <- randomForest(formula, 
                                 data = data,  outscale=TRUE,
                                 importance=TRUE, proximity=TRUE,
                                 keep.forest=TRUE, ntree=5000
      )
      #varImpPlot(subsets.rf, type=1)
      results[j,1] <- paste0(combos[,j], collapse="_")
      results[j,2] <- signif(subsets.rf$rsq[length(subsets.rf$rsq)],3)
      results[j,3] <- signif(rf.predict(data=topdose, descr=paste0(results[j,1]," model on Topdose Data"), rf.model=subsets.rf, title=paste0("RF model trained on ", dataName)),3) 
      results[j,4] <- signif(rf.predict(data=titr, descr=paste0(results[j,1]," model on Titration Data"), rf.model=subsets.rf, title=paste0("RF model trained on ", dataName)),3)
      results[j,5] <- signif(rf.predict(data=del, descr=paste0(results[j,1]," model on Delay Data"), rf.model=subsets.rf, title=paste0("RF model trained on ", dataName)),3)
    }
    resultsAll <- rbind(resultsAll, results)  
  }
  return(resultsAll)
}

topdose.topSubsets <- allcombos(tdOTUs, topdose, "Topdose")
titration.topSubsets <- allcombos(titrOTUs, titr, "Titration")
delay.topSubsets <- allcombos(delOTUs, del, "Delay")

topOTUs <- c("Otu00003", "Otu00039", "Otu00013", "Otu00017", "Otu00020", "Otu00023")
poop2 <- allcombos(topOTUs, toptitdel, "All Data")
poop2[order(poop2$rsq_delay, decreasing = TRUE),]
write.table(poop2, file= "~/Desktop/mothur/abxD01/rf/subsetToptitdel_models_results2.txt", sep="\t", row.names=FALSE)

test3.39 <- toptitdel[,c("Group","nextDayCFU","Otu00003", "Otu00039")]
RF.validate(data = test3.39, iters=100)
# 0.79946012 0.05804047for mean and sd, respectively

toptitdel.rf <- randomForest(nextDayCFU ~ ., 
                           data = toptitdel[,-1],  outscale=TRUE,
                           importance=TRUE, proximity=TRUE,
                           keep.forest=TRUE, ntree=5000
)
print(toptitdel.rf)
plot(toptitdel.rf)
varImpPlot(toptitdel.rf, type=1)
imp<-importance(toptitdel.rf)
imp1<-imp[order(imp[,1], decreasing = TRUE),]
imp2<-imp[order(imp[,2], decreasing = TRUE),]
topimp <- imp1[1:15,]
toppurity <- imp2[1:15,]
write.table(topimp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.1p.importance.txt", sep="\t", row.names=TRUE)
write.table(toppurity, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.1p.importanceNode.txt", sep="\t", row.names=TRUE)


#Plot importance Plot
topimp<-topimp[order(topimp[,1], decreasing=FALSE),]
ids <- read.csv(file = "~/Documents/Github/abxD01/Figure 6/rf.toptitdel.1p.importance.ids.csv", header = TRUE)
#ids <- ids[order]
labels<- paste0(ids$name, " (", ids$otuname, ")")
par(mfrow=c(1, 1)) #+1 to give extra labeling space
par(mar=c(5, 15, 0.5, 2) +0.1, mgp=c(3, 1, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plot(topimp[,1], 1:15, xlab="% Increase in Mean Squared Error", yaxt="n", ylab="", pch=16, cex=1.5, xlim=c(20, 80))
axis(2, at = c(1:15), labels = labels, las=1, cex.axis=.9 )
abline(h=c(1:15), lty="dashed", col='black')


#Plot importance Plot by node purity
toppurity<-toppurity[order(toppurity[,2], decreasing=FALSE),]
ids <- read.csv(file = "~/Documents/Github/abxD01/Figure 6/rf.toptitdel.1p.importanceNode.ids.csv", header = TRUE)
labels<- paste0(ids$name, " (", ids$otuname, ")")
par(mfrow=c(1, 1)) #+1 to give extra labeling space
par(mar=c(5, 15, 0.5, 2) +0.1, mgp=c(3, 1, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plot(toppurity[,2], 1:15, xlab="Increase in Node Purity", yaxt="n", ylab="", pch=16, cex=1.5, xlim=c(0, 600))
axis(2, at = c(1:15), labels = labels, las=1, cex.axis=.9 )
abline(h=c(1:15), lty="dashed", col='black')


# rf.predict(data=topdose, descr=paste0(row.names(results)[i]," model on Topdose Data"), rf.model=toptitdel.rf, title=plot.title)
# signif(rf.predict(data=titr, descr=paste0(row.names(results)[i]," model on Titration Data"), rf.model=toptitdel.rf, title=plot.title),3)
# signif(rf.predict(data=del, descr=paste0(row.names(results)[i]," model on Delay Data"), rf.model=toptitdel.rf, title=plot.title),3)

predictions <- predict(toptitdel.rf, newdata=toptitdel[,-c(1,2)])
actual <- toptitdel$nextDayCFU
expgroup<- as.character(all[,2])
results<-as.data.frame(cbind(actual, predictions))
colnames(results) <- c( "actual", "predict")

ybar <- apply(results, 2, mean)
num<-sum((results$actual-ybar["actual"])*(results$predict-ybar["predict"]))
denA <- sum((results$actual-ybar["actual"])^2)
denB <- sum((results$predict-ybar["predict"])^2)
rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared

results <- cbind(expgroup, results)


#Plot observed vs predicted
#information for plotting
color<-rainbow(7)
colors <- c(control = "black",
          cipro = color[1],
          clinda = "#FFCC00",
          vanc.625 = "#00CC00",
          vanc.3 = "#00CC00",
          vanc.1 = "#00CC00",
          strep5 = "#33FFFF",
          strep.5 = "#33FFFF",
          strep.1 = "#33FFFF",
          cef.5 = color[5],
          cef.3 = color[5],
          cef.1 = color[5],
          amp.5 = color[6],
          amp.5d = color[6],
          metro1 = color[7],
          metro1d = color[7])

# topdose is all circles, middle dose is square, 
# low dose is upward pointing triangle, delayed is starburst
pch <- c(control = 1,
         cipro = 16,
         clinda = 16,
         vanc.625 = 16,
         vanc.3 = 15,
         vanc.1 = 17,
         strep5 = 16,
         strep.5 = 15,
         strep.1 = 17,
         cef.5 = 16,
         cef.3 = 15,
         cef.1 = 17,
         amp.5 = 16,
         amp.5d = 8,
         metro1 = 16,
         metro1d = 8)

#these are needed for the legend
legend.labels <- c(control = "Control",
            cipro = "Ciprofloxacin",
            clinda = "Clindamycin",
            vanc.625 = "Vancomycin 0.625 mg/ml",
            vanc.3 = "Vancomycin 0.3 mg/ml",
            vanc.1 = "Vancomycin 0.1 mg/ml",
            strep5 = "Streptomycin 5 mg/ml",
            strep.5 = "Streptomycin 0.5 mg/ml",
            strep.1 = "Streptomycin 0.1 mg/ml",
            cef.5 = "Cefoperazone 0.5 mg/ml",
            cef.3 = "Cefoperazone 0.3 mg/ml",
            cef.1 = "Cefoperazone 0.1 mg/ml",
            amp.5 = "Ampicillin",
            amp.5d = "Ampicillin +5D",
            metro1 = "Metronidazole",
            metro1d = "Metronidazole +5D")

treatments <- c(control = "Control",
            cipro = "Ciprofloxacin",
            clinda = "Clindamycin",
            vanc.625 = "Vancomycin High",
            vanc.3 = "Vancomycin Medium",
            vanc.1 = "Vancomycin Low",
            strep5 = "Streptomycin High",
            strep.5 = "Streptomycin Medium",
            strep.1 = "Streptomycin Low",
            cef.5 = "Cefoperazone High",
            cef.3 = "Cefoperazone Medium",
            cef.1 = "Cefoperazone Low",
            amp.5 = "Ampicillin",
            amp.5d = "Ampicillin +5D",
            metro1 = "Metronidazole",
            metro1d = "Metronidazole +5D")

par(mfrow=c(1, 1)) 
par(mar=c(5, 5, 4, 2) +0.1, mgp=c(3, 1, 0), las=1) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plot(results[results[,1]=="control",2], 
     results[results[,1]=="control",3], 
     main="", 
     ylab=expression(paste("Predicted Log ", italic("C. difficile"), " CFU/g Feces")), 
     xlab=expression(paste("Actual Log ", italic("C. difficile"), " CFU/g Feces")), 
     ylim=c(0,9), xlim=c(0,9), xaxt='n', yaxt='n', cex=1.5)
#mtext(bquote("r"^"2" ~ " = " ~ .(signif(rsq, 3))), side=3, line=0)
abline(a=0, b=1, lty="dashed", lwd=2, col="black")
axis(1, at = c(0:9), labels=c(0:9))
axis(2, at = c(0:9), labels=c(0:9))


for(i in 2:(length(legend.labels))){
  points(results[results[,1]==names(legend.labels)[i],2], 
         results[results[,1]==names(legend.labels)[i],3], 
         col=colors[i], pch=pch[i], cex=1.5)
}
$legend("bottomright", inset= .07, legend=legend.labels, pch=pch, col=colors, pt.cex=1,  cex=.6, bty="n")




corrs<-read.csv(file="~/Desktop/mothur/abxD01/correlation/toptitdel.correlation.impOTUs.csv", header=TRUE)
corrs <- corrs[1:8,]
otunames <-  row.names(topimp)[15:8]
graphID <- ids[15:8,c(1,3)]

graphOTUxCD(otunames, graphID, corrs)






####################################################################################################
###Correlation for toptitdel against cdiff cfu
c<-1
otu <- c()
cor.spear <- c()
pval.spear <- c()
cor.ken <- c()
pval.ken <- c()
for(i in 3:length(toptitdel)){
  otu[c] <- colnames(toptitdel[i])
  cor.spear[c] <- cor.test(toptitdel[,2],toptitdel[,i], method="spearman")$estimate
  pval.spear[c] <- cor.test(toptitdel[,2],toptitdel[,i], method="spearman")$p.value
  cor.ken[c] <- cor.test(toptitdel[,2],toptitdel[,i], method="kendall")$estimate #good to see because kendall handles ties
  pval.ken[c] <- cor.test(toptitdel[,2],toptitdel[,i], method="kendall")$p.value #but only works if this is tao-b and not tao-a which im not sure about
  c <- c+1
}
pval.spear<-p.adjust(pval.spear, method='BH') #adjust for multiple comparisons
pval.ken<-p.adjust(pval.ken, method='BH') #adjust for multiple comparisons

results = NULL
results <- matrix(c(otu, cor.spear, pval.spear, cor.ken, pval.ken), ncol=5)
colnames(results) <- c('otu','corSpear','pvalSpear', "corKen", "pvalKen")
results <- results[order(results[,3]),]  #order by pvalue column=3
write.table(results[1:dim(results)[1],], file="~/Desktop/mothur/abxD01/correlation/toptidel.1p.correl.txt", sep="\t", row.names=FALSE)
print(results[1:dim(results)[1],])





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


