#Functions for modeling

##########################################################################
# Inputs 
# Assumes the data input has columns "nextDayCFU  Otu1  Otu2 ..."
#
# Returns 
# 
#  
RF.validate <- function(data, iters){
  
  numSamples <- dim(data)[1]
  rsqs <- NULL
  for(i in i:iters){
    trainInd <- sample(1:(numSamples), (2/3*numSamples))
    trainSet <- data[trainInd,]
    testSet <- data[-trainInd,]
    
    
    library(randomForest)
    rf.trainSet <- randomForest(nextDayCFU~., 
                                data = trainSet[,-1],  outscale=TRUE,
                                importance=TRUE, proximity=TRUE,
                                keep.forest=TRUE
    )
    
    #plot(rf.trainSet)
    #varImpPlot(rf.trainSet, type=1)
    
    plot.title <- "RF.model Validation: trained on 2/3 data, tested on 1/3"
    rsqs <- c(rsqs,rf.predict(data=testSet, descr="rf.trainSet on Test Set Data", rf.model=rf.trainSet, title=plot.title, plotgraph = FALSE))
  } #for(i in i:iters)
  
  c(mean(rsqs), sd(rsqs))
  #data.frame(mean=mean(rsqs), sd=sd(rsqs))
  
}






##########################################################################
# Inputs 
#  
# Returns 
# 
#  
RF.analysis <- function(dataFile, otus, plot.title){
  
  library(randomForest)
  
  all <- read.csv(file=dataFile, header=T)
  leng <- dim(all)[2]
  
  cols<- NULL
  for(i in 1:(length(otus))){
   cols <-c(cols, which(names(all)==otus[i]))
  }
  
  cols <- c(which(names(all)=="Group"), which(names(all)=="nextDayCFU"), cols)
  
  topdose <- all[1:99,cols]
  toptit <- all[1:154,cols]
  toptitdel <- all[,cols]
    
  
  #fit the randomforest model
  td.rf <- randomForest(nextDayCFU~., 
                        data = topdose[,-1],  outscale=TRUE,
                        importance=TRUE, proximity=TRUE,
                        keep.forest=TRUE, ntree=5000
  ) #assumes that the nexDayCFU column is right before the first Otu column
  plot(td.rf)
  
  #what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
  varImpPlot(td.rf, type=1)
  #imp<-importance(td.rf)
  #write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.topdose.avg0.01.txt", sep="\t", row.names=T, col.names=T)
  
  # function to give me the prediction results on each data set
 # td <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv"
 # titr <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv"
 # delay <- "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv"
  
  td <- all[1:99,]
  titr <- all[100:154,]
  delay <- all[154:179,]

  results <- data.frame(matrix(ncol=5, nrow=3))
  row.names(results) <- c("topdose", "toptit", "toptitdel")
  colnames(results) <- c("PercVarExpl","td", "titr", "delay", "withinValidate")

  results[1,1] <- signif(td.rf$rsq[length(td.rf$rsq)],3)
  results[1,2] <- signif(rf.predict(data=td, descr="td.rf on Topdose Data", rf.model=td.rf, title=plot.title),3)
  results[1,3] <- signif(rf.predict(data=titr, descr="td.rf on Titration Data", rf.model=td.rf, title=plot.title),3)
  results[1,4] <- signif(rf.predict(data=delay, descr="td.rf on Delay Data", rf.model=td.rf, title=plot.title),3)
  results[1,5] <- signif(RF.validate(topdose, 100)[1],3)
 
 
  ######## Now with the toptit data
  #fit the randomforest model
  rf.toptit <- randomForest(nextDayCFU~., 
                            data = toptit[,-1],  outscale=TRUE,
                            importance=TRUE, proximity=TRUE,
                            keep.forest=TRUE, ntree=5000
  )
  plot(rf.toptit)
 
  #what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
  varImpPlot(rf.toptit, type=1)
  imp<-importance(rf.toptit)
  #write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptit.avg0.01.txt", sep="\t", row.names=T, col.names=T)
    
  results[2,1] <- signif(rf.toptit$rsq[length(rf.toptit$rsq)],3)  
  results[2,2] <- signif(rf.predict(data=td, descr="rf.toptit on Topdose Data", rf.model=rf.toptit, title=plot.title),3)
  results[2,3] <- signif(rf.predict(data=titr, descr="rf.toptit on Titration Data", rf.model=rf.toptit, title=plot.title),3)
  results[2,4] <- signif(rf.predict(data=delay, descr="rf.toptit on Delay Data", rf.model=rf.toptit, title=plot.title), 3)
  results[2,5] <- signif(RF.validate(toptit, 100)[1],3)
 
  
  ######## Now with the toptitdel data
  #fit the randomforest model
  rf.toptitdel <- randomForest(nextDayCFU~., 
                               data = toptitdel[,-1],  outscale=TRUE,
                               importance=TRUE, proximity=TRUE,
                               keep.forest=TRUE, ntree=5000
  )
  plot(rf.toptitdel)
  
  #what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
  varImpPlot(rf.toptitdel, type=1)
  #imp<-importance(rf.toptitdel)
  #write.table(imp, file="~/Desktop/mothur/abxD01/rf/rf.toptitdel.0.5p.txt", sep="\t", row.names=T, col.names=T)
  
  results[3,1] <- signif(rf.toptitdel$rsq[length(rf.toptitdel$rsq)],3)
  results[3,2] <- signif(rf.predict(data=td, descr="rf.toptitdel on Topdose Data", rf.model=rf.toptitdel, title=plot.title), 3)
  results[3,3] <- signif(rf.predict(data=titr, descr="rf.toptitdel on Titration Data", rf.model=rf.toptitdel, title=plot.title), 3)
  results[3,4] <- signif(rf.predict(data=delay, descr="rf.toptitdel on Delay Data", rf.model=rf.toptitdel, title=plot.title), 3)
  results[3,5] <- signif(RF.validate(toptitdel, 100)[1],3)
  
 return(results)
}  



##########################################################################
# Inputs 
#  
# Returns 
# 
#  
makeFormula <- function(x){ #this function will be used to convert combinations of OTUs in formula objects
  as.formula(paste('nextDayCFU', paste(x, collapse=' + '), sep=' ~ '))
}

buildLM <- function(leaps.models, data, actual){ #for models of interest, after going through getModels
  data_results <- NULL
  for(i in 1:(dim(leaps.models)[1])){
    x <- colnames(leaps.models)[which(leaps.models[i,1:(dim(leaps.models)[2]-4)]==1)]
    formula <- makeFormula(x[-1])
    linmodel <- lm(formula, data=data)
    linearModels.results <- lm_Analysis_Tests(linmodel, actual)
    linearModels.rsqs <- RSQcomparisons(linearModels.results, paste(x, collapse="_"))
    data_results <- rbind(data_results, linearModels.rsqs)
    
  }
  return(data_results)
}

buildAllLM <- function(all.leaps.models, data, actual){ #for models of interest, after going through getModels
  data_results <- NULL
  for(i in 1:(dim(all.leaps.models)[1])){
    x <- colnames(all.leaps.models)[which(all.leaps.models[i,1:(dim(all.leaps.models)[2])]==TRUE)]
    formula <- makeFormula(x[-1])
    linmodel <- lm(formula, data=data)
    linearModels.results <- lm_Analysis_Tests(linmodel, actual,plot.graph = FALSE)
    linearModels.rsqs <- RSQcomparisons(linearModels.results, paste(x, collapse="_"))
    data_results <- rbind(data_results, linearModels.rsqs)
  }
  return(data_results)
}

##########################################################################
# Inputs 
# Assumes the data input comes in with "Groups  nextDayCFU  Otu..." headers 
#
#
# Returns 
# 
#  
rf.predict <- function(file=FALSE, data=NULL, descr, rf.model, title, plotgraph=TRUE) {
  if(file){
    file<-read.csv(file=file, header=T)
  } else{
    file <- data 
  }
  file<-file[,-1] #take off groups column
  
  Otu1Col <- min(which(grepl("Otu",names(file))))
  
  #predict the outcome of the testing data
  predictions <- predict(rf.model, newdata=file[ ,(Otu1Col:(dim(file)[2]))])
  
  # what is the proportion variation explained in the outcome of the testing data?
  # i.e., what is 1-(SSerror/SStotal)
  actual <- file$nextDayCFU
  
  results<-as.data.frame(cbind(actual, predictions))
  colnames(results) <- c( "actual", "predict")
  
  #Calculate the rsquared
  ybar <- apply(results, 2, mean)
  num<-sum((results$actual-ybar["actual"])*(results$predict-ybar["predict"]))
  denA <- sum((results$actual-ybar["actual"])^2)
  denB <- sum((results$predict-ybar["predict"])^2)
  rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared
  
  #plot results
  if(plotgraph==TRUE){
    plot(results$actual, results$predict, main=title, ylab="Predicted Values", xlab="Actual Values")
    mtext(paste(descr, "r^2 = ", signif(rsq, 3)), side=3, line=0)
    abline(a=0, b=1, col="red")
  }
  
  return(rsq)
  
}



##########################################################################
# Inputs 
#  
# Returns 
# 
#  

rf.modelNewData <- function(RFmodel, descr, newDataFile ){
  
  newData<-read.csv(file=newDataFile, header=T)
  actual<-as.data.frame(newData[,2]) #save the actual results in new df
  row.names(actual) <- newData[,1] #save the group names
  newData<-newData[,-1] 
  newData<-newData[,-1]
  n <- dim(newData)[1] #number of samples in n
  
  #predict C. difficile colonization using the model on the new (titration) data
  RFmodel.newData.pred <- as.data.frame(predict(RFmodel, newdata=newData))
  
  #create a results data frame and plots results
  results<-cbind(actual, RFmodel.newData.pred)
  names(results) <- c( "actual", "predict")
  
  #Calculate the rsquared
  ybar <- apply(results, 2, mean)
  num<-sum((results$actual-ybar["actual"])*(results$predict-ybar["predict"]))
  denA <- sum((results$actual-ybar["actual"])^2)
  denB <- sum((results$predict-ybar["predict"])^2)
  rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared
  
  #plot results
  plot(results$actual, results$predict, main=RFmodel$call, ylab="Predicted Values", xlab="Actual Values")
  mtext(descr, side=3, line=0)
  mtext(paste("r^2 = ", signif(rsq, 3)), side=1, line=4)
  
  #Calculate the residuals and add them to the results df, ordering them by the magnitude of the residual
  resid <- (results$actual - results$predict)
  results <- cbind(results, resid)
  absResid <- abs(results$actual - results$predict)
  results <- cbind(results, absResid)
  ord_residuals <- order(-results$absResid)
  ord_residuals <- results[ord_residuals,]
  results<-cbind(ord_residuals, rsq) #add rsq column, which will all be one number--the rsq value
  
  #  plot(results$predict, results$resid, main=model$call, xlab="Predicted Values", ylab="Residuals")
  #  mtext(descr, side=3, line=0)
  
  return(results)
}




##########################################################################
# Inputs a regsubsets object and an integer for the number of 
#   variables in the models you're searching for.
# Returns logic matrix showing which OTUs included in the models
# 
#  
getModels <- function(regsubsetsObj, paramNum){
  
  sum.regsubsetsObj <- summary(regsubsetsObj)
  topmodels.regsubsetsObj <- sum.regsubsetsObj$which #matrix of the parameters in each top model for every parameter size
  params <- which(row.names(topmodels.regsubsetsObj) == paramNum) #looking for all the models with 10 variables
  results <- topmodels.regsubsetsObj[params,]
  adjr2 <- sum.regsubsetsObj$adjr2[which(row.names(topmodels.regsubsetsObj) == paramNum)]
  bic <- sum.regsubsetsObj$bic[which(row.names(topmodels.regsubsetsObj) == paramNum)]
  cp <- sum.regsubsetsObj$cp[which(row.names(topmodels.regsubsetsObj) == paramNum)]
  results <- cbind(results, adjr2, bic, cp)
  return(results)
    
}


##########################################################################
# Returns a subset of a matrix, saved to a file
# 
# 
# Input is a file containing the matrix to be searched and a file of ids making up subset
#  
getMatrixSubset <- function(file, ids){
  
  #read in files
  matrix <- read.delim(file=file, header=T, sep="\t", row.names=1)
  ids <- read.delim(file=ids, header=F, sep="\t")
  
  #initialize subset to be filled
  matrixSub <- data.frame(matrix(ncol=dim(ids)[1], nrow=dim(ids)[1]))
  row.names(matrixSub) <- ids[,1]
  colnames(matrixSub) <- ids[,1]
  
  #go through each id in ids file and find 
  for(i in 1:dim(ids)[1]){
    rowIndex <- which(row.names(matrix)==ids[i,1])
    k <- i +1
    #test <- ids[k,1] #see if reached the end of the file
    #stopifnot(!is.na(test))

    for(k in k:dim(ids)[1])
    {
      colIndex <- which(colnames(matrix)==ids[k,1])
      cellValue <- matrix[rowIndex, colIndex]
      
      if(any(is.na(colIndex), is.na(rowIndex))){
        cellValue <- NA
      }
      
      subIndex <- which(colnames(matrixSub) == ids[k,1])
      matrixSub[i, subIndex] <- cellValue
            
    }
    
  }
  
  write.table(matrixSub, file=paste0(file, "_subset.txt"), sep="\t", row.names=T, col.names=T)
  return(matrixSub)
  
}


#getMatrixSubset("~/Desktop/mothur/abxD01/sparCC/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.all.0.03.filter16mintotal.0.03.sparcc_correlation", "~/Desktop/mothur/abxD01/ids.txt")


##########################################################################
# Returns compliments! Plots ranks of models by AdjR^2, BIC, and Mallow's Cp.
# Also plots these models by subset size according to AdjR^2, BIC, and Mallow's Cp.
# 
# Input is a regsubset object from the leaps package and number for the
#   minimum number of parameters in the models to be plotted for the latter 3 plots.
leaps.plots <- function(regsubsetObj, minNumberParameters, maxNumberParameters) {
  
  
  plot(regsubsetObj, scale="adjr2", main=regsubsetObj$call)
  plot(regsubsetObj, scale="bic", main=regsubsetObj$call)
  plot(regsubsetObj, scale="Cp", main=regsubsetObj$call)
  
  library(car)
  
  abbrNames <- substr(regsubsetObj$xnames, 6, 8)
  abbrNames <- abbrNames[-1]
  
  minSubsetSize <- minNumberParameters
  maxSubsetSize <- maxNumberParameters
  subsets(regsubsetObj, names=abbrNames, statistic="adjr2", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, maxSubsetSize), main=regsubsetObj$call)
  subsets(regsubsetObj, names=abbrNames, statistic="bic", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, maxSubsetSize), main=regsubsetObj$call)
  subsets(regsubsetObj, names=abbrNames, statistic="cp", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, maxSubsetSize), main=regsubsetObj$call)
  abline(h=c(1:20), lwd=1)
  abline(v=c(minSubsetSize:11), lwd=1)
  abline(1, 1, col="red") #for mallow's Cp "good" is Cp <= p, parameters
  
  num <- sample(1:10, 1)
  compliments <- cbind( " You're awesome! :D ", " Keep up the amazing work :) ", " Nice analysis ;) ", " You want a PhD? You better work bitch. ", " Just keep writing, just keep writing. ", " OWN IT! ", " Home stretch!!!! You can do it! ", " Be proud of yourself! You are kicking ass! ", " Smart AND beautiful :P ", " Don't let the dragon DRAG ON. ")
  
  return(compliments[num])

}


###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the titration and delay data
# Creates plots of actual values vs fitted or predicted values
# Input is a linear model and df of the actual values
RSQcomparisons <- function(model.results, model.name){
  
  #model.results <- lm_Analysis_Tests(model, actual)
  lm.adjR2 <- model.results[[2]][1, 5]
  td.r2 <- model.results[[3]][1, 5]
  titr.r2 <- model.results[[4]][1, 5]
  delay.r2 <- model.results[[5]][1, 5]
  df <- cbind( model.name, lm.adjR2, td.r2, titr.r2, delay.r2)
  colnames(df) <- c("model", "lm.adjR2", "td.r2", "titr.r2", "delay.r2")
  return(df)
  
}

###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the titration and delay data
# Creates plots of actual values vs fitted or predicted values
lm_Analysis_Tests <- function(model, actualVals, plot.graph=TRUE){
  
  model.results <- lmTests(model, plot.graph = plot.graph) #returns a list
  fitted <- cbind( model.results$lm_fitted_values )
  residuals <- cbind( model.results$lm_residuals )
  rsq <- model.results$rsq_lm
  lm.results <- cbind(actualVals, fitted,  residuals, abs(residuals),  rsq)
  colnames(lm.results) <- c("actual", "predict", "resid", "absResid", "rsq")
  all.results <- list("lm_anova"=model.results$lm_anova, "results_linear"=lm.results, "results_topdose"=model.results$results_topdose, "results_titration"=model.results$results_titration, "results_delay"=model.results$results_delay)
  
  #plot results of linear model
  fitted <- cbind( model.results[[2]] )
  if(plot.graph==TRUE){
    plot(actualVals[,1], fitted[,1], xlab="Actual Values", ylab="Fitted values", main=paste(names(model$coefficients), collapse="+"))
    mtext("Top Dose Data", side=3, line=0)
    mtext(paste("Adj R^2 = ", signif(model.results[[1]], 3)), side=1, line=4)
  }
  
  return(all.results)
}

###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the ttitration and delay data
# Creates plots for new data sets
lmTests <- function(model, plot.graph = TRUE){
  
  model.anova <- anova(model)
  res<-residuals(model)
  lm.fitVals <- predict(model)
  #  plot(lm.fitted, res, xlab="Fitted values", main=model$call, ylab="Residuals", ylim=max(abs(res)) * c(-1, 1))
  model.summary <- summary(model)
  rsq <- model.summary$adj.r.squared
  
  model.results.td <- modelNewData(model, "Topdose Data", "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv", plot.graph = plot.graph)
  
  model.results.titration <- modelNewData(model, "Titration Data", "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv", plot.graph = plot.graph)
  #  colnames(model.results.titration) <- c(paste0(colnames(model.results.titration), "_titration"))
  
  model.results.delay <- modelNewData(model, "Delay Data", "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv", plot.graph = plot.graph)
  #   colnames(model.results.delay) <- c(paste0(colnames(model.results.delay), "_delay"))
  # 
  #   all.results <- cbind(model.results.titration, model.results.delay, rsq)
  #   colnames(all.results)[11] <- "rsq_lmTopdose"
  #   all.results <- all.results[,c(11,1,2,3,4,5,6,7,8,9,10)]
  results_list <- list("rsq_lm"=rsq, "lm_fitted_values"=lm.fitVals, "lm_anova"=model.anova,"lm_residuals"=res, "results_topdose"=model.results.td,
                       "results_titration"=model.results.titration, "results_delay"=model.results.delay)
  return(results_list)
  
}


#######################################################
# newDataFile headings=Group, nextDayCFU, OtuX...etc; descr should be a description of the new data set the anlaysis is performed on
# model should be a linear model, newDataFile = "~/Documents/blah/blah.csv", descr = "Titration Data"
# Returns a results data frame, creates the plot of actual vs predicted values
modelNewData <- function(model, descr, newDataFile, plot.graph=TRUE ){
  
  titr<-read.csv(file=newDataFile, header=T)
  actual<-as.data.frame(titr[,2]) #save the actual results in new df
  row.names(actual) <- titr[,1] #save the group names
  titr<-titr[,-1] 
  titr<-titr[,-1]
  n <- dim(titr)[1] #number of samples in n
  
  #predict C. difficile colonization using the model on the new (titration) data
  predictions <- as.data.frame(predict.lm(model, newdata=titr, se.fit=TRUE))
  
  #create a results data frame and plots results
  results<-cbind(actual, predictions$fit)
  names(results) <- c( "actual", "predict")
  
  #Calculate the rsquared
  ybar <- apply(results, 2, mean)
  num<-sum((results$actual-ybar["actual"])*(results$predict-ybar["predict"]))
  denA <- sum((results$actual-ybar["actual"])^2)
  denB <- sum((results$predict-ybar["predict"])^2)
  rsq <- (num^2)/(denA*denB) #calculated from the square of the sample correlation coefficient, little r squared as opposed to big R squared
  
  #plot results
  if(plot.graph==TRUE){
    plot(results$actual, results$predict, main=paste(names(model$coefficients), collapse="+"), ylab="Predicted Values", xlab="Actual Values")
    mtext(descr, side=3, line=0)
    mtext(paste("r^2 = ", signif(rsq, 3)), side=1, line=4)
  }  

  #Calculate the residuals and add them to the results df, ordering them by the magnitude of the residual
  resid <- (results$actual - results$predict)
  results <- cbind(results, resid)
  absResid <- abs(results$actual - results$predict)
  results <- cbind(results, absResid)
  ord_residuals <- order(-results$absResid)
  ord_residuals <- results[ord_residuals,]
  results<-cbind(ord_residuals, rsq) #add rsq column, which will all be one number--the rsq value
  
  #  plot(results$predict, results$resid, main=model$call, xlab="Predicted Values", ylab="Residuals")
  #  mtext(descr, side=3, line=0)
  
  return(results)
}


