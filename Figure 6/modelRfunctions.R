#Functions for modeling

##########################################################################
# Returns nothing. Plots ranks of models by AdjR^2, BIC, and Mallow's Cp.
# Also plots these models by subset size according to AdjR^2, BIC, and Mallow's Cp.
# 
# Input is a regsubset object from the leaps package and number for the
#   minimum number of parameters in the models to be plotted for the latter 3 plots.
leaps.plots <- function(regsubsetObj, minNumberParameters) {
  
  
  plot(regsubsetObj, scale="adjr2", main=regsubsetObj$call)
  plot(regsubsetObj, scale="bic", main=regsubsetObj$call)
  plot(regsubsetObj, scale="Cp", main=regsubsetObj$call)
  
  library(car)
  
  abbrNames <- substr(regsubsetObj$xnames, 6, 8)
  abbrNames <- abbrNames[-1]
  
  minSubsetSize <- minNumberParameters
  #maxSubsetSize <- 
  subsets(regsubsetObj, names=abbrNames, statistic="adjr2", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12), main=regsubsetObj$call)
  subsets(regsubsetObj, names=abbrNames, statistic="bic", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12), main=regsubsetObj$call)
  subsets(regsubsetObj, names=abbrNames, statistic="cp", legend=FALSE, min.size=minSubsetSize, abbrev=6, cex.subsets=.5, las=1, xlim=c(minSubsetSize, 12), main=regsubsetObj$call)
  abline(h=c(1:11), lwd=1)
  abline(v=c(minSubsetSize:11), lwd=1)
  abline(0, 1, col="red") #for mallow's Cp "good" is Cp <= p, parameters
  
  return(" You're awesome! :D ")

}


###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the titration and delay data
# Creates plots of actual values vs fitted or predicted values
# Input is a linear model and df of the actual values
RSQcomparisons <- function(model.results, model.name){
  
  #model.results <- lm_Analysis_Tests(model, actual)
  lm.adjR2 <- model.results[[2]][1, 5]
  titr.r2 <- model.results[[3]][1, 5]
  delay.r2 <- model.results[[4]][1, 5]
  df <- cbind( model.name, lm.adjR2, titr.r2, delay.r2)
  colnames(df) <- c("model", "lm.adjR2", "titr.r2", "delay.r2")
  return(df)
  
}

###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the titration and delay data
# Creates plots of actual values vs fitted or predicted values
lm_Analysis_Tests <- function(model, actualVals){
  
  #model <- lm_7_3_20_39_15
  # actualVals <- actual
  
  model.results <- lmTests(model) #returns a list
  fitted <- cbind( model.results$lm_fitted_values )
  residuals <- cbind( model.results$lm_residuals )
  rsq <- model.results$rsq_lmTopdose
  lm.results <- cbind(actualVals, fitted,  residuals, abs(residuals),  rsq)
  colnames(lm.results) <- c("actual", "predict", "resid", "absResid", "rsq")
  all.results <- list("lm_anova"=model.results$lm_anova, "results_linear"=lm.results, "results_titration"=model.results$results_titration, "results_delay"=model.results$results_delay)
  
  #plot results of linear model
  fitted <- cbind( model.results[[2]] )
  plot(actualVals[,1], fitted[,1], xlab="Actual Values", ylab="Fitted values", main=model$call)
  mtext("Top Dose Data", side=3, line=0)
  mtext(paste("Adj R^2 = ", signif(model.results[[1]], 3)), side=1, line=4)
  
  return(all.results)
}

###########################################
# Returns a list of results from fitting the linear model to the topdose data 
# and predicting the ttitration and delay data
# Creates plots for new data sets
lmTests <- function(model){
  
  model.anova <- anova(model)
  res<-residuals(model)
  lm.fitVals <- predict(model)
  #  plot(lm.fitted, res, xlab="Fitted values", main=model$call, ylab="Residuals", ylim=max(abs(res)) * c(-1, 1))
  model.summary <- summary(model)
  rsq <- model.summary$adj.r.squared
  
  model.results.titration <- modelNewData(model, "Titration Data", "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.noUntr.csv")
  #  colnames(model.results.titration) <- c(paste0(colnames(model.results.titration), "_titration"))
  
  model.results.delay <- modelNewData(model, "Delay Data", "~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delay.regression.logtrans.filter16mintot.noUntr.csv")
  #   colnames(model.results.delay) <- c(paste0(colnames(model.results.delay), "_delay"))
  # 
  #   all.results <- cbind(model.results.titration, model.results.delay, rsq)
  #   colnames(all.results)[11] <- "rsq_lmTopdose"
  #   all.results <- all.results[,c(11,1,2,3,4,5,6,7,8,9,10)]
  results_list <- list("rsq_lmTopdose"=rsq, "lm_fitted_values"=lm.fitVals, "lm_anova"=model.anova,"lm_residuals"=res, "results_titration"=model.results.titration, "results_delay"=model.results.delay)
  return(results_list)
  
}


#######################################################
# newDataFile headings=Group, nextDayCFU, OtuX...etc; descr should be a description of the new data set the anlaysis is performed on
# model should be a linear model, newDataFile = "~/Documents/blah/blah.csv", descr = "Titration Data"
# Returns a results data frame, creates the plot of actual vs predicted values
modelNewData <- function(model, descr, newDataFile ){
  
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
  plot(results$actual, results$predict, main=model$call, ylab="Predicted Values", xlab="Actual Values")
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


