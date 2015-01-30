library(leaps)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.19otus.rfnegpos.csv", header=T)
actual<-as.data.frame(topdose[,2]) #save the actual results in new df
row.names(actual) <- topdose[,1] #save the group names
td<-topdose[,-1] 
attach(td)
ids<-names(td)
ids = ids[-1] #ids of OTUs in topdose





#Otu00002 + Otu00003 + Otu00006 + Otu00007 + Otu00011 + Otu00013 + Otu00015 + Otu00019 + Otu00020 + Otu00023 + Otu00027 + Otu00029 + Otu00039 + Otu00044 + Otu00065 + Otu00078 + Otu00120 + Otu00053 + Otu00431
#this will make 30 models--the 3 best for each model with 1 to 10 variables
which(ids == "Otu00003")
inGroup <- c(which(ids == "Otu00003"), which(ids == "Otu00007"), which(ids == "Otu000020"), which(ids == "Otu00039"), which(ids == "Otu00013"))
outGroup <- c(which(ids == "Otu00006"), which(ids == "Otu00015"))
leaps.3_7_20_39_13.no6_15<-regsubsets(nextDayCFU ~ ., data=td, nbest=3, nvmax=10, force.in=inGroup, force.out=outGroup)
plot(leaps.3_7_20_39_15.no6_13, scale="adjr2", main="leaps.3_7_20_39_15.no6_13")
plot(leaps.3_7_20_39_15.no6_13, scale="bic", main="leaps.3_7_20_39_15.no6_13")
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

compiled_results
detach(td)

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

  model.results.titration <- modelNewData(model, "Titration Data", "~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.19otus.noUntr.csv")
#  colnames(model.results.titration) <- c(paste0(colnames(model.results.titration), "_titration"))
  
  model.results.delay <- modelNewData(model, "Delay Data", "~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.delay.day0.logtrans.filter16mintotal.19otus.csv")
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

