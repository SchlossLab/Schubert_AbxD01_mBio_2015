#Functions for modeling

##########################################################################
# Inputs 
# Assumes the data input has columns "Group nextDayCFU  Otu1  Otu2 ..."
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
    
   # plot.title <- "RF.model Validation: trained on 2/3 data, tested on 1/3"
    rsqs <- c(rsqs,rf.predict(data=testSet, descr="rf.trainSet on Test Set Data", rf.model=rf.trainSet, title=plot.title, plotgraph = FALSE))
  } #for(i in i:iters)
  
  c(mean(rsqs), sd(rsqs))
  #data.frame(mean=mean(rsqs), sd=sd(rsqs))
  
}

##########################################################################
# Inputs 
# Assumes the data input has columns "nextDayCFU  Otu1  Otu2 ..."
#
# Returns 
# 




# #  
# graphOTUxCD(otunames, otulabels, corrs)
# 
# pretty_otus <- as.data.frame(gsub("Otu0*", "OTU", otus))
# otulabels <- cbind(ids,pretty_otus)
# names(otulabels) <- c("names", "otus")
# 
# 
# otus<-otunames
# ids<-otulabels


graphOTUxCD <- function(otus, ids, corrs){
  
  par(mar=c(0.5,0.5,0.5,0.5))
  design <- matrix(1:8, nrow=4, byrow=T)
  design <- cbind(c(9,9,9,9), design)
  design <- rbind(design, c(0,10,10))
  design <- cbind(design,c(11, 11, 11, 11, 0))
  layout(design, widths=c(0.2,.9,.9, .6), heights=c(1,1, 1, 1,0.3))
  
  for(i in 1:(length(otus))){
    
    otulog <- toptitdel[,otus[i]]
    otulog <- otulog/1625
    zeroes <- otulog == 0
    x_zeroes <- runif(sum(zeroes),.000075,.00016)
    
    #plot controls
    plot(otulog[!zeroes & expgroup=="control"], actual[!zeroes & expgroup=="control"],
         log="x", ylim=c(0,10), xlim=c(0.00008, 1), main="",
         ylab="", 
         xaxt='n', yaxt="n", xlab="")
    
    #plot zeroes in jitter
    points(x=x_zeroes, actual[zeroes], pch=pch[as.character(results[zeroes,1])],
           col=colors[as.character(results[zeroes,1])], cex=1.5)
    
    #add other treatment points
    for(j in 2:(length(legend.labels))){
      points(otulog[!zeroes & names(legend.labels)[j]==expgroup], 
             actual[!zeroes & names(legend.labels)[j]==expgroup], 
             col=colors[j], pch=pch[j], cex=1.5)
    }
    abline(v=c(0.005, 0.05, 0.5), lty="dotted", col="gray")
    abline(v=c(0.001, 0.01, 0.1, 1), lty="longdash", col="gray")
    abline(v=0.0005, lty="dotted", lwd=2, col="black")
    
    test <- grep("aceae", ids[i, 1])
    name<-as.character(ids$name[i])
    number<-as.character(ids$otuname[i])
    corr<-as.character(signif(corrs$corSpear[i], 2))
    library(base)
    graphTitle <- bquote(bolditalic(.(name)) *" " *bold(.(number)) *bold(", ") *bold(rho) *bold(" = ") *bold(.(corr)))
    text(x = 0.00005, y = 9.3, labels = graphTitle, 
         cex = 1, pos=4)
  
    #paste(labels[i], ", \u03c1 = ",signif(corrs[i, 2], 2))
    #bquote(bold(.(labels[i])) ~ bold(", ") ~ bold(rho) ~ bold(" = ") ~ bold(.(signif(corrs[i, 2], 3))))
    #if it's on the bottom row, put a customized axis indicating the % rabund
    if(i == 7 | i==8){
      axis(1, at=c(0.0001,0.001, 0.01, 0.1, 1), labels=c(0,.1, 1, 10, 100), cex.axis=1.5)
    }
    
    #if it's in the first column turn the axis labels to be horizontal
    if(i== 1 | i==3 | i==5 | i==7){
      axis(2, at=c(0, 2, 4, 6, 8, 10), labels =c(0, 2, 4, 6, 8, 10), cex.axis=1.5, las=2)
    }
  }
  
  #for spot 9
  plot.new()
  text(x=0.15, y=0.5, label=expression(paste("Log ", italic("C. difficile"), " CFU/g Feces")), cex=1.5, srt=90)
  
  #for spot 10
  plot.new()
  text(x=0.5, y=0.2, label="% Relative Abundance", cex=1.5)
  
  #for spot 11
  plot.new()
  legend("left", legend=legend.labels, pch=pch, col=colors, pt.cex=2,  cex=1.2, bty="n")
  
}

##########################################################################
# Inputs 
#  
# Returns 
# 
#  
RF.analysis <- function(dataFile, otus, plot.title){
  
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
  trainsets <- list(topdose, titr, del, toptit, titdel, topdel, toptitdel)
  
  #fit the randomforest model
  leng <- length(trainsets)
  results <- data.frame(matrix(ncol=5, nrow=leng))
  row.names(results) <- c("topdose", "titr", "del", "toptit","titdel", "topdel", "toptitdel")
  colnames(results) <- c("PercVarExpl","td", "titr", "delay", "withinValidate")
  library(randomForest)
  for(i in 1:leng){
    trainsets.rf <- randomForest(nextDayCFU~., 
                          data = trainsets[[i]][,-1],  outscale=TRUE,
                          importance=TRUE, proximity=TRUE,
                          keep.forest=TRUE, ntree=5000
    ) #assumes that the nexDayCFU column is right before the first Otu column
    #plot(trainsets.rf)
    varImpPlot(trainsets.rf, type=1)
    # partialPlot(trainsets.rf, trainsets[[i]][,-1], x.var =  "Otu00039") #Otu00003, Otu00013, Otu00023, 27
    
    results[i,1] <- signif(trainsets.rf$rsq[length(trainsets.rf$rsq)],3)
    results[i,2] <- signif(rf.predict(data=topdose, descr=paste0(row.names(results)[i]," model on Topdose Data"), rf.model=trainsets.rf, title=plot.title),3)
    results[i,3] <- signif(rf.predict(data=titr, descr=paste0(row.names(results)[i]," model on Titration Data"), rf.model=trainsets.rf, title=plot.title),3)
    results[i,4] <- signif(rf.predict(data=del, descr=paste0(row.names(results)[i]," model on Delay Data"), rf.model=trainsets.rf, title=plot.title),3)
    results[i,5] <- signif(RF.validate(trainsets[[i]], 100)[1],3)
  }

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
    par(mfrow=c(1, 1)) #+1 to give extra labeling space
    par(mar=c(5, 5, 4, 2) +0.1, mgp=c(3, 1, 0), las=1) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
    plot(results$actual, results$predict, main=title, ylab=expression(paste("Predicted Log ", italic("C. difficile"), " Values")), xlab=expression(paste("Actual Log ", italic("C. difficile"), " Values")), ylim=c(0,9), xlim=c(0,9), xaxt='n', yaxt='n')
    mtext(bquote("r"^"2" ~ " = " ~ .(signif(rsq, 3))), side=3, line=0)
    abline(a=0, b=1, lty="dashed", lwd=2, col="black")
    axis(1, at = c(0:9), labels=c(0:9))
    axis(2, at = c(0:9), labels=c(0:9))
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
file = "~/Desktop/mothur/abxD01/test.txt"
file = "~/Desktop/mothur/abxD01/abxD01.final.an.unique_list.0.03.recovery.vancomycin.thetayc.0.03.lt.ave.dist"

# setRange = vector of 2 numbers signifying the range of days to look over for differences in thetayc or jaccard, etc
getDistDiffs<-function(file, setRange=NULL){
  values <- c()
  values<- scan(file = file, what="character")
  all.ids<-grep("D", values)
  leng <- length(all.ids)
  ltDist <- matrix( nrow=leng, ncol=leng+1)
  idnames <- values[all.ids]
  val.ind<- 2
  for(i in 1:(leng)){
    
    if(i==leng)
    {
      count <- leng
    } else {count <- all.ids[i+1]-all.ids[i]-1}
    
    ltDist[i, 1:(1+count)] <- values[(val.ind):(val.ind+count)]
    val.ind <- val.ind + count +1 
  }
  
  ids <- c()
  days <- c()
  for(i in 1:leng){
    raw.end <- regexpr("D",ltDist[i,1])
    raw.id <- substr(ltDist[i,1], 1, raw.end-1)
    days[i] <- substr(ltDist[i,1], raw.end+1, 100)
    substrInd<-regexpr("[0-9]+[-][0-9]+", raw.id)
    id <- substr(ltDist[i,1], substrInd[1], attr(substrInd, which = "match.length"))
    ids[i] <- id
  }
  uniq.ids<-unique(ids)
  days <- as.numeric(days)

  if (length(setRange)==2){
    start <- setRange[1]
    end <- setRange[2]
  } else {
    start <- range(days)[1]
    end <- range(days)[2] 
  }

  dist.from <- matrix(nrow = (end - start))

  for(i in 1:(length(uniq.ids))){
    id.ind <- ids==uniq.ids[i] #grabs indices in the data frame interested in
    id.days <- days[id.ind] #grabs the days for the first unique id
    id.days <- id.days[order(id.days, decreasing = F)]
    

    start.ind <- c()
    start.ind <- which(days==start & id.ind)
    start.day <- start
    if(length(start.ind)==0){
      local.start <- range(id.days)[1]
      start.ind <- which(days==local.start & id.ind)
      start.day <- local.start
    }
    
    # initialize a "distance from" matrix with NAs based on the entire range
    # This way if we are missing any samples we have an NA instead
    id.from <- matrix(nrow=(end-start))
    rownames(id.from) <- c((start+1):end)
    colnames(id.from) <- idnames[which(ids==uniq.ids[i] & days==start.day)]
    
    if(length(id.days)>1){
      for(j in 2:(length(id.days))){          
        # index to search ltDist:
        index <- which(days==id.days[j] & id.ind)
        
        # find matching index in id.from to put in right spot:
        id.from.index <- which(rownames(id.from)==id.days[j])
        
        # Note: the cols of ltDist are off by 1 because the ids are in col 1
        # That's why there's "+1" for the col index
        if( is.na( ltDist[start.ind, (index+1)] ) ) {
          id.from[id.from.index, 1] <- ltDist[index, (start.ind+1)]
        } else( id.from[id.from.index, 1] <- ltDist[start.ind, (index+1)] )
        
      } # end for(j in 2:(length(id.days)))
    } # end if(length(id.days)>1)
    
    dist.from <- cbind(dist.from,id.from)
    
  } # end for(i in 1:(length(uniq.ids)))
  
  dist.from <- dist.from[,-1]
  
  return(dist.from)
}

plot.dist.from <- function(dist.from){
  matchPos <- match(colnames(dist.from), metadata$sample)
  abx <- metadata$abx[matchPos]
  time <- metadata$time[matchPos]
  day <- metadata$day[matchPos]
  group <- metadata$group[matchPos]
  delay.ids<- colnames(dist.from)[time=="delayed" & day == "-11"]
  delay.abx <- abx[time=="delayed" & day == "-11"]
  delay.group <- group[time=="delayed" & day == "-11"]
  ids <- cbind(delay.ids, delay.abx, delay.group)
  
#   dist.from <- dist.from[ids, ]
#   na.ind <- which(is.na(dist.from))
#   max <- max(c(dist.from)[-na.ind])
#   
  colors<-rainbow(length(unique(delay.abx)))
  colors <- c(amp="#6600FF",
              metro="#FF00DBFF")
  
  plot(x="", y="",xlab="Day", ylab=paste("ThetaYC distance from Day", start), 
           ylim=c(.2,1), xlim=c(-5, 0))
  # Fill remaining points
  for(i in 1:(dim(ids)[1])){
    notNApts <- !is.na(dist.from[6:11,ids[i,1]])
    
    points(row.names(dist.from)[which(notNApts)+5], dist.from[which(notNApts)+5,ids[i]], 
           col=colors[delay.abx[i]], pch=which(unique(delay.group)==delay.group[i]), cex=1)
    lines(row.names(dist.from)[which(notNApts)+5], dist.from[which(notNApts)+5,ids[i]], 
          col=colors[delay.abx[i]], pch=which(unique(delay.group)==delay.group[i]), cex=1)
  }
  legend("bottom", legend=names(colors), col=colors, pch=1, bty="n")
  
  subdist <- dist.from[6:11, ids[,1]]
  ampdist <- subdist[,ids[,2]=="amp"]
  metrodist <- subdist[,ids[,2]=="metro"]
  avg<-function(numbers){
    #print(numbers)
    index <- which(!is.na(numbers))
   # print(numbers[index])
    sumTot <- sum(as.numeric(numbers[index]))
  # print(sumTot)
    total<-sum(!is.na(numbers))
    result <-sumTot/total
    #print(result)
    return(result)
  }
  apply(ampdist, 1, avg)
  plot(x="", y="",xlab="Day", ylab=paste("ThetaYC distance from Day", start), 
       ylim=c(.2,1), xlim=c(-5, 0))
  points(-5:0, apply(ampdist, 1, avg), pch=16, col=colors[1])
  lines(-5:0, apply(ampdist, 1, avg), pch=16, col=colors[1])
  points(-5:0, apply(metrodist, 1, avg), pch=16, col=colors[2])
  lines(-5:0, apply(metrodist, 1, avg), pch=16, col=colors[2])
  legend("bottom", legend=c("Ampicillin", "Metronidazole"), col=colors, pch=16, bty="n", horiz=TRUE, text.width=1)
  
}


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


