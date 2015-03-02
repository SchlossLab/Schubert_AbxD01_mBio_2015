#################################################
#Purpose: Convert csv matrix of columns you want converted into barcharts
#input: csv file following this heading format format:  sampleID, expgroup, (nextDayCFU,) Otu001, Otu002, Otu003, etc
#ex:file<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.5.subsample.relabund.topdose.forlogscale2.csv
#Note: It is assumed that anything after the expgroup column are vectors containing data to be divide by the groups in 'expgroup' 
#      and averaged with SD measured. 
#Output: matrices avg[],  std[], bp[y] with the barplot coordinates info stored for each subgraph
#################################################
# Parameters to change, an example:
# CSV file: Group  expgroup  Otu001... (limited by most abund, end with 'Other', OTUs normalized +0.001, expgroups #'d by graph order & sorted by first graph)
# 


# file <- "~/Documents/Github/abxD01/Figure 5/abxD01.final.tx.2.subsample.alldelay.forlogscale.csv"
# fileIDS <- "~/Documents/Github/abxD01/Figure 5/alldelay_tx2_barchart_ids.csv"
# #graphLabels<-c("Cefoperazone", "Streptomycin", "Vancomycin")
# sortbyphyl<-TRUE
# graphbyphyl<-FALSE
# divide<-TRUE


# Highlight all and run!
#################################################
stackedbarcharts <- function(file, fileIDS, graphLabels=NULL, sortbyphyl=TRUE, graphbyphyl=FALSE, excludeGroup=NA, divide=FALSE) {
  file<-read.csv(file=file, header=T)
  
  fileIDS<-read.csv(file=fileIDS, header=T)
  abx <- graphLabels
  
  if(graphbyphyl==TRUE)
  {
    sortbyphyl<-TRUE
  }
  
  id.info <- NULL
  if(divide == TRUE) {
    id.info <- file[,1:4]
    file<-file[,-1] # get rid of the extra divide column in this filetype
    
  }

  file<-file[,-1] #delete the first col of group names
  
  
  avgs<-NULL
  avgs <-data.frame(levels(file$expgroup))
  colnames(avgs) <- c("expgroup")
  stds<-NULL
  stds <-data.frame(levels(file$expgroup))
  colnames(stds) <- c("expgroup")
  SEs<-NULL
  SEs <-data.frame(levels(file$expgroup))
  colnames(SEs) <- c("expgroup")
  
  #Calculate average cdiff
  avgCD<-NULL
  sdCD<-NULL
  avgCD<-tapply(file$nextDayCFU, file$expgroup, mean)
  avgCD<-format(avgCD, scientific=TRUE, digits=2)
  sdCD<-tapply(file$nextDayCFU, file$expgroup, sd)
  file<-file[,-2] #delete the nexDayCFU column
  
  library(plotrix)
  
  for(i in 1:ncol(file)){
    
    if(i==1){
      ids<-names(file)
    }
    else{
      columni <- data.frame(file[,i])
      colnames(columni) <- c("columndata")
      
      #calculates mean for each column in file by expgroup and stores it in avgs data frame
      poop<- aggregate(columni$columndata~file$expgroup, FUN=mean) #poop stores the aggregated means for each group for this specific column
      colnames(poop) <- c(ids[1], ids[i])
      
      avgs<- merge(avgs, poop, by ="expgroup")
      
      #calculates st dev for each column in file by expgroup and stores it in stds data frame
      poop2<- aggregate(columni$columndata~file$expgroup, FUN=sd) #poop stores the aggregated means for each group for this specific column
      colnames(poop2) <- c(ids[1], ids[i])
      
      stds<- merge(stds, poop2, by ="expgroup")
      
      
      #calculates lengthfor each column in file by expgroup and stores it in lengths data frame
      poop3<- aggregate(columni$columndata~file$expgroup, FUN=std.error) #poop stores the aggregated means for each group for this specific column
      colnames(poop3) <- c(ids[1], ids[i])
      
      SEs<- merge(SEs, poop3, by ="expgroup")
    }
    
  }
  
  
  
  ordIDS <- fileIDS[order(fileIDS[,2]), ] #sort by the phyla
  totalOTU <- dim(ordIDS)[1]
  
  ########for doing sort by phylum then within phylum (by relabund)
  if(sortbyphyl == TRUE){
    ordIDS<-cbind(ordIDS, 0)
    
    for (i in 1:totalOTU){ #loop to fill out the final column for the phyla groups, while phyla is sorted!!
      if(i==1){
        ordIDS[i, 4] <- 1
      }
      else{
        if(ordIDS[i, 2]==ordIDS[i-1, 2]){
          ordIDS[i, 4] <- ordIDS[i-1, 4]
        }
        else{
          ordIDS[i, 4] <- ordIDS[i-1, 4] + 1
        }
        
      }
    }
    
    names(ordIDS)[4]<-paste("phynum")
    
    ### WILL GIVE AN ERROR MESSAGEWarning message:
    #In `[<-.factor`(`*tmp*`, ri, value = 0) :
    #  invalid factor level, NA generated
    avgs<-rbind(0,avgs) 
    
    numphyla<-max(ordIDS[,4])
    
    #test_avg <- avgs[, rev(order(avgs[1,2:length(avgs)-1]))] #sort all but other column (last)
    
    physums <- data.frame(c(1:numphyla))
    physums <- cbind(physums, 0)
    colnames(physums) <- c("Num", "Sum")
    
    leng<-length(avgs)
    
    #jtime<-NULL
    #itime<-NULL
    #inside<-FALSE
    #count<-0
    j<-1
    ##This loop is to fill the physums data frame with the sum of the relabund for each phylum represented
    for(i in 2:(leng-1)){ #searching through avgs, so leng is number of otus, -1 because i dont want to deal with "Other" column
      found<-FALSE
      #itime<-rbind(itime, i)
      j<-1
      while(found==FALSE){ #searching through ordIDS
        test<-ordIDS[j, 1] == names(avgs)[i]
        if(test) { 
          pnum<-ordIDS[j, 4]
          found<-TRUE
          inside<-TRUE
          avgs[1,i]<-pnum
        }
        else{
          j <- j + 1
          #count<-count+1
          #jtime<-rbind(jtime, j)    
        }
      } #WARNING: could cause an error if dont find the ID, theres no check for the end of the list 
      physums[pnum, 2] <- physums[pnum, 2] + avgs[2, i] ##add the relabund of the otu   
    }
    
    
    y<-rev(order(physums$Sum))
    x<-avgs[1, 2:(leng-1)]
    
    
    ord_physums<-physums[rev(order(physums[,2])),]
    sort_avgs<-data.frame(avgs$expgroup)
    names(sort_avgs)[1]<-c("expgroup")
    
    #loop is to sort by decreasing phyla and perform internal sort within phylum, return the new sorted avgs to go into mavgs
    for(i in 1:numphyla){ 
      val<-NULL
      phy<-which(x %in% c(ord_physums[i,1])) #phy returns the indices for where that phyNum is
      for(j in 1:length(phy)){ #using locations of phylum i to find in avgs
        val[j]<-avgs[2,phy[j]+1]
      }
      phy<-phy[rev(order(val))]
      
      for(j in 1:length(phy)){ #now loop through the sorted sub cols and add them to the sort_avgs
        sort_avgs<-cbind(sort_avgs, avgs[,phy[j]+1])
        names(sort_avgs)[length(sort_avgs)]<-names(avgs)[phy[j]+1]
      }
    }
    
    attach(avgs)
    sort_avgs<-cbind(sort_avgs, Other)
    detach(avgs)
    
    
    col <- c(1)
    for( i in 2:(length(SEs)) ) {
      col <- c(col, which(names(SEs)==names(sort_avgs)[i]))
    }
    
    sort_stderr <- SEs[,col]
    
    sort_phy <- sort_avgs[1,]
    sort_phy<-sort_phy[, -1]
    sort_avgs <- sort_avgs[-1,] #get rid of phylum first row from earlier sorting
    
    #loop for getting the right ids in the right order as sort
    ids<-data.frame(1:(totalOTU))
    
    row.names(sort_avgs)<-sort_avgs$expgroup
    sort_avgs<-sort_avgs[,-1]
    
    row.names(sort_stderr)<-sort_stderr$expgroup
    sort_stderr<-sort_stderr[,-1]
    
    
    for (i in 1:(length(sort_avgs)-1)){
      found<-FALSE
      #itime<-rbind(itime, i)
      j<-1
      while(found==FALSE){ #searching through ordIDS
        test<-ordIDS[j, 1] == names(sort_avgs)[i]
        if(test) { 
          ids[i,2]<-ordIDS[j, 3]
          ids[i,3]<-ordIDS[j,2]        
          
          found<-TRUE
        }
        else{
          j <- j + 1
          #count<-count+1
          #jtime<-rbind(jtime, j)
        }
      }
      
    }
    
    
    names(ids)<-c("x", "name", "phyname")
    ids<-rbind(ids, data.frame(x=totalOTU+1, name= "Other", phyname="Other"))
    ids<-cbind(ids, t(sort_phy))
    names(ids)<-c("x", "name", "phyname", "phynum")
    
    mavgs<-as.matrix(sort_avgs)
    
    leng<-dim(mavgs)[2]
    
    numgr <- length(unique(file$expgroup))
  } ##end if(sortbyphyl == TRUE){
  
  
  ####order the results by the first group alpha numerically
  if(sortbyphyl == FALSE){
    
    row.names(avgs)<-avgs$expgroup
    avgs<-avgs[,-1]
    ordered_avgs<- avgs[, rev(order(avgs[1,1:length(avgs)-1]))] #sort all but other column (last)
    attach(avgs)
    ordered_avgs <- cbind(ordered_avgs, Other) #put other back on
    detach(avgs)
    
    mavgs<-as.matrix(ordered_avgs)
    leng<-dim(mavgs)[2]
    
    numgr <- length(unique(file$expgroup))
    
    ids<-data.frame(1:(totalOTU))
    
    for (i in 1:(length(ordered_avgs)-1)){
      found<-FALSE
      #itime<-rbind(itime, i)
      j<-1
      while(found==FALSE){ #searching through ordIDS
        test<-ordIDS[j, 1] == names(ordered_avgs)[i]
        if(test) { 
          ids[i,2]<-ordIDS[j, 3]
          found<-TRUE
        }
        else{
          j <- j + 1
          #count<-count+1
          #jtime<-rbind(jtime, j)
        }
      }
    }
    
    names(ids)<-c("x", "name")
    ids<-rbind(ids, data.frame(x=totalOTU+1, name= "Other"))
    
  }  #end if(sortbyphyl == FALSE){
  
  
  ###started to do the stds table... didnt complete
  matrix.se<-as.matrix(sort_stderr)
  
  
  
  
  #####################################################################################################################
  #####################################################################################################################
  #####################################################################################################################
  ###PLOT PARAMETERS
  library(plotrix)
  if(graphbyphyl==FALSE && divide==FALSE){
    par(mfrow=c(numgr+1, 1)) #+1 to give extra labeling space
    par(mar=c(0.3, 8, 0.5, 2) +0.1, mgp=c(4.5, 1, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
    color_transparent <- adjustcolor("black", alpha.f = 0.1) 
    color <- gray.colors((numphyla+1), start=0, end=1, alpha=NULL)
    
    count <- NULL
    for (n in 1:(numphyla)){
      count[n] <- length(which(ids$phynum==n))
    }
    
    colors <- NULL
    for( n in 1:(length(color)-1) ) {
      phyNumOrder <- unique(ids$phynum)
      times <- count[phyNumOrder[n]] 
      temp <- rep( color[n], times )
      colors <- c(colors, temp)
    }
    colors <- c(colors, "#FFFFFF")
    
    for(j in 1:numgr){
      
      if(j != numgr){
        bp <- barplot(mavgs[j, 1:leng], ylab=NULL, col=colors,  yaxt="n", xaxt="n", ylim=c(0.001, 1), log="y", cex.names=3)
        #error.bar(bp[k], mavgs[j, 1:leng[2]], mstds[j, 1:leng[2]])  #the bp[k] was for storing the barplot locations to use for these errors
        # k <- k+1
        d <- as.data.frame( cbind(bp, mavgs[j,], matrix.se[j,] ))
        names(d) <- c("x", "y", "se")
        
        with (
          data = d
          , expr = errbar(x, y, y+se, y-se, add=T, pch=".", cap=.01)
        )
        
        axis(2, las=1, at=c(.001, .01, .1, 1), labels=c(0, .01, .1, 1), cex.axis=1.1)
        mtext(abx[j], side=2, line=6, cex=.8)
        mtext(avgCD[j], side=2, line=4.5, cex=.8)
        abline(h=c(0.001, 1), lwd=3) #min/max
        abline(h=c(0.01, 0.1), col=color_transparent, lty="longdash", lwd=2)
        abline(h=c(0.0025, 0.005, 0.0075, 0.025, 0.05, 0.075, 0.25, 0.5, 0.75), col=color_transparent, lty="dashed")
      }
      
      else{ #the last graph needs different margins    
        label<-barplot(mavgs[j, 1:leng], col=colors,ylab=NULL,  yaxt="n", xaxt="n", ylim=c(0.001, 1), log="y", cex.names=3)
        
        d <- as.data.frame( cbind(label, mavgs[j,], matrix.se[j,] ))
        names(d) <- c("x", "y", "se")
        
        with (
          data = d
          , expr = errbar(x, y, y+se, y-se, add=T, pch=".", cap=.01)
        )
        
        #error.bar(bp[k], mavgs[j, 1:leng[2]], mstds[j, 1:leng[2]])  #the bp[k] was for storing the barplot locations to use for these errors
        axis(2, las=1, at=c(.001, .01, .1, 1), labels=c(0, .01, .1, 1), cex.axis=1.1)
        mtext(abx[j], side=2, line=6, cex=.8)
        mtext(avgCD[j], side=2, line=4.5, cex=.8)
        axis(1, at=(label[,1]), labels=FALSE)
        text(label[,1]+.13, .0005, label=ids[,2], xpd=NA, pos=2, srt=45, cex=1.2)
        #text(-3.9,.0001, label=expression(paste(italic("C.d."), " CFU/g Feces:")), xpd=NA, pos=2, srt=90, cex=1.1)
        abline(h=c(0.001,1), lwd=3) #min/max
        abline(h=c(0.01, 0.1), col=color_transparent, lty="longdash", lwd=2)
        abline(h=c(0.0025, 0.005, 0.0075, 0.025, 0.05, 0.075, 0.25, 0.5, 0.75), col=color_transparent, lty="dashed") 
      }  
      
    }
    par(mfrow=c(1, 1))
  } #if(graphbyphyl=FALSE)
  #####END 1ST PLOT PARAMETERS
  
  
  
  ################################
  ###2nd PLOT PARAMETERS
  if(graphbyphyl==TRUE){  
    par(mfrow=c(numphyla+1, 1)) #+2 to give extra labeling space--numphyl doesn't include the Other group
    #  par(mfrow=c(3, 1)) #temporarily for testing
    par(mar=c(2.5, 8, 0.5, 2) +0.1, mgp=c(4.5, 1, 0)) #default is  par(mar=c(5, 4, 4, 2 ) +0.1, mgp=c(3, 1, 0)) bot/left/top/right, also default mgp is c(3,1,0)
    colors= gray.colors(numgr, start=0, end=1, alpha=NULL)
    currphy=ids[1, 4]
    k <- 1
    idleng<-dim(mavgs)[2]
    j <- 1
    
    while(j < (idleng+1)){
      leng <- 0
      
      
      while(currphy==ids[k, 4])
      {
        leng <- leng+1
        k <- k+1
        if( k > idleng ){
          break}
      }
      
      label<-barplot(mavgs[,j:(leng+j-1)],, beside=TRUE, ylab=ids[j,3], col=colors, yaxt="n", xaxt="n", ylim=c(0.001, 1), log="y", cex.names=5)
      
      #CTR+SHIFT+C=comment block of code out    
      
      #Will perform the pairwise wilcox test for each OTU and fill "statLetter" with the lettering scheme for the graph to show statistical signifiance
      #This will work for the titration which has 3 different comparisons
      statLetter <- NULL
      statLetter <- matrix(c("NA"),nrow=numgr, ncol=leng)
      otus <- dimnames(mavgs)[2][[1]][j:(leng+j-1)]
      m <- j
      for (i in 1:leng){      
        results.wilcox <- pairwise.wilcox.test(file[,which(names(file)==otus[i])], file$expgroup, p.adj="BH")
        #results.wilcox <- pairwise.t.test(file[,which(names(file)==otus[i])], file$expgroup, p.adj="BH")
        
        #double check that all values were calculated, if not put as n.s.
        if(results.wilcox$p.value[1] == "NaN"){
          results.wilcox$p.value[1] <- 10 #a value greater than 0.05
        }
        if(results.wilcox$p.value[2] == "NaN"){
          results.wilcox$p.value[2] <- 10
        }
        if(results.wilcox$p.value[4] == "NaN"){
          results.wilcox$p.value[4] <- 10
        }
        
        if( results.wilcox$p.value[1] >= 0.05 ){
          if( results.wilcox$p.value[2] >= 0.05 ){
            if( results.wilcox$p.value[4] >= 0.05 ){
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "a"
              statLetter[3, i] <- "a"
            } else{ 
              statLetter[1, i] <- "ab"
              statLetter[2, i] <- "a"
              statLetter[3, i] <- "b"            
            }
          } else{
            if( results.wilcox$p.value[4] >= 0.05 ){
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "ab"
              statLetter[3, i] <- "b"
            } else{
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "a"
              statLetter[3, i] <- "b"            
            }
          }
        } else{
          if( results.wilcox$p.value[2] >= 0.05 ){
            if( results.wilcox$p.value[4] >= 0.05 ){
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "b"
              statLetter[3, i] <- "ab"
            } else{
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "b"
              statLetter[3, i] <- "a"
            }
          } else{
            if( results.wilcox$p.value[4] >= 0.05 ){
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "b"
              statLetter[3, i] <- "b"
            } else{
              statLetter[1, i] <- "a"
              statLetter[2, i] <- "b"
              statLetter[3, i] <- "c"
            }
          }
        }  
      }
      
      #For labeling letters above each bar
      m <- j
      for(n in 1:leng){
        text( label[,n], mavgs[,m], labels=statLetter[,n], cex=1, col="red", pos=3, offset=0.15, xpd=TRUE)
        m <- m+1
      }
      
      
      axis(2, las=1, at=c(.001, .01, .1, 1), labels=c(0, .01, .1, 1), cex.axis=1.1)
      # mtext("Relative Abundance", side=2, line=6, cex=.8)
      #mtext(avgCD[j], side=2, line=4.5, cex=.8)
      abline(h=c(0.001, 1), lwd=3) #min/max
      
      color_transparent <- adjustcolor("black", alpha.f = 0.1) 
      abline(h=c(0.01, 0.1), col=color_transparent, lty="longdash", lwd=2)
      abline(h=c(0.0025, 0.005, 0.0075, 0.025, 0.05, 0.075, 0.25, 0.5, 0.75), col=color_transparent, lty="dashed")
      
      labelAVG=apply(label, 2, mean)
      axis(1, at=(labelAVG), labels=FALSE)
      text(labelAVG+.13, .0005, label=ids[j:(leng+j-1),2], xpd=NA, pos=1, srt=15, cex=.8)
      
      j<-j+leng
      currphy <- ids[j, 4]
    }
    par(mfrow=c(1, 1))
  } #if(graphbyphyl=TRUE) 
  
  
  
  ################################
  ###3rd PLOT PARAMETERS
  if(divide == TRUE){ 

    if(!is.na(excludeGroup)){
      mavgs <- mavgs[-which(grepl(excludeGroup, row.names(mavgs))),]
      matrix.se <- matrix.se[-which(grepl(excludeGroup, row.names(matrix.se))),]
      id.info <- id.info[-which(grepl(excludeGroup, id.info$expgroup)),]
    }
    
    barplotBeside(id.info = id.info, matrix.of.avgs = mavgs, matrix.of.se = matrix.se, ids = ids, file=file)
    
  } # if(divide=TRUE){
  
}


