# The most number of side by side groups should be 3. Not made for others. 

#run this alyx!! to fix missing stat letter
barplotBeside <- function(id.info, matrix.of.avgs, matrix.of.se, ids, file){

  library(Hmisc)
  groups <- unique(id.info$abx)
  numgr <- length(groups)
  leng<-dim(mavgs)[2]
  
  par(mfrow=c(numgr+1, 1)) #+1 to give extra labeling space
  par(mar=c(0.3, 8, 0.5, 2) +0.1, mgp=c(4.5, 1, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
  color_transparent <- adjustcolor("black", alpha.f = 0.2) 
  color <- gray.colors(3, start=0.2, end=1, alpha=NULL)

  for(j in 1:numgr){
    
    graphGroups <- as.character(unique(id.info[id.info$abx == groups[j], "expgroup"]))
    
    ind <- NULL
    for(a in 1:(length(graphGroups))){
      ind <- c(ind, which(row.names(matrix.of.avgs)==graphGroups[a]))
    }
    mavgs <- matrix.of.avgs[ind,]
    matrix.se <- matrix.of.se[ind,]
    otus <- dimnames(mavgs)[2][[1]]
    
    if(j != numgr){
      bp <- barplot(mavgs, ylab=NULL, col=color,  beside=TRUE, yaxt="n", xaxt="n", ylim=c(0.001, 1), log="y", cex.names=.5)

      # add standard error bars
      d <- as.data.frame( cbind( c(bp), c(mavgs), c(matrix.se) ) )
      names(d) <- c("x", "y", "se")
      library(Hmisc)
      with (
        data = d
        , expr = errbar(x, y, y+se, y-se, add=T, pch=".", cap=.01)
      )
      
      # Statistical differences between groups
      statLetter <- findStatLetters(file, otus, graphGroups)
      
      statLetter <- as.data.frame(statLetter)
      bp <- as.data.frame(bp)
      mavgs <- as.data.frame(mavgs)
      matrix.se <- as.data.frame(matrix.se)
      
      #For labeling letters above each bar
      for(n in 1:leng){
        text( bp[,n], (mavgs[,n]+matrix.se[,n]), labels=statLetter[,n], cex=1, col="red", pos=3, offset=0.15, xpd=TRUE)
      }
      
      # add labels and pretty the graph
      axis(2, las=1, at=c(.001, .01, .1, 1), labels=c(0, .01, .1, 1), cex.axis=1.1)
      mtext(groups[j], side=2, line=6, cex=.8)
      abline(h=c(0.001, 1), lwd=3) #min/max
      abline(h=c(0.01, 0.1), col=color_transparent, lty="longdash", lwd=2)
      abline(h=c(0.0025, 0.005, 0.0075, 0.025, 0.05, 0.075, 0.25, 0.5, 0.75), col=color_transparent, lty="dashed")
    
      } else{ #the last graph needs different margins    
     
        bp<-barplot(mavgs, col=color,beside=TRUE, ylab=NULL,  yaxt="n", xaxt="n", ylim=c(0.001, 1), log="y", cex.names=3)
        
        # add standard error bars
        d <- as.data.frame( cbind( c(bp), c(mavgs), c(matrix.se) ) )
        names(d) <- c("x", "y", "se")
        library(Hmisc)
        with (
          data = d
          , expr = errbar(x, y, y+se, y-se, add=T, pch=".", cap=.01)
        )
        
        # Statistical differences between groups
        statLetter <- findStatLetters(file, otus, graphGroups)
        statLetter <- as.data.frame(statLetter)
        
        bp <- as.data.frame(bp)
        mavgs <- as.data.frame(mavgs)
        matrix.se <- as.data.frame(matrix.se)
        
        #For labeling letters above each bar
        for(n in 1:leng){
          text( bp[,n], (mavgs[,n]+matrix.se[,n]), labels=statLetter[,n], cex=1, col="red", pos=3, offset=0.15, xpd=TRUE)
        }
        
        # add labels and pretty the graph
        labelAVG <- apply(bp, 2, mean)
        axis(1, at=(labelAVG), labels=FALSE)
        text(labelAVG+.13, .0005, label=ids[,2], xpd=NA, pos=2, srt=45, cex=1.2)
        axis(2, las=1, at=c(.001, .01, .1, 1), labels=c(0, .01, .1, 1), cex.axis=1.1)
        mtext(groups[j], side=2, line=6, cex=.8)
        abline(h=c(0.001,1), lwd=3) #min/max
        abline(h=c(0.01, 0.1), col=color_transparent, lty="longdash", lwd=2)
        abline(h=c(0.0025, 0.005, 0.0075, 0.025, 0.05, 0.075, 0.25, 0.5, 0.75), col=color_transparent, lty="dashed")     
    } # else last graph
        
  } # for(j in 1:numgr){
 
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
  leg<-c("Untreated", "On Time","+5d Recovery" )
  #pch<-c(16, 16, 16)
  colors<-color
  colors[3] <- "black"
  legend("bottom", legend=leg, pch=c(15, 15, 22), col=colors, pt.cex=2.5,  cex=1.5, bty="n", horiz=TRUE, text.width=.2)
  
  par(mfrow=c(1, 1)) #default
  par(mar=c(5, 4, 4, 2) +0.1, mgp=c(3, 1, 0)) #default
}