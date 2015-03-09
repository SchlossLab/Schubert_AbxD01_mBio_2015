#function to plot the error bars, either SEM or SD
##################################################
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

###Deconstructed Stacked bar charts
###################################
###Deconstructed Stacked bar charts

c<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.2.subsample.2.pick.relabund.D0.csv", header=T)
b<-read.csv("~/Desktop/mothur/abxD01/correlation/abxD01.final.tx.2.subsample.2.pick.relabund.D01.csv", header=T)

avgCD<-tapply(b$nextDayCFU, b$exp, mean)
sdCD<-tapply(b$nextDayCFU, b$exp, sd)

avg1<-tapply(c$Otu001, c$exp, mean)
sd1<-tapply(c$Otu001, c$exp, sd)

avg14<-tapply(c$Otu014, c$exp, mean)
sd14<-tapply(c$Otu014, c$exp, sd)

avg27<-tapply(c$Otu027, c$exp, mean)
sd27<-tapply(c$Otu027, c$exp, sd)

avg7<-tapply(c$Otu007, c$exp, mean)
sd7<-tapply(c$Otu007, c$exp, sd)

par(mfrow=c(5, 1)) 

plotCD<-barplot(avgCD, col="black", cex.names=.6) 
error.bar(plotCD, avgCD, sdCD) #plot SD

leg<-c("cdiff CFU", "porphyromonadaceae -0.65", "alistipes -0.68", "clostridia -0.58", "enterobacteriaceae +0.55")
col<-c("black", "blue","darkgreen", "darkorange", "red")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=1, bty="n")

plot1<-barplot(avg1, col="blue", ylim=c(0, .8), cex.names=.6) #how can i change this to log scale
error.bar(plot1, avg1, sd1) #plot SD

plot14<-barplot(avg14, col="darkgreen", ylim=c(0, .8), cex.names=.6)
error.bar(plot14, avg14, sd14) #plot SD

plot27<-barplot(avg27, col="darkorange", ylim=c(0, .8), cex.names=.6)
error.bar(plot27, avg27, sd27) #plot SD

plot7<-barplot(avg7, col="red", ylim=c(0, .8), cex.names=.6)
error.bar(plot7, avg7, sd7) #plot SD

par(mfrow=c(1, 1)) 

###no set scale

par(mfrow=c(5, 1)) 

plotCD<-barplot(avgCD, col="black", cex.names=.6) 
error.bar(plotCD, avgCD, sdCD) #plot SD



plot1<-barplot(avg1, col="blue", cex.names=.6) 
error.bar(plot1, avg1, sd1) #plot SD

leg<-c("cdiff CFU", "porphyromonadaceae -0.65", "alistipes -0.68", "clostridia -0.58", "enterobacteriaceae +0.55")
col<-c("black", "blue","darkgreen", "darkorange", "red")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=.7, bty="n")

plot14<-barplot(avg14, col="darkgreen", cex.names=.6)
error.bar(plot14, avg14, sd14) #plot SD

plot27<-barplot(avg27, col="darkorange", cex.names=.6)
error.bar(plot27, avg27, sd27) #plot SD

plot7<-barplot(avg7, col="red",  cex.names=.6)
error.bar(plot7, avg7, sd7) #plot SD

par(mfrow=c(1, 1)) 

#to call error.bar, call the function at the top of this file
#error.bar(plot1, avg1, 1.96*sd1/10)  #SEM=SD/sqqrt(N)    so need to have N's for each





##loop that will make all the avg/std tables for all otus in the file (contains only otus that want to make chart for)
##the loop will then make all the barcharts, and put them together using par(mfrow=c(#OTUSincluded, 1))
#meta<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.2.subsample.2.pick.relabund.D0.csv", header=T)
#meta<-meta[1:56]
#for(i in 3:length(meta)){
#  otu[c] <- colnames(meta[i])
#  cor[c] <- cor.test(meta[,2],meta[,i], method="spearman")$estimate
#  pval[c] <- cor.test(meta[,2],meta[,i], method="spearman")$p.value
#  c <- c+1
#}



###phylum level barcharts across topdoses
###################################
###phylum level barcharts across topdoses

#################################################
#Purpose: Convert csv matrix of columns you want converted into barcharts
#input: csv file following this heading format format:  sampleID, expgroup, (nextDayCFU,) Otu001, Otu002, Otu003, etc
#ex:file<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.5.subsample.relabund.topdose.forlogscale2.csv
#Note: It is assumed that anything after the expgroup column are vectors containing data to be split by the groups in 'expgroup' 
#      and averaged with SD measured. 
#Output: matrices avg[],  std[], bp[y] with the barplot coordinates info stored for each subgraph
#################################################

file<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.5.subsample.relabund.allceftitr.forlogscale.csv", header=T)
file<-file[,-1] #delete the first col of group names
#attributes(file)
fileIDS<-read.csv("~/Desktop/mothur/abxD01/barcharts/alltitr_tx5_barchart_ids.csv", header=T)

avgs=NULL
avgs <-data.frame(levels(file$expgroup))
colnames(avgs) <- c("expgroup")
stds=NULL
stds <-data.frame(levels(file$expgroup))
colnames(stds) <- c("expgroup")


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
  }
}
  
ordIDS <- fileIDS[order(fileIDS[,2]), ] #sort by the phyla
totalOTU <- dim(ordIDS)[1]

########for doing sort by phylum then within phylum
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

##add row for phylum number to the avgs?
#avgs<-rbind(avgs, 0)[c(9, 1, 2, 3, 4, 5, 6, 7, 8),] #this is hard coded for the top doses... also dont rerun this multiple times!!! might need to change
avgs<-rbind(0,avgs) #this is hard coded for the top doses... also dont rerun this multiple times!!! might need to change, will give error

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

for(i in 1:numphyla){ #loop is to sort by decreasing phyla and perform internal sort within phylum, return the new sorted avgs to go into mavgs
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

sort_avgs <- sort_avgs[-1,] #get rid of phylum first row from earlier sorting

#loop for getting the right ids in the right order as sort
ids<-data.frame(1:(totalOTU))

row.names(sort_avgs)<-sort_avgs$expgroup
sort_avgs<-sort_avgs[,-1]

for (i in 1:(length(sort_avgs)-1)){
  found<-FALSE
  #itime<-rbind(itime, i)
  j<-1
  while(found==FALSE){ #searching through ordIDS
    test<-ordIDS[j, 1] == names(sort_avgs)[i]
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


mavgs<-as.matrix(sort_avgs)

leng<-dim(mavgs)[2]

numgr <- nlevels(file$expgroup)


####order the results by the first group alpha numerically

#row.names(avgs)<-avgs$expgroup
#avgs<-avgs[,-1]


####order the results by the first group alpha numerically
#ordered_avgs<- avgs[, rev(order(avgs[1,1:length(avgs)-1]))] #sort all but other column (last)
#attach(avgs)
#ordered_avgs <- cbind(ordered_avgs, Other) #put other back on
#detach(avgs)
#ids<-names(ordered_avgs) #change order for the ids too
#mavgs<-as.matrix(ordered_avgs)

#stds<-stds[,-1]
#mstds<-as.matrix(stds)



par(mfrow=c(numgr+1, 1)) #+1 to give extra labeling space
par(mar=c(0.1, 7, 0.5, 2) +0.1, mgp=c(4.5, 1, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)

abx<-c("Untreated", "Strep 5.0", "Strep 0.5", "Strep 0.1")

for(j in 1:numgr){
 
  if(j != numgr){
    barplot(mavgs[j, 1:leng], ylab=abx[j], col="black", yaxt="n", xaxt="n", ylim=c(0.0001, 1), log="y", cex.names=3)
    #error.bar(bp[k], mavgs[j, 1:leng[2]], mstds[j, 1:leng[2]])  #the bp[k] was for storing the barplot locations to use for these errors
   # k <- k+1
    axis(2, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1), cex.axis=1.1)
    abline(h=c(0.0001, 1), lwd=3)
    abline(h=c(0.001, 0.01, 0.1), col="black", lty="longdash", lwd=1.5)
    abline(h=c(.25,.5, .75), col="black", lty="dashed")
  }
  
  else{ #the last graph needs different margins    
    label<-barplot(mavgs[j, 1:leng], col="black",ylab=abx[j],  yaxt="n", xaxt="n", ylim=c(0.0001, 1), log="y", cex.names=3)
    #error.bar(bp[k], mavgs[j, 1:leng[2]], mstds[j, 1:leng[2]])  #the bp[k] was for storing the barplot locations to use for these errors
    axis(2, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1), cex.axis=1.1)
    axis(1, at=(label[,1]), labels=FALSE)
    text(label[,1]+.13, .0001-.00008, label=ids[,2], xpd=NA, pos=2, srt=45, cex=1.2)
    abline(h=c(0.0001,1), lwd=3)
    abline(h=c(0.001, 0.01, 0.1), col="black", lty="longdash", lwd=1.5)
    abline(h=c(.25,.5, .75), col="black", lty="dashed")
    
  }  
  
}
par(mfrow=c(1, 1))


##############END CODE###########################
#################################################


c<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.5.subsample.relabund.topdose.csv", header=T)



avgCD<-tapply(c$nextDayCFU, c$expgroup, mean)
sdCD<-tapply(c$nextDayCFU, c$expgroup, sd)

avg1<-tapply(c$Otu001, c$expgroup, mean)
sd1<-tapply(c$Otu001, c$expgroup, sd)

avg14<-tapply(c$Otu014, c$expgroup, mean)
sd14<-tapply(c$Otu014, c$expgroup, sd)

avg27<-tapply(c$Otu027, c$expgroup, mean)
sd27<-tapply(c$Otu027, c$expgroup, sd)

avg7<-tapply(c$Otu007, c$expgroup, mean)
sd7<-tapply(c$Otu007, c$expgroup, sd)

par(mfrow=c(5, 1)) 

plotCD<-barplot(avgCD, col="black", cex.names=.6) 
error.bar(plotCD, avgCD, sdCD) #plot SD

leg<-c("cdiff CFU", "porphyromonadaceae -0.65", "alistipes -0.68", "clostridia -0.58", "enterobacteriaceae +0.55")
col<-c("black", "blue","darkgreen", "darkorange", "red")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=1, bty="n")

plot1<-barplot(avg1, col="blue", ylim=c(0, .8), cex.names=.6) #how can i change this to log scale
error.bar(plot1, avg1, sd1) #plot SD

plot14<-barplot(avg14, col="darkgreen", ylim=c(0, .8), cex.names=.6)
error.bar(plot14, avg14, sd14) #plot SD

plot27<-barplot(avg27, col="darkorange", ylim=c(0, .8), cex.names=.6)
error.bar(plot27, avg27, sd27) #plot SD

plot7<-barplot(avg7, col="red", ylim=c(0, .8), cex.names=.6)
error.bar(plot7, avg7, sd7) #plot SD

par(mfrow=c(1, 1)) 

##paper figure: cdiff topdose CFU
c<-read.csv("~/Desktop/mothur/abxD01/barcharts/topdose_nextDayCFU.csv", header=T)

avgCD<-tapply(c$tnextDayCFU, c$expgroup, mean)
sdCD<-tapply(c$tnextDayCFU, c$expgroup, sd)
abx<-c("Untreated", "Ciprofloxacin", "Clindamycin", "Cefoperazone", "Streptomycin", "Vancomycin", "Metronidazole", "Ampicillin")

par(mar=c(8, 8, 2, 2) +0.1, mgp=c(5, 2, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plotCD<-barplot(avgCD, col="black", log="y", yaxt="n", xaxt="n", ylim=c(1, 1e9), cex.lab=2, ylab=expression(paste(italic("C. difficile"), " CFU/g Feces"))) 
error.bar(plotCD, avgCD, sdCD, lwd=7) #plot SD
ticksy<-seq(0, 9, by=1)
ticksx<-seq(1,8, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, lwd=2, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), cex.axis=1.7, labels=labelsy)
axis(1, lwd=2, at=c(plotCD), labels=FALSE)
text(plotCD+.05, 0.5, label=abx, xpd=NA, pos=2, srt=45, cex=1.2)

#text(plotCD[2], 8*10^8, label="TITLE", cex=2.2, pos=1)



#####order level barcharts for day 0/1 microbiome/CFU for the primary set

c<-read.csv("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.3.subsample.3.pick.relabund.vanc.D01.csv", header=T)

avgCD<-tapply(c$nextDayCFU, c$exp, mean)
sdCD<-tapply(c$nextDayCFU, c$exp, sd)

avg1<-tapply(c$Otu001, c$exp, mean)
sd1<-tapply(c$Otu001, c$exp, sd)

avg05<-tapply(c$Otu005, c$exp, mean)
sd05<-tapply(c$Otu005, c$exp, sd)

avg04<-tapply(c$Otu004, c$exp, mean)
sd04<-tapply(c$Otu004, c$exp, sd)

avg3<-tapply(c$Otu003, c$exp, mean)
sd3<-tapply(c$Otu003, c$exp, sd)

avg10<-tapply(c$Otu010, c$exp, mean)
sd10<-tapply(c$Otu010, c$exp, sd)

avg07<-tapply(c$Otu007, c$exp, mean)
sd07<-tapply(c$Otu007, c$exp, sd)

avg02<-tapply(c$Otu002, c$exp, mean)
sd02<-tapply(c$Otu002, c$exp, sd)

avgOther<-tapply(c$Other, c$exp, mean)
sdOther<-tapply(c$Other, c$exp, sd)
  

#paper figure: cdiff titration barplots
par(mar=c(5, 8, 3, 4) +0.1, mgp=c(5, 2, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plotCD<-barplot(avgCD, col="black", log="y", yaxt="n", xaxt="n", ylim=c(1, 1e9), cex.lab=2, ylab=expression(paste(italic("C. difficile"), " CFU/g Feces"))) 
error.bar(plotCD, avgCD, sdCD, lwd=7) #plot SD
ticksy<-seq(0, 9, by=1)
ticksx<-seq(1,7, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, lwd=2, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), cex.axis=1.7, labels=labelsy)
axis(1, las=1, lwd=2, at=c(plotCD), labels=c("0.1 mg/ml", "0.3 mg/ml", "0.625 mg/ml"), cex.axis=1.7)
text(plotCD[2], 8*10^8, label="Vancomycin Titrations", cex=2.2, pos=1)


par(mfrow=c(5, 2)) 

plot1<-barplot(avg1, col="darkblue", cex.names=.6) #how can i change this to log scale
error.bar(plot1, avg1, sd1) #plot SD

plot05<-barplot(avg05, col="darkgreen", cex.names=.6)
error.bar(plot05, avg05, sd05) #plot SD

plot04<-barplot(avg04, col="darkorange", cex.names=.6)
error.bar(plot04, avg04, sd04) #plot SD

plot3<-barplot(avg3, col="darkred",  cex.names=.6)
error.bar(plot3, avg3, sd3) #plot SD

plot10<-barplot(avg10, col="blue", cex.names=.6)
error.bar(plot10, avg10, sd10) #plot SD

plot07<-barplot(avg07, col="green", cex.names=.6)
error.bar(plot07, avg07, sd07) #plot SD

plot02<-barplot(avg02, col="orange", cex.names=.6)
error.bar(plot02, avg02, sd02) #plot SD

plotOther<-barplot(avgOther, col="red", cex.names=.6)
error.bar(plotOther, avgOther, sdOther) #plot SD

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Bacteroidales", "Enterobacteriaceae", "Akkermansia", "Lactobacillales")
col<-c("darkblue","darkgreen", "darkorange", "darkred")
#pch<-c(16, 16, 16)
legend("center", legend=leg, pch=16, col=col, cex=2, bty="n", ncol=2)

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Bacillales", "Erysipelotricaceae", "Clostridiales", "Other")
col<-c("blue", "green", "orange", "red")
#pch<-c(16, 16, 16)
legend("center", legend=leg, pch=16, col=col, cex=2, bty="n", ncol=2)

par(mfrow=c(1, 1)) 


###metro delayed barcharts (day zeros)-reasoning for doing experiment
b<-read.delim("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.2.subsample.2.pick.metro.relabund.d0s.txt", header=T)

avg1<-tapply(b$Otu001, b$order, mean)
sd1<-tapply(b$Otu001, b$order, sd)

avg11<-tapply(b$Otu011, b$order, mean)
sd11<-tapply(b$Otu011, b$order, sd)

avg02<-tapply(b$Otu002, b$order, mean)
sd02<-tapply(b$Otu002, b$order, sd)

avg03<-tapply(b$Otu003, b$order, mean)
sd03<-tapply(b$Otu003, b$order, sd)

avg7<-tapply(b$Otu007, b$order, mean)
sd7<-tapply(b$Otu007, b$order, sd)

avg05<-tapply(b$Otu005, b$order, mean)
sd05<-tapply(b$Otu005, b$order, sd)

avg16<-tapply(b$Otu016, b$order, mean)
sd16<-tapply(b$Otu016, b$order, sd)

avg06<-tapply(b$Otu006, b$order, mean)
sd06<-tapply(b$Otu006, b$order, sd)

avg04<-tapply(b$Otu004, b$order, mean)
sd04<-tapply(b$Otu004, b$order, sd)

#avg14<-tapply(b$Otu014, b$order, mean)
#sd14<-tapply(b$Otu014, b$order, sd)

avgCD<-tapply(b$nextDayCFU, b$order, mean)
sdCD<-tapply(b$nextDayCFU, b$order, sd)

#inv simpson for the metro d0s, do a test of significance
t.test(b$invsimpson~b$order)
t.test(b$Otu001~b$order)
t.test(b$Otu011~b$order)
t.test(b$Otu002~b$order)
t.test(b$Otu003~b$order) #not significantly diff
t.test(b$Otu007~b$order) 
t.test(b$Otu005~b$order)
t.test(b$Otu016~b$order)
t.test(b$Otu006~b$order) #not significantly diff
t.test(b$Otu004~b$order)
t.test(b$nextDayCFU~b$order)

avgInvSimpson<-tapply(b$invsimpson, b$order, mean)
sdInvSimpson<-tapply(b$invsimpson, b$order, sd)

plotInvSimpson<-barplot(avgInvSimpson, col="black", cex.axis=1.3, xaxt="n", ylim=c(0,7))
error.bar(plotInvSimpson, avgInvSimpson, sdInvSimpson) #plot SD
axis(1, at=plotInvSimpson, labels=c("D0", "D0 Recovered"), cex.axis=1.7)

par(mar=c(4, 9, 3, 4) +0.1, mgp=c(6, 2, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)
plotCD<-barplot(avgCD, col="black", log="y", yaxt="n", xaxt="n", ylim=c(1, 1e9), cex.lab=2.5, ylab="C. difficile CFU/g Feces") #how can i change this to log scale
error.bar(plotCD, avgCD, sdCD, lwd=7) #plot SD
ticksy<-seq(0, 9, by=1)
ticksx<-seq(1,7, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, lwd=2, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), cex.axis=2, labels=labelsy)
axis(1, las=1, lwd=2, at=c(plotCD), labels=c("D1", "D1 Recovered"), cex.axis=2)


par(mfrow=c(4, 3)) 

par(mar=c(2, 4, 2, 3) +0.1)
plot1<-barplot(avg1, col="darkblue", cex.axis=1.3, xaxt="n") 
error.bar(plot1, avg1, sd1) #plot SD

par(mar=c(2, 3, 2, 3) +0.1)
plot03<-barplot(avg03, col="steelblue", cex.axis=1.3, xaxt="n")
error.bar(plot03, avg03, sd03) #plot SD

plot11<-barplot(avg11, col="darkgreen", cex.axis=1.3, xaxt="n")
error.bar(plot11, avg11, sd11) #plot SD

par(mar=c(2, 4, 2, 3) +0.1)
plot02<-barplot(avg02, col="darkred", cex.axis=1.3, xaxt="n")
error.bar(plot02, avg02, sd02) #plot SD

par(mar=c(2, 3, 2, 3) +0.1)
plot05<-barplot(avg05, col="orange2", cex.axis=1.3, xaxt="n")
error.bar(plot05, avg05, sd05) #plot SD

plot04<-barplot(avg04, col="gold1", cex.axis=1.3, xaxt="n")
error.bar(plot04, avg04, sd04) #plot SD

par(mar=c(2, 4, 2, 3) +0.1)
plot16<-barplot(avg16, col="deeppink3", cex.axis=1.3, xaxt="n")
error.bar(plot16, avg16, sd16) #plot SD
axis(1, at=plot16, labels=c("D0", "D0 Recovered"), cex.axis=1.7)

par(mar=c(2, 3, 2, 3) +0.1)
plot7<-barplot(avg7, col="blueviolet", cex.axis=1.3, xaxt="n")
error.bar(plot7, avg7, sd7) #plot SD
axis(1, at=plot7, labels=c("D0", "D0 Recovered"), cex.axis=1.7)

plot06<-barplot(avg06, col="black", cex.axis=1.3, xaxt="n")
error.bar(plot06, avg06, sd06) #plot SD
axis(1, at=plot06, labels=c("D0", "D0 Recovered"), cex.axis=1.7)

#plot14<-barplot(avg14, col="gold1", cex.names=1)
#error.bar(plot14, avg14, sd14) #plot SD

#plotOther<-barplot(avgOther, col="darkorange1", cex.names=1)
#error.bar(plotOther, avgOther, sdOther) #plot SD

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Porphyromonadaceae", "Bacteroides","Erysipelotrichaceae" )
col<-c("darkblue",  "steelblue", "darkgreen" )
#pch<-c(16, 16, 16)
legend("left", legend=leg, pch=15, col=col,  cex=1.7, bty="n", ncol=1)

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Lachnospiraceae",  "Lactobacillus", "Ruminococcaceae" )
col<-c("darkred", "orange2", "gold1" )
#pch<-c(16, 16, 16)
legend("left", legend=leg, pch=15, col=col,  cex=1.7, bty="n", ncol=1)

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Bifidobacterium",  "Enterobacteriaceae","Akkermansia" )
col<-c("deeppink3", "blueviolet", "black" )
#pch<-c(16, 16, 16)
legend("left", legend=leg, pch=15, col=col,  cex=1.7, bty="n", ncol=1)


par(mfrow=c(1, 1)) 




###amp delayed barcharts (over time)--shows recovery days and reasoning for doing delays in first place
b<-read.delim("~/Desktop/mothur/abxD01/barcharts/abxD01.final.tx.2.subsample.2.pick.ampD.relabund.top3.txt", header=T)

avg1<-tapply(b$Otu001, b$day, mean)
sd1<-tapply(b$Otu001, b$day, sd)

avg11<-tapply(b$Otu011, b$day, mean)
sd11<-tapply(b$Otu011, b$day, sd)

avg02<-tapply(b$Otu002, b$day, mean)
sd02<-tapply(b$Otu002, b$day, sd)

avg03<-tapply(b$Otu003, b$day, mean)
sd03<-tapply(b$Otu003, b$day, sd)

avg7<-tapply(b$Otu007, b$day, mean)
sd7<-tapply(b$Otu007, b$day, sd)

avg05<-tapply(b$Otu005, b$day, mean)
sd05<-tapply(b$Otu005, b$day, sd)

avg19<-tapply(b$Otu019, b$day, mean)
sd19<-tapply(b$Otu019, b$day, sd)

avg06<-tapply(b$Otu006, b$day, mean)
sd06<-tapply(b$Otu006, b$day, sd)

avg28<-tapply(b$Otu028, b$day, mean)
sd28<-tapply(b$Otu028, b$day, sd)

avg14<-tapply(b$Otu014, b$day, mean)
sd14<-tapply(b$Otu014, b$day, sd)

avg24<-tapply(b$Otu024, b$day, mean)
sd24<-tapply(b$Otu024, b$day, sd)

avg18<-tapply(b$Otu018, b$day, mean)
sd18<-tapply(b$Otu018, b$day, sd)

avgOther<-tapply(b$Other, b$day, mean)
sdOther<-tapply(b$Other, b$day, sd)


par(mfrow=c(5, 3)) 

plot1<-barplot(avg1, col="black", cex.names=1) 
error.bar(plot1, avg1, sd1) #plot SD

plot11<-barplot(avg11, col="darkblue", cex.names=1) 
error.bar(plot11, avg11, sd11) #plot SD

plot02<-barplot(avg02, col="darkgreen", cex.names=1)
error.bar(plot02, avg02, sd02) #plot SD

plot03<-barplot(avg03, col="darkgoldenrod3", cex.names=1)
error.bar(plot03, avg03, sd03) #plot SD

plot7<-barplot(avg7, col="darkorange4", cex.names=1)
error.bar(plot7, avg7, sd7) #plot SD

plot05<-barplot(avg05, col="firebrick3", cex.names=1)
error.bar(plot05, avg05, sd05) #plot SD

plot19<-barplot(avg19, col="grey", cex.names=1)
error.bar(plot19, avg19, sd19) #plot SD

plot06<-barplot(avg06, col="blue", cex.names=1)
error.bar(plot06, avg06, sd06) #plot SD

plot28<-barplot(avg28, col="green2", cex.names=1)
error.bar(plot28, avg28, sd28) #plot SD

plot14<-barplot(avg14, col="gold1", cex.names=1)
error.bar(plot14, avg14, sd14) #plot SD

plot24<-barplot(avg24, col="darkorange1", cex.names=1)
error.bar(plot24, avg24, sd24) #plot SD

plot18<-barplot(avg18, col="red", cex.names=1)
error.bar(plot18, avg18, sd18) #plot SD

plotOther<-barplot(avgOther, col="blueviolet", cex.names=1)
error.bar(plotOther, avgOther, sdOther) #plot SD

plot(1, type="n", axes=F, xlab="", ylab="")


leg<-c("Porphyromonadaceae", "Erysipelotrichaceae", "Lachnospiraceae", "Bacteroides", "Enterobacteriaceae", "Lactobacillus")
col<-c("black", "darkblue","darkgreen", "darkgoldenrod3", "darkorange4", "firebrick3")
#pch<-c(16, 16, 16)
legend("right", legend=leg, pch=15, col=col,  cex=1.1, bty="n", ncol=2)

plot(1, type="n", axes=F, xlab="", ylab="")

leg<-c("Clostridium sensu stricto", "Akkermansia", "Enterococcus", "Alistipes","Pseudomonas", "Staphylococcus", "Other")
col<-c("grey", "blue", "green2", "gold1", "darkorange1", "red", "blueviolet")
#pch<-c(16, 16, 16)
legend("left", legend=leg, pch=15, col=col,  cex=1.1, bty="n", ncol=2)

par(mfrow=c(1, 1)) 



###correlation plots--relabund vs Cdiff CFU
##########################################
###correlation plots--relabund vs Cdiff CFU
c<-read.csv("~/Desktop/mothur/abxD01/correlation/graphs/abxD01.final.tx.2.subsample.2.pick.relabund.D01.csv", header=T)

#porphyro
amp.5d<-c[c$exp=="amp.5d", c(10,34)] # 10 is oty, second number is the cdiff
amp.5<-c[c$exp=="amp.5", c(10,34)]
cef.1<-c[c$exp=="cef.1", c(10,34)]
cef.3<-c[c$exp=="cef.3", c(10,34)]
cef.5<-c[c$exp=="cef.5", c(10,34)]
cipro<-c[c$exp=="cipro", c(10,34)]
clinda<-c[c$exp=="clinda", c(10,34)]
metro1d<-c[c$exp=="metro1d", c(10,34)]
metro1<-c[c$exp=="metro1", c(10,34)]
strep.1<-c[c$exp=="strep.1", c(10,34)]
strep.5<-c[c$exp=="strep.5", c(10,34)]
strep5<-c[c$exp=="strep5", c(10,34)]
vanc.1<-c[c$exp=="vanc.1", c(10,34)]
vanc.3<-c[c$exp=="vanc.3", c(10,34)]
vanc1.625<-c[c$exp=="vanc1.625", c(10,34)]


plot(cipro, pch=16, col="red", xlim=c(0, 1), log="y", ylim=c(1, 1e+10), cex=3) #want to change this to include every log
points(clinda,  cex=3, col="dark orange", pch=16)
points(vanc1.625,  cex=3, col="gold", pch=16)
points(cef.5,  cex=3, col="forestgreen", pch=16)
points(strep5,  cex=3, col="darkblue", pch=16)
points(amp.5,  cex=3, col="deeppink3", pch=16)
points(metro1,  cex=3, col="green2", pch=16)

points(metro1d,  cex=3, col="darkturquoise", pch=16)
points(amp.5d,  cex=3, col="blueviolet", pch=16)

points(vanc.3,  cex=2, col="gold", pch=16)
points(cef.3,  cex=2, col="forestgreen", pch=16)
points(strep.5,  cex=2, col="darkblue", pch=16)
points(strep.1,  cex=1, col="darkblue", pch=16)
points(cef.1,  cex=1, col="forestgreen", pch=16)
points(vanc.1,  cex=1, col="gold", pch=16)

leg<-c("Cipro", "Clinda", "Vanc", "Cef", "Strep", "Amp", "Amp delay", "Metro", "Metro delay")
col<-c("red","darkorange", "gold", "forestgreen", "darkblue", "deeppink3", "blueviolet", "green2", "darkturquoise")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=1, bty="n")

##want to make a loop that will make all the graphs im interested in and save them all

c<-read.csv("~/Desktop/mothur/abxD01/correlation/graphs/abxD01.final.tx.2.subsample.2.pick.relabund.D01.csv", header=T)

#alistipes
amp.5d<-c[c$exp=="amp.5d", c(13,34)] # 13 is CFU, second number is the OTU column
amp.5<-c[c$exp=="amp.5", c(13,34)]
cef.1<-c[c$exp=="cef.1", c(13,34)]
cef.3<-c[c$exp=="cef.3", c(13,34)]
cef.5<-c[c$exp=="cef.5", c(13,34)]
cipro<-c[c$exp=="cipro", c(13,34)]
clinda<-c[c$exp=="clinda", c(13,34)]
metro1d<-c[c$exp=="metro1d", c(13,34)]
metro1<-c[c$exp=="metro1", c(13,34)]
strep.1<-c[c$exp=="strep.1", c(13,34)]
strep.5<-c[c$exp=="strep.5", c(13,34)]
strep5<-c[c$exp=="strep5", c(13,34)]
vanc.1<-c[c$exp=="vanc.1", c(13,34)]
vanc.3<-c[c$exp=="vanc.3", c(13,34)]
vanc1.625<-c[c$exp=="vanc1.625", c(13,34)]


plot(cipro, pch=16, col="red", xlim=c(0, .09), log="y", ylim=c(1, 1e+10), cex=3) #want to change this to include every log
points(clinda,  cex=3, col="dark orange", pch=16)
points(vanc1.625,  cex=3, col="gold", pch=16)
points(cef.5,  cex=3, col="forestgreen", pch=16)
points(strep5,  cex=3, col="darkblue", pch=16)
points(amp.5d,  cex=3, col="blueviolet", pch=16)
points(amp.5,  cex=3, col="deeppink3", pch=16)
points(metro1d,  cex=3, col="darkturquoise", pch=16)
points(metro1,  cex=3, col="green2", pch=16)
points(vanc.3,  cex=2, col="gold", pch=16)
points(cef.3,  cex=2, col="forestgreen", pch=16)
points(strep.5,  cex=2, col="darkblue", pch=16)
points(strep.1,  cex=1, col="darkblue", pch=16)
points(cef.1,  cex=1, col="forestgreen", pch=16)
points(vanc.1,  cex=1, col="gold", pch=16)

leg<-c("Cipro", "Clinda", "Vanc", "Cef", "Strep", "Amp", "Amp delay", "Metro", "Metro delay")
col<-c("red","darkorange", "gold", "forestgreen", "darkblue", "deeppink3", "blueviolet", "green2", "darkturquoise")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=1, bty="n")


#ecoli
c<-read.csv("~/Desktop/mothur/abxD01/correlation/graphs/abxD01.final.tx.2.subsample.2.pick.relabund.D01.csv", header=T)

amp.5d<-c[c$exp=="amp.5d", c(15,34)] # 15 is CFU, second number is the OTU column
amp.5<-c[c$exp=="amp.5", c(15,34)]
cef.1<-c[c$exp=="cef.1", c(15,34)]
cef.3<-c[c$exp=="cef.3", c(15,34)]
cef.5<-c[c$exp=="cef.5", c(15,34)]
cipro<-c[c$exp=="cipro", c(15,34)]
clinda<-c[c$exp=="clinda", c(15,34)]
metro1d<-c[c$exp=="metro1d", c(15,34)]
metro1<-c[c$exp=="metro1", c(15,34)]
strep.1<-c[c$exp=="strep.1", c(15,34)]
strep.5<-c[c$exp=="strep.5", c(15,34)]
strep5<-c[c$exp=="strep5", c(15,34)]
vanc.1<-c[c$exp=="vanc.1", c(15,34)]
vanc.3<-c[c$exp=="vanc.3", c(15,34)]
vanc1.625<-c[c$exp=="vanc1.625", c(15,34)]


plot(cipro, pch=16, col="red", xlim=c(0, 1), log="y", ylim=c(1, 1e+10), cex=3) #want to change this to include every log
points(clinda,  cex=3, col="dark orange", pch=16)
points(vanc1.625,  cex=3, col="gold", pch=16)
points(cef.5,  cex=3, col="forestgreen", pch=16)
points(strep5,  cex=3, col="darkblue", pch=16)
points(amp.5d,  cex=3, col="blueviolet", pch=16)
points(amp.5,  cex=3, col="deeppink3", pch=16)
points(metro1d,  cex=3, col="darkturquoise", pch=16)
points(metro1,  cex=3, col="green2", pch=16)
points(vanc.3,  cex=2, col="gold", pch=16)
points(cef.3,  cex=2, col="forestgreen", pch=16)
points(strep.5,  cex=2, col="darkblue", pch=16)
points(strep.1,  cex=1, col="darkblue", pch=16)
points(cef.1,  cex=1, col="forestgreen", pch=16)
points(vanc.1,  cex=1, col="gold", pch=16)

leg<-c("Cipro", "Clinda", "Vanc", "Cef", "Strep", "Amp", "Amp delay", "Metro", "Metro delay")
col<-c("red","darkorange", "gold", "forestgreen", "darkblue", "deeppink3", "blueviolet", "green2", "darkturquoise")
#pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=16, col=col,  cex=1, bty="n")

###OTUS
c<-read.csv("~/Desktop/mothur/abxD01/correlation/graphs/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.relabund.D01.csv", header=T)



#data breakdown
amp.5d<-c[c$exp=="amp.5d"&c$day=="0", c(39,2)] # 2 is cdiff CFU, first number is the OTU
amp.5<-c[c$exp=="amp.5"&c$day=="0", c(39,2)]
cef.1<-c[c$exp=="cef.1"&c$day=="0", c(39,2)]
cef.3<-c[c$exp=="cef.3"&c$day=="0", c(39,2)]
cef.5<-c[c$exp=="cef.5"&c$day=="0", c(39,2)]
cipro<-c[c$exp=="cipro"&c$day=="0", c(39,2)]
clinda<-c[c$exp=="clinda"&c$day=="0", c(39,2)]
metro1d<-c[c$exp=="metro1d"&c$day=="0", c(39,2)]
metro1<-c[c$exp=="metro1"&c$day=="0", c(39,2)]
strep.1<-c[c$exp=="strep.1"&c$day=="0", c(39,2)]
strep.5<-c[c$exp=="strep.5"&c$day=="0", c(39,2)]
strep5<-c[c$exp=="strep5"&c$day=="0", c(39,2)]
vanc.1<-c[c$exp=="vanc.1"&c$day=="0", c(39,2)]
vanc.3<-c[c$exp=="vanc.3"&c$day=="0", c(39,2)]
vanc1.625<-c[c$exp=="vanc.625"&c$day=="1", c(39,2)]
untreated<-c[c$preABX=="preAbx", c(39,2)]



#correlation of top dose
plot(untreated, pch=16, col="black", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=2, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(clinda,  cex=2, col="darkorange", pch=16)
points(vanc1.625,  cex=2, col="darkgoldenrod1", pch=16)
points(cef.5,  cex=2, col="forestgreen", pch=16)
points(strep5,  cex=2, col="darkblue", pch=16)
points(amp.5,  cex=2, col="deeppink3", pch=16)
points(metro1,  cex=2, col="blueviolet", pch=16)
points(cipro, cex=2, col="red", pch=16)
leg<-c("Ciprofloxacin 10 mg/kg", "Clindamycin 10 mg/kg", "Vancomycin 0.625 mg/ml", "Cefoperazone 0.5 mg/ml", "Streptomycin 0.5 mg/ml", "Ampicillin 0.5 mg/ml", "Metronidazole 1 mg/ml", "Untreated")
col<-c("red","darkorange", "darkgoldenrod1", "forestgreen", "darkblue", "deeppink3", "blueviolet", "black")
legend(.08,10000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend(.00007,10000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("right", legend=leg, pch=16, col=col,  cex=1, bty="n")



#correlation by vancomycin treatment
plot(vanc1.625, pch=16, col="darkgoldenrod1", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(vanc.3,  cex=2, col="darkgoldenrod3", pch=16)
points(vanc.1,  cex=1, col="darkgoldenrod4", pch=16)
leg<-c("Vancomycin 0.625 mg/ml", "Vancomycin 0.3 mg/ml", "Vancomycin 0.1 mg/ml")
col<-c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")



#correlation by cef treatment
plot(cef.5, pch=16, col="green", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(cef.3,  cex=2, col="green4", pch=16)
points(cef.1,  cex=1, col="darkgreen", pch=16)
leg<-c("Cefoperazone 0.5 mg/ml", "Cefoperazone 0.3 mg/ml", "Cefoperazone 0.1 mg/ml")
col<-c("green", "green4", "darkgreen")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(.08,10000000, legend=leg, pch=16, col=col,  cex=1, bty="n")


#correlation by strep treatment
plot(strep5, pch=16, col="darkturquoise", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(strep.5,  cex=2, col="dodgerblue2", pch=16)
points(strep.1,  cex=1, col="darkblue", pch=16)
leg<-c("Streptomycin 5 mg/ml", "Streptomycin 0.5 mg/ml", "Streptomycin 0.1 mg/ml")
col<-c("darkturquoise", "dodgerblue2", "darkblue")
legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")


#correlation by amp treatment
plot(amp.5, pch=16, col="deeppink3", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(amp.5d,  cex=3, col="lightpink3", pch=16)
leg<-c("Ampicillin 0.5 mg/ml", "Ampicillin 0.5 mg/ml delayed")
col<-c("deeppink3", "lightpink3")
legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")


#correlation by metro treatment
plot(metro1, pch=16, col="purple4", log="xy", xlim=c(0.0001, 1), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Staphylococcus OTU25") #want to change this to include every log
axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
#axis(1, las=1, at=seq(0,1, by=.1))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(metro1d,  cex=3, col="mediumpurple1", pch=16)
leg<-c("Metronidazole 1 mg/ml", "Metronidazole 1 mg/ml delayed")
col<-c("purple4", "mediumpurple1")
legend(.055,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")







#####INVERSE SIMPSON correlation graphs
####cfu by inverse simpson
c<-read.csv("~/Desktop/mothur/abxD01/correlation/graphs/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.relabund.D01.csv", header=T)

#data breakdown
amp.5d<-c[c$exp=="amp.5d"&c$day=="0", c(14,2)] # 2 is cdiff CFU, first number is the OTU
amp.5<-c[c$exp=="amp.5"&c$day=="0", c(14,2)]
cef.1<-c[c$exp=="cef.1"&c$day=="0", c(14,2)]
cef.3<-c[c$exp=="cef.3"&c$day=="0", c(14,2)]
cef.5<-c[c$exp=="cef.5"&c$day=="0", c(14,2)]
cipro<-c[c$exp=="cipro"&c$day=="0", c(14,2)]
clinda<-c[c$exp=="clinda"&c$day=="0", c(14,2)]
metro1d<-c[c$exp=="metro1d"&c$day=="0", c(14,2)]
metro1<-c[c$exp=="metro1"&c$day=="0", c(14,2)]
strep.1<-c[c$exp=="strep.1"&c$day=="0", c(14,2)]
strep.5<-c[c$exp=="strep.5"&c$day=="0", c(14,2)]
strep5<-c[c$exp=="strep5"&c$day=="0", c(14,2)]
vanc.1<-c[c$exp=="vanc.1"&c$day=="0", c(14,2)]
vanc.3<-c[c$exp=="vanc.3"&c$day=="0", c(14,2)]
vanc1.625<-c[c$exp=="vanc.625"&c$day=="1", c(14,2)]
untreated<-c[c$preABX=="preAbx", c(14,2)]



#correlation of top dose
par(mar=c(9, 9, 3, 3) +0.1, mgp=c(6, 2, 0)) #default is 5.1 4.1 4.1 2.1, bot/left/top/right, also default mgp is c(3,1,0)

plot(untreated, pch=16, col="black", log="y", xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3,cex.lab=2.3, yaxt="n",  xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, cex.axis=2, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1,cex.axis=2, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(clinda,  cex=3, col="darkorange", pch=16)
points(vanc1.625,  cex=3, col="darkgoldenrod1", pch=16)
points(cef.5,  cex=3, col="forestgreen", pch=16)
points(strep5,  cex=3, col="darkblue", pch=16)
points(amp.5,  cex=3, col="deeppink3", pch=16)
points(metro1,  cex=3, col="blueviolet", pch=16)
points(cipro, cex=3, col="red", pch=16)
leg<-c("Ciprofloxacin 10 mg/kg", "Clindamycin 10 mg/kg", "Vancomycin 0.625 mg/ml", "Cefoperazone 0.5 mg/ml", "Streptomycin 0.5 mg/ml", "Ampicillin 0.5 mg/ml", "Metronidazole 1 mg/ml", "Untreated")
col<-c("red","darkorange", "darkgoldenrod1", "forestgreen", "darkblue", "deeppink3", "blueviolet", "black")
legend(17,1000000000, legend=leg, pch=16, col=col,  cex=1.3, bty="n")
#legend(.00007,10000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("top", legend=leg, pch=16, col=col,  cex=1, bty="n")



#correlation by vancomycin treatment
plot(vanc1.625, pch=16, col="darkgoldenrod1", log="y",xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(vanc.3,  cex=2, col="darkgoldenrod3", pch=16)
points(vanc.1,  cex=1, col="darkgoldenrod4", pch=16)
leg<-c("Vancomycin 0.625 mg/ml", "Vancomycin 0.3 mg/ml", "Vancomycin 0.1 mg/ml")
col<-c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(18,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")



#correlation by cef treatment
plot(cef.5, pch=16, col="green", log="y", xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(cef.3,  cex=2, col="green4", pch=16)
points(cef.1,  cex=1, col="darkgreen", pch=16)
leg<-c("Cefoperazone 0.5 mg/ml", "Cefoperazone 0.3 mg/ml", "Cefoperazone 0.1 mg/ml")
col<-c("green", "green4", "darkgreen")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(18,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")




#correlation by strep treatment
plot(strep5, pch=16, col="darkturquoise", log="y",xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(strep.5,  cex=2, col="dodgerblue2", pch=16)
points(strep.1,  cex=1, col="darkblue", pch=16)
leg<-c("Streptomycin 5 mg/ml", "Streptomycin 0.5 mg/ml", "Streptomycin 0.1 mg/ml")
col<-c("darkturquoise", "dodgerblue2", "darkblue")
#legend(.07,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(18,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")




#correlation by amp treatment
plot(amp.5, pch=16, col="deeppink3", log="y",xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(amp.5d,  cex=3, col="lightpink3", pch=16)
leg<-c("Ampicillin 0.5 mg/ml", "Ampicillin 0.5 mg/ml delayed")
col<-c("deeppink3", "lightpink3")
#legend(.05,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(18,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")


#correlation by metro treatment
plot(metro1, pch=16, col="purple4", log="y",xlim=c(0, 30), ylim=c(1, 1e+9),  cex=3, yaxt="n", xaxt="n", ylab="C. difficile CFU/g feces", xlab="Inverse Simpson Index") #want to change this to include every log
#axis(1, las=1, at=c(.0001, .001, .01, .1, 1), labels=c(0, .001, .01, .1, 1))
axis(1, las=1, at=seq(0,30, by=5))
ticksy<-seq(0, 9, by=1)
labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)
points(metro1d,  cex=3, col="mediumpurple1", pch=16)
leg<-c("Metronidazole 1 mg/ml", "Metronidazole 1 mg/ml delayed")
col<-c("purple4", "mediumpurple1")
#legend(.04,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")
#legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")
legend(18,1000000000, legend=leg, pch=16, col=col,  cex=1, bty="n")



#################################
##paper figure: Heatmap for correlation values
library(gplots)
correl<- read.csv("~/Desktop/mothur/abxD01/correlation/correl_heatmap_0.01rel_topdose2_newtitr.csv", header=T)

row.names(correl)<-correl$OTU
correl_matrix<-data.matrix(correl)
correl_matrix<-correl_matrix[,-1]
correl_matrix<-data.matrix(correl_matrix)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
breaks = c(seq(-1,-.33,length=100),seq(-.33,.33,length=100),seq(.33,1,length=100))#changed to equal numbers

#par(mar=c(4, 6, 3, 4.5) +0.1) #default is 5.1 4.1 4.1 2.1

side<-scan("~/Desktop/mothur/abxD01/correlation/correl_heatmapSIDE_0.01rel_topdose2_newtitr.csv", what=",", sep=",")
side2<-read.csv("~/Desktop/mothur/abxD01/correlation/correl_heatmapSIDE_0.01rel_topdose2_newtitr.csv", header=T)
side3<-as.matrix(side2)
#lmat=rbind(c(4,3), c(2,1))
#lmat=rbind(c(4,3), c(1,2))

#lmat order: 1)row dendrogram, 2) heatmap, 3)col?, 4)key?--key disappeared
lmat = rbind(3:4, 1:2)
lmat = rbind(c(2,4), c(1,3))
lmat = rbind(c(0,3), c(2,1))
lmat=NULL
lwid=NULL
lhei=NULL
lwid = c(.5, 2, 4)
lhei = c(.7, 4)
  
correl_heatmap<-heatmap.2(correl_matrix, scale="none",  col=my_palette, breaks=breaks, density.info="none", lhei=lhei,lwid=lwid, RowSideColors=side3[,2], labRow=side3[,1], trace="none",margins=c(8, 25), dendrogram="none",  Rowv=FALSE, Colv=FALSE,  na.color="black", cexCol=1.5,  key=TRUE, keysize=1, cexRow=1) #then plot was 7x8in, portrait
#correl_heatmap<-heatmap.2(t(correl_matrix), scale="none", col=my_palette, density.info=c("none"), trace=c("none"), breaks=colors, Rowv=NA, Colv=NA,  na.color="light gray", cexCol=0.8, key=TRUE, keysize=1, cexRow=.7) #then plot was 7x8in, portrait

#combine with a histogram along the y axis (OTU list) 
c<- read.csv("~/Desktop/mothur/abxD01/correlation/correl_heatmap_topdose_titrations.csv", header=T)
library(sm)
library(vioplot)
porph<-c$cor[c$grouping=="porphyro"]
vioplot(porph)



#par(mfrow=c(1, 1)) 



###box and whisker correlation values
#####################################
###box and whisker correlation values
c<-read.csv("~/Desktop/mothur/abxD01/correlation/corr_pos_topdose_boxwhisker.csv", header=T)
library(gplots)
boxplot.n(c$cor~c$last,  shrink=2, xaxt="n", cex.axis=2, ylim=c(.3, .55))
axis(1, at=1:5)

library(vioplot)
porph<-c$cor[c$name=="Porphyromonadaceae"]
vioplot(porph)


###stripchart correlation values
#####################################
###stripchart correlation values
par(mfrow=c(1, 1))
c<-read.csv("~/Desktop/mothur/abxD01/correlation/corr_pos_topdose_boxwhisker.csv", header=T)
par(mar=c(4, 6, 3, 4.5) +0.1) #default is 5.1 4.1 4.1 2.1
stripchart(c$cor~c$last, vertical=TRUE, ylab="Correlation", cex=3, pch=15, col="dark red", cex.lab=2.5, xaxt="n", yaxt="n")
axis(1, cex.axis=1.3, at=(1:5), labels=c("E. coli", "Firmicutes", "Pseudomonas", "Rhodobacter", "Streptococcus") )
axis(2, cex.axis=2, at=seq(0.3,0.55, by=.05))


c<-read.csv("~/Desktop/mothur/abxD01/correlation/corr_neg_topdose_boxwhisker.csv", header=T)
par(mar=c(14, 6, 2, 2) +0.1) #default is 5.1 4.1 4.1 2.1
stripchart(c$cor~c$name, vertical=TRUE, ylab="Correlation", cex=3, pch=15, col="darkblue", xaxt="n", method="jitter", cex.lab=2.5,yaxt="n")
axis(1, cex.axis=1, at=(1:nlevels(c$name)), labels=FALSE)
text(1:8, par("usr")[3]-0.03, label=levels(c$name), xpd=NA, pos=2, srt=45, cex=1.5)
text(1:8, par("usr")[3]-0.005, label=c("n=3", "n=1", "n=1", "n=7", "n=16", "n=2", "n=14", "n=7"), xpd=NA, pos=1, cex=1.25)
axis(2, cex.axis=1.5, at=seq(-0.3,-0.8, by=-.05))

##Paper Figure: correlations for topdose 
c<-read.csv("~/Desktop/mothur/abxD01/correlation/corr_allSig_topdose2_stripchart.csv", header=T)
labels=c("Lachnospiraceae", "Ruminococcaceae", "Clostridia", "Lactobacillales", "Firmicutes", "Bacillales", "Porphyromonadaceae", "Bacteroidales", "Bacteroidetes", "Actinobacteria", "Proteobacteria", "Anaeroplasma", "Deinococcus", "Unclassified")
ns=c("n=18", "n=9", "n=8", "n=4", "n=2", "n=1", "n=15", "n=4", "n=1", "n=4", "n=3", "n=1", "n=1", "n=3")

par(mar=c(12, 7, 2, 2) +0.1, mgp=c(5, 1, 0)) #default is 5.1 4.1 4.1 2.1 [bottom, left, top, right space], mgp=c(3, 1, 0) [label line location for x/y location labels, tick mark labels location, tick mark locations]
stripchart(c$cor~c$order, vertical=TRUE, ylab="Correlation", ylim=c(.85, -.85), cex.axis=1.1, cex=2, pch=21, lwd=3, col="black", xaxt="n", method="jitter",  cex.lab=1.7, yaxt="n")
axis(1, cex.axis=1, at=(1:nlevels(c$name)), labels=FALSE)
text(1:nlevels(c$name)+.1, par("usr")[3]-0.16, label=labels, xpd=NA, pos=2, srt=45, cex=1.2)
text(1:nlevels(c$name), par("usr")[3]-0.01, label=ns, xpd=NA, pos=1, cex=.9)
axis(2, cex.axis=1.1, at=seq(-0.85,+0.85, by=.1), las=1)
abline(v=c(1:nlevels(c$name)), col="dark gray")
#abline(h=c(.3, -.3), col="dark gray", lty="dashed", lwd=4) #change to something else
abline(h=0, col="black", lwd=2)


###2d NMDS all abx, d0
#######################
###2d NMDS all abx, d0

c<-read.delim("~/Desktop/mothur/abxD01/abxD01.final.pick.tx.1.subsample.thetayc.1.lt.nmds", header=T)

clinda<-c[c$abx=="clinda"&c$day=="0", c(2:3)]
vanc6<-c[c$abx=="vanc"&c$day=="0"&c$dose=="0.625", c(2:3)]
vanc3<-c[c$abx=="vanc"&c$day=="0"&c$dose=="0.3", c(2:3)]
vanc1<-c[c$abx=="vanc"&c$day=="0"&c$dose=="0.1", c(2:3)]
cipro<-c[c$abx=="cipro"&c$day=="0", c(2:3)]
amp6d<-c[c$abx=="amp"&c$day=="0"&c$time=="delayed", c(2:3)]
amp1d<-c[c$abx=="amp"&c$day=="0"&c$time=="on-time", c(2:3)]
cef5<-c[c$abx=="cef"&c$day=="0"&c$dose=="0.5", c(2:3)]
cef3<-c[c$abx=="cef"&c$day=="0"&c$dose=="0.3", c(2:3)]
cef1<-c[c$abx=="cef"&c$day=="0"&c$dose=="0.1", c(2:3)]
metro6d<-c[c$abx=="metro"&c$day=="0"&c$time=="delayed", c(2:3)]
metro1d<-c[c$abx=="metro"&c$day=="0"&c$time=="on-time", c(2:3)]
strep50<-c[c$abx=="strep"&c$day=="0"&c$dose=="5", c(2:3)]
strep5<-c[c$abx=="strep"&c$day=="0"&c$dose=="0.5", c(2:3)]
strep1<-c[c$abx=="strep"&c$day=="0"&c$dose=="0.1", c(2:3)]
untr<-c[c$preABX=="preAbx", c(2:3)]

plot(cipro,  cex=2, xlim=c(-.9,.9), ylim=c(-.9, .9), pch=16, col="red", ylab="Axis 2", xlab="Axis 1")
points(clinda,  cex=2, col="dark orange", pch=16)
points(vanc6,  cex=2, col="gold", pch=16)
points(vanc3,  cex=1.5, col="gold", pch=16)
points(vanc1,  cex=1, col="gold", pch=16)
points(cef5,  cex=2, col="forestgreen", pch=16)
points(cef3,  cex=1.5, col="forestgreen", pch=16)
points(cef1,  cex=1, col="forestgreen", pch=16)
points(strep50,  cex=2, col="darkblue", pch=16)
points(strep5,  cex=1.5, col="darkblue", pch=16)
points(strep1,  cex=1, col="darkblue", pch=16)
points(amp6d,  cex=2, col="blueviolet", pch=16)
points(amp1d,  cex=2, col="deeppink3", pch=16)
points(metro6d,  cex=2, col="darkturquoise", pch=16)
points(metro1d,  cex=2, col="green2", pch=16)
points(untr,  cex=2, col="black", pch=16)

leg<-c("Cipro", "Clinda", "Vanc", "Cef", "Strep", "Amp 1d", "Amp 6d", "Metro 1d", "Metro 6d")
col<-c("red","darkorange", "gold", "forestgreen", "darkblue", "deeppink3", "blueviolet", "green2", "darkturquoise")
#pch<-c(16, 16, 16)
legend("topleft", legend=leg, pch=16, col=col,  cex=.7, bty="n")

##NMDS for matt, start from above
plot(untr,  cex=2, xlim=c(-.9,.9), ylim=c(-.9, .9), pch=16, col="gray", ylab="Axis 2", xlab="Axis 1")
points(vanc6,  cex=2, col="dark orange", pch=16)
points(cef5,  cex=2, col="forestgreen", pch=16)
points(strep50,  cex=2, col="darkblue", pch=16)
points(amp1d,  cex=2, col="blueviolet", pch=16)
points(metro1d,  cex=2, col="darkturquoise", pch=16)
points(clinda,  cex=2, col="red", pch=16)

leg<-c( "Clindamycin", "Vancomycin", "Cefoperazone", "Streptomycin", "Ampicillin", "Metronidazole", "Untreated")
col<-c("red","darkorange", "forestgreen", "darkblue", "blueviolet", "darkturquoise", "gray")
pch<-c(16, 16, 16, 16, 16, 16, 1)
legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")



###correlation analysis
########################
###correlation analysis

meta <- read.table('~/Desktop/mothur/abxD01/correlation/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.correl.topdose2.txt',header=T)
#meta<-meta[1:96] #change based on number of OTUs above .05%, then add 2 for first two cols
c<-1
otu <- c()
cor.spear <- c()
pval.spear <- c()
cor.ken = c()
pval.ken = c()
for(i in 3:length(meta)){
  otu[c] <- colnames(meta[i])
  cor.spear[c] <- cor.test(meta[,2],meta[,i], method="spearman")$estimate
  pval.spear[c] <- cor.test(meta[,2],meta[,i], method="spearman")$p.value
  cor.ken[c] <- cor.test(meta[,2],meta[,i], method="kendall")$estimate #good to see because kendall handles ties
  pval.ken[c] <- cor.test(meta[,2],meta[,i], method="kendall")$p.value #but only works if this is tao-b and not tao-a which im not sure about
  c <- c+1
}

pval.spear<-p.adjust(pval.spear, method='BH') #adjust for multiple comparisons
pval.ken<-p.adjust(pval.ken, method='BH') #adjust for multiple comparisons

results = NULL
results <- matrix(c(otu, cor.spear, pval.spear, cor.ken, pval.ken), ncol=5)
colnames(results) <- c('otu','corSpear','pvalSpear', "corKen", "pvalKen")
results <- results[order(results[,3]),]  #order by pvalue column=3
#head(results, n=100) #shows first n lines
write.table(results[1:dim(results)[1],], file="~/Desktop/mothur/abxD01/correlation/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.correl.topdose2.results.txt", sep="\t", row.names=FALSE)

#combo<-merge(results_top2, results_less, by="otu")
write.table(combo, file="~/Desktop/mothur/abxD01/correlation/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.correl.topdose2.compare.results.txt", sep="\t")


###Mock Error analysis
#######################

s <- read.table(file="~/Desktop/mothur/abxD01/mock23.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.summary", header=T)
ct <- read.table(file="~/Desktop/mothur/abxD01/mock23.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.pick.count_table", header=T)
rownames(s) <- s[,1]
rownames(ct) <- ct[,1]
no.chim <- s$numparents==1
s.good <- s[no.chim,]
query <- rownames(s.good)
ct.good <- ct[as.character(query),]
s.good[,1]==ct.good[,1]
sum(ct.good$mock3 * s.good$mismatches)/sum(ct.good$mock3 * s.good$total)



##RF analysis
library(randomForest)
a<-read.delim("~/Desktop/mothur/abxD01/rf/abxD01.final.tx.2.subsample.2.pick.shared.rf.topdose2.category.txt", header=T)
#nextDayCFU is continuous therefore want randomForest to use regression and not classification method
a500.rf<- randomForest(cdiff ~ ., data=a, importance=TRUE, proximity=TRUE, ntree=500)
a1000.rf<- randomForest(cdiff ~ ., data=a, importance=TRUE, proximity=TRUE, ntree=1000)

print(a500.rf)
print(a1000.rf)
varImpPlot(a1000.rf) #Gini and mean decrease accuracy
title(main="4category trees=1000")

#rf.topdose.regressionCD<-a1000.rf
#rf.topdose.binomialCD<-a1000.rf
rf.topdose.categoryCD<-a1000.rf


importance(a.rf) #list of all importances by otu, the plot in table form
head(importance(a1000.rf))
MDSplot(a500.rf, a$treatment)
MDSplot(a1000.rf, a$nextDayCFU)
partialPlot(a.rf, a, Otu0003, "IP")
a$cluster<-factor(a$cluster) #for if get this message:
#Error in randomForest.default(m, y, ...) : 
#NA/NaN/Inf in foreign function call (arg 2)
#In addition: Warning message:
#  In randomForest.default(m, y, ...) : NAs introduced by coercion

leg<-c("Case", "Diarrheal Control", "Non-diarrheal Control")
col<-c("red2", "royalblue2", "green4")
pch<-c(16, 16, 16)
legend("topright", legend=leg, pch=pch, col=col, cex=.7)
print(a1000.rf)

##RF model training and testing data
######################################
##RF model training and testing data

topdose<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose.regression.csv", header=T)
training<-topdose[,-1]
newtit<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.csv", header=T)
testing2<-newtit[,-1] #new titration data only

#make data set called training=topdose original, testing=titration data, then also test on ampicillin, and metronidazole
#cdiff in first column, shared in subsequent columns... first do the rf like before then save it and put in the titration
library(randomForest)


#fit the randomforest model
model <- randomForest(nextDayCFU~., 
                      data = training, 
                      importance=TRUE, proximity=TRUE,outscale=TRUE,
                      keep.forest=TRUE, ntree=1000
)
print(model)

####
#Call:
#  randomForest(formula = nextDayCFU ~ ., data = training, importance = TRUE,      proximity = TRUE, outscale = TRUE, keep.forest = TRUE, ntree = 1000) 
#Type of random forest: regression
#Number of trees: 1000
#No. of variables tried at each split: 441

#Mean of squared residuals: 6.306925e+14
#% Var explained: 50.36
#####


#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(model, type=1)

#predict the outcome of the testing data
predicted <- predict(model, newdata=testing2[ ,-1])

# what is the proportion variation explained in the outcome of the testing data?
# i.e., what is 1-(SSerror/SStotal)
actual <- testing2$nextDayCFU
rsq <- 1-sum((actual-predicted)^2)/sum((actual-mean(actual))^2)
print(rsq)


##testing to see how including less untreated samples affects the %incMSE results compared to file with n=+100 for untreated
######################
####testing to see how including less untreated samples affects the %incMSE results compared to file with n=+100 for untreated

topdose_less<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose.regression.lessuntreated.csv", header=T)
training_less<-topdose_less[,-1]

library(randomForest)

#fit the randomforest model
model <- randomForest(nextDayCFU~., 
                      data = training_less,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=1000
)
print(model)

####
#Call:
#  randomForest(formula = nextDayCFU ~ ., data = training_less,      outscale = TRUE, importance = TRUE, proximity = TRUE, keep.forest = TRUE,      ntree = 1000) 
#Type of random forest: regression
#Number of trees: 1000
#No. of variables tried at each split: 441

#Mean of squared residuals: 1.870489e+15
#% Var explained: 22.49

####

#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(model, type=1)

#saved as rf.topdose.regression.CD.nodiversity.lessuntreated.pdf
#didnt see that much difference in the most important OTUs

###Testing the best linear model with the new titration data
######################
###Testing the best linear model with the new titration data
#Top Candidate model includes OTUs 6, 7, 20, 39
#Second candidate model includes OTUs 6, 7, 15, 20, 39

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.15otus.rfnegpos.csv", header=T)
newtit2<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.lm.newtitration.logtrans.filter16mintotal.15otus.csv", header=T)
newtit2<-newtit2[,-1]
actual<-as.data.frame(newtit2[,1])
newtit<-newtit2[,-1]
#newtit1<-newtit2[,-1]

m1<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039, data=topdose)
m2<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00015, data=topdose)
m3<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00013 + Otu00015, data=topdose)
m4<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015, data=topdose)
m5<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00015 + Otu00027, data=topdose)
m6<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00013, data=topdose)
m7<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00015 + Otu00023 + Otu00027, data=topdose)
m8<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00013 + Otu00015 + Otu00027, data=topdose)
m9<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00027, data=topdose)
m10<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00027 + Otu00023, data=topdose)
m11<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00013 + Otu00027, data=topdose)
m12<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00013 + Otu00015 + Otu00023 + Otu00027, data=topdose)
m13<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00003 + Otu00015 + Otu00023 + Otu00027, data=topdose)
m14<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00013 + Otu00023, data=topdose)
m15<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00013 + Otu00023 + Otu00027, data=topdose)
m16<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00003 + Otu00011 + Otu00015 + Otu00023 + Otu00027, data=topdose)
m17<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00003 + Otu00011 + Otu00015 + Otu00013 + Otu00023 + Otu00027, data=topdose)
m18<-lm(nextDayCFU ~ Otu00006 + Otu00007 + Otu00020 + Otu00039 + Otu00011 + Otu00015 + Otu00013 + Otu00023 + Otu00027 + Otu00078, data=topdose)


predictm1 <- as.data.frame(predict.lm(m1, newdata=newtit))
predictm2 <- as.data.frame(predict.lm(m2, newdata=newtit))
predictm3 <- as.data.frame(predict.lm(m3, newdata=newtit))
predictm4 <- as.data.frame(predict.lm(m4, newdata=newtit))
predictm5 <- as.data.frame(predict.lm(m5, newdata=newtit))
predictm6 <- as.data.frame(predict.lm(m6, newdata=newtit))
predictm7 <- as.data.frame(predict.lm(m7, newdata=newtit))
predictm8 <- as.data.frame(predict.lm(m8, newdata=newtit))
predictm9 <- as.data.frame(predict.lm(m9, newdata=newtit))
predictm10 <- as.data.frame(predict.lm(m10, newdata=newtit))
predictm11 <- as.data.frame(predict.lm(m11, newdata=newtit))
predictm12 <- as.data.frame(predict.lm(m12, newdata=newtit))
predictm13 <- as.data.frame(predict.lm(m13, newdata=newtit))
predictm14 <- as.data.frame(predict.lm(m14, newdata=newtit))
predictm15 <- as.data.frame(predict.lm(m15, newdata=newtit))
predictm16 <- as.data.frame(predict.lm(m16, newdata=newtit))
predictm17 <- as.data.frame(predict.lm(m17, newdata=newtit))
predictm18 <- as.data.frame(predict.lm(m18, newdata=newtit))




#rsq calc for both models
ybar = colMeans(actual)[1]
SStot = sum((actual-ybar)^2)

SSres1 = sum((actual-predictm1)^2)
rsq1 = 1-(SSres1/SStot)
rsq1

SSres2 = sum((actual-predictm2)^2)
rsq2 = 1-(SSres2/SStot)
rsq2

SSres3 = sum((actual-predictm3)^2)
rsq3 = 1-(SSres3/SStot)
rsq3

SSres4 = sum((actual-predictm4)^2)
rsq4 = 1-(SSres4/SStot)
rsq4

SSres5 = sum((actual-predictm5)^2)
rsq5 = 1-(SSres5/SStot)
rsq5

SSres6 = sum((actual-predictm6)^2)
rsq6 = 1-(SSres6/SStot)
rsq6

SSres7 = sum((actual-predictm7)^2)
rsq7 = 1-(SSres7/SStot)
rsq7

SSres8 = sum((actual-predictm8)^2)
rsq8 = 1-(SSres8/SStot)
rsq8

SSres9 = sum((actual-predictm9)^2)
rsq9 = 1-(SSres9/SStot)
rsq9

SSres10 = sum((actual-predictm10)^2)
rsq10 = 1-(SSres10/SStot)
rsq10

SSres11 = sum((actual-predictm11)^2)
rsq11 = 1-(SSres11/SStot)
rsq11

SSres12 = sum((actual-predictm12)^2)
rsq12 = 1-(SSres12/SStot)
rsq12

SSres13 = sum((actual-predictm13)^2)
rsq13 = 1-(SSres13/SStot)
rsq13

SSres14 = sum((actual-predictm14)^2)
rsq14 = 1-(SSres14/SStot)
rsq14

SSres15 = sum((actual-predictm15)^2)
rsq15 = 1-(SSres15/SStot)
rsq15

SSres16 = sum((actual-predictm16)^2)
rsq16 = 1-(SSres16/SStot)
rsq16

SSres17 = sum((actual-predictm17)^2)
rsq17 = 1-(SSres17/SStot)
rsq17

SSres18 = sum((actual-predictm18)^2)
rsq18 = 1-(SSres18/SStot)
rsq18



####testing rf regression using log(CFU+1) transformation
######################
####testing rf regression using log(CFU+1) transformation

topdose_logtrans<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.topdose2.regression.logtrans.filter16mintot.csv", header=T)
training_logtrans<-topdose_logtrans[,-1]
newtit_logtrans<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.newtitration.regression.logtrans.filter16mintot.csv", header=T)
testing_logtrans<-newtit_logtrans[,-1] #new titration data only
delay_logtrans<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.rf.delays.regression.logtrans.csv", header=T)
testing2_logtrans<-delay_logtrans[,-1]

library(randomForest)

#fit the randomforest model
model_logtrans <- randomForest(nextDayCFU~., 
                      data = training_logtrans,  outscale=TRUE,
                      importance=TRUE, proximity=TRUE,
                      keep.forest=TRUE, ntree=5000
)
plot(model_logtrans)

print(model_logtrans)


#what are the important variables (via permutation) #type 1 is mean decrease in accuracy, type 2 is mean decrease in node impurity
varImpPlot(model_logtrans, type=1)
imp<-importance(model_logtrans)
most.imp<-rownames(imp[order(imp[, "MeanDecreaseAccuracy"], decreasing=T)[1:5],])

#predict the outcome of the testing data
predicted_logtrans <- predict(model_logtrans, newdata=testing_logtrans[ ,-1])
predicted2_logtrans <- predict(model_logtrans, newdata=testing2_logtrans[ ,-1])

# what is the proportion variation explained in the outcome of the testing data?
# i.e., what is 1-(SSerror/SStotal)
actual_logtrans <- testing_logtrans$nextDayCFU
actual2_logtrans <- testing2_logtrans$nextDayCFU

ybar = mean(actual_logtrans)
SSres = sum((actual_logtrans-predicted_logtrans)^2)
SStot = sum((actual_logtrans-ybar)^2)
rsq = 1-(SSres/SStot)

#rsq_logtrans <- 1-sum((actual_logtrans-predicted_logtrans)^2)/sum((actual_logtrans-mean(actual_logtrans))^2)
print(rsq)

ybar2 = mean(actual2_logtrans)
SSres2 = sum((actual2_logtrans-predicted2_logtrans)^2)
SStot2 = sum((actual2_logtrans-ybar2)^2)
rsq2 = 1-(SSres2/SStot2)

print(rsq2)

write.table(cbind(actual_logtrans, predicted_logtrans), file="~/Desktop/mothur/abxD01/rf/RF.regression.topTitr.txt", sep="\t", row.names=FALSE)


par(mar=c(5, 4, 4, 2) +0.1, mgp=c(3, 1, 0)) #default is 5.1 4.1 4.1 2.1 [bottom, left, top, right space], mgp=c(3, 1, 0) [label line location for x/y location labels, tick mark labels location, tick mark locations]
plot(actual_logtrans,  pch=16, col="black",ylab="log(C. difficile CFU + 1)", ylim=c(0, 8), cex=2, xlab="Samples in New Data Set")
points(predicted_logtrans, cex=2, pch=16, col="magenta")
leg<-c("Actual", "Predicted", "Rsquared=+0.789")
col<-c("black", "magenta", "white")
#pch<-c(16, 16, 16)
legend("topleft", legend=leg, pch=16, col=col,  cex=1, bty="n")

for(i in 1:(length(predicted_logtrans))){
  lines(c(i, i), c(predicted_logtrans[i], actual_logtrans[i]))
}
    

#ticksy<-seq(0, 9, by=1)
#labelsy<-sapply(ticksy, function(i) as.expression(bquote(10^ .(i))))
#axis(2, las=1, at=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), labels=labelsy)


##  negative binomial model
#######################################
##  negative binomial model
require(foreign)
require(ggplot2)
require(MASS)

topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.test1.csv", header=T)
td<-topdose[,-1]
model<-glm.nb(nextDayCFU ~ Otu00027 + Otu00039 + Otu00023 + Otu00006 + Otu00020 + Otu00013 + Otu00015, data=td)
model<-glm.nb(td$nextDayCFU~., data= td)

test1<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.shared.newtitration.test1.csv", header=T)
test1<-test1[,-1]

predicted<-predict(model, newdata=test1[,-1], type="response")
actual <- as.data.frame(test1$nextDayCFU)

predict2<-as.data.frame(predict(model, newdata=td[,-1], type="response"))
actual2<- td$nextDayCFU

ybar = colMeans(actual)[1]
SSres = sum((actual-predicted)^2)
SStot = sum((actual-ybar)^2)
rsq = 1-(SSres/SStot)

print(rsq)


ybar = mean(actual2)
SSres = sum((actual2-predict2)^2)
SStot = sum((actual2-ybar)^2)
rsq = 1-(SSres/SStot)

print(rsq)


##zero inflated negative binomial model
#######################################
##zero inflated negative binomial model

topdose<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.topdose.regression.csv", header=T)
training<-topdose[,-1]

newtit<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.newtitration.regression.csv", header=T)
testing2<-newtit[,-1] #new titration data only



test<-read.csv("~/Desktop/mothur/abxD01/rf/test_topdose.csv", header=T)
trainingtop12<-test[,-1]

newtittest2<-read.csv("~/Desktop/mothur/abxD01/rf/test_newtitration.csv", header=T)
testingtop12<-newtittest2[,-1] #new titration data only




library(pscl)
test_topdose.zinb <- zeroinfl(nextDayCFU ~ ., data = trainingtop12, dist = "negbin", EM = TRUE)
summary(test_topdose.zinb)
AIC(test_topdose.zinb) #to compare models
#1600.852 - this model has 12 OTUs chosen based on RF regression of topdose data

##AIC for just otu3 model is 1653.16

predicted<-predict(test_topdose.zinb, newdata=testingtop12[,-1])
actual <- testing2$nextDayCFU

rsq <- 1-sum((actual-predicted)^2)/sum((actual-mean(actual))^2)##this might be wrong
print(rsq)


#> test<-read.csv("~/Desktop/mothur/abxD01/rf/test_topdose.csv", header=T)
#> trainingtop12<-test[,-1]
#> topdose.zinb <- zeroinfl(nextDayCFU ~ ., data = trainingtop12, dist = "negbin", EM = TRUE)
#Warning messages:
#  1: glm.fit: fitted rates numerically 0 occurred 
#2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
#> summary(topdose.zinb)

#Call:
#  zeroinfl(formula = nextDayCFU ~ ., data = trainingtop12, dist = "negbin", EM = TRUE)

#Pearson residuals:
#  Min         1Q     Median         3Q        Max 
#-1.475e+00 -1.022e-01 -1.092e-02 -9.809e-06  4.201e+00 

#Count model coefficients (negbin with log link):
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) 19.2543592  1.0244923  18.794  < 2e-16 ***
#  Otu00003    -0.0638724  0.0484852  -1.317  0.18772    
#Otu00006    -0.0175728  0.0066645  -2.637  0.00837 ** 
#  Otu00015    -0.0126167  0.0236002  -0.535  0.59293    
#Otu00005    -0.0004825  0.0007037  -0.686  0.49290    
#Otu00007     0.0128220  0.0051479   2.491  0.01275 *  
#  Otu00012     0.0078286  0.0180617   0.433  0.66470    
#Otu00011    -0.0010224  0.0006965  -1.468  0.14212    
#Otu00020    -0.0446065  0.0327500  -1.362  0.17319    
#Otu00019    -0.0086725  0.0132595  -0.654  0.51307    
#Otu00092    -0.4197526  0.2317090  -1.812  0.07006 .  
#Otu00010    -0.0009715  0.0010073  -0.964  0.33480    
#Otu00027    -4.7702157  0.4955347  -9.626  < 2e-16 ***
#  Log(theta)   1.1095995  0.4715603   2.353  0.01862 *  
#  
#  Zero-inflation model coefficients (binomial with logit link):
#  Estimate Std. Error z value Pr(>|z|)  
#(Intercept) -7.211e+00  4.762e+00  -1.514   0.1300  
#Otu00003     6.118e-04  3.862e-02   0.016   0.9874  
#Otu00006     1.820e-02  1.421e-02   1.281   0.2002  
#Otu00015     1.428e-01  8.357e-02   1.709   0.0874 .
#Otu00005     3.755e-03  4.847e-03   0.775   0.4385  
#Otu00007     1.159e-02  8.113e-03   1.429   0.1531  
#Otu00012     1.685e-02  3.704e-02   0.455   0.6492  
#Otu00011     2.638e-03  3.112e-03   0.848   0.3966  
#Otu00020     5.461e-02  3.520e-02   1.551   0.1208  
#Otu00019    -9.982e-05  3.081e-02  -0.003   0.9974  
#Otu00092     1.032e+00  6.859e-01   1.504   0.1325  
#Otu00010     1.419e-03  5.890e-03   0.241   0.8097  
#Otu00027    -1.059e+00  9.779e-01  -1.083   0.2790  
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Theta = 3.0331 
#Number of iterations in BFGS optimization: 1 
#Log-likelihood: -773.4 on 27 Df
#> AIC(topdose.zinb)
#[1] 1600.852




##zero inflated negative binomial model, currently editing/selecting OTUs to include
#######################################
##zero inflated negative binomial model, currently editing/selecting OTUs to include

library(pscl)

topdose<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.topdose.regression.csv", header=T)
training<-topdose[,-1]

#newtit<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.newtitration.regression.csv", header=T)
#testing2<-newtit[,-1] #new titration data only

topdose.zinb <- zeroinfl(nextDayCFU ~ ., data = training, dist = "negbin", EM = TRUE)
summary(topdose.zinb)
AIC(topdose.zinb)

predicted<-predict(topdose.zinb, newdata=testing2[,-1])
actual <- testing2$nextDayCFU

rsq <- 1-sum((actual-predicted)^2)/sum((actual-mean(actual))^2)##this might be wrong
print(rsq)





##zero inflated poisson regression model
#######################################
##zero inflated poisson regression model

topdose<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.topdose.regression.csv", header=T)
training<-topdose[,-1]

newtit<-read.csv("~/Desktop/mothur/abxD01/rf/abxD01.final.an.unique_list.0.03.subsample.0.03.pick.shared.zinb.newtitration.regression.csv", header=T)
testing2<-newtit[,-1] #new titration data only



test<-read.csv("~/Desktop/mothur/abxD01/rf/test_topdose.csv", header=T)
trainingtop12<-test[,-1]

newtittest2<-read.csv("~/Desktop/mothur/abxD01/rf/test_newtitration.csv", header=T)
testingtop12<-newtittest2[,-1] #new titration data only




library(pscl)
test_topdose.zinb <- zeroinfl(nextDayCFU ~ ., data = trainingtop12)
summary(test_topdose.zinb)
AIC(test_topdose.zinb) #to compare models



##Linear model, as normal CFU and log(CFU+1)
#######################################
##Linear model, as normal CFU and log(CFU+1)
library(pscl)

test<-read.csv("~/Desktop/mothur/abxD01/rf/test_topdose_logtrans.csv", header=T)
trainingtop12<-test[,-1]

linear<-lm(nextDayCFU ~ ., data = trainingtop12)
summary(linear)
AIC(linear)


#for normal CFU
#> linear<-lm(nextDayCFU ~ ., data = trainingtop12)
#> summary(linear)

#Call:
#  lm(formula = nextDayCFU ~ ., data = trainingtop12)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-55469893 -10195252  -1387191   7193409 156131990 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 78885977   11160076   7.069 4.66e-11 ***
#  Otu00003      -97451     110286  -0.884 0.378233    
#Otu00006     -175970      34632  -5.081 1.04e-06 ***
#  Otu00015     -615762     206779  -2.978 0.003357 ** 
#  Otu00005      -28187      14790  -1.906 0.058476 .  
#Otu00007      -87149      24134  -3.611 0.000408 ***
#  Otu00012       26859     154246   0.174 0.861983    
#Otu00011      -19689       9643  -2.042 0.042816 *  
#  Otu00020      -98397      39637  -2.482 0.014085 *  
#  Otu00019       48691      74453   0.654 0.514068    
#Otu00092    -1502490     764649  -1.965 0.051164 .  
#Otu00010        5087      15222   0.334 0.738684    
#Otu00027      -62804     160579  -0.391 0.696240    
#---
 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 24700000 on 159 degrees of freedom
#Multiple R-squared:  0.5561,  Adjusted R-squared:  0.5227 
#F-statistic:  16.6 on 12 and 159 DF,  p-value: < 2.2e-16

#> AIC(linear)
#[1] 6358.245


#For log(CFU+1)

#Call:
#  lm(formula = nextDayCFU ~ ., data = trainingtop12)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-7.2392 -0.7050  0.0337  0.6847  5.2012 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.887e+00  6.391e-01   9.212  < 2e-16 ***
#  Otu00003    -1.556e-02  6.316e-03  -2.463  0.01483 *  
#  Otu00006    -1.180e-02  1.983e-03  -5.948 1.67e-08 ***
#  Otu00015    -5.482e-02  1.184e-02  -4.629 7.57e-06 ***
#  Otu00005    -6.431e-05  8.469e-04  -0.076  0.93957    
#Otu00007    -7.627e-03  1.382e-03  -5.519 1.36e-07 ***
#  Otu00012     2.299e-02  8.833e-03   2.603  0.01012 *  
#  Otu00011     1.018e-03  5.522e-04   1.843  0.06718 .  
#Otu00020    -9.931e-03  2.270e-03  -4.375 2.19e-05 ***
#  Otu00019     3.277e-03  4.264e-03   0.769  0.44332    
#Otu00092    -1.270e-01  4.379e-02  -2.901  0.00424 ** 
#  Otu00010     1.938e-03  8.717e-04   2.223  0.02763 *  
#  Otu00027    -4.034e-03  9.196e-03  -0.439  0.66149    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.414 on 159 degrees of freedom
#Multiple R-squared:  0.8234,  Adjusted R-squared:  0.8101 
#F-statistic:  61.8 on 12 and 159 DF,  p-value: < 2.2e-16

#> AIC(linear)
#[1] 621.8525





