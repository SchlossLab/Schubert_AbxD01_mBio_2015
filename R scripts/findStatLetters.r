

findStatLetters <- function(file, otus, graphGroups){

  statLetter <- NULL
  statLetter <- matrix(c("NA"),nrow=3, ncol=length(otus))
  
  ind <- NULL
  for(a in 1:(length(graphGroups))){
    ind <- c(ind, which(file$expgroup==graphGroups[a]))
  }
  file <- file[ind,]
  
  for (i in 1:(length(otus))) {      
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
  } # for (i in 1:(length(otus))) {      

  return(statLetter)
}