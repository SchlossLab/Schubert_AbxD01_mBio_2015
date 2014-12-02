require(foreign)
require(ggplot2)
require(MASS)


## linear model
#######################################
## linear model
## File headings: SampleName nextDayCFU OTU1 OTU2 ... total OTU candidates being tested, candidate number should be at least 10
## nextDayCFU should be log transformed
## IN file type is *.csv




topdose<-read.csv("~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.19otus.rfnegpos.csv", header=T)
td<-topdose[,-1]
ids<-names(td)
ids = ids[-1]

numOTU=length(td)-1

results = data.frame(matrix(0, ncol=numOTU+6)) #table as follows: model# AIC numVar intercept OTU1B OTU2B ... etc all included
colnames(results)[1:6] = c("Modelnum", "AIC", "BIC", "rsq", "vars", "intCoeff")
colnames(results)[7:(numOTU+6)] = ids
#make results headings

#results = NULL
#results = data.frame(NULL)

######################################################
#all combos for a model with 1 variable
for(i in 1:numOTU){ 
  m1 = NULL
  err = NULL
  err<-tryCatch(
{
  m1<-lm(td[,1] ~ td[,(i+1)])
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
  )

if(is.null(err)) {
  results[i, 1] = i
  results[i, 2] = c("NA")
  results[i, 3] = c("NA")
  results[i, 4] = c("NA")
  results[i, 5] = c(1)
  results[i, 6] = names(td)[(i+1)] #add more names as in formula
}
else{
  results[i, 1] = i #model number
  results[i, 2] = AIC(m1) #AIC 
  results[i, 3] = BIC(m1) #BIC 
  coeff=m1$coefficients
  names(coeff)[2] = names(td)[i+1]
  vars = length(coeff)-1
  results[i, 4] = summary(m1)$r.squared
  results[i, 5] = vars #number of variables in model
  results[i, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    for(k in 2:(vars+1)){
      if(names(coeff)[k]==names(results)[j]){
        results[i, j] =   coeff[k]  #print coefficient value
      }
      else{
        if(k>vars){
          results[i, j] = c("NA")
        }
      }  
    }
  } 
}
}

######################################################
#all combos for a model with 2 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU)){ 
  for(b in (i+1):(numOTU+1)){
    m1 = NULL
    err = NULL
    err<-tryCatch(
{
  m1<-lm(td[,1] ~ td[,(i)]+td[,(b)])
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
    )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(2) #change based on vars#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
}
else{
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
        
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
    }
  } 
}
rowNum = rowNum+1
  }
}


######################################################
#all combos for a model with 3 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-1)){ #change subtracted number according to number of variables in the model: vars-1=2
  for(b in (i+1):(numOTU)){ #change subtracted number according to number of variables in the model: vars-2=1
    for(c in (b+1):(numOTU+1)){
      m1 = NULL
      err = NULL
      err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] ) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
      )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(3) #change based on var num
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} 
else{
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else

rowNum = rowNum+1
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}


######################################################
#all combos for a model with 4 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-2)){ #change for num varibales in model
  for(b in (i+1):(numOTU-1)){ #change for num varibales in model
    for(c in (b+1):(numOTU)){ #change for num varibales in model
      for(d in (c+1):(numOTU+1)){ #change for num varibales in model
        m1 = NULL
        err = NULL
        err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] ) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
        )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(4) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} 
else{
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}

#########################################################
#all combos for a model with 5 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-3)){ #change for num varibales in model
  for(b in (i+1):(numOTU-2)){ #change for num varibales in model
    for(c in (b+1):(numOTU-1)){ #change for num varibales in model
      for(d in (c+1):(numOTU)){ #change for num varibales in model
        for(e in (d+1):(numOTU+1)){
          m1 = NULL
          err = NULL
          err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] ) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
          )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(5) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}


######################################################
#all combos for a model with 6 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-4)){ #change for num varibales in model
  for(b in (i+1):(numOTU-3)){ #change for num varibales in model
    for(c in (b+1):(numOTU-2)){ #change for num varibales in model
      for(d in (c+1):(numOTU-1)){ #change for num varibales in model
        for(e in (d+1):(numOTU)){
          for(f in (e+1):(numOTU+1)){
            
            m1 = NULL
            err = NULL
            err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] + td[,f]) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
            )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA") 
  results[rowNum, 4] = c("NA") 
  results[rowNum, 5] = c(6) #changed based on var#
  results[rowNum, 6] = c("NA") 
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[f]==names(results)[j]){
      results[rowNum, j] =   names(td)[f]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  names(coeff)[7] = names(td)[(f)] #add more names as in formula
  
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1

          } #for(f in (e+1):(numOTU+1))
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}


############################################################
#all combos for a model with 7 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-5)){ #change for num varibales in model
  for(b in (i+1):(numOTU-4)){ #change for num varibales in model
    for(c in (b+1):(numOTU-3)){ #change for num varibales in model
      for(d in (c+1):(numOTU-2)){ #change for num varibales in model
        for(e in (d+1):(numOTU-1)){
          for(f in (e+1):(numOTU)){
            for(g in (f+1):(numOTU+1)){ 
              m1 = NULL
              err = NULL
              err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] + td[,f] + td[,g]) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
              )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(7) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[f]==names(results)[j]){
      results[rowNum, j] =   names(td)[f]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[g]==names(results)[j]){
      results[rowNum, j] =   names(td)[g]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  names(coeff)[7] = names(td)[(f)] #add more names as in formula
  names(coeff)[8] = names(td)[(g)] #add more names as in formula
  
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1

            } #for(g in (f+1):(numOTU+1))
          } #for(f in (e+1):(numOTU+1))
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}


###############################################################
#all combos for a model with 8 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-6)){ #change for num varibales in model
  for(b in (i+1):(numOTU-5)){ #change for num varibales in model
    for(c in (b+1):(numOTU-4)){ #change for num varibales in model
      for(d in (c+1):(numOTU-3)){ #change for num varibales in model
        for(e in (d+1):(numOTU-2)){
          for(f in (e+1):(numOTU-1)){
            for(g in (f+1):(numOTU)){ 
              for(h in (g+1):(numOTU+1)){ 
                
                m1 = NULL
                err = NULL
                err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] + td[,f] + td[,g] + td[,h]) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
                )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(8) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[f]==names(results)[j]){
      results[rowNum, j] =   names(td)[f]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[g]==names(results)[j]){
      results[rowNum, j] =   names(td)[g]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[h]==names(results)[j]){
      results[rowNum, j] =   names(td)[h]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  names(coeff)[7] = names(td)[(f)] #add more names as in formula
  names(coeff)[8] = names(td)[(g)] #add more names as in formula
  names(coeff)[9] = names(td)[(h)] #add more names as in formula
  
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1

              } #for(h in (g+1):(numOTU+1))
            } #for(g in (f+1):(numOTU+1))
          } #for(f in (e+1):(numOTU+1))
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}



############################################################
#all combos for a model with 9 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-7)){ #change for num varibales in model
  for(b in (i+1):(numOTU-6)){ #change for num varibales in model
    for(c in (b+1):(numOTU-5)){ #change for num varibales in model
      for(d in (c+1):(numOTU-4)){ #change for num varibales in model
        for(e in (d+1):(numOTU-3)){
          for(f in (e+1):(numOTU-2)){
            for(g in (f+1):(numOTU-1)){ 
              for(h in (g+1):(numOTU)){ 
                for(p in (h+1): (numOTU+1)){
                  
                  m1 = NULL
                  err = NULL
                  err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] + td[,f] + td[,g] + td[,h] + td[,p]) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
                  )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(9) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[f]==names(results)[j]){
      results[rowNum, j] =   names(td)[f]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[g]==names(results)[j]){
      results[rowNum, j] =   names(td)[g]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[h]==names(results)[j]){
      results[rowNum, j] =   names(td)[h]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[p]==names(results)[j]){
      results[rowNum, j] =   names(td)[p]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  names(coeff)[7] = names(td)[(f)] #add more names as in formula
  names(coeff)[8] = names(td)[(g)] #add more names as in formula
  names(coeff)[9] = names(td)[(h)] #add more names as in formula
  names(coeff)[10] = names(td)[(p)] #add more names as in formula
  
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1

                } # for(p in (h+1): (numOTU+1))
              } #for(h in (g+1):(numOTU+1))
            } #for(g in (f+1):(numOTU+1))
          } #for(f in (e+1):(numOTU+1))
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}


#########################################################
#all combos for a model with 10 variables
totmodels = dim(results)[1]
rowNum = totmodels +1
for(i in 2:(numOTU-8)){ #change for num varibales in model
  for(b in (i+1):(numOTU-7)){ #change for num varibales in model
    for(c in (b+1):(numOTU-6)){ #change for num varibales in model
      for(d in (c+1):(numOTU-5)){ #change for num varibales in model
        for(e in (d+1):(numOTU-4)){
          for(f in (e+1):(numOTU-3)){
            for(g in (f+1):(numOTU-2)){ 
              for(h in (g+1):(numOTU-1)){ 
                for(p in (h+1): (numOTU)){
                  for(q in (p+1): (numOTU+1)){
                    
                    m1 = NULL
                    err = NULL
                    err<-tryCatch(
{
  m1<-lm( td[,1] ~ td[,(i)] + td[,(b)] + td[,(c)] + td[,d] + td[,e] + td[,f] + td[,g] + td[,h] + td[,p] + td[,q]) #add more vars if needed
},
error = function(cond) {
  message(cond)
  return(NULL)     
},
finally = {}
                    )

if(is.null(err)) {
  results[rowNum, 1] = rowNum
  results[rowNum, 2] = c("NA")
  results[rowNum, 3] = c("NA")
  results[rowNum, 4] = c("NA")
  results[rowNum, 5] = c(10) #changed based on var#
  results[rowNum, 6] = c("NA")
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    
    if(names(td)[i]==names(results)[j]){
      results[rowNum, j] =   names(td)[i]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[b]==names(results)[j]){
      results[rowNum, j] =   names(td)[b]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[c]==names(results)[j]){
      results[rowNum, j] =   names(td)[c]  #print coefficient value
      found = TRUE
    }   
    else if(names(td)[d]==names(results)[j]){
      results[rowNum, j] =   names(td)[d]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[e]==names(results)[j]){
      results[rowNum, j] =   names(td)[e]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[f]==names(results)[j]){
      results[rowNum, j] =   names(td)[f]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[g]==names(results)[j]){
      results[rowNum, j] =   names(td)[g]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[h]==names(results)[j]){
      results[rowNum, j] =   names(td)[h]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[p]==names(results)[j]){
      results[rowNum, j] =   names(td)[p]  #print coefficient value
      found = TRUE
    }
    else if(names(td)[q]==names(results)[j]){
      results[rowNum, j] =   names(td)[q]  #print coefficient value
      found = TRUE
    }
    else if(found){
    }
    else{
      results[rowNum, j] = c("NA")  
    }  
    
  } #for(j in 5:(numOTU+4))
  
} #if(is.null(err))
else{ #no error in model building
  results[rowNum, 1] = rowNum #model number
  results[rowNum, 2] = AIC(m1) #AIC 
  results[rowNum, 3] = BIC(m1) #BIC 
  
  coeff=m1$coefficients      
  names(coeff)[2] = names(td)[(i)]
  names(coeff)[3] = names(td)[(b)]
  names(coeff)[4] = names(td)[(c)] #add more names as in formula
  names(coeff)[5] = names(td)[(d)] #add more names as in formula
  names(coeff)[6] = names(td)[(e)] #add more names as in formula
  names(coeff)[7] = names(td)[(f)] #add more names as in formula
  names(coeff)[8] = names(td)[(g)] #add more names as in formula
  names(coeff)[9] = names(td)[(h)] #add more names as in formula
  names(coeff)[10] = names(td)[(p)] #add more names as in formula
  names(coeff)[11] = names(td)[(q)] #add more names as in formula
  
  
  vars = length(coeff)-1
  results[rowNum, 4] = summary(m1)$r.squared
  results[rowNum, 5] = vars #number of variables in model
  results[rowNum, 6] = coeff[1] #intercept coefficient
  
  for(j in 7:(numOTU+6)){  #fill results columns for each OTU whether NA or a coefficient value
    found = FALSE
    for(k in 2:(vars+1)){
      
      if(names(coeff)[k]==names(results)[j]){
        results[rowNum, j] =   coeff[k]  #print coefficient value
        found = TRUE
      }
      else if(found){
      }
      else{
        if(k>vars){
          results[rowNum, j] = c("NA")
        }
      }  
      
    } #for(k in 2:(vars+1))
  } #for(j in 5:(numOTU+4))
  
} #else
rowNum = rowNum+1

                  } # for(q in (p+1): (numOTU+1))
                } # for(p in (h+1): (numOTU+1))
              } #for(h in (g+1):(numOTU+1))
            } #for(g in (f+1):(numOTU+1))
          } #for(f in (e+1):(numOTU+1))
        } #for(e in (d+1):(numOTU+1))
      } #for(d in (c+1):(numOTU+1))
    } #for(c in cstart:numOTU)
  } #for(b in bstart:numOTU)
}

write.table(results[1:dim(results)[1],], file="~/Desktop/mothur/abxD01/model/abxD01.final.an.unique_list.0.03.subsample.filter16mintotal.shared.topdose2.logtrans.19otus.rfnegpos.results.txt", sep="\t", row.names=FALSE)
