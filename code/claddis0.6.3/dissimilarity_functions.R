###################Hopkins and St. John 2018#####################
##function for applying proposed family of metrics. Input: character matrix in Claddis object format; table showing hierarchical relationships between characters; alpha value.
##Input required:
##matr = character matrix in Claddis object format (see "ReadMorphNexus")
##type = table with three columns with headings "Char", "Type", and "Sub"
##      "Char" is the character number
##      "Type" is the type of character, either "S" for secondary or "T" for tertiary.  In our experience, quaternary characters are rarely to never used
##      "Sub" is the character numeber that the secondary or tertiary character is contingent on
##alpha = desired alpha value.  Should vary between 0 and 1. 
##      If alpha = 0, only primary characters will contribute to dissimilarity estimate
##      If alpha = 1, then shared primary characters will be weighted by the secondary characters shared; same for secondary characters with tertiary characters
##      See text for further explanation.

##note that delta{ijk} = 1 when characters are comparable, consistent with Gower's initial formulation. It may be preferable to replace the total number of 
#comparable characters with the maximum character total possible (Lloyd's MORB metric) especially if there are ordered characters. This is not implemented here
#but could be done by modifying lines 90 and 212 for the new metric and Gower's, respectively 

## Updated June 2022 to score unlimited nesting of character dependencies.  
##      Dependencies can be entered as 'P', 'S', and 'T' if only primary, secondary and tertiary dependencies, 
##          or as 1 for primary, 2 for secondary, ... for more levels.
##      Also updated weightings to handle more levels.
#################


alpha.coefficient<-function(matr,type,alpha){
  #Update for new versions of claddis (extra level, matrix_1 in object & lower case attribute names):
  matr <- m$matrix_1
  pairs<-combn(nrow(matr$matrix),2)
  HSJ<-matrix(NA,nrow=nrow(matr$matrix),ncol=nrow(matr$matrix))
  alpha<-alpha
  
  #Update to use numbers for levels, instead of P, S, and T:
  type$Type[which(type$Type =='P')] <-1
  type$Type[which(type$Type =='S')] <-2
  type$Type[which(type$Type =='T')] <-3
  type$Type <- as.numeric(type$Type)
  max_levels <- max(type$Type)
  
  
  for (i in 1:ncol(pairs)){
    sim.temp<-matrix(abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],]))))
    #correct values for ordered multistate characters
    if (any(matr$ordering=='ord')){
      m.o<-which(matr$ordering=='ord')
      sim.temp[m.o]<-abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],][m.o]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],][m.o])))#/matr$maximum_values[m.o]
    }
    #correct values for unordered multistate characters
    if (any(matr$ordering=='unord'&matr$maximum_values>1)){
      m.u<-which(matr$ordering=='unord'&matr$maximum_values>1)
      sim.temp[m.u]<-replace(sim.temp[m.u],which(sim.temp[m.u]>1),1)
    }
    #correct polymorphic characters
    if (length(grep("&", unique(c(matr$matrix[pairs[1,i],], matr$matrix[pairs[2,i],])))) > 0){
      polym.states <- sort(c(grep("&", matr$matrix[pairs[1,i],]), grep("&", matr$matrix[pairs[2,i],])))
      for (v in 1:length(polym.states)) {
        char.state1<-strsplit(matr$matrix[pairs[1,i],][polym.states[v]],"&")[[1]]
        char.state2<-strsplit(matr$matrix[pairs[2,i],][polym.states[v]],"&")[[1]]
        int.value <- intersect(char.state1, char.state2)
        if (length(int.value)>0){
          sim.temp[polym.states[v]]<-0
        }
        if (length(int.value)==0 & anyNA(c(char.state1,char.state2))==FALSE){
          if (matr$ordering[polym.states[v]]=='unord'){
            sim.temp[polym.states[v]]<-1
          }
          if (matr$ordering[polym.states[v]]=='ord'){
            pairs.poly<-matrix(0,nrow=length(char.state1),ncol=length(char.state2))
            for (m in 1:length(char.state1)){
              for (n in 1:length(char.state2)){
                pairs.poly[m,n]<-abs(as.numeric(char.state1[m])-as.numeric(char.state2[n]))
              }
            }
            sim.temp[polym.states[v]]<-min(pairs.poly)
          }
        }
      }
    }


    #Weight all characters:
    sim.temp <- sim.temp*matr$character_weights
    
    #Start with the max_level, compute it's contribution to the score, and then work down until 
    #   level 1 (primary) characters:
    for (m in max_levels:1){
      if (any(type$Type==m)){
        te<-which(type$Type==m)
        for (j in 1:length(unique(type$Sub[te]))){
          te.sub<-which(type$Sub==as.character(unique(type$Sub[te])[j]))
          te.sim<-sim.temp[te.sub]
          s.te<-as.numeric(as.character(unique(na.omit(type$Sub[te])))[j])
          if (length(na.omit(te.sim))>0) {
            #Need to add in the weighted count for the denominator:
            wts <- matr$character_weights[te.sub]
            wts[which(is.na(te.sim))] <- NA
            sim.temp[s.te]<-1-(alpha*(1-(sum(na.omit(te.sim))/sum(na.omit(wts))))+(1-alpha))  #alpha is applied as if it were a similarity measure, therefore the distance is converted to a similarity by substracted from one, then the whole thing is substracted from one to convert back to dissimilarity
            sim.temp[te.sub]<-NA
          }
        }
      }
      
    }
    
    #calculate total dissimilarity
    wt.comp.char<-sum(na.omit(cbind(sim.temp,matr$character_weights))[,2])
    HSJ[pairs[1,i],pairs[2,i]]<-HSJ[pairs[2,i],pairs[1,i]]<-sum(na.omit(sim.temp*matr$character_weights))/wt.comp.char
  }
  diag(HSJ)<-0
  return(HSJ)
}





###############
##Implementation of Will's GED
##matr = character matrix in Claddis object format (see "ReadMorphNexus")

GED<-function(matr){
  pairs<-combn(nrow(matr$matrix),2)
  
  GED<-matrix(NA,nrow=nrow(matr$matrix),ncol=nrow(matr$matrix))
  for (i in 1:ncol(pairs)){
    sim.temp<-matrix(abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],]))))
    #correct values for ordered multistate characters
    if (any(matr$ordering=='ord')){
      m.o<-which(matr$ordering=='ord')
      sim.temp[m.o]<-abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],][m.o]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],][m.o])))
    }
    #correct values for unordered multistate characters
    if (any(matr$ordering=='unord'&matr$maximum_values>1)){
      m.u<-which(matr$ordering=='unord'&matr$maximum_values>1)
      sim.temp[m.u]<-replace(sim.temp[m.u],which(sim.temp[m.u]>1),1)
    }
    #correct polymorphic characters
    if (length(grep("&", unique(c(matr$matrix[pairs[1,i],], matr$matrix[pairs[2,i],])))) > 0){
      polym.states <- sort(c(grep("&", matr$matrix[pairs[1,i],]), grep("&", matr$matrix[pairs[2,i],])))
      for (v in 1:length(polym.states)) {
        char.state1<-strsplit(matr$matrix[pairs[1,i],][polym.states[v]],"&")[[1]]
        char.state2<-strsplit(matr$matrix[pairs[2,i],][polym.states[v]],"&")[[1]]
        int.value <- intersect(char.state1, char.state2)
        if (length(int.value)>0){
          sim.temp[polym.states[v]]<-0
        }
        if (length(int.value)==0 & anyNA(c(char.state1,char.state2))==FALSE){
          if (matr$ordering[polym.states[v]]=='unord'){
            sim.temp[polym.states[v]]<-1
          }
          if (matr$ordering[polym.states[v]]=='ord'){
            pairs.poly<-matrix(0,nrow=length(char.state1),ncol=length(char.state2))
            for (m in 1:length(char.state1)){
              for (n in 1:length(char.state2)){
                pairs.poly[m,n]<-abs(as.numeric(char.state1[m])-as.numeric(char.state2[n]))
              }
            }
            sim.temp[polym.states[v]]<-min(pairs.poly)
          }
        }
      }
    }
    
    #replace missing values
    Sij<-mean(na.omit(sim.temp))
    sim.temp[which(is.na(sim.temp)==TRUE)]<-Sij
    
    #calculate total dissimilarity
    #wt.comp.char<-sum(na.omit(cbind(sim.temp,matr$character_weights))[,2])
    GED[pairs[1,i],pairs[2,i]]<-GED[pairs[2,i],pairs[1,i]]<-sqrt(sum((sim.temp*matr$character_weights)^2))
  }
  diag(GED)<-0
  return(GED)
}



###############
##Implementation of Gower's coefficient of dissimilarity
##matr = character matrix in Claddis object format (see "ReadMorphNexus")
#Note: this script and functions in Claddis do not standardize ordered multistate by range (as if they were quantitative; this follows Wills 1998)
#If weights are according to P&K (part of the Claddis object format), this returns P&K weighting scheme

gower<-function(matr){
  pairs<-combn(nrow(matr$matrix),2)
  gower<-matrix(NA,nrow=nrow(matr$matrix),ncol=nrow(matr$matrix))
  for (i in 1:ncol(pairs)){
    sim.temp<-matrix(abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],]))))
    #correct values for ordered multistate characters
    if (any(matr$ordering=='ord')){
      m.o<-which(matr$ordering=='ord')
      sim.temp[m.o]<-abs(suppressWarnings(as.numeric(matr$matrix[pairs[1,i],][m.o]))-
                           suppressWarnings(as.numeric(matr$matrix[pairs[2,i],][m.o])))
    }
    #correct values for unordered multistate characters
    if (any(matr$ordering=='unord'&matr$maximum_values>1)){
      m.u<-which(matr$ordering=='unord'&matr$maximum_values>1)
      sim.temp[m.u]<-replace(sim.temp[m.u],which(sim.temp[m.u]>1),1)
    }
    #correct polymorphic characters
    if (length(grep("&", unique(c(matr$matrix[pairs[1,i],], matr$matrix[pairs[2,i],])))) > 0){
      polym.states <- sort(c(grep("&", matr$matrix[pairs[1,i],]), grep("&", matr$matrix[pairs[2,i],])))
      for (v in 1:length(polym.states)) {
        char.state1<-strsplit(matr$matrix[pairs[1,i],][polym.states[v]],"&")[[1]]
        char.state2<-strsplit(matr$matrix[pairs[2,i],][polym.states[v]],"&")[[1]]
        int.value <- intersect(char.state1, char.state2)
        if (length(int.value)>0){
          sim.temp[polym.states[v]]<-0
        }
        if (length(int.value)==0 & anyNA(c(char.state1,char.state2))==FALSE){
          if (matr$ordering[polym.states[v]]=='unord'){
            sim.temp[polym.states[v]]<-1
          }
          if (matr$ordering[polym.states[v]]=='ord'){
            pairs.poly<-matrix(0,nrow=length(char.state1),ncol=length(char.state2))
            for (m in 1:length(char.state1)){
              for (n in 1:length(char.state2)){
                pairs.poly[m,n]<-abs(as.numeric(char.state1[m])-as.numeric(char.state2[n]))
              }
            }
            sim.temp[polym.states[v]]<-min(pairs.poly)
          }
        }
      }
    }
    #calculate total dissimilarity
    wt.comp.char<-sum(na.omit(cbind(sim.temp,matr$character_weights))[,2])
    gower[pairs[1,i],pairs[2,i]]<-gower[pairs[2,i],pairs[1,i]]<-sum(na.omit(sim.temp*matr$character_weights))/wt.comp.char
  }
  diag(gower)<-0
  return(gower)
}