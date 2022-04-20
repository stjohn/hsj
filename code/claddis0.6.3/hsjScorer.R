# Scorer for phylogenetic trees using HSJ disparity
#   Katherine St. John & Melanie Hopkins, 2020
#   Updated, June 2021 to use Claddis 0.6.3 
#     (function calls and variable names changed in Claddis)

library(ape)
library(phangorn)
library(Claddis)
library(TreeTools)
library(TreeSearch)

#To score a single tree:
#   Expects a phylo tree that agrees with the phyDat object,
#     as well as a morphy (datObj), a table of character types.  
#     alpha can range from 0 (no secondaries included) to 1.0.
hsjTS <- function(tree, phyObj, datObj, charTypes, alpha = 0.5) {
  tree <- RenumberTips(tree, names(phyObj))
  nTaxa <- length(tree$tip.label)
  score <- hsjScorer(tree$edge[,1], tree$edge[,2], phyObj, datObj, charTypes, alpha,
                     inPostorder = FALSE, nTaxa)
  score
}

#Wrapper to run TreeSearch's ratchet with HSJ computing the tree scores:
#   Need to call a very slightly modified hsjBootstrap to allow own EdgeListSearch to be used.
hsjRatchet <- function(startTree, phyObj, datObj, charTypes, alpha = 0.5, verbosity=5)
{
  t <- Ratchet(startTree, 
               phyObj, 
               TreeScorer = function(parent, child, morphyObj, inPostorder = FALSE, nTaxa) 
                    hsjScorer(parent, child, phyObj, datObj, charTypes, alpha, inPostorder, nTaxa),
               Bootstrapper = function (edgeList, morphyObj, InitializeData = PhyDat2Morphy, 
                                        CleanUpData = UnloadMorphy, TreeScorer, 
                                        EdgeSwapper = SPRSwap, ratchIter = 20, searchHits = 10, verbosity = 1L, 
                                        stopAtPeak = FALSE, stopAtPlateau = 0L, forestSize = 1L, ...) 
                              hsjBootstrap(edgeList, morphyObj, InitializeData, CleanUpData, 
                                 TreeScorer = function(parent, child, morphyObj, inPostorder = FALSE, nTaxa) 
                                      hsjScorer(parent, child, phyObj, datObj, charTypes, alpha, inPostorder, nTaxa), 
                                 EdgeSwapper = NNISwap, maxIter=20, maxHits=10, verbosity = 1L, 
                                 stopAtPeak = FALSE, stopAtPlateau = 0L, forestSize = 1L),
               verbosity = verbosity)
  t
}

#Wrapper to run TreeSearch's search with HSJ computing the tree scores:

#Where all the real work gets done:
hsjScorer <- function(parent, child, phyDat,clad, charTypes, alpha=0.5, 
                      inPostorder = FALSE, nTaxa)
{
    #Set up a tree from the parent & child lists passed in, 
    #   and then reorder the tree edges for post-order traversal:
    tree <- RandomTree(phyDat, root = TRUE)
    #Search expects the tips to be in the same order as in the phyDat object:
    tree <- RenumberTips(tree,names(phyDat))
    #Update to use with new version of TreeSearch and TreeTools:
    edgeList <- PostorderEdges(cbind(parent, child))
    tree$edge <- edgeList
    #Set up extended matrix:
    numNodes = length(parent)+1
    nTaxa = length(parent)/2 +1  
    
    #Make a first pass with Fitch's parsimony to get internal labels for the tree:
    clad2 <- initialInternalLabels(tree, clad, nTaxa, numNodes)
  
    #Make sure charTypes has a heading:
    colnames(charTypes)<-c('Char','Type','Sub')
    #Identify the controlling characters & primaries:   
    controlling <-unique(charTypes[which(is.na(charTypes$Sub)=='FALSE'),3])    
    primaries <- which(charTypes$Type=='P')
    controllingPrimaries <- controlling[ controlling %in% primaries ]
    
    #Get the "special" characters (i.e. those that are secondary and tertiary)    
    #ADD IN:  allow type levels to also be specified by numbers, i.e. 1 = P, 2 = S, 3 = T
    # Or inside the loop?
    #secondaries <- which(charTypes$Type=='S')
    #tertiaries <- which(charTypes$Type == 'T')

    if (length(controllingPrimaries) > 0){
        if (length(controllingPrimaries) == length(primaries)) {
            nonC <- NULL     #There's no non-controlling primaries
        }
        else {
            #Set non-controlling to be the primaries - controlling primaries:
            nonC <- primaries[! primaries %in% controllingPrimaries]
        }
    }
    else { #There's no controlling primaries:
        nonC <- primaries
    }
  
    #For each controlling primary, compute the contribution of it:
    #   First check its secondaries if any controlling.
    #       If so, compute tertiary contribution & save best labels 
    #   Next, using secondaries & precomputed tertiary labels, compute
    #       contribution for each primary in one pass and add to the sum.

    overallScore = 0
    for (p in controllingPrimaries)
    {
        #Find the secondaries for p:
        sec <- which(charTypes$Sub == p)
    
        #For each secondary, check if it's controlling:
        controllingSec <- controlling[ controlling %in% sec ]
        
        #If so, compute the best labeling for the tertiary characters:
        for (s in controllingSec) {
            tert <- which( charTypes$Sub %in% sec)
            restricted <- restrict_clad(clad2,c(s,tert))
            restrictedTypes <- renumber_types(charTypes,c(s,tert))
            newClad = hsjHelper(tree,restricted, nTaxa, numNodes, alpha, cPrim = TRUE, 
                      charTypes = restrictedTypes,preset=0)$clad
            #Add newly fixed tertiary labels to clad object:
            clad2$matrix_1$matrix[,tert] <- newClad$matrix_1$matrix[,c(2:(length(tert)+1))]
        }
        
        #Find all the tertiary characters that are controlled, eventually, by p:
        tert <- which(charTypes$Sub %in% controllingSec)
    
        #Restrict the character matrix to only the controlling primary and its secondaries & tertiaries:
        restricted <- restrict_clad(clad2,c(p,sec,tert))
        #Restrict the types matrix to only the controlling primary and dependents, renumbered:
        restrictedTypes <- renumber_types(charTypes,c(p,sec,tert))

        #Compute the score for this restricted matrix and add to the total
        overallScore <- overallScore + hsjHelper(tree,restricted, nTaxa, numNodes, alpha, cPrim = TRUE, 
                                             charTypes = restrictedTypes, preset=length(tert))$score
  }
  
  #Add in the score for primaries that are not controlling (pure Fitch):
  if (length(nonC) > 0)
  {
    restricted <- restrict_clad(clad2,nonC)
    #Compute the score for this restricted matrix and add to the total
    nonCScore <- hsjHelper(tree,restricted, nTaxa, numNodes, alpha, cPrim=FALSE)$score
    overallScore <- overallScore + nonCScore
  }
  #Return overall score: 
  overallScore
}


#Restrict a claddis object to a subset of characters:
restrict_clad <- function(clad,chars)
{
    restricted <- clad
    restricted$matrix_1$matrix <- restricted$matrix_1$matrix[,chars]
    restricted$matrix_1$ordering <- restricted$matrix_1$ordering[chars]
    restricted$matrix_1$character_weights <- restricted$matrix_1$character_weights[chars]
    restricted$matrix_1$minimum_values <- restricted$matrix_1$minimum_values[chars]
    restricted$matrix_1$maximum_values <- restricted$matrix_1$maximum_values[chars]
    restricted$matrix_1$CharChanges <- restricted$matrix_1$CharChanges[chars]
  
    restricted
}


#Restrict & renumber the types matrix to only the controlling primary and dependents:
# For example:
#   1 P NA
#   2 S 1
#   4 S 1
#   7 T 4
# becomes:
#   1 P NA
#   2 S 1
#   3 S 1
#   4 T 3   
renumber_types <- function(charTypes, chars){

    rt <- charTypes[chars,]
    #First adjust the sub to reflect the new index
    rt$Sub[1] <- NA
    for (i in 2:length(rt$Sub))
    {
        rt$Sub[i] <- which(rt$Char == rt$Sub[i] )
    }
    #Then adjust the char to reflect it:
    for (i in 1:length(rt$Sub))
    {
        rt$Char[i] <- i
    }

print(rt)    
    as.data.frame(rt)
}


hsjHelper <- function(tree, clad, nTaxa, numNodes, alpha,cPrim=FALSE,charTypes=NULL,preset=0)
{
  #Set up running counters:
  numChanges <- 0
  score <- 0
  
  #If no controlling primaries, in the case where Fitch scoring works:
  if (cPrim == FALSE){
    #Then count differences using the character changes that were stored when the internal
    # nodes were labeled:
    s <- sum(clad$matrix_1$CharChanges)
    return( list(score=sum(clad$matrix_1$CharChanges),clad=clad))  
  }
  
  #How to use this to prescore for tertiaries and to score in general?  Sometimes have secondaries in 
  #   charge-- what will alpha function do if there's no primary?
  
  #Otherwise, we have a controlling primary (in first position) and secondaries:
  #   update to also handle tertiary characters:
  #   By construction, the types are organized as c(p,sec,tert) so the indices  
  #     for the secondaries range from 2 to length of secondaries, and the 
  #     possibly preset columns are at the end.

  #Ignore the last preset number of columns since they've been set already (default is 0 have been preset):
  secondaries <- c(2:(nrow(charTypes) - preset))
  
  #Make a matrix for bookkeeping for each controlling primary (best score if node 
  # is score 0 or scored 1)
  scores <- matrix(-1, nrow = length(clad$matrix_1$matrix[,1]),ncol = 2)
  #Set up values for tips:
  for (i in (1:nTaxa)) {
    if (is.na(clad$matrix_1$matrix[i,1])) {
      #Set to a large value that won't be picked up:
      scores[i,]<-10
    }
    else if (clad$matrix_1$matrix[i,1] == '0') {
      scores[i,1]<-0
      scores[i,2]<-10
    }
    else{
      scores[i,1]<-10
      scores[i,2]<-0
    }    
  }
  
  colnames(scores) <- c("Best for Label 0","Best for Label 1")
  score <- 0
  
  #Make a pass through the tree and score via Fitch:
  #Since in post-order, can assume edges come in pairs to the same parent:
  for (i in seq(1,nrow(tree$edge),2))
  {
    par = tree$edge[i,1]
    c1 = tree$edge[i,2]
    c2 = tree$edge[i+1,2]
    
    #Check first to see if the controlling primary is missing for either child,
    #  If that's the case, the score is that of the child:
    if ( is.na(clad$matrix_1$matrix[c1,1]) )
    {
      if ( is.na(clad$matrix_1$matrix[c2,1]) )
      {
        #Both c1 and c2 have missing data as their state for the controlling primary.
        # set scores to 0, since nothing can be contributed if both missing:
        scores[par,1] = 0
        scores[par,2] = 0
      }
      else
      {
        #c1 is missing, but c2 has a state, then use the values for absence and presence of c2:
        scores[par,1] = min(scores[c2,1],scores[c2,2]+1) 
        scores[par,2] = min(scores[c2,2]+1,scores[c2,2])
      }
    }
    else
    {
      #c1 is not missing, but check if c2 is:
      if ( is.na(clad$matrix_1$matrix[c2,1]))
      {
        #c2 is missing, but c1 has a state, then use the values for absence and presence of c1:
        scores[par,1] = min(scores[c1,1],scores[c1,2]+1) 
        scores[par,2] = min(scores[c1,1]+1,scores[c1,2])
      }
      
      ####  In the case where neither c1 or c2 is missing:
      else 
      {
        c1M1 <- 0
        c1M2 <- 0
        c2M1 <- 0
        c2M2 <- 0
        #Special cases for leaves:
        if ( c1 <= nTaxa) {
          if ( clad$matrix_1$matrix[c1,1] == "0" ) {
            c1M1 = 0  #leaf has 0 as its primary character
            c1M2 = 10 #something big, so, it won't get chosen
          }
          else {
            c1M1 = 10 #something big, so, it won't get chosen
            c1M2 = 0  #leaf has 1 as its primary character
          }
        }
        else {
          c1M1 = max(0, scores[c1,1])
          c1M2 = max(0, scores[c1,2])
        }
        if ( c2 <= nTaxa) {
          if ( clad$matrix_1$matrix[c2,1] == "0" ) {
            c2M1 = 0  #leaf has 0 as its primary character
            c2M2 = 10 #something big, so, it won't get chosen
          }
          else {
            c2M1 = 10 #something big, so, it won't get chosen
            c2M2 = 0  #leaf has 1 as its primary character
          }
        }
        else {
          c2M1 = max(0, scores[c2,1])
          c2M2 = max(0, scores[c2,2])
        }
        
        #By definition, the parent is not a leaf node, so, can have either value:
        #When it's 0, we can compute it directly (first making sure that we're not 
        #picking up the -1 that are placeholders)
        scores[par,1] = min(c1M1+c2M1, c1M1+c2M2+1, c1M2+c2M1+1, c1M2+c2M2+2)
        
        
        #To compute when the primary character is present, 3 out of 4 cases can be computed 
        # directly, but we need to compute the 4th one (all three are present) using the 
        # alpha coefficient:
        #Special cases for leaves:
        if ( (c1 <= nTaxa) & (clad$matrix_1$matrix[c1,1] == "0") ) {
          #It's a leaf that has primary state is 0, so, no need to compute the alpha distance:
          scores[par,2] = min(c2M1+2,c2M2+1)
        }
        else if ( (c2 <= nTaxa) & (clad$matrix_1$matrix[c2,1] == "0") ) {
          #It's a leaf that has primary state is 0, so, no need to compute the alpha distance:
          scores[par,2] = min(c1M1+2,c1M2+1)
        }
        else {
          #The children both have computed values for the primary being present:
          tmp <- clad
          tmp$matrix_1$matrix[c1, 1] <- "1"
          tmp$matrix_1$matrix[c2, 1] <- "1"
          tmp$matrix_1$matrix[par, 1] <- "1"
          for (s in secondaries) {
            #Look for any overlap in the labels
            if ( (!is.na(tmp$matrix_1$matrix[[par,s]])) & (tmp$matrix_1$matrix[[par, s]] != "")) {
              #The parent label is non-empty:
              paru <- sort(unlist(strsplit(tmp$matrix_1$matrix[par, s], '/')))
              if ( (!is.na(tmp$matrix_1$matrix[[c1,s]])) & (tmp$matrix_1$matrix[[c1, s]] != "")) {
                #First child's label is non-empty:
                c1u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c1, s], '/')))
                overlap <- c1u[(c1u %in% paru)]
                if (length(overlap) == 0) {
                  #Nothing in common, so, set c1 to any of its labels
                  tmp$matrix_1$matrix[[c1,s]] <- c1u[[1]]
                  if ( (!is.na(tmp$matrix_1$matrix[[c2,s]])) & (tmp$matrix_1$matrix[[c2, s]] != "")) {
                    #c2 has a label:
                    c2u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c2, s], '/')))
                    overlap2 <- c2u[(c2u %in% paru)]    
                    if (length(overlap2) == 0){
                      #None have anything in common, so, set to the first of each:
                      tmp$matrix_1$matrix[[c2,s]] <- c2u[[1]]
                      tmp$matrix_1$matrix[[par,s]] <- paru[[1]]
                    }
                    else {
                      #No overlap with c1
                      #And to compute distance, the last two to something in common:
                      tmp$matrix_1$matrix[[c2,s]] <- overlap2[[1]]
                      tmp$matrix_1$matrix[[par,s]] <- overlap2[[1]]
                    }
                  }
                  else {
                    #c2 doesn't have a label, so, no need to set it, just set par to one of its labels:
                    tmp$matrix_1$matrix[[par,s]] <- paru[[1]]
                  }
                }
                else {
                  #There's something in common, so, check against c2:
                  if ( (!is.na(tmp$matrix_1$matrix[[c2,s]])) & (tmp$matrix_1$matrix[[c2, s]] != "")) {
                    #c2 has a label:
                    c2u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c2, s], '/')))
                    overlap2 <- c2u[(c2u %in% overlap)]
                    if (length(overlap2) == 0){
                      #No overlap with c2
                      #And to compute distance, the first two to something in common:
                      tmp$matrix_1$matrix[[c1,s]] <- overlap[[1]]
                      tmp$matrix_1$matrix[[par,s]] <- overlap[[1]]
                      tmp$matrix_1$matrix[[c2,s]] <- c2u[[1]]
                    }
                    else {
                      #All three have something in common
                      tmp$matrix_1$matrix[[c1,s]] <- overlap2[[1]]
                      tmp$matrix_1$matrix[[c2,s]] <- overlap2[[1]]
                      tmp$matrix_1$matrix[[par,s]] <- overlap2[[1]]
                    }
                  }
                  else {
                    #c2 doesn't have a label, so, only set c1 and par:
                    tmp$matrix_1$matrix[[c1,s]] <- overlap[[1]]
                    tmp$matrix_1$matrix[[par,s]] <- overlap[[1]]
                  }
                }
              }            
              else {
                #c1 is empty, check c2:
                if ( (!is.na(tmp$matrix_1$matrix[[c2,s]])) & (tmp$matrix_1$matrix[[c2, s]] != "")) {
                  #c2 has a label:
                  c2u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c2, s], '/')))
                  overlap2 <- c2u[(c2u %in% paru)]
                  if (length(overlap2) == 0){
                    #No overlap with c2:
                    tmp$matrix_1$matrix[[par,s]] <- paru[[1]]
                    tmp$matrix_1$matrix[[c2,s]] <- c2u[[1]]
                  }
                  else {
                    #c2 and par have something in common:
                    tmp$matrix_1$matrix[[c2,s]] <- overlap2[[1]]
                    tmp$matrix_1$matrix[[par,s]] <- overlap2[[1]]
                  }
                }
                else {
                  #c1 and c2 are empty, so, set par to be one of its labels:
                  tmp$matrix_1$matrix[[par,s]] <- paru[[1]]
                }
              }
            }
            else {
              #The parent has an empty label, so, the kids can be set to the first label for each:
              if ( (!is.na(tmp$matrix_1$matrix[[c1,s]])) & (tmp$matrix_1$matrix[[c1, s]] != "")) {
                c1u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c1, s], '/')))
                tmp$matrix_1$matrix[[c1u,s]] <- c1u[[1]]
              }
              if ( (!is.na(tmp$matrix_1$matrix[[c2,s]])) & (tmp$matrix_1$matrix[[c2, s]] != "")) {
                c2u <- sort(unlist(strsplit(tmp$matrix_1$matrix[c2, s], '/')))
                tmp$matrix_1$matrix[[c2u,s]] <- c2u[[1]]                
              }
            }
          }  
          #Restrict matrix to the 3 nodes:
          tmp$matrix_1$matrix <- tmp$matrix_1$matrix[c(c1, c2, par), ]
          
          d = alpha.coefficient(tmp$matrix_1, type = charTypes, alpha = alpha)
          #d = MorphDistMatrix(CladisticMatrix = tmp, Distance = "GC",
          #                    TransformDistances = "none", InapplicableBehaviour = "HSJ",
          #                    CharacterDependencies = charTypes, Alpha = alpha)
          
          #R starts counting at 1, so, add 1 to access the "0" and "1" columns in score
          #The last term adds the difference if l1 and p have different labels and l2 and p are different
          scores[par, 2] = min(c1M1 + c2M1 + 2, c1M1 + c2M2 + 1, c1M2 + c2M1 +
                                 1, c1M2 + c2M2 + d[1, 3] + d[2, 3])
          if ( (c1M1 + c2M2 + 1 <= scores[par,2]) & (c1M1 + c2M2 + 1 < c1M2 + c2M2 + d[1, 3] + d[2, 3])) {
            #### If the lowest score is due to a child being present and other absent, copy 
            #### the present child's secondary character states to the parent, since that 
            #### needs to be fixed to propagate up the tree:
            clad$matrix_1$matrix[par,secondaries] = clad$matrix_1$matrix[c2,secondaries]
          }
          if ( (c1M2 + c2M1 + 1 <= scores[par,2]) & (c1M2 + c2M1 + 1 < c1M2 + c2M2 + d[1, 3] + d[2, 3])) {
            #### If the lowest score is due to a child being present and other absent, copy 
            #### the present child's secondary character states to the parent, since that 
            #### needs to be fixed to propagate up the tree:
            clad$matrix_1$matrix[par,secondaries] = clad$matrix_1$matrix[c1,secondaries]
          }          
        }
      }
      score = min(scores[par,1],scores[par,2])
    }
  }
  #print(paste('score for par',par,'is',scores[par,1],scores[par,2],sep=" "))
  #ADD TO THE RETURN:  the clad with the labels computed for the secondaries
  return( list(score=score,clad=clad) )  
}

#Does a first pass of the tree to get initial labels for the internal nodes:
initialInternalLabels <- function(tree, clad, nTaxa, numNodes)
{
  newClad <- clad
  
  #Num chars:
  k <- length(clad$matrix_1$matrix[1,])
  #Chnages for each character:
  newClad$matrix_1$CharChanges <- c(rep(0,k))
  #Empty row:
  emptyElt <- c(rep("",k))
  addedNames <- paste("t",seq(nTaxa+1,numNodes),sep="")
  extendedNames <- c(row.names(newClad$matrix_1$matrix),addedNames)
  #Look up how to do this without a loop-- should be able to rep the rows...
  for (i in 1:tree$Nnode)
  {
    newClad$matrix_1$matrix <- rbind(newClad$matrix_1$matrix,emptyElt)
  }
  #Add in the names for the new rows:
  row.names(newClad$matrix_1$matrix) <- extendedNames
  #For each node with children, assign a label:
  #Since in post-order, can assume edges come in pairs to the same parent:
  for (i in seq(1,nrow(tree$edge),2))
  {
    par = tree$edge[i,1]
    c1 = tree$edge[i,2]
    c2 = tree$edge[i+1,2]
    #Compute new label and amount of change it contributes:
    for (j in 1:k) {
      #If c1 has empty value (i.e. missing data), set parent to c2's value:  
      if (is.na(newClad$matrix_1$matrix[c1,j])) {
        if (is.na(newClad$matrix_1$matrix[c2,j])) {
          newClad$matrix_1$matrix[par,j] <- NA
        }
        else {
          newClad$matrix_1$matrix[par,j] <-newClad$matrix_1$matrix[c2,j]
        }
      }
      else if (is.na(newClad$matrix_1$matrix[c2,j])) {
        newClad$matrix_1$matrix[par,j] <- newClad$matrix_1$matrix[c1,j]
      }
      else {
        #Both are non-empty, so, can split up entries:
        #   Split on both '/' and '&' (used by PhyDat for polymorphisms)
        c1u <-  sort(unlist(strsplit(newClad$matrix_1$matrix[c1,j],'/|&')))
        c2u <-  sort(unlist(strsplit(newClad$matrix_1$matrix[c2,j],'/|&')))
        overlap <- sort(c1u[c1u %in% c2u])
        if (length(overlap) > 1) {
          newClad$matrix_1$matrix[par,j] <- paste(overlap,collapse='/')
        }
        else if (length(overlap) == 1) {
          newClad$matrix_1$matrix[par,j] <-overlap
        }
        else {
          #No overlap, so use the union of the two:
          newClad$matrix_1$matrix[par,j] <- paste(sort(c(c1u,c2u)),collapse='/')
          newClad$matrix_1$CharChanges[j] <- newClad$matrix_1$CharChanges[j] + 1
        }
      }
    }
  }
  #print(newClad$matrix_1$CharChanges)  
  newClad 
}




#The EdgeListSearch() in MorphyBootstrap() hardcodes to Morphy.  This function is 
# almost identical to MorphyBootstrap() except it allows the EdgeListSearch() to
# be specified.
hsjBootstrap <- function (edgeList, morphyObj, InitializeData = PhyDat2Morphy, 
                          CleanUpData = UnloadMorphy, TreeScorer, 
                          EdgeSwapper = NNISwap, maxIter = 10, maxHits = 20, verbosity = 1L, 
                          stopAtPeak = FALSE, stopAtPlateau = 0L, forestSize = 1L, ...) 
{
  if (verbosity > 1L) {
    print('Inside hsj boot strap')
  }
  startWeights <- MorphyWeights(morphyObj)[1, ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace = TRUE), 
                         length(startWeights))
  errors <- vapply(eachChar, function(i) mpl_set_charac_weight(i, 
                                                               resampling[i], morphyObj), integer(1))
  if (any(errors)) {
    stop("Error resampling morphy object: ", mpl_translate_error(unique(errors[errors < 
                                                                                 0L])))
  }
  if (error <- mpl_apply_tipdata(morphyObj)) {
    stop("Error applying tip data: ", mpl_translate_error(error))
  }
  res <- #EdgeListSearch(edgeList[1:2], morphyObj, EdgeSwapper = EdgeSwapper, 
    #                maxIter = maxIter, maxHits = maxHits, stopAtPeak = stopAtPeak, 
    #                stopAtPlateau = stopAtPlateau, verbosity = verbosity - 
    #                  1L, ...)
    EdgeListSearch(edgeList[1:2], morphyObj, TreeScorer = TreeScorer, 
                   EdgeSwapper = EdgeSwapper, maxIter = maxIter, maxHits = maxHits,
                   forestSize = forestSize, stopAtPeak = stopAtPeak, stopAtPlateau = stopAtPlateau, 
                   verbosity = verbosity, ...)
  if (verbosity > 1) {
    print('ran EdgeListSearch!')
    print(res)
  }
  errors <- vapply(eachChar, function(i) mpl_set_charac_weight(i, 
                                                               startWeights[i], morphyObj), integer(1))
  if (any(errors)) 
    stop("Error resampling morphy object: ", mpl_translate_error(unique(errors[errors < 
                                                                                 0L])))
  if (error <- mpl_apply_tipdata(morphyObj)) 
    stop("Error applying tip data: ", mpl_translate_error(error))
  res[1:2]
}

