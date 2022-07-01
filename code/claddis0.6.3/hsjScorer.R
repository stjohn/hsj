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
    #edgeList <- PostorderEdges(cbind(parent, child))
    edgeList <- tree$edge[order(tree$edge[,1],decreasing = TRUE),]
    tree$edge <- edgeList
    #Set up extended matrix:
    numNodes = length(parent)+1
    nTaxa = length(parent)/2 +1  
    
    #Make a first pass with Fitch's parsimony to get internal labels for the tree:
    clad2 <- initialInternalLabels(tree, clad, nTaxa, numNodes)
  
    #Make sure charTypes has a heading:
    colnames(charTypes)<-c('Char','Type','Sub')
    #And update to use numbers for levels, instead of P, S, and T:
    charTypes$Type[which(charTypes$Type =='P')] <-1
    charTypes$Type[which(charTypes$Type =='S')] <-2
    charTypes$Type[which(charTypes$Type =='T')] <-3
    charTypes$Type <- as.numeric(charTypes$Type)
    maxLevel <- max(charTypes$Type)
    
    #Make a matrix for bookkeeping for each possible label for the controlling primary:
    num_symbols = length(clad$matrix_1$characters$symbols)
    scoring <- array(-1, dim = c( length(charTypes$Char),numNodes,num_symbols,3) )
    dimnames(scoring)[[3]] <- c(clad$matrix_1$characters$symbols)
    dimnames(scoring)[[4]] <- c('score','contrib1','contrib2')
    #Fill in values for tips:
    scoring[,1:nTaxa,,'score'] <- 10
    for (m in (1:length(charTypes$Char))){
      for (i in (1:nTaxa)) {
        c = clad$matrix_1$matrix[i,m]
        if (c != ""){
          scoring[m,i,clad$matrix_1$matrix[i,m],'score'] <- 0
        }
      }      
    }

        
    #Working from the most nested level (maxLevel), 
    #   compute the possible scores for each label and store it... 
    #   then use that stored value(s) in the next level
    if (maxLevel > 1) {
      for (lev in maxLevel:2){
        if (any(charTypes$Type==lev)){
          #Identify all characters at that level        
          te<-which(charTypes$Type==lev)
          for (j in 1:length(unique(charTypes$Sub[te]))){
            #Find the controlling character
            te.sub<-which(charTypes$Char==as.character(unique(charTypes$Sub[te])[j]))
            #Restrict to the controlling and dependents & compute values to be used in next level:
            restricted <- restrict_clad(clad2,c(j,te.sub))
            restrictedTypes <- renumber_types(charTypes,c(j,te.sub))
            restrictedScores <- scoring[c(j,te.sub),,,]
            #Add newly computed scores for the controlling character to the scores:
            #Want to set the whole row:
            scoring[j,,,] <- hsjHelper(tree,restricted, nTaxa, numNodes, alpha, restrictedScores, cPrim = TRUE, 
                                charTypes = restrictedTypes)
          }
        }
      }
    }
    #With all the contributions of nested characters computed, sum up the scores of the primaries:
    prims<- which(charTypes$Type == 1)
    #Get the root of the tree:
    rootIndex <- tree$edge[numNodes-1,1]
    overallScore <- 0
    for (p in prims){
      #If nested contributions computed, use those:
      if (min(scoring[p,rootIndex,,'score']) > -1){
        overallScore <- overallScore + min(scoring[p,rootIndex,,'score'])
      }
      else {#Else use the initial labels score (Fitch):
        overallScore <- overallScore + sum(clad2$matrix_1$CharChanges[p])
      }
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
renumber_types <- function(charTypes, chars)
{
    rt <- charTypes[chars,]
    #First adjust the sub to reflect the new index
    rt$Sub[1] <- NA
    if (length(rt$Sub) > 1) {
      for (i in 2:length(rt$Sub))
      {
        rt$Sub[i] <- which(rt$Char == rt$Sub[i] )
      }
      #Then adjust the char to reflect it:
      for (i in 1:length(rt$Sub))
      {
        rt$Char[i] <- i
      }
      as.data.frame(rt)
    }
}


hsjHelper <- function(tree, clad, nTaxa, numNodes, alpha, scoring, cPrim=FALSE,charTypes=NULL)
{
  #For each possible value of the label for the parent, compute it's best score.
  #Pre-computed the contributions of each clade, so, only computing additional 
  #from "crossing branches" of the tree:
  compute_parent_label_score <- function(p_label, base)
  {
    #Keep track of the currentMin:
    currentMin <- base[1,1] +2 #placeholder
    contrib1 <- -1 #placeholder
    contrib2 <- -1
    #Check how much the edges add.  The easy cases are when p_label != c1 or c2 label (since it's just +1 for each)
    #Hard case:  if equal, and defined, clean up and send to dissimilarity function to get answer.
    for (c1_lab in clad$matrix_1$characters$symbols) {
      for (c2_lab in clad$matrix_1$characters$symbols) {
        #For each combination of labels
        if ( p_label == c1_lab ) {
          if ( p_label == c2_lab ) {
            #Set up a temporary claddis object to pass to the dissimilarity function:
            tmp = clad
            #All have the same label-- check if all defined and resolve overlaps:
            for (s in secondaries) {
              #Check to see if the secondary value has been pre-computed:
              if ( !is.null(scoring) & scoring[s,par,p_label,"contrib1"] > -1 ) {
                #Set the type of s to ordered, the parent to 0, and the kids to the precomputed values:
                tmp$matrix_1$matrix[[c1,s]] <- scoring[s,par,p_label,"contrib1"]
                tmp$matrix_1$matrix[[c2,s]] <- scoring[s,par,p_label,"contrib2"]
                tmp$matrix_1$matrix[[par,s]] <- 0
                tmp$matrix_1$ordering[[s]] <-"ord"
              }
              #Look for any overlap in the labels
              else if ( (!is.na(tmp$matrix_1$matrix[[par,s]])) & (tmp$matrix_1$matrix[[par, s]] != "")) {
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
                        #No overlap with c1 & to compute distance, the last two to something in common:
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
                      ### SHOULD WE ALSO SET C2?? ###  No, it should stay ""...              
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
            #Set the parent label to p_label to run through the dissimilarity function:
            tmp$matrix_1$matrix[par,1] <- p_label
            tmp$matrix_1$matrix[c1,1] <- c1_lab
            tmp$matrix_1$matrix[c2,1] <- c2_lab
            #Restrict matrix to the 2 nodes:
            save <- tmp$matrix_1$matrix
            tmp$matrix_1$matrix <- tmp$matrix_1$matrix[c(par,c1), ]
            distances <- alpha.coefficient(tmp,charTypes,alpha)
            d1 <- distances[1,2]
            tmp$matrix_1$matrix <- save[c(par,c2), ]
            distances <- alpha.coefficient(tmp,charTypes,alpha)            
            d2 <- distances[1,2]           
            tmp <- base[c1_lab,c2_lab] + d1 + d2
            if (tmp <= currentMin)
            {
                currentMin <- min(currentMin, tmp)
                contrib1 <- d1
                contrib2 <- d2
            }

            #Save the labels that gave the best score:
            #min_labels = clad$matrix_1$matrix[c2,secondaries]
            
            #d = alpha.coefficient(tmp$matrix_1, type = charTypes, alpha = alpha)
            #d = MorphDistMatrix(CladisticMatrix = tmp, Distance = "GC",
            #                    TransformDistances = "none", InapplicableBehaviour = "HSJ",
            #                    CharacterDependencies = charTypes, Alpha = alpha)
            
            #R starts counting at 1, so, add 1 to access the "0" and "1" columns in score
            #The last term adds the difference if l1 and p have different labels and l2 and p are different
            
            #scores[par, 2] = min(c1M1 + c2M1 + 2, c1M1 + c2M2 + 1, c1M2 + c2M1 +
            #                       1, c1M2 + c2M2 + d1 + d2)
            #if ( (c1M1 + c2M2 + 1 <= scores[par,2]) & (c1M1 + c2M2 + 1 < c1M2 + c2M2 + d1 + d2)) {
              #### If the lowest score is due to a child being present and other absent, copy 
              #### the present child's secondary character states to the parent, since that 
              #### needs to be fixed to propagate up the tree:
           #   clad$matrix_1$matrix[par,secondaries] = clad$matrix_1$matrix[c2,secondaries]
           # }
           # if ( (c1M2 + c2M1 + 1 <= scores[par,2]) & (c1M2 + c2M1 + 1 < c1M2 + c2M2 + d1 + d2)) {
              #### If the lowest score is due to a child being present and other absent, copy 
              #### the present child's secondary character states to the parent, since that 
              #### needs to be fixed to propagate up the tree:
           #   clad$matrix_1$matrix[par,secondaries] = clad$matrix_1$matrix[c1,secondaries]
           # }          
          }
          else { #p_label == c1_label but c2_label are different
            tmp <- base[c1_lab,c2_lab] + 1
            if (tmp <= currentMin)
            {
              currentMin <- min(currentMin, tmp)
              contrib1 <- -1
              contrib2 <- -1
            }            
          
          }
        }
        else {
          tmp <- base[c1_lab,c2_lab] + 1
          if (tmp <= currentMin)
          {
            currentMin <- min(currentMin, tmp)
            contrib1 <- -1
            contrib2 <- -1
          }          
          #No need for the else where all p_label != c1_label and p_label != c2_label 
          #since captured by initialization of currentMin.
        }
        
      }
    }
    #Return smallest found & their contributions for nesting:
    return( c(currentMin,contrib1,contrib2) )
  }    
  
  #Set up running counters:
  running_score <- 0  #Best score so far-- 

  #Make a pass through the tree and update scores via Fitch:
  #Since in post-order, can assume edges come in pairs to the same parent:
  for (i in seq(1,nrow(tree$edge),2))
  {
    par = tree$edge[i,1]
    c1 = tree$edge[i,2]
    c2 = tree$edge[i+1,2]
    
    #Set up a matrix of the scores contributed by the children clades.
    #Those contribute the same no matter what label the parent has, so, 
    #can be computed once and added to each parent label:
    num_symbols = length(clad$matrix_1$characters$symbols)
    base <- matrix(-1, nrow = num_symbols,ncol = num_symbols)
    colnames(base) <- c(clad$matrix_1$characters$symbols)
    rownames(base) <- c(clad$matrix_1$characters$symbols)
    for (c1s in 1:num_symbols) {
      for (c2s in 1:num_symbols) {
        base[c1s,c2s] = scoring[1,c1,c1s,"score"] + scoring[1,c2,c2s,"score"]
      }
    }
    
    #Compute the best score for each possible label for the parent:
    secondaries <- c(2:(nrow(charTypes)))
    for (p_label in clad$matrix_1$characters$symbols) {
       s <- compute_parent_label_score(p_label, base)
       scoring[1,par,p_label,"score"] <- s[1]      
       scoring[1,par,p_label,"contrib1"] <- s[2]
       scoring[1,par,p_label,"contrib2"] <- s[3]       
    }
    running_score <- min( scoring[1,par,,])
  }

  #Return the score computed as well as the contribution if needed for nested computations:
  return( scoring[1,,,] )  
}


#Does a first pass of the tree to get initial labels for the internal nodes:
initialInternalLabels <- function(tree, clad, nTaxa, numNodes)
{
  newClad <- clad
  k <- length(clad$matrix_1$matrix[1,])
  #Changes for each character:
  newClad$matrix_1$CharChanges <- c(rep(0,k))
  emptyElt <- c(rep("",k))
  addedNames <- paste("t",seq(nTaxa+1,numNodes),sep="")
  extendedNames <- c(rownames(newClad$matrix_1$matrix),addedNames)
  #Look up how to do this without a loop-- should be able to rep the rows...
  for (i in 1:tree$Nnode)
  {
    newClad$matrix_1$matrix <- rbind(newClad$matrix_1$matrix,emptyElt)
  }
  #Add in the names for the new rows:
  rownames(newClad$matrix_1$matrix) <- extendedNames
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

