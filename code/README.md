# HSJ:  Code #

This directory contains the following R scripts:

+ `dissimilarity_functions.R`:  functions from [Hopkins and St. John, 2018](https://doi.org/10.1098/rspb.2018.1784) that compute the dissimilarity between two taxa, where the contribution of secondary characters (if present) can be scaled by a parameter alpha.  These functions have been modified from those published in [Hopkins and St. John, 2018](https://doi.org/10.1098/rspb.2018.1784) to accommodate changes in the hierarchical structure of lists expected by the Claddis R package (informally a "Claddis object") [Lloyd 2016](https://doi.org/10.1111/bij.12746).

+ `hsjScorer.R`: builds on the extendible framework of [Brazeau et al, 2019.](https://doi.org/10.1093/sysbio/syy083) and contains methods for scoring trees (`hsjTS()`) as well wrappers to use the scorer in the framework's ratchet search function.

## Dependencies: ##

The code relies on the following packages:

```
library(ape)
library(phangorn)
library(TreeSearch)
library(TreeTools)
library(Claddis)
```

## Sample Code: ##

Assuming the above packages have been loaded, we can source the code files:
```
source('dissimilarity_functions.R')
source('hsjScorer.R')
```
To run the tree scorer, you need a tree, characters for that tree, and how the character are related.   Using the character matrix [matrix2-4.nex](../examples/matrix2-4.nex) and types [type2-4.txt](../examples/type2-4.txt) from the [examples](../examples), we can set up the data structures.  First, we use the nexus file with the character matrix to create PhyDat and Claddis objects:
```
madPhy <- ReadAsPhyDat('matrix2-4.nex')
madDat <- ReadMorphNexus('matrix2-4.nex')
```
We also need to read in the description of the types (which are primary, which are secondary, and the dependencies between them):
```
madTyp <- read.table('type2-4.txt',header=TRUE)
```
Let's generate a random tree using the phyDat object:
```
firstTree <- firstTree <- RandomTree(madPhy)
```
To score the tree with &alpha; = 0 (no contribution of secondaries to the score):
```
hsjTS(firstTree, madPhy, madDat, madTyp, alpha = 0)
```
To score the tree with &alpha; = 0.5 (50% contribution of secondaries to the score):
```
hsjTS(firstTree, madPhy, madDat, madTyp, alpha = 0.5)
```
To run a ratchet search starting from that tree with &alpha; = 0.5:
```
hsjRatchet(firstTree, madPhy, madDat, madTyp, alpha = 0.5)
```
