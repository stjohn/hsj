# HSJ:  Code #
This directory contains the following R scripts and uses Claddis 0.3.4 (which was the most recent compiled binary when the associated paper was published and of June 2021).  Claddis has gone through multiple revisions that include renaming functions and variables-- for a version of our scripts that work with Claddis 0.6.3 (requires R 3.6+), see [here](claddis0.6.3/README.md).

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

## Data Assumptions: ##

The code assumes that the hierarchical characters are coded using reductive coding sensu Strong & Lipscomb, 1999 and expects as input:
+ tree: a tree in the phylo format (ape)
+ phyObj: a phyDat object with the same taxa as the tree
+ morObj: a Claddis (morph) object on the same taxa as the tree
+ typObj: a table with a row for each character, containing the character number, the character type ('P' for primary, 'S' for secondary), and the dependency (NA for primaries; and the number of the controlling primary for secondaries)
+ alpha: a number between 0 and 1.0.  Default is 0.5.

We assume that the "controlling primaries" are present/absent characters and coded as either "0" (absent) or "1" (present).




## Sample Code: ##

Assuming the above packages have been loaded, we can source the code files:
```
source('dissimilarity_functions.R')
source('hsjScorer.R')
```
To run the tree scorer, you need a tree, characters for that tree, and how the character are related.   Using the character matrix [matrix2-4.nex](../examples/matrix2-4.nex) and types [type2-4.txt](../examples/type2-4.txt) from the [examples](../examples), we can set up the data structures.  First, we use the nexus file with the character matrix to create PhyDat and Claddis objects:
```
madPhy <- ReadAsPhyDat('../examples/matrix2-4.nex')
madDat <- ReadMorphNexus('../examples/matrix2-4.nex')
```
We also need to read in the description of the types (which are primary, which are secondary, and the dependencies between them):
```
madTyp <- read.table('../examples/type2-4.txt', header=TRUE)
```
Let's generate a random rooted tree using the phyDat object:
```
firstTree <- RandomTree(madPhy, root=TRUE)
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
(The ratchet will take some time to run.)
