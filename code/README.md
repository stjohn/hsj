# HSJ:  Code #

This directory contains the following R scripts:

+ `dissimilarity_functions.R`:  functions from [Hopkins and St. John 2018](https://doi.org/10.1098/rspb.2018.1784) that compute the dissimilarity between two taxa, where the contribution of secondary characters (if present) can be scaled by a parameter alpha.  These functions have been modified from those published in [Hopkins and St. John 2018](https://doi.org/10.1098/rspb.2018.1784) to accommodate changes in the hierarchical structure of lists expected by the Claddis R package (informally a "Claddis object") [Lloyd 2016](https://doi.org/10.1111/bij.12746).

+ `hsjScorer.R`: builds on the extendible framework of [Brazeau et al, 2019.](https://doi.org/10.1093/sysbio/syy083) and contains methods for scoring trees (`hsjTS()`) as well wrappers (`hsjSearch()` and `hsjRatchet()`) to use the scorer in the framework's search and ratchet functions.

## Dependencies: ##

The code relies on the following packages:

'''library(ape)
library(phangorn)
library(TreeSearch)
library(TreeTools)
library(Claddis)
'''

## Sample Code: ##

Create phy, dat, and typ objects and call score, search, ratchet.
