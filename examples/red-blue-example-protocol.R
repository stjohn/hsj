###script for red/blue tail example matrix

library(TreeSearch)
library(phangorn)
library(Claddis)
library(TreeTools)

source('../code/hsjScorer.R')
source('../code/dissimilarity_functions.R')

##### First example:
#Read in the data files:
mor2.4 <- ReadMorphNexus('matrix2-4.nex')
typ2.4 <- read.table('type2-4.txt')
phy2.4 <- ReadAsPhyDat('matrix2-4.nex')
# Traditional Fitch algorithm with inapplicable characters treated as missing:
phy2.4fitch <- ReadAsPhyDat('matrix2-4_prepped_for_fitch.nex')
# Some sample trees:
trees2.4 <- read.nexus('matrix2-4.nex')
t1 <- tree2.4[[1]]

#Running traditional Fitch:
fitch(t1, phy2.4fitch)
#Running morphy:
Fitch(t1, phy2.4)
#Running hsj scorer with alpha = 0.5:
hsjTS(t1,phy2.4,mor2.4,typ2.4,alpha=0.5)
#Scoring for different values of alpha:
hsjTS(t1,phy2.4,mor2.4,typ2.4,alpha=0.1)
hsjTS(t1,phy2.4,mor2.4,typ2.4,alpha=1.0)

#Repeating for second example tree:
t2 <- trees2.4[[2]]
fitch(t2, phy2.4fitch)
Fitch(t2, phy2.4)
hsjTS(t2,phy2.4,mor2.4,typ2.4,alpha=0.5)


#Repeating for third example tree:
t3 <- trees2.4[[3]]
fitch(t3, phy2.4fitch)
Fitch(t3, phy2.4)
hsjTS(t3,phy2.4,mor2.4,typ2.4,alpha=0.5)



##### Second example:
#Read in the data files:
mor4.3 <- ReadMorphNexus('matrix4-3.nex')
typ4.3 <- read.table('type4-3.txt')
phy4.3 <- ReadAsPhyDat('matrix4-3.nex')
phy4.3fitch <- ReadAsPhyDat('matrix4-3_prepped_for_fitch.nex')
trees4.3 <- read.nexus('matrix4-3.nex')

#Running methods on first tree:
t1 <- trees4.3[[1]]
fitch(t1, phy4.3fitch)
Fitch(t1, phy4.3)
hsjTS(t1, phy4.3, mor4.3, typ4.3, alpha = 0.5)

#Repeating for second example tree:
t2 <- trees4.3[[2]]
fitch(t2, phy4.3fitch)
Fitch(t2, phy4.3)
hsjTS(t2, phy4.3, mor4.3, typ4.3, alpha = 0.5)

#Repeating for third example tree:
t3 <- trees4.3[[3]]
fitch(t3, phy4.3fitch)
Fitch(t3, phy4.3)
hsjTS(t3, phy4.3, mor4.3, typ4.3, alpha = 0.5)

