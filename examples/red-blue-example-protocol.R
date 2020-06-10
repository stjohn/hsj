###script for red/blue tail example matrix

library(TreeSearch)
library(phangorn)
library(Claddis)
library(TreeTools)

setwd("~/Desktop/maddison-example")
source('~/github/phangs/code/phangsScorer7.R')
source('~/github/phangs/code/dissimilarity_functions.R')

# ####
claddisObj<-ReadMorphNexus('matrix2-4_fixed.nex')
typesObj<-read.table('type2-4.txt')
colnames(typesObj)<-c('Char','Type','Sub')

phyDatObj<-ReadAsPhyDat('matrix2-4_fixed.nex')
phyDatObj.fitch<-ReadAsPhyDat('matrix2-4_fitch_fixed.nex')

mydata<-phangs(phyDatObj = phyDatObj,claddisObj = claddisObj,typesObj = typesObj, alpha = 0.5)
tree<-read.nexus('matrix2-4_fixed.nex')
tree.1<-tree[[1]]
fitch(tree = tree.1, data = phyDatObj.fitch)
Fitch(tree = tree.1, dataset = phyDatObj)
hsjScorer(parent = tree.1$edge[,1],child = tree.1$edge[,2],phangsObj = mydata, nTaxa = 14)

tree.2<-tree[[2]]
fitch(tree = tree.2, data = phyDatObj.fitch)
Fitch(tree = tree.2, dataset = phyDatObj)
hsjScorer(parent = tree.2$edge[,1],child = tree.2$edge[,2],phangsObj = mydata, nTaxa = 14)

tree.3<-tree[[3]]
fitch(tree = tree.3, data = phyDatObj.fitch)
Fitch(tree = tree.3, dataset = phyDatObj)
hsjScorer(parent = tree.3$edge[,1],child = tree.3$edge[,2],phangsObj = mydata, nTaxa = 14)

####
claddisObj<-ReadMorphNexus('matrix4-3_fixed.nex')
typesObj<-read.table('type4-3.txt')
colnames(typesObj)<-c('Char','Type','Sub')

phyDatObj<-ReadAsPhyDat('matrix4-3_fixed.nex')
phyDatObj.fitch<-ReadAsPhyDat('matrix4-3_fitch_fixed.nex')

mydata<-phangs(phyDatObj = phyDatObj,claddisObj = claddisObj,typesObj = typesObj, alpha = 0.5)
tree<-read.nexus('matrix4-3_fixed.nex')
tree.1<-tree[[1]]
fitch(tree = tree.1, data = phyDatObj.fitch)
Fitch(tree = tree.1, dataset = phyDatObj)
hsjScorer(parent = tree.1$edge[,1],child = tree.1$edge[,2],phangsObj = mydata, nTaxa = 14)

tree.2<-tree[[2]]
fitch(tree = tree.2, data = phyDatObj.fitch)
Fitch(tree = tree.2, dataset = phyDatObj)
hsjScorer(parent = tree.2$edge[,1],child = tree.2$edge[,2],phangsObj = mydata, nTaxa = 14)

tree.3<-tree[[3]]
fitch(tree = tree.3, data = phyDatObj.fitch)
Fitch(tree = tree.3, dataset = phyDatObj)
hsjScorer(parent = tree.3$edge[,1],child = tree.3$edge[,2],phangsObj = mydata, nTaxa = 14)
