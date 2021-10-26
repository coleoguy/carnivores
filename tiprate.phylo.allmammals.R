#set working directory as figures/heatmap
#load in packages that are needed
library(phytools)

#load in tree data
trees <- read.nexus("carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance 
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method="extend")
}

#load in chromosome data
chroms <- read.csv("mammal_chroms_update.csv")
chroms <- chroms[,5:6]
chroms[,2] <- as.numeric(chroms[,2])/2
chroms <- na.omit(chroms)

###PRUNE DATA###----------------------------------------------------------------

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% chroms$tree.name]
#empty list to store pruned trees
trees.drop <- list()

#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  trees.drop[[i]] <- drop.tip(trees[[i]], tip = missing)
}

#remove tree 52 which is not rooted or fully dichotomous
trees.pruned <- list()
trees.pruned <- trees.drop[-c(52)]

#select a single tree to use
pruned.tree <- trees.pruned[[sample(1:99, 1)]]
rm(trees, trees.pruned)

#load in tip rate data
tips <- read.csv("tip.rates_allmammals.csv", row.names = 1)
#subset out columns that we need
tips <- tips[100]

#this loop will make sure the tip labels and the chromosome data are in the
#same order
foo2 <- tips
sp2 <- c()
for(i in 1:nrow(tips)){
  hit2 <- which(row.names(tips)==pruned.tree$tip.label[i])
  foo2[i, ] <- tips[hit2, ]
  sp2[i] <- row.names(tips)[hit2]
}
row.names(foo2) <- sp2

tips <- foo2$Average
names(tips) <- row.names(foo2)

#plot tree with bars
plotTree.wBars(tree = pruned.tree,
               x = abs(tips))

#export as pdf 7" x 7"