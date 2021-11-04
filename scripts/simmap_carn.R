#load in necessary packages
library(phytools)
library(geiger)

#read in carnivore phylogenies
trees <- read.nexus("carnivora.nex")
for(i in 1:100){
  trees[[i]] <- force.ultrametric(trees[[i]], method="extend")
}

#read in range size data 
dat <- read.csv("carn.discrete.csv", as.is = T)
rownames(dat) <- dat[,1]

#make range size data a names vector of discrete traits
carn_dat <- dat[,-c(1:2)]
names(carn_dat) <- row.names(dat)

#make a vector to store the trees
trees.pruned <- list()
#loop through to drop tips that are not in our dat from the trees
for(i in 1:100){
  trees.pruned[[i]] <- treedata(phy = trees[[i]], data = carn_dat)[[1]]
}
rm(trees, i)

#make a vector to store the simmaps
histories <- list()

#loop through making simmaps for each of the 100 trees
for(i in 1:100){
  print(paste("working on tree", i))
  histories[[i]] <- make.simmap(trees.pruned[[i]], 
                                carn_dat, 
                                model="ARD", 
                                pi="estimated")
}

#make the class of the simmaps of type Phylo
class(histories) <- "Phylo"

plot <- ReorderData(trees.pruned[[1]], dat, taxa.names = row.names(dat))
barplot(dat$x)
abline(h=29301531.71, col = c("blue"))
sp <- c()
for(i in 1:nrow(carn_dat)){
  hit <- which(names(carn_dat)==trees.pruned[[i]]$tip.label)
  foo[i, ] <- carn_dat[hit, ]
  sp[i] <- row.names(carn_dat)[hit]
}
row.names(foo) <- sp