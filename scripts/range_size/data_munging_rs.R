#Michelle Jonika
#October 25, 2022

#scrip that matches up tree data and chromosome number data, as well as 
#discretizes our binary trait

###LOAD IN PACKAGES###----------------------------------------------------------

library(phytools)
library(chromePlus)
library(viridis)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method = "extend")
}
#load in chromosome data
chroms <- read.csv("../../data/chroms.csv")

#load in range size
range <- read.csv("../../data/range_size/range_size.csv")
#change column names to be informative
colnames(range) <- c("species", "range.size")

###PRUNE DATA###----------------------------------------------------------------

#prune chromosome number and combnine with range size
dat.pruned <- range
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
# Columns are out of order lets fix them immediately
dat.pruned <- dat.pruned[, c(1, 3, 2)]

#this loop creates 100 datasets, sampling a chromosome number for each species 
#when there is more than one
datalist <- list()
for(j in 1:100){
  for(i in 1:nrow(range)){
    hit <- which(chroms$species == range$species[i])
    if(length(hit) > 1){
      hit <- sample(hit, 1)
    }
    dat.pruned$hap.chrom[i] <- chroms[hit, 2]
  }
  datalist[[j]] <- dat.pruned
}

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% datalist[[1]]$species]
#empty list to store pruned trees
trees.pruned <- list()
#empty vector to store tree depths
#keep tree depths to correct for depth when analysing rates
tree.depths <- c()

#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  tree.depths[i] <- max(branching.times(cur.tree))
  # a rounding error also makes tree 52 fail due to a multitomy in the genus
  # pusa
  if(i == 52){
    cur.tree <- multi2di(cur.tree)
    
  }
  cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
  trees.pruned[[i]] <- cur.tree
}

#rm old data and clean up environment
rm(trees, cur.tree, i, j, missing, chroms, dat.pruned, range, hit)

#write out our pruned trees
write.nexus(trees.pruned, file = "../../data/range_size/carnivora_rs_pruned.nex")

#write out our prunned tree's tree depths
write.csv(tree.depths, file = "../../data/range_size/rs_treedepths.csv")

###DISCRETIZE RANGE SIZES###----------------------------------------------------
#discretize range size based on the median
for(i in 1:100){
  x <- median(datalist[[1]]$range.size)
  datalist[[i]]$range.size <- as.numeric(datalist[[1]]$range.size >= x)
}

#0 = small; 1 = large pop size
rm(i, x)

###ORDER TREE AND DATA###-------------------------------------------------------
#loop through to make sure all of the datalists are in the same order as each
#of the trees
for(i in 1:length(datalist)){
  #empty vector to store correct order in 
  neworder <- c()
  #loop through to store the correct order of the data to match the tree tips
  for(j in 1:length(trees.pruned[[i]]$tip.label)){
    neworder[j] <- which(trees.pruned[[i]]$tip.label[j] == datalist[[i]]$species)
  }
  #reorder the data into a new data frame
  datalist[[i]] <- datalist[[i]][neworder,]
}

#remove variables that aren't needed anymore
rm(i, j, neworder, trees.pruned, tree.depths)

#write out the environment to store the data lists 
save.image("../../data/range_size/datalists_rangesize.RData")


