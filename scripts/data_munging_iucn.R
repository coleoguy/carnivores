#Michelle Jonika
#November 8, 2022

#scrip that matches up tree data and chromosome number data, as well as 
#discretizes our binary trait

###LOAD IN PACKAGES###----------------------------------------------------------

library(phytools)
library(chromePlus)
library(viridis)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in chromosome data
chroms <- read.csv("../data/chroms.csv")

#load in IUCN data
iucn <- read.csv("../data/iucndata.csv")

###PRUNE DATA###----------------------------------------------------------------

#prune chromosome number and combnine with range size
dat.pruned <- iucn
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
# Columns are out of order lets fix them immediately
dat.pruned <- dat.pruned[, c(1, 3, 2)]
colnames(dat.pruned) <- c("species", "hap.chrom", "status")

#this loop creates 100 datasets, sampling a chromosome number for each species 
#when there is more than one
datalist <- list()
for(j in 1:100){
  for(i in 1:nrow(iucn)){
    hit <- which(chroms$species == iucn$name[i])
    if(length(hit) > 1){
      hit <- sample(hit, 1)
    }
    dat.pruned$hap.chrom[i] <- chroms[hit, 2]
  }
  datalist[[j]] <- dat.pruned
}

#rm old data and clean up environment
rm(i, j, chroms, dat.pruned, iucn, hit)

###DISCRETIZE RANGE SIZES###----------------------------------------------------
#discretize range size based on the median
for(i in 1:100){
  datalist[[i]]$status[which(!datalist[[i]]$status == "Least Concern")] <- 1
  datalist[[i]]$status[which(datalist[[i]]$status == "Least Concern")] <- 0
  datalist[[i]]$status <- as.numeric(datalist[[i]]$status)
}

#0 = least concern; 1 = all other iucn statuses
rm(i)

###ORDER TREE AND DATA###-------------------------------------------------------
#loop through to make sure all of the datalists are in the same order as each
#of the trees
for(i in 1:length(datalist)){
  #empty vector to store correct order in 
  neworder <- c()
  #loop through to store the correct order of the data to match the tree tips
  for(j in 1:length(trees[[i]]$tip.label)){
    neworder[j] <- which(trees[[i]]$tip.label[j] == datalist[[i]]$species)
  }
  #reorder the data into a new data frame
  datalist[[i]] <- datalist[[i]][neworder,]
}

#remove variables that aren't needed anymore
rm(i, j, neworder, trees)

#write out the environment to store the data lists 
save.image("~/Documents/GitHub/carnivores/data/datalists_iucn.RData")


