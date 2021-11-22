#Michelle Jonika
#september 7, 2021

#creates a phylogenetic plot visualizing chromosome number and range size as a
#discrete trait

###LOAD IN PACKAGES###----------------------------------------------------------

library(phytools)
library(chromePlus)
library(viridis)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method = "extend")
}
#load in chromosome data
chroms <- read.csv("../data/chroms.csv")

#load in range size
range <- read.csv("../data/calc.carn.range.sizes.csv")
#change column names to be informative
colnames(range) <- c("species", "range.size")

###PRUNE DATA###----------------------------------------------------------------

#prune chromosome number and combnine with range size
dat.pruned <- range
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
# TODO columns are out of order lets fix them immediately
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

###DISCRETIZE RANGE SIZES###----------------------------------------------------
#discretize range size based on the median
for(i in 1:100){
  x <- median(datalist[[1]]$range.size)
  datalist[[i]]$range.size <- as.numeric(datalist[[1]]$range.size >= x)
}
#0 = small; 1 = large pop size
rm(i, x)

###CONTINUOUS TRAIT MAP###------------------------------------------------------

#create a vector of chromosome numbers for barplot
chroms <- datalist[[1]]$hap.chrom
names(chroms) <- datalist[[1]]$species

#use plot tree with bars to create a phylogenetic tree with barplots for the 
#chromosome number
plotTree.barplot(tree = trees.pruned[[1]],
                 x = chroms,
                 lwd=4,
                 args.plotTree = 
                   list(ftype = "off"),
                 args.barplot = 
                   list(col = viridis(2, option= "G", 
                                      end = 0.6)[datalist[[1]]$range.size + 1],
                        xlab = "Chromosome Number"),
                 args.axis = list(at = seq(0,40, by = 5)))
legend(x = "right", 
       legend = c("Large Range Size", "Small Range Size"), 
       pch = 22, 
       pt.cex = 2, 
       box.col = "transparent", 
       pt.bg = viridis(2, option= "G", end = 0.6))

#export 9"x9"






