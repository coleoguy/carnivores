#Michelle Jonika
#November 2, 2021

#calculates the tip rate and average tip rate for chromosome number
#across 100 trees
library(phytools)
library(geiger)


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


###CALCULATE TIP RATE###--------------------------------------------------------
#pull out data needed for calculating tip rates
chrom <- datalist[[1]]$hap.chrom
names(chrom) <- datalist[[1]]$species

chrom2 <- chrom
#loop through calculating
for(k in 1:length(trees.pruned)){
  print(k)
  for(i in 1:length(chrom)){
    #store current trees tip labels
    current <- trees.pruned[[k]]$tip.label[i]
    #which tip in the tree matches the chromosome number 
    hit <- which(names(chrom) == current)
    #store the chromosome number match in the chrom 
    chrom2[i] <- chrom[hit]
  }
  # estimate ancestral states
  foo <- ace(x = chrom2, phy = trees.pruned[[k]], model = "ML")
  # get tip branches
  tip.branch <- c()
  for(i in 1:nrow(trees.pruned[[k]]$edge)){
    #store the tip number
    val <- trees.pruned[[k]]$edge[i, 2]
    #if the tip value is not in the edge value then it is a terminal tip
    if(!val %in% trees.pruned[[k]]$edge[, 1]){
      #if a terminal tip then store that value in the vector
      tip.branch <- c(tip.branch, i)
    }
  }
  #create an empty vector for the ancestral state
  anc.state <- c()
  for(j in 1:length(tip.branch)){
    # get node chromosome number estimate
    trees.pruned[[k]]$edge.length[tip.branch]
    #get the node number in the tree
    wanted.node <- trees.pruned[[k]]$edge[tip.branch[j], 1]
    #store the ancestral state estimate for the given node
    anc.state[j] <- as.numeric(foo$ace[names(foo$ace) == wanted.node])
  }
  #create an empty vector for the current tip state
  curr.state <- c()
  for(j in 1:length(tip.branch)){
    #get the current species we are working on
    curr.sp <- trees.pruned[[k]]$tip.label[j]
    #store the current species' current chromosome state
    curr.state <- chrom[names(chrom) == curr.sp]
  }
  #print the ancestral state - the current state
  print(anc.state - curr.state)
  #calculate the tip rates by substracting the current state from the ancestral
  #state and dividing by the branch length
  tip.rates <- (anc.state - curr.state) /
    trees.pruned[[k]]$edge.length[tip.branch]
  #name the tip rates based on the tip label
  names(tip.rates) <- trees.pruned[[k]]$tip.label
  if(k == 1){
    tipp.rates <- tip.rates
  }else{
    tipp.rates <- cbind(tipp.rates, tip.rates)
  }
}
#create informative column names for each trees values
colnames(tipp.rates) <- paste("tree", 1:100)

#average the sum of the tip rates across all trees
Average <- rowSums(tipp.rates)/100
#bind the value of the averages to the results with all tip rates
tipp.rates <- cbind(tipp.rates, Average)
#write out a csv storing all the tip rate and average tip rate values
write.csv(tipp.rates, file = "../results/tip.rates.csv")
