#load in necessary packages
library(phytools)
library(geiger)

#load in tree data
trees <- read.nexus("carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method="extend")
}

#load in chromosome data
chroms <- read.csv("chroms.csv")

#load in range size
range <- read.csv("calc.carn.range.sizes.csv")
#change column names to be informative
colnames(range) <- c("species", "range.size")

###PRUNE DATA###----------------------------------------------------------------

#prune chromosome number and combnine with range size
dat.pruned <- range
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
#this loop samples a chromosome number for each species when there is more than
#one
for(i in 1:nrow(range)){
  hit <- which(chroms$species == range$species[i])
  if(length(hit)>1)  hit <- sample(hit, 1)
  dat.pruned[i, 3] <- chroms[hit, 2]
}
#rm old data and clean up environment
rm(chroms, range, hit, i)

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% dat.pruned$species]
#empty list to store pruned trees
trees.drop <- list()

#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  trees.drop[[i]] <- drop.tip(trees[[i]], tip = missing)
}


#pull out data needed for calculating tip rates
chrom <- dat.pruned$hap.chrom
names(chrom) <- dat.pruned$species
rm(dat.pruned)

#remove tree 52 which is not rooted or fully dichotomous
trees.pruned <- list()
trees.pruned <- trees.drop[-c(52)]

rm(trees, trees.drop, i , missing)


chrom2 <- chrom

#loop through calculating
for(k in 1:99){
  print(k)
  for(i in 1:length(chrom)){
    current <- trees.pruned[[k]]$tip.label[i]
    hit <- which(names(chrom) == current)
    chrom2[i] <- chrom[hit]
  }
  # estimate ancestral states
  foo <- ace(x = chrom2, phy = trees.pruned[[k]], model = "ML")
  # get tip branches
  tip.branch <- c()
  for(i in 1:nrow(trees.pruned[[k]]$edge)){
    val <- trees.pruned[[k]]$edge[i, 2]
    if(!val %in% trees.pruned[[k]]$edge[, 1]){
      tip.branch <- c(tip.branch, i)
    }
  }

  anc.state <- c()
  for(j in 1:length(tip.branch)){
    # get node microsatellite estimate
    trees.pruned[[k]]$edge.length[tip.branch]
    wanted.node <- trees.pruned[[k]]$edge[tip.branch[j], 1]
    anc.state[j] <- as.numeric(foo$ace[names(foo$ace) == wanted.node])
  }
  curr.state <- c()
  for(j in 1:length(tip.branch)){
    curr.sp <- trees.pruned[[k]]$tip.label[j]
    curr.state <- chrom[names(chrom) == curr.sp]
  }
  print(anc.state - curr.state)
  tip.rates <- (anc.state - curr.state) /
    trees.pruned[[k]]$edge.length[tip.branch]
  names(tip.rates) <- trees.pruned[[k]]$tip.label
  if(k == 1){
    tipp.rates <- tip.rates
  }else{
    tipp.rates <- cbind(tipp.rates, tip.rates)
  }
}
colnames(tipp.rates) <- paste("tree", 1:99)

#tipp.rates[,] <- abs(tipp.rates)
Average <- rowSums(tipp.rates)/99
tipp.rates <- cbind(tipp.rates, Average)
write.csv(tipp.rates, file = "tip.rates.csv")
