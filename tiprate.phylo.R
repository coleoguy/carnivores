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

#remove tree 52 which is not rooted or fully dichotomous
trees.pruned <- list()
trees.pruned <- trees.drop[-c(52)]

#select a single tree to use
pruned.tree <- trees.pruned[[sample(1:99, 1)]]
rm(trees, trees.pruned)

#load in tip rate data
tips <- read.csv("tip.rates.csv", row.names = 1)
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
x <- c()
for(i in 1:110){
  x[i] <- tips[which(names(tips) == pruned.tree$tip.label[i])]
}
names(x) <- pruned.tree$tip.label
#plot tree with bars
plotTree.wBars(tree = pruned.tree,
               x = x)

#export as pdf 7" x 7"
plot(pruned.tree, cex=.5)
