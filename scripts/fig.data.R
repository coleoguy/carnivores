#create a phylogenetic plot visualizing chromosome number and range size
#as a discrete trait

###LOAD IN PACKAGES###----------------------------------------------------------

library(phytools)
library(chromePlus)
library(viridis)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivora.nex")
for(i in 1:100){
  trees[[i]] <- force.ultrametric(trees[[i]], method="extend")
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
#add emtpy third column for chromosome number
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
trees.pruned <- list()
#empty vector to store tree depths
#keep tree depths to correct for depth when analysing rates
tree.depths <- c()
#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  tree.depths[i] <- max(branching.times(cur.tree))
  cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
  trees.pruned[[i]] <- cur.tree
}
#rm old data and clean up environment
rm(trees, cur.tree, i,  missing, tree.depths)

###DISCRETIZE RANGE SIZES###----------------------------------------------------

#stores as a variable the 75th quantile cutoff value
quant <- quantile(dat.pruned$range.size, 0.75)
#assigns the range values into two groups for the model.
dat.pruned$range.size <- dat.pruned$range.size < quant
#changes T/F to 1/0
dat.pruned$range.size[dat.pruned$range.size == "TRUE"] <- 1
#quant is the the threshold cutoff for discretization
#quant: ***6.7E7*** [Replace if modified]
#rm quant variable to clean up environment
rm(quant)

###DISCRETE TRAIT MAP###--------------------------------------------------------

#create a vector of tip states for make.simmap function
tips <- dat.pruned$range.size
names(tips) <- dat.pruned$species

#only use the first tree
tree <- trees.pruned[[1]]
#make.simmap requires the tree to be of class "phylo"
as.phylo(tree)

#make simmaps of the tree using the tree and tip states
maps.est <- make.simmap(tree = tree,
                        x = tips,
                        model = "ARD",
                        pi = "estimated",
                        nsim=100)
#create a density map using the simmap trees
smap.est <- densityMap(trees = maps.est,
                       res=50, lwd=.75,
                       check=FALSE,
                       type="fan",
                       direction="rightwards",
                       fsize=.25)
#set the colors for the 2 states
smap.est2 <- setMap(smap.est, viridis(2, end = 0.95, 
                                      option = "G"))
#plot the discrete trait map
plot(smap.est2, fsize=c(.00001, .4),
     lwd=5,type="fan")

#export 6"x6"

###CONTINUOUS TRAIT MAP###------------------------------------------------------
tree <- trees.pruned[1]

#paint subtrees
ss <- getStates(tree, "tips")

#create a vector of chromosome numbers for barplot
chroms <- dat.pruned$hap.chrom
names(chroms) <- dat.pruned$species

#use plot tree with bars to create a phylogenetic tree with barplots for the 
#chromosome number
plotTree.barplot(tree = smap.est2$tree,
               x = chroms,
               args.plotTree = list(ftype = "off"),
               col = (viridis(3, option = "magma")[dat.pruned$range.size + 1]),
               args.axis = list(at = seq(0,50, by = 10)))
add.color.bar(cols = viridis(2, end = 0.95, option = "G"), title = "Range Size")

#export 6"x6"






