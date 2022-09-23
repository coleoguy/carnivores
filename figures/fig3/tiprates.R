#plotting carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(castor)
library(ape)
library(geiger)
library(chromePlus)
library(diversitree)
library(viridis)

###LOAD IN AND PRUNE DATA###----------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method = "extend")
}
#only keep the first tree
tree <- trees[[1]]
#clean up the environment
rm(trees, i)

#load in chromosome data
chroms <- read.csv("../data/chroms.csv")

#load in range size
range <- read.csv("../data/calc.carn.range.sizes.csv")
#change column names to be informative
colnames(range) <- c("species", "range.size")

#prune chromosome number and combnine with range size
dat.pruned <- range
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
#Columns are out of order, lets fix them 
dat.pruned <- dat.pruned[, c(1, 3, 2)]

#this loop creates 100 datasets, sampling a chromosome number for each species
#when there is more than one
datalist <- list()
for(i in 1:nrow(range)){
  hit <- which(chroms$species == range$species[i])
  if(length(hit) > 1){
    hit <- sample(hit, 1)
  }
  dat.pruned$hap.chrom[i] <- chroms[hit, 2]
  datalist <- dat.pruned
}

#prune and scale trees
#find tips that are missing from the dataset
missing <- tree$tip.label[!tree$tip.label %in% datalist$species]

#drop any missing tips between tree and data
tree <- drop.tip(tree, tip = missing)

#discretize range size based on the median
x <- median(datalist$range.size)
datalist$range.size <- as.numeric(datalist$range.size >= x)
#0 = small; 1 = large pop size

#rm old data and clean up environment
rm(i, missing, chroms, dat.pruned, range, hit, x)

#load in the tip rate data
tips <- read.csv("../results/tiprates.csv", row.names = 1)


#loop to prune the chromsome number to those in the final chromosome number
#dataset
m <- matrix(nrow = 110, ncol = 3)
datalist.ordered <- data.frame(m)
for(i in 1:nrow(datalist)){
  hit <- which(datalist$species == rownames(tips)[i])
  datalist.ordered[i, ] <- datalist[hit, ]
}

#bind in the tip rates in the appropriate order
all_data_order <- cbind(datalist.ordered, tips)
colnames(all_data_order) <- c("species", "hap.chrom", "range.size", "tip.rate")


#plot the tip rate data
barplot(height = all_data_order$tip.rate + 0.025,
        col = c("#FDE725FF", "#39568CFF")[all_data_order$range.size + 1],
        horiz = T,
        xlim = c(0,3),
        xlab = "Tip Rates")

#add a legend to the plot
legend(x = "topright", 
       legend = c("Small Range Size", "Large Range Size"), 
       pch = 22, 
       pt.cex = 2, 
       box.col = "transparent", 
       pt.bg = c("#FDE725FF", "#39568CFF"))

#export as PDF 8.5"x11"





