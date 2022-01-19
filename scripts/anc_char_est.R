#testing carnivore chromosome data and range size in chromePlus

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
rm(tree, i, missing, chroms, dat.pruned, range, hit, x)

#read in chromeplus data
chromplus <- read.csv("../results/carn.med.hb.csv")

###MAKE Q MATRIX###-------------------------------------------------------------
#store chrom.range
chrom.range <- range(datalist$hap.chrom) + c(-1, 1)

p.mat <- list()
#run datatoMatrix function necessary for diversitree likelihood function
p.mat <- datatoMatrix(x=datalist, range = chrom.range, hyper = T)

#get the ancestral states
source("helperfunctions.R")

#create Q matrix from given data
Q <- getQ(data = p.mat,
          hyper = T,
          polyploidy = F)

#fill in the Q matrix with data from chromplus
Qfilled <- fillQ(data = p.mat,
                 Q = Q,
                 hyper = T,
                 polyploidy = F)

#calculate tip rates given the probability and Q matrix
tiprates <- GetTipRates(tree = tree,
                        Q = Qfilled,
                        tip.states = NULL,
                        tip.probability = p.mat,
                        hyper = T)
all_data <- cbind(datalist, tiprates)
barplot(height = all_data$tiprates + 0.025,
        col = viridis(2, option = "G", end = 0.6)[all_data$range.size + 1],
        horiz = T,
        xlim = c(0,3),
        main = "Range Size")
