#this script goes thorugh posterior predictive modelling for 
#the carnivore data

###write a loop that stores each tree's data from carn.all in its own list---###

#load in data
carn.all <- read.csv("results_75small25large/carn_data_all.csv")

#empty list to store
dat_sep <- list()
#vector that stores the beginning row value for each tree
x  <- seq(1, 245001, by=2500)

counter <- 1
#loop that stores each tree's rate data into a separate list
for(i in x){
  dat_sep[[counter]] <- carn.all[i:(i+2499), ]
  counter <- counter + 1
}

###POSTERIOR PREDICTIVE MODELS FOR 25 SMALL/75 LARGE-------------------------###

#LOAD IN DATA

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

#PRUNE DATA

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

#DISCRETIZE RANGE SIZES

#look at a histagram of the data
hist(dat.pruned$range.size, breaks =200)
#adds a line to the histogram where the 75th quantile is
abline(v=quantile(dat.pruned$range.size, 0.75), col = "blue")
#stores as a variable the 75th quantile cutoff value
quant <- quantile(dat.pruned$range.size, 0.75)
#assigns the range values into two groups for the model.
#0 = small; 1 = large pop size
dat.pruned$range.size <- as.numeric(dat.pruned$range.size > quant)
#quant is the the threshold cutoff for discretization
#quant: ***6.7E7*** [Replace if modified]
#rm quant variable to clean up environment
rm(quant)


#ANCESTRAL CHARACTER ESTIMATION FOR ROOT HYPERSTATE
library(phytools)

#make named character vector with binary trait states
range_bin <- dat.pruned$range.size 
names(range_bin) <- dat.pruned$species

#run ancestral character estimation
ace <- ace(range_bin, trees.pruned[[1]], type = "discrete", method = "ML")

#remove tree 52 which couldn't produce results for the model
trees <- trees.pruned[-c(52)]

#store median chromosome number
median <- 19

hyp <- c(0,1)
hyperstate_samp <- sample(hyp, size = 99, replace = T, 
                          prob = c(0.7455696, 0.2544304))
sim <- list()
for(i in 1:99){
  dat <- dat_sep[[i]][dat_sep[[i]]$i >= 2000, ]
  samp <- dat[sample(nrow(dat), 1), ]
  samp <- samp[c(2,4,3,5,6:7)] %>%
    add_column(demi1 = "0",
               demi2 = "0",
               poly1 = "0",
               poly2 = "0",
               .before = 5) %>%
    add_column(chrom_root = median,
               hyper_root = hyperstate_samp[i],
               .after = 10)
  rates <- as.numeric(samp[1, ])
  sim[[i]] <- simChrom(tree = trees[[i]], pars = rates, 
                       limits = c(14, 41), model = "ChromPlus")
}

med <- avg <- ranges <- c()

for(i in 1:99){
  med[i] <- median(sim[[i]]$chrom.num)
  avg[i] <- mean(sim[[i]]$chrom.num)
  temp <- range(sim[[i]]$chrom.num)
  ranges[i] <- (temp[2] - temp[1])
}

hist(med)
abline(v=19)
hist(avg)
abline(v=21.72727)
hist(ranges)
abline(v=23)

