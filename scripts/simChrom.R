#Michelle Jonika
#12 September
#This code runs the code for SimChrom analysis of the carnivore data


###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

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
#columns are out of order lets fix them immediately
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

#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  if(i == 52){
    cur.tree <- multi2di(cur.tree)
    
  }
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

###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------

#store chrom.range
chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1) #(15, 40)
chrom.range <- as.numeric(chrom.range)
median <- median(datalist[[52]]$hap.chrom) #19

#load in the results from running the chromplus model
results <- read.csv("../results/carn.med.hb.csv")

#create a list to store simChrom results in
results.sim <- list()
#loop through to create the rates needed as input into simChrom and run
#the chromosome simulations across the carnivore data
for(i in 1:100){
  #store sampling a row from the ChromPlus results to give rates needed
  #for the simChrom simulations
  hit <- sample(1:nrow(results), 1)
  #create the parameters vector needed for simChrom simulations
  pars <- unlist(c(results[hit, c(2,4,3,5)],
                 0,0,0,0,
                 results[hit, 6:7],
                 median,sample(0:1,1,.5)))
  #runs simChrom simulations and stores in results.sim
  results.sim[[i]] <- simChrom(tree=trees.pruned[[i]], 
                               pars=pars, 
                               limits=c(1,100), 
                               model="ChromPlus" )
}
#create vector to store the variance of simChrom results
chrom.variance <- c()
#loop through to calculate variance in simChrom results
for(i in 1:100){
  chrom.variance[i] <- var(results.sim[[i]]$chrom.num)
}

#load in chromosome data
chroms.num <- read.csv("../data/chroms.csv")
#load in range size
range.num <- read.csv("../data/calc.carn.range.sizes.csv")







hap.chrom <- c()
for(i in 1:nrow(range.num)){
    hit <- which(chroms.num$species == range.num$X[i])
    if(length(hit) > 1){
      hit <- sample(hit, 1)
    }
    hap.chrom[i] <- chroms.num[hit, 2]
  }


emp.var <- var(hap.chrom)

plot(y=rep(1,))





hist(chrom.variance)
plot(density(chrom.variance),
     xlim = c(-0.01, 40))
abline(v = emp.var, col = "red")






