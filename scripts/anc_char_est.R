#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(castor)
library(plotrix)
library(coda)
library(ggplot2)
library(viridis)
library(ape)

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


###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------

#store chrom.range
chrom.range <- range(datalist[[50]]$hap.chrom) + c(-1, 1)

p.mat <- list()
#run datatoMatrix function necessary for diversitree likelihood function
for(i in 1:100){
  p.mat[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
}

#make a likelihood function
#states = named vector of character states
#k = number of states to model; 1 to k for make.mkn number of columns in pmat
#strict = T allows for missing states
#control (ode) = uses an ODE based approach to compute only the k variables over
#time; more efficient when k is large
lk.mk <- make.mkn(trees.pruned[[50]], states = p.mat[[50]],
                  k = ncol(p.mat[[50]]), strict = F,
                  control = list(method = "ode"))

###CONSTRAIN LIKELIHOOD FUNCTION###---------------------------------------------

#constrain the likelihood function to remove states that are not biologically
#realistic
con.lik <- constrainMkn(data = p.mat[[50]],
                        lik = lk.mk,
                        polyploidy = F,
                        hyper = T,
                        verbose = T, 
                        oneway = F,
                        constrain = list(drop.demi = T,
                                         drop.poly= T))

###READ IN CHROMPLUS DATA###----------------------------------------------------
chromplus <- read.csv("../results/carn.med.hb.csv")

###GET Q MATRIX MAT###----------------------------------------------------------
#get q matrix matrix
parMat <- con.lik[[2]]
#not interested in all of the rates, discard those not needed from table
parMat[!(parMat %in% c(1,2,3,4,8,9))] <- 0

# rate1 ascending aneuploidy - state1 
# rate2 descending aneuploidy - state1 
# rate3 ascending aneuploidy - state2 
# rate4 descending aneuploidy - state2 
# rate8 transitions from 1 to 2 for hyperstate 
# rate9 transitions from 2 to 1 for hyperstate

#fill the q matrix
parMat[parMat == 1] <- mean(chromplus$asc1)
parMat[parMat == 2] <- mean(chromplus$desc1)
parMat[parMat == 3] <- mean(chromplus$asc2)
parMat[parMat == 4] <- mean(chromplus$desc2)
parMat[parMat == 8] <- mean(chromplus$tran12)
parMat[parMat == 9] <- mean(chromplus$tran21)

#rename the columns so that states now start from 1
states <- 1:ncol(p.mat[[50]])
names(states) <- colnames(p.mat[[50]])
colnames(parMat) <- rownames(parMat) <- colnames(p.mat[[50]]) <- states

#make row sum to 0
for(i in 1:nrow(parMat)){
  parMat[i,i] <- 0 -sum(parMat[i,])
}

#choose a single tree to use in ancestral character estimation
tree <- trees.pruned[[50]]

#sort p.mat to match with the tree tip order
p.mat[[50]] <- p.mat[[50]][tree$tip.label,]
#sort MCMC data to match with the tree tip order
rownames(datalist[[50]]) <- datalist[[50]]$species
datalist[[50]] <- datalist[[50]][tree$tip.label, ]

#get the ancestral states
asr <- asr_mk_model(tree = tree, 
                    tip_states = NULL, 
                    tip_priors = p.mat[[50]],
                    transition_matrix = parMat,
                    Nstates = ncol(parMat),
                    include_ancestral_likelihoods = T,
                    root_prior = "empirical")

# make states
rng <- range(datalist[[50]]$hap.chrom) + c(-2,2)
chrom.states <- as.data.frame(matrix(data = NA,
                                     nrow = length(rng[1]:rng[2]) * 2,
                                     ncol = 3))
colnames(chrom.states) <- c("chroms", "binary", "state")
# fill chrom states
chrom.states$chroms <- rep(rng[1]:rng[2], 2)
chrom.states$binary <- rep(c(1,0), each = length(rng[1]:rng[2]))

for(i in 1:nrow(chrom.states)){
  if (nrow(chrom.states) < 100) 
    pad <- 2
  if (nrow(chrom.states) >= 100) 
    pad <- 3
  if (nrow(chrom.states) < 10) 
    pad <- 1
  chrom.states$state[i] <- sprintf(paste("%0", pad, "d", sep = ""),i)
}
# get the most likely chromosome number state for each node
node.chroms <- as.data.frame(matrix(data = NA,
                                    nrow = Nnode(tree),
                                    ncol = 3))
colnames(node.chroms) <- c("node", "state", "chromNum")
# fill the table
for(i in 1:nrow(node.chroms)){
  node.chroms$node[i] <- i + Ntip(tree)
  poss.state <- which.max(asr$ancestral_likelihoods[i,])
  if(length(poss.state) > 1){
    node.chroms$state[i] <- sample(poss.state,1) 
  }else{
    node.chroms$state[i] <- poss.state     
  }
}

for(i in 1:nrow(node.chroms)){
  hit.state <- which(as.numeric(chrom.states$state) == node.chroms$state[i])
  node.chroms$chromNum[i] <- chrom.states$chroms[hit.state]
}

# get chromosome number at the tip
dat.tips <- as.data.frame(matrix(data = NA,
                                 nrow = Ntip(tree),
                                 ncol = 2))
colnames(dat.tips) <- c("tip", "chrom")
for(i in 1:Ntip(tree)){
  dat.tips$tip[i] <- i
  dat.tips$chrom[i] <- datalist[[50]]$hap.chrom[datalist[[50]]$species == tree$tip.label[i]]
}
# get the branch specific rates
branch.rates <- as.data.frame(matrix(data = NA,
                                     nrow = Nedge(tree),
                                     ncol = 3))
colnames(branch.rates) <- c('branch', "rate", "rateClass")
counter <- 1
for(i in 1:nrow(node.chroms)){
  hit.edge <- which(tree$edge[,1] == node.chroms$node[i])
  for(j in 1:length(hit.edge)){
    end.node <- tree$edge[hit.edge[j],2]
    if(end.node %in% c(1:Ntip(tree))){
      br.rate <- abs(dat.tips$chrom[end.node] - node.chroms$chromNum[i]) / tree$edge.length[hit.edge[j]]
      branch.rates$branch[counter] <- counter
      branch.rates$rate[counter] <- br.rate
      counter <- counter + 1
    }
    else{
      br.rate <- abs(node.chroms$chromNum[node.chroms$node == end.node] - node.chroms$chromNum[i]) / tree$edge.length[hit.edge[j]]
      branch.rates$branch[counter] <- counter
      branch.rates$rate[counter] <- br.rate
      counter <- counter + 1
    }
  }
}

# get quantile
low <- 0.25
high <- 0.90

# process branch rates
# Assign a rate class for each branch depending on the branch specific rates
for(i in 1:nrow(branch.rates)){
  if(branch.rates$rate[i] <= quantile(x = branch.rates$rate, low)){
    branch.rates$rateClass[i] <- 1
  }
  if(branch.rates$rate[i] > quantile(x = branch.rates$rate, low) & branch.rates$rate[i] < quantile(x = branch.rates$rate, high)){
    branch.rates$rateClass[i] <- 2
  }
  if(branch.rates$rate[i] >= quantile(x = branch.rates$rate, high)){
    branch.rates$rateClass[i] <- 3
  }
}

# colour the tip branches only
branch.rates$rateClass[!(tree$edge[,2] %in% (1:110))] <- 0

# make a maps object and add it to tree in order to colour branches
maps <- vector(mode = "list", length = Nedge(tree))
for(i in 1:length(maps)){
  maps[[i]] <- tree$edge.length[i]
  names(maps[[i]]) <- branch.rates$rateClass[i]
}
tree$maps <- maps

# replace chromosome number with branch rates
# isolate tip branches
tip.branches <- matrix(data = NA,
                       nrow = Ntip(tree),
                       ncol = 3)
counter <- 1
for(i in 1:nrow(tree$edge)){
  if(tree$edge[i,2] %in% (1:Ntip(tree))){
    tip.branches[counter,1] <- counter
    tip.branches[counter,2] <- i
    tip.branches[counter,3] <- branch.rates$rate[branch.rates$branch == i]
    counter <- counter + 1
  }
}
datalist[[50]]$tipRates <- tip.branches[,3]
datalist[[50]]$tipRates
write <- datalist[[50]]

#save 
write.csv(write, "../results/rates.csv")
