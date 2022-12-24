#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(castor)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/datalists_range.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in tree depths
tree.depths <- read.csv("../data/treedepths.csv")
colnames(tree.depths) <- c("tree", "tree.depth")

#store chrom.range
chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)

#empty vector to store diversitree likelihood  
p.mat <- list()
#run datatoMatrix function necessary for diversitree likelihood function
for(i in 1:100){
  p.mat[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
}

#read in chromeplus data
chromplus <- read.csv("../results/rangesize.csv")

#get the ancestral states
source("getQ.function.R")

tiprates <- list()
for(i in 1:length(p.mat)){
  #create Q matrix from given data
  Q <- getQ(data = p.mat[[1]],
            hyper = T,
            polyploidy = F)
  
  # here there are rates that we are not interested at
  # discard them from the table
  Q[!(Q %in% c(1,2,3,4,8,9))] <- 0
  
  # fill the qmat
  Q[Q == 1] <- mean(chromplus$asc1)
  Q[Q == 2] <- mean(chromplus$desc1)
  Q[Q == 3] <- mean(chromplus$asc2)
  Q[Q == 4] <- mean(chromplus$desc2)
  Q[Q == 8] <- mean(chromplus$tran12)
  Q[Q == 9] <- mean(chromplus$tran21)
  
  # fill the diagonal so that row sums are zero
  diag(Q) <- -rowSums(Q)
  
  #calculate tip rates given the probability and Q matrix
  tiprates[[i]] <- GetTipRates(tree = trees[[i]],
                          Q = Qfilled,
                          tip.states = NULL,
                          tip.probability = p.mat[[i]],
                          hyper = T)
  
}

#need to back transform the rates by the tree depths
tiprate_final <- list()
for(i in 1:length(tiprates)){
  tiprate_final[[i]] <- tiprates[[i]]/tree.depths[i,2]
}

#combine all estimates for tip rates across the posterior distribution
tiprate_all <- do.call("cbind", tiprate_final)

#calculate the average tip rate across the posterior distribution
for(i in 1:nrow(tiprate_all)){
  average <- rowMeans(tiprate_all[,])
}
tiprates_complete <- cbind(tiprate_all, average)
labels <- c()
for(i in 1:100){
  labels[i] <- paste("Tree", i)
  labels[101] <- "Average"
}
colnames(tiprates_complete) <- labels
#write out the results for plotting later
write.csv(tiprates_complete, "../results/tiprates_new.csv")
