library(ape)
library(chromePlus)
library(castor)
library(phytools)
#set the seed
set.seed(1)
#make a tree
tree <- rcoal(10)

#make the data; we will need chromosome number and a binary trait
species <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")
chrom <- c(19, 44, 19, 44, 19, 19, 19, 19, 19, 44)
trait <- c("A", "A", "B", "B", "A", "A", "B", "B", "A", "B")
data <- data.frame(species, chrom, trait)

#make the tree unit length
tree_depth <- max(branching.times(tree))

tree$edge.length <-  tree$edge.length/max(branching.times(tree))
tree_unit <- tree

#calculate the range of chromosome number
chrom_range <- range(data$chrom) + c(-1, 1)

#change the binary trait to the likelihood of being in state 1/2
data$trait <- c(0,0,1,1,0,0,1,1,0,1)

#empty vector to store diversitree likelihood  
p.mat <- list()
#run datatoMatrix function necessary for diversitree likelihood function
p.mat <- datatoMatrix(x=data, range = chrom_range, hyper = T)

#read in a function to create a Q matrix
source("getQ.function.R")
#use function to create a Q matrix for the given data
Q <- getQ(data = p.mat,
          hyper = T,
          polyploidy = F)

# here there are rates that we are not interested at
# discard them from the table
Q[!(Q %in% c(1,2,3,4,8,9))] <- 0

# fill the qmat
Q[Q == 1] <- 0
Q[Q == 2] <- 0
Q[Q == 3] <- 0.003
Q[Q == 4] <- 0.003
Q[Q == 8] <- 0
Q[Q == 9] <- 0.003

# fill the diagonal so that row sums are zero
diag(Q) <- -rowSums(Q)

rm(list=ls()[-c(5,6,9)])

#read in function to calculate tip rates
source("tip.rates.function.R")

tip.rates <- c()
tip.rates <- GetTipRates(tree = tree, 
                         Q = Q, 
                         tip.states = NULL, 
                         hyper = TRUE, 
                         p.mat = p.mat)
