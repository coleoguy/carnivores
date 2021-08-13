# Heath, Michelle, and Sarah testing this old holocentric
# monocentric idea on a tree of insects

# load packages
library(ape) # basic phylo tools
library(chromePlus) # chromosome models
library(diversitree) # basic likelihood functions
library(doMC) # this allows multicore runs on a mac

# load custom functions
source("functions.R")

# read in the data and the tree
trees <- read.nexus("../data/misof.backbone.nex")
dat <- read.csv("../data/data.invert.csv", as.is = T)[, -c(6:9)]

# set the MCMC chain length
iter <- 50

# we scale our trees so lets store that info for transforming
# rates back into millions of years
rate.scalers <- c()

#############################
#                           #
# TEST RUN TO GET W         #
#                           #
#############################

# slim the data to include only the desired data
# and generate the format table needed by chromPlus
foo <- getData(trees, dat)
trees.pruned <- foo[[1]]
chroms <- foo[[2]]
ncol(chroms)
rm(foo)
lk.mk <- make.mkn(trees.pruned[[1]], states = chroms,
                  k = ncol(chroms), strict = F,
                  control = list(method = "ode"))
# we will do this next bit twice to get w with and without polyploidy
# first with polyploidy
con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                        polyploidy = F, verbose = F,
                        constrain = list(drop.demi = T, drop.poly = T))
prior <- make.prior.exponential(.5)
temp.wop <- mcmc(con.lk.mk,
                 x.init = runif(6, 0, 1),
                 prior = prior,
                 w = 1,
                 nsteps = 20,
                 upper = 50,
                 lower = 0)
temp.wop <- temp.wop[-c(1:10), ]
w.wop <- diff(sapply(temp.wop[2:7],
                    quantile, c(.05, .95)))



# now with polyploidy
con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                        polyploidy = F, verbose = F,
                        constrain = list(drop.demi = T, drop.poly = F))
prior <- make.prior.exponential(.5)
temp.wp <- mcmc(con.lk.mk,
                 x.init = runif(8, 0, 1),
                 prior = prior,
                 w = 1,
                 nsteps = 20,
                 upper = 50,
                 lower = 0)
temp.wp <- temp.wp[-c(1:10), ]
w.wp <- diff(sapply(temp.wp[2:9],
                     quantile, c(.05, .95)))

#################################
#                               #
#  Now we can do our full run   #
#                               #
#################################
#iter <- 1 # just for testing
# this will allow us to run on 14 cores
registerDoMC(14)
result <- list()

# we will loop through all 100 trees
# fitting model without polyploidy
x <- foreach(i=1:100) %dopar%{
  # slim the data to include only the desired data
  # and generate the format table needed by chromPlus
  foo <- getData(trees, dat)
  ncol(foo[[2]])
  trees.pruned <- foo[[1]]
  chroms <- foo[[2]]
  rm(foo)
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees.pruned[[i]], states = chroms,
                    k = ncol(chroms), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 1),
                      prior = prior,
                      w = w.wop,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
  # just in case we have a crash lets write results for each tree
  write.csv(result[[i]], file=paste("/results/tree.nop",i,".csv", sep=""))
}
results <- list()
for(i in 1:100){
 scaler <- max(branching.times(trees.pruned[[i]]))
 x[,2:7] <- x[,2:7]/scaler
 results[[i]] <- x
}
save(results, file="rates.centromeres.wop")

# we will loop through all 100 trees
# fitting model with polyploidy
result <- list()
# we will loop through all 100 trees
x <- foreach(i=1:100) %dopar%{
  # slim the data to include only the desired data
  # and generate the format table needed by chromPlus
  foo <- getData(trees, dat)
  ncol(foo[[2]])
  trees.pruned <- foo[[1]]
  chroms <- foo[[2]]
  rm(foo)
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees.pruned[[i]], states = chroms,
                    k = ncol(chroms), strict = F,
                    control = list(method = "ode"))

  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = F))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(8, 0, 1),
                      prior = prior,
                      w = w.wp,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
  # just in case we have a crash lets write results for each tree
  #write.csv(result[[i]], file=paste("/results/tree.p",i,".csv", sep=""))
}
results <- list()
depths <- getData(trees, dat)[[3]]
for(i in 1:100){
  #scaler <- max(branching.times(trees.pruned[[i]]))
  x[[i]][,2:9] <- x[[i]][,2:9]/depths
  results[[i]] <- x[[i]]
}
save(results, file="rates.centromeres.wp")
# NOTE monocentric is state 2 and holocentric is state 1 in results



