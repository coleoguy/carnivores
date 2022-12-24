# Michelle Jonika
# 10 November
#This code munges the data from the each of the ChromPlus replicates to gather 
#the final posterior dataset for visualization

#load in chromosome number and binary trait data
load("../data/datalists_iucn.RData")

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#assign prior from exponential distribution
prior <- make.prdatalist[[i]]ior.exponential(2)

# from primary analysis we can get our w
w <- c(14.44759, 12.13979, 16.06275, 14.23117, 2.477789, 3.194005)

#set iter to 500 for the number of steps to take in the model
iter <- 500

#load in each of the replicates
load("trial1_iucn.Rdata")
load("trial2_iucn.Rdata")
load("trial3_iucn.Rdata")
load("trial4_iucn.Rdata")
load("trial5_iucn.Rdata")
load("trial6_iucn.Rdata")
load("trial7_iucn.Rdata")
load("trial8_iucn.Rdata")
load("trial9_iucn.Rdata")

#create an empty 
results <- list()

#make a list of the replicates
trials <- list(x1,x2,x3,x4,x5,x6,x7,x8,x9)

#loop that retains the MCMC run with the highest final posterior probability
for(i in 1:length(x1)){
  #loops through each tree and stores the final posterior probability
  liks <- c(trials[[1]][[i]]$p[500],
            trials[[2]][[i]]$p[500],
            trials[[3]][[i]]$p[500],
            trials[[4]][[i]]$p[500],
            trials[[5]][[i]]$p[500],
            trials[[6]][[i]]$p[500],
            trials[[7]][[i]]$p[500],
            trials[[8]][[i]]$p[500],
            trials[[9]][[i]]$p[500])
  #store the tree with the highest posterior probabiliyt for the ifinal dataset
  results[[i]] <- trials[[which.max(liks)]][[i]]
}


#vector to store the likelihood differences
rate_swap <- c()
#loop that swaps rates between IUCN binary categories to see if that impacts 
#convergence
for(i in 1:length(results)){
  #store chrom.range
  chrom.range <- range(datalist[[i]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  #make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  
  #constrain our model to be biologically realistic for chromosomes
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  
  #likelihood value for the original rates
  found <- con.lk.mk(pars = as.numeric(results[[3]][500, 2:7]))
  
  #likelihood value or swapped rates
  tried <- con.lk.mk(pars = as.numeric(results[[3]][500, c(4,5,2,3,6,7)]))
  
  #print out the result of the rate swapping
  print(tried-found)
  
  #store the difference between the likelihood ratio of those found in the model
  #versus those tried when you swap rates between the two states
  rate_swap[i] <- tried-found
}

#check for trees that are potentially stuck in a local optimum for the MCMC model
local_op <- which(rate_swap > 0)

#vector to store new results for trees that need to be run again
results.uncon <- c()
#load in the datalists for running the MCMC models
load("../data/datalists_iucn.RData")
#loop that reruns MCMC model on those trees that are stuck in local optimums
for(i in 1:length(local_op)){
  j <- local_op[i]
  #store chrom.range
  chrom.range <- range(datalist[[j]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[j]] <- datatoMatrix(x=datalist[[j]], range = chrom.range, hyper = T)
  
  #make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[j]], states = datalist[[j]],
                    k = ncol(datalist[[j]]), strict = F,
                    control = list(method = "ode"))
  
  #constrain our model to be biologically realistic for chromosomes
  con.lk.mk<-constrainMkn(data = datalist[[j]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  
  # now we are ready to run our inference run
  results.uncon[[i]] <- mcmc(con.lk.mk,
                              x.init =  runif(6, 0, 10),
                              prior = prior,
                              w = w,
                              nsteps = iter,
                              upper = 50,
                              lower = 0)
}

#vector to store the likelihood differences
rate_swap2 <- c()
#loop that swaps rates between IUCN binary categories to see if that impacts 
#convergence
for(i in 1:length(results.uncon)){
  #store chrom.range
  chrom.range <- range(datalist[[i]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  #make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  #store the difference between the likelihood ratio of those found in the model
  #versus those tried when you swap rates between the two states
  #constrain our model to be biologically realistic for chromosomes
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  
  #likelihood value for the original rates
  found <- con.lk.mk(pars = as.numeric(results.uncon[[3]][500, 2:7]))
  
  #likelihood value or swapped rates
  tried <- con.lk.mk(pars = as.numeric(results.uncon[[3]][500, c(4,5,2,3,6,7)]))
  
  #print out the result of the rate swapping
  print(tried-found)
  
  #store the difference between the likelihood ratio of those found in the model
  #versus those tried when you swap rates between the two states
  rate_swap2[i] <- tried-found
}

#check for trees that are potentially stuck in a local optimum for the MCMC model
local_op2 <- which(rate_swap2 > 0)

#save the results after rerunning trees that did not meet convergence
save(results.uncon, file="../results/results_uncon.Rdata")

#only keep the results from teh original analysis that were not stuck in a local
#optima
results_orig <- results[-c(local_op)]

#combine the results from the original analysis with the results from the rerun
#analysis for those trees that did not meet convergence
results_total <- c(results_orig, results.uncon)

#store the names of the trees that are contained within each list

#create a vector that contains 1 to 100
labels <- c(1:100)
#remove the trees that needed to be rerun
labels <- names[-c(local_op)]
#add back in the trees that needed to be rerun in the correct order
labels <- c(names, local_op)

#store the names of each tree that corresponds to each list in the correct order 
#in the total results
names(results_total) <- labels

#save the results after rerunning trees that did not meet convergence
save(results_total, file="../results/results_total.Rdata")

#load in the tree depths to transform back into units of millions of years (MY)
tree.depths <- read.csv("../data/treedepths.csv")

#reorder the tree depths to correspond to the order they are in the lists
tree_depths_reorder <- tree.depths[labels,]

#only process post burn-in results
post.burn <- results_total[[1]][451:500, 2:8]

#transform the post burn in results back into MY from the tree depths
post.burn[,1:6] <- post.burn[,1:6]/tree_depths_reorder[1,2]

#loop that organizes rates into a pretty table
for(i in 2:100){
  #pulls results for specific tree
  temp <- results_total[[i]]
  #transforms the rates back by their tree depths
  temp[,2:7] <- temp[,2:7]/tree_depths_reorder[i,2]
  ##bind in post-burn-in samples from each tree after they have been back 
  #transformed
  post.burn <- rbind(post.burn, temp[451:500,2:8])
}

#save the results output
write.csv(post.burn, file="../results/iucn.csv")
