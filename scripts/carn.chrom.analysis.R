#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/datalists_range.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in tree depths
tree.depths <- read.csv("../data/treedepths.csv")
colnames(tree.depths) <- c("tree", "tree.depth")


###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------
#store chrom.range
chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)

#run datatoMatrix function necessary for diversitree likelihood function
for(i in 1:100){
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
}

#make a likelihood function
#states = named vector of character states
#k = number of states to model; 1 to k for make.mkn number of columns in pmat
#strict = T allows for missing states
#control (ode) = uses an ODE based approach to compute only the k variables over
#time; more efficient when k is large
lk.mk <- make.mkn(trees[[52]], states = datalist[[52]],
                  k = ncol(datalist[[52]]), strict = F,
                  control = list(method = "ode"))

###CONSTRAIN LIKELIHOOD FUNCTION###---------------------------------------------

#constrain the likelihood function to remove states that are not biologically
#realistic
con.lik <- constrainMkn(data = datalist[[52]],
                        lik = lk.mk,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly= T))

#check to make sure we have the parameters we are expecting
argnames(con.lik)
rm(chrom.range, i, lk.mk)
######### FIT 1 TREE AS TEST##########
#assign prior from exponential distribution
prior <- make.prior.exponential(2)

#fit a biologically realistic model using diversitree's mcmc
#con.lik2 = likelihood function for mcmc to run on
#x.init = initial parameter location
#prior = an optional prior probability distribution
#w = tuning parameter for the  sampler
#nsteps = number of  mcmc steps to  take
#upper = upperr bounds on  parameter space
#lower = lower bounds on parameter space
temp <- mcmc(lik = con.lik,
             x.init = runif(6, 0, 1),
             prior = prior,
             w = 10,
             nsteps = 500,
             upper = 20,
             lower = 0)

###FULL PARALLEL RUN OF TREES###------------------------------------------------

#tuning parameters for w for full run
tune <- temp[-c(1:400), ]
w <- diff(sapply(tune[2:7],
                 quantile, c(.05, .95)))

# w <- c(14.44759, 12.13979, 16.06275, 14.23117, 2.477789, 3.194005)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 25)

#create empty list to store results
result <- list()
#set iter to 500 for the number of steps to take in the model
iter <- 500


# we will loop through all 100 trees
# fitting model
x1 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x1, file="trial1.Rdata")

x2 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x2, file="trial2.Rdata")

x3 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x3, file="trial3.Rdata")

x4 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x4, file="trial4.Rdata")

x5 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x5, file="trial5.Rdata")

x6 <- foreach(i = 1:100) %dopar%{
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x6, file="trial6.Rdata")

##### Checking for convergence ###########
# After checking runs I found that some runs stayend in a low prob
# area with high rates in clades that have high pop sizes. To see if these
# really are global optimum I first identified which runs these were

#assign variable to store those runs that don't converge
uncon <- c()
#loops through to identify runs that have a low probability and don't reach 
#convergence
for(i in 1:100){
  #looks at the mean rates of the two ascending rates from the model to see if 
  #they are less than 0
  if(mean(x[[i]]$asc2[200:300])-mean(x[[i]]$asc1[200:300]) < 0){
    #stores the current run that didn't meet convergence
    uncon <- c(uncon, i)
    #runs that didn't meet convergence (29)
    #4,14,16,19,20,27,28,29,30,34,38,41,42,49,50,53,55,58,61,66,67,69,75,76,79,81,88,96,100
  }
}
#loop that swaps rates between high and low pop to see if that impacts convergence
for(i in uncon){
  #make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  #constrain our model to be biologically realistic for chromosomes
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  #likelihood value for original rates
  found <- con.lk.mk(pars = as.numeric(x[[i]][300, 2:7]))
  #likelihood value for swapped rates
  tried <- con.lk.mk(pars = as.numeric(x[[i]][300, c(4,5,2,3,6,7)]))
  #print out the difference between the rate swapping
  print(tried - found)
  #1: 10.65769
  #2: 11.85909
  #3: 12.84475
  #4: 10.31322
  #5: 11.93795
  #6: 16.56795
  #7: 20.91868
  #8: 13.26349
  #9: 9.032648
  #10: 9.601858
  #11: 17.70944
  #12: 7.862768
  #13: 13.93501
  #14: 23.22894
  #15: 8.381647
  #16: 14.04385
  #17: 20.81376
  #18: 16.70518
  #19: 14.529
  #20: 13.19418
  #21: 16.0233
  #22: 13.11986
  #23: 11.78891
  #24: 13.78249
  #25: 13.07813
  #26: 13.55449
  #27: 14.48879
  #28: -10.61469
  #29: 11.03922
}

##### Processing results #########
#only process post burn-in results
post.burn <- x[[1]][451:500, 2:8]
#transform the post burn in results back into MY from the tree depths
post.burn[,1:6] <- post.burn[,1:6]/tree.depths[1,2]

#loop that organizes rates into a pretty table
for(i in 2:100){
  #pulls results for specific tree
  temp <- x[[i]]
  #transforms the rates back by their tree depths
  temp[,2:7] <- temp[,2:7]/tree.depths[i,2]
  ##bind in post-burn-in samples from each ttree after they have been back 
  #transformed
  post.burn <- rbind(post.burn, temp[451:500,2:8])
}

#save the results output
write.csv(post.burn,file="../results/rangesize.csv")


