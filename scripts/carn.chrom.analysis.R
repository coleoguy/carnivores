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
uncon3 <- c()
#loops through to identify runs that have a low probability and don't reach 
#convergence
for(i in 1:100){
  #looks at the mean rates of the two ascending rates from the model to see if 
  #they are less than 0
  if(mean(x3[[i]]$asc2[451:500])-mean(x3[[i]]$asc1[451:500]) < 0){
    #stores the current run that didn't meet convergence
    uncon3 <- c(uncon3, i)
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
  found <- con.lk.mk(pars = as.numeric(x1[[i]][500, 2:7]))
  #likelihood value for swapped rates
  tried <- con.lk.mk(pars = as.numeric(x1[[i]][500, c(4,5,2,3,6,7)]))
  
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  as.numeric(x1[[i]][500, c(4,5,2,3,6,7)]),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
  #print out the difference between the rate swapping
  print(tried - found)
}


