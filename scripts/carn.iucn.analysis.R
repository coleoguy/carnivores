#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/datalists_iucn.RData")

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in tree depths
tree.depths <- read.csv("../data/treedepths.csv")

#assign prior from exponential distribution
prior <- make.prior.exponential(2)

# from primary analysis we can get our w
w <- c(0.6598417, 0.7470406, 4.915778, 4.435216, 2.065825, 2.17056)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 3)

#create empty list to store results
result <- list()
#set iter to 500 for the number of steps to take in the model
iter <- 500

x <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
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
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}

##### Processing results #########
#only process post burn-in results
iucn <- x[[1]][451:500,]
iucn[,2:7] <- iucn[,2:7]/tree.depths[1,2] 
for(i in 2:100){
  foo <- x[[i]][451:500,]
  foo[,2:7] <- foo[,2:7]/tree.depths[i,2] 
  iucn <- rbind(iucn, foo)
}
hist(iucn$asc1-iucn$asc2)
hist(iucn$desc1-iucn$desc2)

#save the results output
write.csv(iucn,file="../results/carn.iucn.csv")


