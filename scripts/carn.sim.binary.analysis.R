# Heath Blackmon and Michelle Jonika
# December 27, 2022
# Script that simulates neutrally evolving binary traits under chromePlus
# model to estimate a false positive rate

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(coda)
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/range_size/carnivora_rs_pruned.nex")

#load in binary trait data
load("../data/range_size/datalists_rangesize.RData")

#load in tree depths
tree.depths <- read.csv("../data/range_size/rs_treedepths.csv")
colnames(tree.depths) <- c("tree", "tree.depth")

#assign prior from exponential distribution
prior <- make.prior.exponential(2)

# from primary analysis we can get our w
w <- c(14.44759, 12.13979, 16.06275, 14.23117, 2.477789, 3.194005)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 1)

#set iter to 200  for the number of steps to take in the model
iter <- 200
pbrows <- 151:200
treesamp <- sample(1:100)
x <- foreach(j = 1:100) %dopar%{
  # our true trait (range size) evolves at a rate of 2 ish
  # simulate a bunch of null traits and makes sure that the delta R we see is 
  # significant relative to a random trait.
  sim.bin <- sim.character(trees[[treesamp[j]]], pars = c(2, 2),
                           x0 = sample(0:1, 1), model = "mk2")
  for(k in 1:100){
    datalist[[k]]$range.size <- sim.bin
  }
  for(i in 1:100){
    chrom.range <- range(datalist[[i]]$hap.chrom) + c(-1, 1)
    datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
    lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                      k = ncol(datalist[[i]]), strict = F,
                      control = list(method = "ode"))
    con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                            polyploidy = F, verbose = F,
                            constrain = list(drop.demi = T, drop.poly = T))
    result <- mcmc(con.lk.mk, x.init =  runif(6, 0, 10),
                   prior = prior, w = w, nsteps = iter,
                   upper = 50, lower = 0)
    result[,2:7] <- result[,2:7]/tree.depths$tree.depth[i]
    if(i == 1){
      postburn <- result[pbrows,]
    }else{
      postburn <- rbind(postburn, result[pbrows,])
    }
  }
  postburn
}
#
sig.asc <- sig.desc <- c()
for(i in 1:100){
  curdat <- x[[i]]
  dr.asc <- dr.desc <- c()
  dr.asc <- curdat$asc1 - curdat$asc2
  dr.desc <- curdat$desc1 - curdat$desc2
  cint <- HPDinterval(as.mcmc(dr.asc))[1:2]
  if(cint[1] < 0 & cint[2] > 0){
    sig.asc[i] <- F
  }else{
    sig.asc[i] <- T
  }
  cint <- HPDinterval(as.mcmc(dr.desc))[1:2]
  if(cint[1] < 0 & cint[2] > 0){
    sig.desc[i] <- F
  }else{
    sig.desc[i] <- T
  }
}
sum(sig.desc)/100
sum(sig.asc)/100
#
#
#
#

  #store chrom.range
  # make the basic likelihood function for the data
 # now we constrain our model to be biologically realistic for
  # chromosomes.
  # now we are ready to run our inference run
}


##### Processing results #########
#only process post burn-in results
delta.r.sim <- as.data.frame(matrix(NA,100,4))
colnames(delta.r.sim) <-c("tree","fission","fusions", "mean")
for(i in 1:100){
  temp <- x[[i]]
  temp[,2:7] <- temp[,2:7]/tree.depths$tree.depth[i]
  temp <- temp[451:500,]
  delta.r.sim$fission[i] <- mean(temp$asc2 - temp$asc1)
  delta.r.sim$fusions[i] <- mean(temp$desc2 - temp$desc1)

}
delta.r.sim$tree <- c(1:100)

obs <- read.csv("../results/range_size/rs.csv")
obs.fis <- mean(obs$asc2 - obs$asc1)
obs.fus <- mean(obs$desc2 - obs$desc1)

hist(abs(delta.r.sim$fission))
abline(v=0.163)
hist(abs(delta.r.sim$fusions))
abline(v=0.101)
hist(abs(delta.r.sim$mean))
abline(v=0.264)
#save the results output
write.csv(post.burn,file="../results/carn_delta_sim.csv")


