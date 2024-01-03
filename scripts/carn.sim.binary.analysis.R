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
w <-  c(14.44759, 12.13979, 16.06275, 14.23117, 2.477789, 3.194005)
# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 35)

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

for(i in 1:100){
  temp <- x[i]
  write.csv(temp, paste0("../results/trait_sim/trait", i, ".csv"))
}

sig.asc <- sig.desc <- c()
mn.asc <- mn.desc <- c()
for(i in 1:100){
  curdat <- x[[i]]
  dr.asc <- dr.desc <- c()
  dr.asc <- curdat$asc1 - curdat$asc2
  dr.desc <- curdat$desc1 - curdat$desc2
  mn.asc[i] <- mean(dr.asc)
  mn.desc[i] <- mean(dr.desc)
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


mn.mn.rates <- c()
for(i in 1:100){
  mn.mn.rates[i] <- abs(mn.asc[i] + mn.desc[i]) / 2
}
emp <- (.101+.163)/2
sum(mn.mn.rates >= emp) / 100
