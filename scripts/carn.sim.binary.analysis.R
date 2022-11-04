#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

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


#assign prior from exponential distribution
prior <- make.prior.exponential(2)

# from primary analysis we can get our w
w <- c(0.6598417, 0.7470406, 4.915778, 4.435216, 2.065825, 2.17056)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 3)

#create empty list to store results
result <- list()
#set iter to 200  for the number of steps to take in the model
iter <- 200

x <- foreach(i = 1:100) %dopar%{
  
  
  # our true trait (range size) evolves at a rate of 2 ish
  # simulate a bunch of null traits and makes sure that the delta R we see is 
  # significant relative to a random trait.
  datalist[[i]]$range.size <-  sim.character(trees.pruned[[i]], 
                                             pars=c(2,2), 
                                             x0=sample(0:1,1), 
                                             model="mk2")
    #store chrom.range
    chrom.range <- range(datalist[[i]]$hap.chrom) + c(-1, 1)
    datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees.pruned[[i]], states = datalist[[i]],
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



##### Processing results #########
#only process post burn-in results
delta.r.sim <- as.data.frame(matrix(NA,100,4))
colnames(delta.r.sim) <-c("tree","fission","fusions", "mean")
for(i in 1:100){
  temp <- x[[i]]
  temp[,2:7] <- temp[,2:7]/tree.depths[i]
  temp <- temp[51:100,]
  delta.r.sim$fission[i] <- mean(temp$asc1 - temp$asc2)
  delta.r.sim$fusions[i] <- mean(temp$desc1 - temp$desc2)
  delta.r.sim$mean[i] <- mean(temp$asc1 + temp$desc1) - mean(temp$asc2 + temp$desc2)
  
}

obs.fis <- mean(carn_med_hb$asc1 - carn_med_hb$asc2)
obs.fus <- mean(carn_med_hb$desc1 - carn_med_hb$desc2)
obs.mean <- mean(carn_med_hb$asc1 + carn_med_hb$desc1) - mean(carn_med_hb$asc2 + carn_med_hb$desc2)

hist(abs(delta.r.sim$fission))
abline(v=0.164)
hist(abs(delta.r.sim$fusion))
abline(v=0.101)
hist(abs(delta.r.sim$mean))
abline(v=0.26)
#save the results output
write.csv(post.burn,file="../results/carn.delta.sim.hb.csv")


