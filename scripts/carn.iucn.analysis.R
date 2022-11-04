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

#load in iucn
iucn <- read.csv("../data/iucndata.csv")
#change column names to be informative

###PRUNE DATA###----------------------------------------------------------------

#prune chromosome number and combine with range size
dat.pruned <- iucn
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
  for(i in 1:nrow(dat.pruned)){
    hit <- which(chroms$species == dat.pruned$name[i])
    if(length(hit) > 1){
      hit <- sample(hit, 1)
    }
      dat.pruned$hap.chrom[i] <- chroms[hit, 2]
  }
  datalist[[j]] <- dat.pruned
}

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% datalist[[1]]$name]
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
rm(trees, cur.tree, i, j, missing, chroms, dat.pruned, iucn, hit)


#assign prior from exponential distribution
prior <- make.prior.exponential(4)

# from primary analysis we can get our w
w <- c(0.6598417, 0.7470406, 4.915778, 4.435216, 2.065825, 2.17056)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 3)

#create empty list to store results
result <- list()
#set iter to 200  for the number of steps to take in the model
iter <- 300

x <- foreach(i = 1:100) %dopar%{
    iuc <- rep(0, 110)
    iuc[datalist[[i]]$status != "Least Concern"] <- 1
    datalist[[i]]$status <- iuc
    # with this coding state 1 should rare animals and state 2 should be
    # common
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
                      x.init =  runif(6, 0, 1),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}



##### Processing results #########
#only process post burn-in results
iucn <- x[[1]][251:300,]
iucn[,2:7] <- iucn[,2:7]/tree.depths[1] 
for(i in 2:100){
  foo <- x[[i]][251:300,]
  foo[,2:7] <- foo[,2:7]/tree.depths[i] 
  iucn <- rbind(iucn, foo)
}
hist(iucn$asc1-iucn$asc2)
hist(iucn$desc1-iucn$desc2)

#save the results output
write.csv(iucn,file="../results/carn.iucn.csv")


