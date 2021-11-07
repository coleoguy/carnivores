#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
# I switched to doMC it is a lot easier on a mac it takes environment with it
# library(parallel)
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


#this loop samples a chromosome number for each species when there is more than
#one
# TODO here I changed the code to create 100 datasets
# i did this because I didn't want to mess with reading
# and editing the list in each of the parallel runs it will
# lets us simplify the code down in the parallel run part
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

# TODO previously you had written out the tree depth file
# for analysing results a better approach is to use them
# in this script at the end and convert rates back to MY
# units prior to saving results

#rm old data and clean up environment
rm(trees, cur.tree, i, j, missing, chroms, dat.pruned, range, hit)

###DISCRETIZE RANGE SIZES###----------------------------------------------------
# TODO after a lot of playing I think you can simplify and just say
# cut it at the median has less of a feel of cherry picking
for(i in 1:100){
  x <- median(datalist[[1]]$range.size)
  datalist[[i]]$range.size <- as.numeric(datalist[[1]]$range.size >= x)
}
#0 = small; 1 = large pop size
rm(i, x)

# # look at a histagram of the data
# hist(dat.pruned$range.size, breaks =200)
# # adds a line to the histogram where the 75th quantile is
# abline(v=quantile(dat.pruned$range.size, 0.25), col = "blue")
# # stores as a variable the 75th quantile cutoff value
# quant <- quantile(dat.pruned$range.size, 0.25)
#
# #assigns the range values into two groups for the model.
#
# dat.pruned$range.size <- as.numeric(dat.pruned$range.size > quant)
# dat.pruned$range.size <- as.numeric(dat.pruned$range.size >= median(dat.pruned$range.size))
# # quant is the the threshold cutoff for discretization
# # quant: ***6.7E7*** [Replace if modified]
# # rm quant variable to clean up environment


###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------

# TODO already put in correct order above
# convert chromosome number and binary trait state to format for diversitree
# d.data <- as.data.frame(dat.pruned[,c(1,3,2)])

#store chrom.range
chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)

#run datatoMatrix function necessary for diversitree likelihood function
# TODO here I replicate what you had but for all datasets
for(i in 1:100){
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)

}

#make a likelihood function
#states = named vector of character states
#k = number of states to model; 1 to k for make.mkn number of columns in pmat
#strict = T allows for missing states
#control (ode) = uses an ODE based approach to compute only the k variables over
#time; more efficient when k is large
lk.mk <- make.mkn(trees.pruned[[52]], states = datalist[[52]],
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
# TODO I spent a few hours playing with the
# prior on many different trees and basically
# what I found was that when the prior is low
# it can wander around a lot with high rates in
# large pops before eventually finding a better
# parameter space with high rates associated with
# small pops so I went with a relatively high prior
# as this tended to get us honed in on this high
# probability parameter space much more quickly.
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
             nsteps = 300,
             upper = 20,
             lower = 0)
######### FIT 1 TREE AS TEST##########

###EVALUATE RESULTS FROM 1 TREE#######
# TODO lets get rid of this eventually/now
#evaluate convergence
plot(temp$p)
#evaluate asc1
plot(temp$asc1)
#evaluate desc1
plot(temp$desc1)
#evaluate asc2
plot(temp$asc2)
#evaluate desc2
plot(temp$desc2)
#evaluate tran12
plot(temp$tran12)
#evaluate tran21
plot(temp$tran21)

###FULL PARALLEL RUN OF TREES###------------------------------------------------

#tuning parameters for w for full run
tune <- temp[-c(1:100), ]
w <- diff(sapply(tune[2:7],
                 quantile, c(.05, .95)))


# TODO I switch doMC because it takes the environment with it into
# the parallel processing so greatly eases the coding
# this is the number of cores to use
registerDoMC(detectCores(all.tests = T) - 3)

#create empty list to store results
result <- list()
#set iter to 100  for the number of steps to take in the model
iter <- 300


# we will loop through all 100 trees
# fitting model

x <- foreach(i = 1:100) %dopar%{
  # library(phytools)
  # library(chromePlus)
  # library(diversitree)
  # library(parallel)
  # library(doSNOW)
  # TODO becuase I made two lists ready to go (trees and data)
  # we can get rid of a bunch of stuff in here
  ###Sampling of chromosome dataset ###
  #load in tree data
  #trees <- read.nexus("../data/carnivora.nex")
  #for(j in 1:100){
    #trees are ultrametric, this line corrects for the fact that the tolerance
    #for being ultrametric is not met by some trees
  #  trees[[j]] <- force.ultrametric(trees[[j]], method="extend")
  #}
  #load in chromosome data
  #chroms <- read.csv("../data/chroms.csv")
  #this loop samples a chromosome number for each species when there is more than
  #one
  #for(k in 1:nrow(range)){
  #  hit <- which(chroms$species == range$species[k])
  #  if(length(hit)>1)  hit <- sample(hit, 1)
  #  dat.pruned[k, 3] <- chroms[hit, 2]
  #}
  #rm old data and clean up environment
  #rm(chroms, range, hit, k, j)
  #prune and scale trees
  ###END CHROMOSOME SAMPLING###

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
  # just in case we have a crash lets write results for each tree
  # TODO lets not I think I fixed it so we just dont have any crashes that seems better
  # write.csv(result[[i]], file=paste("../result/tree.carn",i,".csv", sep=""))
}
##### Checking for convergence ###########
# TODO after checking runs I found that some runs stayend in a low prob
# area with high rates in clades that have high pop sizes. To see if these
# really are global optimum I first identified which runs these were
uncon <- c()
for(i in 1:100){
  if(mean(x[[i]]$asc2[200:300])-mean(x[[i]]$asc1[200:300]) < 0){
    uncon <- c(uncon, i)
  }

}
# TODO Then I just tried swapping rates between high and
# low pop to see what happened
for(i in uncon){
  lk.mk <- make.mkn(trees.pruned[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  found <- con.lk.mk(pars = as.numeric(x[[i]][300, 2:7]))
  tried <- con.lk.mk(pars = as.numeric(x[[i]][300, c(4,5,2,3,6,7)]))
  print(tried - found)
}


##### Checking for convergence ###########

##### Processing results #########
post.burn <- x[[1]][201:300, 2:8]
post.burn[,1:6] <- post.burn[,1:6]/tree.depths[1]

for(i in 2:100){
  temp <- x[[i]]
  temp[,2:7] <- temp[,2:7]/tree.depths[i]
  post.burn <- rbind(post.burn, temp[201:300,2:8])
}
# now save this
write.csv(post.burn,file="../results/carn.med.hb.csv")
##### Processing results #########

