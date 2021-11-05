#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(parallel)
library(doSNOW)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance 
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method="extend")
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
#this loop samples a chromosome number for each species when there is more than 
#one
for(i in 1:nrow(range)){
  hit <- which(chroms$species == range$species[i])
  if(length(hit)>1)  hit <- sample(hit, 1)
    dat.pruned[i, 3] <- chroms[hit, 2]
}
#rm old data and clean up environment
rm(chroms, range, hit, i)

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% dat.pruned$species]
#empty list to store pruned trees
trees.pruned <- list()
#empty vector to store tree depths
#keep tree depths to correct for depth when analysing rates
tree.depths <- c()

#loop that drops missing tips and stores tree depths for further use
for(i in 1:100){
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  tree.depths[i] <- max(branching.times(cur.tree))
  cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
  trees.pruned[[i]] <- cur.tree
}

#write out the tree depth file for analysing results 
write.csv(tree.depths, "../results/25small75larg_tree_depths.csv")

#rm old data and clean up environment
rm(trees, cur.tree, i,  missing, tree.depths)

###DISCRETIZE RANGE SIZES###----------------------------------------------------

#look at a histagram of the data
hist(dat.pruned$range.size, breaks =200)
#adds a line to the histogram where the 75th quantile is
abline(v=quantile(dat.pruned$range.size, 0.25), col = "blue")
#stores as a variable the 75th quantile cutoff value
quant <- quantile(dat.pruned$range.size, 0.25)
#assigns the range values into two groups for the model.
#0 = small; 1 = large pop size
dat.pruned$range.size <- as.numeric(dat.pruned$range.size > quant)
#quant is the the threshold cutoff for discretization
#quant: ***6.7E7*** [Replace if modified]
#rm quant variable to clean up environment
rm(quant)

###VISUALIZE RANGE/CHROM DATA###------------------------------------------------

#visualize the range of the chromosome data
range(dat.pruned$hap.chrom) #16-39
#view a histogram of the haploid chromosome data
hist(dat.pruned$hap.chrom)

#view a histogram of the discretized range size data
hist(dat.pruned$range.size)

#view a plot with both the continuous and discrete data
boxplot(dat.pruned$hap.chrom ~ dat.pruned$range.size, 
        boxwex = 0.5,
        outpch = 16,
        outline = F, 
        xlab = "Discretized Range Size",
        ylab = "Chromosome Number",
        main="Chromosome Number by Discretized Range Size")
stripchart(dat.pruned$hap.chrom ~ dat.pruned$range.size, 
           vertical = TRUE, 
           method = "jitter", 
           add = TRUE, 
           pch = 20, 
           col = 'blue')

###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------

#convert chromosome number and binary trait state to format for diversitree
d.data <- as.data.frame(dat.pruned[,c(1,3,2)])

#store chrom.range for each tree 
chrom.range <- range(d.data$hap.chrom) + c(-2,2)

#run datatoMatrix function necessary for diversitree likelihood function
p.mat <- datatoMatrix(x=d.data, range=chrom.range, hyper=T)

#make a likelihood function
#states = named vector of character states
#k = number of states to model; 1 to k for make.mkn number of columns in pmat
#strict = T allows for missing states
#control (ode) = uses an ODE based approach to compute only the k variables over
#time; more efficient when k is large
lk.mk <- make.mkn(trees.pruned[[1]], states = p.mat,
                  k = ncol(p.mat), strict = F,
                  control = list(method = "ode"))

###CONSTRAIN LIKELIHOOD FUNCTION###---------------------------------------------

#constrain the likelihood function to remove states that are not biologically
#realistic
con.lik <- constrainMkn(data = p.mat,
                        lik = lk.mk,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly= T))

#check to make sure we have the parameters we are expecting
argnames(con.lik)

###FIT 1 TREE AS TEST###--------------------------------------------------------
prior <- make.prior.exponential(0.5)

#fit a biologically realistic model using diversitree's mcmc
#con.lik2 = likelihood function for mcmc to run on
#x.init = initial parameter location
#prior = an optional prior probability distribution
#w = tuning parameter for the  sampler
#nsteps = number of  mcmc steps to  take
#upper = upperr bounds on  parameter space
#lower = lower bounds on parameter space
temp <- mcmc(lik = con.lik, 
             x.init = runif(6,0,1), 
             prior = prior,
             w = 1, 
             nsteps = 2500, 
             upper = 50,
             lower = 0)

###EVALUATE RESULTS FROM 1 TREE###----------------------------------------------

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
tune <- temp[-c(1:550), ]
w <- diff(sapply(tune[2:7],
                 quantile, c(.05, .95)))

# this will allow us to run on more cores
cores <- detectCores(all.tests = T)
#take out  two cores for other basic functions
cores <- cores-2

#make a cluster
cl <- makeCluster(cores, outfile = " ")
#register the cluster
registerDoSNOW(cl)

#create empty list to store results
result <- list()
#set iter to 100  for the number of steps to take in the model
iter <- 2500


# we will loop through all 100 trees
# fitting model 

x <- foreach(i=1:100) %dopar%{
  library(phytools)
  library(chromePlus)
  library(diversitree)
  library(parallel)
  library(doSNOW)
  
  ###Sampling of chromosome dataset ###
  #load in tree data
  trees <- read.nexus("../data/carnivora.nex")
  for(j in 1:100){
    #trees are ultrametric, this line corrects for the fact that the tolerance 
    #for being ultrametric is not met by some trees
    trees[[j]] <- force.ultrametric(trees[[j]], method="extend")
  }
  #load in chromosome data
  chroms <- read.csv("../data/chroms.csv")
  #load in range size
  range <- read.csv("../data/calc.carn.range.sizes.csv")
  #change column names to be informative
  colnames(range) <- c("species", "range.size")
  #prune chromosome number and combnine with range size
  dat.pruned <- range
  #add empty third column for chromosome number
  dat.pruned[, 3] <- NA
  #name the third column
  colnames(dat.pruned)[3] <- "hap.chrom"
  #this loop samples a chromosome number for each species when there is more than 
  #one
  for(k in 1:nrow(range)){
    hit <- which(chroms$species == range$species[k])
    if(length(hit)>1)  hit <- sample(hit, 1)
    dat.pruned[k, 3] <- chroms[hit, 2]
  }
  #rm old data and clean up environment
  rm(chroms, range, hit, k, j)
  #prune and scale trees
  #find tips that are missing from the dataset
  missing <- trees[[i]]$tip.label[!trees[[i]]$tip.label %in% dat.pruned$species]
  #empty list to store pruned trees
  trees.pruned <- list()
  #empty vector to store tree depths
  #keep tree depths to correct for depth when analysing rates
  tree.depths <- c()
  #loop that drops missing tips and stores tree depths for further use
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  tree.depths[i] <- max(branching.times(cur.tree))
  cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
  trees.pruned[[i]] <- cur.tree
  #rm old data and clean up environment
  rm(trees, cur.tree,  missing)
  ###END CHROMOSOME SAMPLING###
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees.pruned[[i]], states = p.mat,
                    k = ncol(p.mat), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = p.mat, lik = lk.mk, hyper = T,
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
  # just in case we have a crash lets write results for each tree
  write.csv(result[[i]], file=paste("../result/tree.carn",i,".csv", sep=""))
}

#stop the cluster
stopCluster(cl)
