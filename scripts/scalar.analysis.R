#### PACKAGES ####
library(phytools)
library(doSNOW)
library(ape)
library(chromePlus)
library(diversitree)
library(expm)
library(castor)
library(viridis)
library(ggplot2)
source("functions.R")

#### LOAD DATA ####
load("../data/range_size/datalists_rangesize.RData")

fam.dat <- read.csv("../data/fam_data.csv",
                    as.is=T)[,-1]

trees <- read.nexus("../data/range_size/carnivora_rs_pruned.nex")

mat <- as.matrix(read.csv("../data/Q_matrix.csv"))[,-1]

#### CHANGE RED PANDA and LINSANG FAMILY ####

fam.dat$Family[which(fam.dat$species == "Ailurus_fulgens")] <- "Ailuridae"
fam.dat$Family[which(fam.dat$species == "Prionodon_linsang")] <- "Prionodontidae"

#### GET FAMILY DATA ####

#convert to factor
fam.dat$Family <- as.factor(fam.dat$Family)

#get counts
fam.count <- summary(fam.dat$Family)
fam.interest <- names(which(fam.count >= 5))

#### MATRIX MANIPULATION ####

#get matrices to numeric
mat <- matrix(as.numeric(mat),ncol=ncol(mat),nrow=nrow(mat))
mat <- mat[1:26,1:26]
#### OBJECTS TO STORE ####
#df for scalar means
scalar.means <- as.data.frame(matrix(NA, nrow=100,ncol=length(fam.interest)))
colnames(scalar.means) <- fam.interest

clade.edge.lengths <- as.data.frame(matrix(NA, nrow=100,ncol=length(fam.interest)))
colnames(clade.edge.lengths) <- fam.interest

scalar.means.multiplied <- as.data.frame(matrix(NA, nrow=100,ncol=length(fam.interest)))
colnames(scalar.means.multiplied) <- fam.interest

scalars.above <- as.data.frame(matrix(NA, nrow=100,ncol=length(fam.interest)))
colnames(scalars.above) <- fam.interest

scalars.below <- as.data.frame(matrix(NA, nrow=100,ncol=length(fam.interest)))
colnames(scalars.below) <- fam.interest

#list to store
#scaled.trees <- list()
load("../results/scaled.trees.RData")

#### SET UP PARALLELIZATION ####

# Define number of clusters
nClust <- 100

# Set up clusters, print intermediates to 
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

#### LOOP THROUGH POSTERIOR ####

scaled.trees <- foreach(i=1:100,
                   .verbose = T,
                   .packages = c("phytools","chromePlus","ape",
                                 "diversitree","expm")) %dopar% {
  
  print(paste0("TREE ", i))
  
  #### GET DATA FOR cURRENT POSTERIOR ####
  dat <- datalist[[i]]
  tree <- trees[[i]]
  
  #### ESTIMATE RATES ####
  
  #make datamatrix
  dat.matrix <- datatoMatrix(dat,
                              c(15,40),
                              hyper = F)
  
  #make model
  model <- make.mkn(tree,
                    dat.matrix,
                    ncol(dat.matrix),
                    strict=F,
                    control=list(method="ode"))
  
  #constrain model
  model.con <- constrainMkn(dat.matrix,
                            model,
                            hyper = F,
                            polyploidy = F,
                            verbose = T,
                            constrain = list(drop.poly=T,
                                             drop.demi=T))
  
  argnames(model.con$`likelihood function`)
  
  #Check parMat
  parMat <- model.con$`parameter matrix`
  
  #Run mcmc
  print("ESTIMATING RATES")
  
  model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                  x.init=runif(2,0,1),
                                  prior=make.prior.exponential(r=2),
                                  upper=c(50,50),
                                  lower=0,
                                  nsteps = 500,
                                  w=1)
  
  #### BUILD QMATRIX ####
  
  #Extract post burn portion
  model.mcmc.postburn <- model.mcmc[450:500,]
  
  #Get mean params
  params <- c(mean(model.mcmc.postburn$asc1),
              mean(model.mcmc.postburn$desc1))
  
  #sub into matrix
  print("BUILDING Q-MATRIX")
  
  parMat[parMat == "asc1"] <- params[1]
  parMat[parMat == "desc1"] <- params[2]
  
  Qmat <- matrix(as.numeric(parMat),ncol=ncol(parMat),nrow=nrow(parMat))
  diag(Qmat) <- -rowSums(Qmat)
  
  #### DATA MANIPULATION ####
  
  #get simulation states for taxa
  for(j in 1:nrow(dat)){
    
    dat$sim.state[j] <- dat$hap.chrom[j] -  14
    
  }
  
  #Get tip states
  tip.states <- dat$sim.state
  names(tip.states) <- dat$species

  #### RUN SCALAR ANALYSIS ####
  
  #full model
  result <- scaleTreeRates(tree=tree,
                                tip.states=tip.states,
                                model=mat,
                                fixedQ = Qmat,
                                max.ratio = 2,
                           nbins = 10)
                                 }

#Close cluster connection
stopCluster(cl)
rm(cl)

#### LOOP THROUGH FAMILIES OF INTEREST ####

for(i in 1:100){
  for(j in 1:length(fam.interest)){
    
    #assign tree
    tree <- scaled.trees[[i]]
    
    #get species in family
    species <- fam.dat$species[which(fam.dat$Family == fam.interest[j])]
    
    #get tip indices
    tip.indices <- which(tree$tip.label %in% species == TRUE)
    
    #get lower and upper bounds of clade
    lower.tip <- tree$edge[min(which(tree$edge[,2] %in% tip.indices)),2]
    upper.tip <- tree$edge[max(which(tree$edge[,2] %in% tip.indices)),2]
    
    #last common ancestral node and edge leading to node
    anc.node <- get_pairwise_mrcas(tree,lower.tip,upper.tip)
    anc.edge <- which(tree$edge[,2] == anc.node)
    
    #get edges associated with clade
    clade.edges <- anc.edge:max(which(tree$edge[,2] %in% tip.indices))
    
    #get total edge length in clade
    clade.edge.length <- sum(tree$edge.length[clade.edges])
    
    #get mean of clade scalars
    scalars.above[i,j] <- length(which(tree$scalar[clade.edges] >= 1))/length(clade.edges)
    scalars.below[i,j] <- length(which(tree$scalar[clade.edges] < 1))/length(clade.edges)
    scalar.means[i,j] <- mean(scaled.trees[[i]]$scalar[clade.edges])
    clade.edge.lengths[i,j] <- clade.edge.length
    scalar.means.multiplied[i,j] <- mean(scaled.trees[[i]]$scalar[clade.edges]) * clade.edge.length
  }
}

#combine data
scalars.summarized <- as.data.frame(matrix(NA,ncol = 3,
                                           nrow=14))
scalars.summarized[,1] <- fam.interest
scalars.summarized[,2] <- c(rep("Depressed",7),rep("Elevated",7))
colnames(scalars.summarized) <- c("family","scalar.state","proportion")

for(i in 1:nrow(scalars.summarized)){
  
  if(i <= 7){
    scalars.summarized[i,3] <- mean(scalars.below[,i])
  } else {
    
    scalars.summarized[i,3] <- mean(scalars.above[,i-7])
  }
  
}

#### SAVE OUTPUTS ####
save(scaled.trees,file="../results/scaled_trees.RData")
save(list = c("scalar.means",
              "clade.edge.lengths",
              "scalar.means.multiplied",
              "scalars.summarized"),file = "../results/summarized_scalar_outputs.RData")





