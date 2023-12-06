#rerunning chromplus analysis with canidae removed

#load in tree data
trees <- read.nexus("../../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method = "extend")
}

tree <- trees[[1]]
#load in chromosome data
chroms <- read.csv("../../data/chroms.csv")

#load in range size
range <- read.csv("../../data/range_size/range_size.csv")
#change column names to be informative
colnames(range) <- c("species", "range.size")

###PRUNE DATA###----------------------------------------------------------------

drop.species <- c("Canis_latrans",
                  "Canis_lupus",
                  "Chrysocyon_brachyurus",
                  "Cuon_alpinus",
                  "Lycaon_pictus",
                  "Nyctereutes_procyonoides",
                  "Otocyon_megalotis",
                  "Speothos_venaticus",
                  "Urocyon_cinereoargenteus",
                  "Vulpes_corsac",
                  "Vulpes_macrotis",
                  "Vulpes_vulpes")

store.drop <- c()
for(i in 1:length(drop.species)){
  store.drop[i] <- which(range$species == drop.species[i])
}

range_new <- range[-c(store.drop),]


#prune chromosome number and combnine with range size
dat.pruned <- range_new
#add empty third column for chromosome number
dat.pruned[, 3]  <- NA
#name the third column
colnames(dat.pruned)[3] <- "hap.chrom"
# Columns are out of order lets fix them immediately
dat.pruned <- dat.pruned[, c(1, 3, 2)]


#this loop creates 100 datasets, sampling a chromosome number for each species 
#when there is more than one
for(i in 1:nrow(range_new)){
  hit <- which(chroms$species == range_new$species[i])
  if(length(hit) > 1){
    hit <- sample(hit, 1)
  }
  dat.pruned$hap.chrom[i] <- chroms[hit, 2]
}


#prune and scale trees
#find tips that are missing from the dataset
missing <- tree$tip.label[!tree$tip.label %in% dat.pruned$species]

#drops missing tips and stores tree depths for further use

cur.tree <- drop.tip(tree, tip = missing)
tree.depth <- max(branching.times(cur.tree))
cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
tree.pruned <- cur.tree

rm(chroms, cur.tree, range, range_new, tree, trees, drop.species, i, missing, 
   store.drop, hit)

###DISCRETIZE RANGE SIZES###----------------------------------------------------
#discretize range size based on the median
x <- median(dat.pruned$range.size)
dat.pruned$range.size <- as.numeric(dat.pruned$range.size >= x)
#0 = small; 1 = large pop size


###ORDER TREE AND DATA###-------------------------------------------------------

#empty vector to store correct order in 
neworder <- c()
#loop through to store the correct order of the data to match the tree tips
for(j in 1:length(tree.pruned$tip.label)){
  neworder[j] <- which(tree.pruned$tip.label[j] == dat.pruned$species)
}
#reorder the data into a new data frame
dat.pruned <- dat.pruned[neworder,]

rm(j, neworder, x)

###MAKE LIKELIHOOD FUNCTION###--------------------------------------------------
#store chrom.range
chrom.range <- range(dat.pruned$hap.chrom) + c(-1, 1)

#run datatoMatrix function necessary for diversitree likelihood function
datalist <- datatoMatrix(x=dat.pruned, range = chrom.range, hyper = T)


#make a likelihood function
#states = named vector of character states
#k = number of states to model; 1 to k for make.mkn number of columns in pmat
#strict = T allows for missing states
#control (ode) = uses an ODE based approach to compute only the k variables over
#time; more efficient when k is large
lk.mk <- make.mkn(tree.pruned, states = datalist,
                  k = ncol(datalist), strict = F,
                  control = list(method = "ode"))

###CONSTRAIN LIKELIHOOD FUNCTION###---------------------------------------------

#constrain the likelihood function to remove states that are not biologically
#realistic
con.lik <- constrainMkn(data = datalist,
                        lik = lk.mk,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly= T))

#check to make sure we have the parameters we are expecting
argnames(con.lik)
rm(chrom.range, lk.mk)
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

write.csv(temp, "../../results/range_size/rs_wocanid.csv")


#read in chromplus data for range size
fission <- temp$asc2 - temp$asc1
fusion <- temp$desc2 - temp$desc1
data_munge <- data.frame(fission, fusion)

plotChromeplus(data = data_munge,
               colors = c("#FDE725FF", "#39568CFF"),
               main_title = "",
               x_title = "rate difference (per MY)\n small - large range size",
               alpha_geom = 0.75,
               alpha_line = 0.45)



