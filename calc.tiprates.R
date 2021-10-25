#load in necessary packages
library(phytools)
library(geiger)

#read in insect phylogeny
trees <- read.nexus("../data/post.nex")

#read in the microsatellite data
dat.mic <- read.csv("../results/micRocounter_results_TII.csv",
                    as.is = T, row.names = 4)


# match up tree and data
trees.pruned <- list()

#loop through to drop any unmatching data or tree tips
for(i in 1:100){
  trees.pruned[[i]] <- treedata(phy = trees[[i]], data = dat.mic)[[1]]
}
rm(trees, i)

#pull out data needed for calculating tip rates
bpMbp <- dat.mic$all#/(dat.mic$gsz/1000)
names(bpMbp) <- row.names(dat.mic)
rm(dat.mic)


#loop through calculating 
for(k in 1:100){
  print(k)
  # estimate ancestral states
  foo <- ace(x = bpMbp, phy = trees.pruned[[k]], model = "BM")
  # get tip branches
  tip.branch <- c()
  for(i in 1:nrow(trees.pruned[[k]]$edge)){
    val <- trees.pruned[[k]]$edge[i, 2]
    if(!val %in% trees.pruned[[k]]$edge[, 1]){
      tip.branch <- c(tip.branch, i)
    }
  }
  
  anc.state <- c()
  for(j in 1:length(tip.branch)){
    # get node microsatellite estimate
    trees.pruned[[k]]$edge.length[tip.branch]
    wanted.node <- trees.pruned[[k]]$edge[tip.branch[j], 1]
    anc.state[j] <- as.numeric(foo$ace[names(foo$ace) == wanted.node])
  }
  curr.state <- c()
  for(j in 1:length(tip.branch)){
    curr.sp <- trees.pruned[[k]]$tip.label[j]
    curr.state <- bpMbp[names(bpMbp) == curr.sp]
  }
  print(anc.state - curr.state)
  tip.rates <- (anc.state - curr.state) /
    trees.pruned[[k]]$edge.length[tip.branch]
  names(tip.rates) <- trees.pruned[[k]]$tip.label
  if(k == 1){
    tipp.rates <- tip.rates
  }else{
    tipp.rates <- cbind(tipp.rates, tip.rates)
  }
}
colnames(tipp.rates) <- paste("tree", 1:100)

#tipp.rates[,] <- abs(tipp.rates)
Average <- rowSums(tipp.rates)/100
tipp.rates <- cbind(tipp.rates, Average)
write.csv(tipp.rates, file = "../results/tip.rates.csv")
