#Michelle Jonika
#september 7, 2021

#creates a phylogenetic plot visualizing chromosome number and range size as a
#discrete trait

###LOAD IN PACKAGES###----------------------------------------------------------
library(viridis)

###LOAD IN DATA###--------------------------------------------------------------

#load in tree data
trees <- read.nexus("../data/carnivora.nex")
for(i in 1:100){
  #trees are ultrametric, this line corrects for the fact that the tolerance
  #for being ultrametric is not met by some trees
  trees[[i]] <- force.ultrametric(trees[[i]], method = "extend")
}

#load in tip rate data
tips <- read.csv("../results/rates.csv", row.names = 1)


###PRUNE DATA###----------------------------------------------------------------

#prune and scale trees
#find tips that are missing from the dataset
missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% tips$species]
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
rm(trees, cur.tree, i, missing, tree.depths)

###ORDER TIP LABELS AND CHROMOSOME DATA###--------------------------------------
#this loop will make sure the tip labels and the chromosome data are in the
#same order

#create an empty tip rate vector
tiprates <- c()
#store the tip rates in a vector
rates <- tips$tipRates
#name the tip rate vector
names(rates) <- tips$species
#loop through tip rates and tip labels to make sure they are in the same order 
for(i in 1:110){
  tiprates[i] <- rates[which(names(rates) == trees.pruned[[50]]$tip.label[i])]
}
#name the newly ordered vector
names(tiprates) <- trees.pruned[[50]]$tip.label
#clean up environment
rm(i, rates)

barplot(height = tips$ + 1, 
        col=viridis(2, option= "G", end = 0.6)[tips$range.size + 1],
        horiz = T,
        xlim = c(0,400),
        main = "Chromosome Number Tip Rates")
legend(x = "topright", 
       legend = c("Small Range Size", "Large Range Size"), 
       pch = 22, 
       pt.cex = 2, 
       box.col = "transparent", 
       pt.bg = viridis(2, option= "G", end = 0.6))

#export to PDF 8.5"x11"






