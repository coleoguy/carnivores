#Michelle Jonika
#December 27, 2022

#creates a phylogenetic plot visualizing tip rates and range size as a
#discrete trait

#load in libraries needed
library(phytools)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../../data/range_size/datalists_rangesize.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../../data/carnivora_rs_pruned.nex")
#reduce to one tree from dataset that is needed
tree <- trees[[1]]
#save the order of the tree tip labels in order from 1:length of data rather
#than in phylogenetic order
test <- untangle(tree, "read.tree")
rm(trees, tree)

#load in tip rate data
tip_rates <- read.csv("../../results/tiprates.csv")

#put the data together
all_data <- data.frame(datalist[[1]], tip_rates$Average)

#loop through to make sure all of the datalists are in the same order as each
#of the trees
data_all <- c()
for(i in 1:length(all_data)){
  #empty vector to store correct order in 
  neworder <- c()
  #loop through to store the correct order of the data to match the tree tips
  for(j in 1:length(test$tip.label)){
    neworder[j] <- which(test$tip.label[j] == all_data$species)
  }
  #reorder the data into a new data frame
  data_all <- all_data[neworder,]
}
rm(datalist, i, j, neworder, tip_rates, all_data)

#make a plot of the chromosome nuber data
plot(y = 1:110, x = data_all$tip_rates.Average, 
     xlab = "Tip Rates",
     xlim = c(0,5.5), 
     pch = 21,
     cex = 0.9,
     col = "black",
     bg = c("#FDE725FF", "#39568CFF")[data_all$range.size + 1])

#export as PDF 8.5"x11"

