#Michelle Jonika
#January 16, 2023

#creates a plot showing the tip rates as a function of range size category

#load in libraries needed
library(phytools)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/datalists_range.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")
#reduce to one tree from dataset that is needed
tree <- trees[[1]]
#save the order of the tree tip labels in order from 1:length of data rather
#than in phylogenetic order
test <- untangle(tree, "read.tree")
rm(trees, tree)

#load in tip rate data
tip_rates <- read.csv("../results/tiprates_new.csv")

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
rm(datalist, i, j, neworder, tip_rates, all_data, test)

#subset out species with small range size
carn_small <- data_all[data_all$range.size == "0",]
#calculate the mean population size in species with small range size
mean_smallRS <- mean(carn_small$tip_rates.Average, na.rm = T)

#subset out species with large range size
carn_large <- data_all[data_all$range.size == "1",]
#calculate the mean population size in species with large range size
mean_largeRS <- mean(carn_large$tip_rates.Average, na.rm = T)

#store the small range sizes in a vector
hitsmall <- which(data_all$range.size == "0")
#change all the small range sizes to "name "Small Range Size" instead of 0
data_all[hitsmall, 3] <- "Small Range Size"

#store the large range sizes in a vector
hitlarge <- which(data_all$range.size == "1")
#change all the small range sizes to "name "Small Range Size" instead of 0
data_all[hitlarge, 3] <- "Large Range Size"

#create a beeswarm plot of the data
beeswarm(data_all$tip_rates.Average[which(data_all$tip_rates.Average > 0.5)] ~ 
           data_all$range.size[which(data_all$tip_rates.Average > 0.5)],
         pch = 16,
         xlab = "",
         ylab = "Tip Rate",
         ylim = c(0,5.5))
lines(x = c(0.9,1.1), 
      y = rep(mean_largeRS, 2),
      lwd = 4,
      col = "red")
lines(x = c(1.9,2.1), 
      y = rep(mean_smallRS, 2), 
      lwd = 4,
      col = "red")

#save as pdf 6"x6"
