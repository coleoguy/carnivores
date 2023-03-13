# Michelle Jonika
# 10 November
#This code creates a figure that visualizes the estimated population sizes for
#species with small and large range sizes

#load in libraries
library(beeswarm)

#load in data needed
carn <- read.csv("../../data/Carn_data_incomplete.csv")
#prune data to only informative columns for this plot
carn_prune <- carn[,-c(2:3,6:9,12:36)]

#subset out species with small range size
carn_small <- carn_prune[carn_prune$rangesize == "0",]
#calculate the mean population size in species with small range size
mean_smallRS <- mean(carn_small$final.pop, na.rm = T)

#subset out species with large range size
carn_large <- carn_prune[carn_prune$rangesize == "1",]
#calculate the mean population size in species with large range size
mean_largeRS <- mean(carn_large$final.pop, na.rm = T)

#store the small range sizes in a vector
hitsmall <- which(carn_prune$rangesize == "0")
#change all the small range sizes to "name "Small Range Size" instead of 0
carn_prune[hitsmall, 5] <- "Small Range Size"

#store the large range sizes in a vector
hitlarge <- which(carn_prune$rangesize == "1")
#change all the small range sizes to "name "Small Range Size" instead of 0
carn_prune[hitlarge, 5] <- "Large Range Size"

#create a beeswarm plot of the data
beeswarm(log(carn_prune$final.pop)~carn_prune$rangesize,
         pch = 16,
         xlab = "",
         ylab = "Estimated Population Size")
segments(x0 = 0.9,
         y0 = log(mean_largeRS),
         x1 = 1.1,
         y1 = log(mean_largeRS),
         lwd = 3)
segments(x0 = 1.9,
         y0 = log(mean_smallRS),
         x1 = 2.1,
         y1 = log(mean_smallRS),
         lwd = 3)

#save as pdf 6"x6"