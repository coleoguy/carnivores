#plotting carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)

library(viridis)

library(castor)
library(ape)
library(geiger)
library(chromePlus)
library(diversitree)

###LOAD IN DATA###--------------------------------------------------------------
#load in chromosome number and binary trait data
load("../data/datalists_range.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in tiprates
tip_rates <- read.csv("../results/tiprates_new.csv")

#bind together the tip rate results and the range size to color by
plotting <- data.frame(tip_rates$X, tip_rates$Average, datalist[[1]]$hap.chrom, datalist[[1]]$range.size)

#plot the tip rate data
barplot(height = plotting$tip_rates.Average + 0.025,
        col = c("#FDE725FF", "#39568CFF")[plotting$datalist..1...range.size + 1],
        horiz = T,
        xlim = c(0,6),
        xlab = "Tip Rates")

#add a legend to the plot
legend(x = "topright", 
       legend = c("Small Range Size", "Large Range Size"), 
       pch = 22, 
       pt.cex = 2, 
       box.col = "transparent", 
       pt.bg = c("#FDE725FF", "#39568CFF"))

#export as PDF 8.5"x11"





