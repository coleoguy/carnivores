#read in the libraries needed to plot
library(tidyr)
library(ggplot2)
library(stringr)

#read in the simulation of chromosme number results
chroms <- read.csv("../../results/carn_delta_sim.csv")
chroms <- chroms[,-1]

#load in chromosome data from epirical dataset
chroms.num <- read.csv("../../data/chroms.csv")
#load in range size data from empirical dataset
range.num <- read.csv("../../results/calc_carn_range_sizes.csv")

#vector to store the final empirical chromosome number dataset
hap.chrom <- c()
#loop to prune the chromsome number to those in the final chromosome number
#dataset
for(i in 1:nrow(range.num)){
  hit <- which(chroms.num$species == range.num$X[i])
  if(length(hit) > 1){
    hit <- sample(hit, 1)
  }
  hap.chrom[i] <- chroms.num[hit, 2]
}

#bind the empirical chromosome number data into the simulatio chromosome 
#number data
chrom_all <- cbind(chroms, hap.chrom)

#add column names to the chromosome number data
colnames(chrom_all) <- c(paste("Tree", rep(1:100)), "Emp")


#plotting

#plot the first trees density
plot(density(chrom_all[,1]), 
     lwd=.2, 
     ylim=c(0,.7),
     xlim=c(7,44),
     xlab = "Chromosome Number",
     ylab = "Density",
     main = "")
#plot the other 100 trees of density
for(i in 2:101){
  if(i!=101){
    lines(density(chrom_all[,i]), lwd=.2)
  }else{
    #makes the empirical dataset (column 101) stand out with a thick red line
    lines(density(chrom_all[,i]), lwd=2,col="red")
  }
}

#export as PDF 6"x6"


