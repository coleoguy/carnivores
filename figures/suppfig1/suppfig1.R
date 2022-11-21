# Michelle Jonika
# 10 November
#This code makes the supplemental lfigure showing our MCMC parameter space
#and the corresponding posterior proability associated with each parameter
#space

#load in libraries needed 
library(viridis)

#load in each of the replicates
load("../results/iucn/trial1_iucn.Rdata")
load("../results/iucn/trial2_iucn.Rdata")
load("../results/iucn/trial3_iucn.Rdata")
load("../results/iucn/trial4_iucn.Rdata")
load("../results/iucn/trial5_iucn.Rdata")
load("../results/iucn/trial6_iucn.Rdata")

#create an empty 
results <- list()

#make a list of the replicates
trials <- list(x1,x2,x3,x4,x5,x6)

#loop that retains the MCMC run with the highest final posterior probability
for(i in 1:length(x1)){
  #loops through each tree and stores the final posterior probability
  liks <- c(trials[[1]][[i]]$p[500],
            trials[[2]][[i]]$p[500],
            trials[[3]][[i]]$p[500],
            trials[[4]][[i]]$p[500],
            trials[[5]][[i]]$p[500],
            trials[[6]][[i]]$p[500])
  #store the tree with the highest posterior probabiliyt for the ifinal dataset
  results[[i]] <- trials[[which.max(liks)]][[i]]
}

#create an empty vector
drdesc <- drasc <- c()
#store the results of rates of fusion in small range size minus rates of fusion
#in small range size in the first tree
drdesc <- results[[1]]$desc2 - results[[1]]$desc1
#store the results of the posterior probability in the first tree
drp <- results[[1]]$p

#loop through to store each of the posterior probabilities and rates for the 
#rest of the 99 trees
for(i in 2:100){
  #store the results of rates of fusion in small range size minus rates of fusion
  #in small range size in the first tree
  drdesc <- c(drdesc, results[[i]]$desc2 - results[[i]]$desc1)
  #store the results of the posterior probability in the first tree
  drp <- c(drp, results[[i]]$p)
}


runs <- rep(1:500, 100)
names(drdesc) <- runs
#plot the parameter space explored 
op <- par(mar = c(5,6,4,2) + 0.1)
plot(drdesc~drp, 
     pch = 19,
     cex = 0.9,
     col = c(viridis(500, alpha = 0.5))[factor(names(drdesc))],
     xlab = "Posterior Probability",
     ylab = "Fusion Parameter Space in Species with Small \n Range Sizes - Species with Large Rnage Sizes")
par(op)

#save as PDF 8"x6" 

#initially plot showing the fusion parameters explored throughout the 500 
#generations of the MCMC model in the first tree
op <- par(mar = c(5,6,4,2) + 0.1)
plot(results[[1]]$desc2 - results[[1]]$desc1, 
     type="l", 
     cex=.4, 
     col=rgb(1,0,0,.2),
     ylim=c(-20,30),
     xlab = "Generation",
     ylab = "Fusion Parameter Space in Species with Small \n Range Sizes - Species with Large Rnage Sizes")
par(op)
#plot showing the fusion parameters explored throughout the 500 generations of 
#the MCMC model in the other 99 trees
for(i in 2:100){
  lines(results[[i]]$desc2 - results[[i]]$desc1, pch=16,cex=.4, col=rainbow(100, alpha=.2)[i])
}

#save as PDF 8"x6" 