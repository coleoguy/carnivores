# Michelle Jonika
# 10 November
#This code munges the data from the each of the ChromPlus replicates to gather 
#the final posterior dataset for visualization

#load in each of the replicates
load("../../results/body_mass/trial1_bm.Rdata")
load("../../results/body_mass/trial2_bm.Rdata")
load("../../results/body_mass/trial3_bm.Rdata")
load("../../results/body_mass/trial4_bm.Rdata")
load("../../results/body_mass/trial5_bm.Rdata")
load("../../results/body_mass/trial6_bm.Rdata")

#create an empty 
results <- list()

#make a list of the replicates
trials <- list(x1,x2,x3,x4,x5,x6)

#loop that retains the MCMC run with the highest final posterior probability
for(i in 1:length(x1)){
  #loops through each tree and stores the final posterior probability
  liks <- c(trials[[1]][[i]]$p[1000],
            trials[[2]][[i]]$p[1000],
            trials[[3]][[i]]$p[1000],
            trials[[4]][[i]]$p[1000],
            trials[[5]][[i]]$p[1000],
            trials[[6]][[i]]$p[1000])
  #store the tree with the highest posterior probabiliyt for the final dataset
  results[[i]] <- trials[[which.max(liks)]][[i]]
}

#load in the ttree depths to transform back into units of millions of years (MY)
tree.depths <- read.csv("../../data/body_mass/bm_treedepths.csv")

#only process post burn-in results
post.burn <- results[[1]][951:1000, 2:8]

#transform the post burn in results back into MY from the tree depths
post.burn[,1:6] <- post.burn[,1:6]/tree.depths[1,2]

#loop that organizes rates into a pretty table
for(i in 2:100){
  #pulls results for specific tree
  temp <- results[[i]]
  #transforms the rates back by their tree depths
  temp[,2:7] <- temp[,2:7]/tree.depths[i,2]
  ##bind in post-burn-in samples from each ttree after they have been back 
  #transformed
  post.burn <- rbind(post.burn, temp[951:1000,2:8])
}

#save the results output
write.csv(post.burn, file="../../results/body_mass/bm.csv")

#assign variable to store those runs that don't converge
uncon1 <- c()
#loops through to identify runs that have a low probability and don't reach 
#convergence



for(i in 1:100){
  #looks at the mean rates of the two ascending rates from the model to see if 
  #they are less than 0
  if(mean(results[[i]]$asc2[950:1000])-mean(results[[i]]$asc1[950:1000]) < 0){
    #stores the current run that didn't meet convergence
    uncon1 <- c(uncon1, i)
  }
}
