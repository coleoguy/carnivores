# Michelle Jonika
# 10 November
#This code munges the data from the each of the ChromPlus replicates to gather 
#the final posterior dataset for visualization

#load in each of the replicates
load("../scripts/trial1_bm.Rdata")

#load in the ttree depths to transform back into units of millions of years (MY)
tree.depths <- read.csv("../data/body_mass/bm_treedepths.csv")

#only process post burn-in results
post.burn <- results[[1]][451:500, 2:8]

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
  post.burn <- rbind(post.burn, temp[451:500,2:8])
}

#save the results output
write.csv(post.burn, file="../results/body_mass/initial_bm.csv")

