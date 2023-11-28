# Michelle Jonika
# 10 November
#This code munges the data from the each of the ChromPlus replicates to gather 
#the final posterior dataset for visualization

#load in each of the replicates
load("../results/range_size/trial1_rs.Rdata")
load("../results/range_size/trial2_rs.Rdata")
load("../results/range_size/trial3_rs.Rdata")
load("../results/range_size/trial4_rs.Rdata")
load("../results/range_size/trial5_rs.Rdata")
load("../results/range_size/trial6_rs.Rdata")

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

#load in the ttree depths to transform back into units of millions of years (MY)
tree.depths <- read.csv("../data/range_size/rs_treedepths.csv")

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
write.csv(post.burn, file="../results/range_size/rs.csv")

