# Michelle Jonika
# 10 November
#This code munges the data from the each of the ChromPlus replicates to gather 
#the final posterior dataset for visualization

#load in each of the replicates
load("../results/iucn/trial1_iucn.Rdata")
load("../results/iucn/trial2_iucn.Rdata")
load("../results/iucn/trial3_iucn.Rdata")
load("../results/iucn/trial4_iucn.Rdata")
load("../results/iucn/trial5_iucn.Rdata")
load("../results/iucn/trial6_iucn.Rdata")
load("../results/iucn/trial7_iucn.Rdata")
load("../results/iucn/trial8_iucn.Rdata")
load("../results/iucn/trial9_iucn.Rdata")

#create an empty 
results <- list()

#make a list of the replicates
trials <- list(x1,x2,x3,x4,x5,x6,x7,x8,x9)

#loop that retains the MCMC run with the highest final posterior probability
for(i in 1:length(x1)){
  #loops through each tree and stores the final posterior probability
  liks <- c(trials[[1]][[i]]$p[500],
            trials[[2]][[i]]$p[500],
            trials[[3]][[i]]$p[500],
            trials[[4]][[i]]$p[500],
            trials[[5]][[i]]$p[500],
            trials[[6]][[i]]$p[500],
            trials[[7]][[i]]$p[500],
            trials[[8]][[i]]$p[500],
            trials[[9]][[i]]$p[500])
  #store the tree with the highest posterior probabiliyt for the ifinal dataset
  results[[i]] <- trials[[which.max(liks)]][[i]]
}

#assign variable to store those runs that don't converge
uncon2 <- c()
#loops through to identify runs that have a low probability and don't reach 
#convergence
for(i in 1:100){
  #looks at the mean rates of the two ascending rates from the model to see if 
  #they are less than 0
  if(mean(x1[[i]]$asc2[451:500])-mean(x1[[i]]$asc1[451:500]) < 0){
    #stores the current run that didn't meet convergence
    uncon <- c(uncon, i)
  }
}

#loop that swaps rates between high and low pop to see if that impacts convergence
result.uncon <- c()
for(i in uncon){
  j <- uncon[i]
  
  #store chrom.range
  chrom.range <- range(datalist[[j]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[j]] <- datatoMatrix(x=datalist[[j]], range = chrom.range, hyper = T)
  
  #make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[j]], states = datalist[[j]],
                    k = ncol(datalist[[j]]), strict = F,
                    control = list(method = "ode"))
  
  #constrain our model to be biologically realistic for chromosomes
  con.lk.mk<-constrainMkn(data = datalist[[j]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  
  # now we are ready to run our inference run
  result.uncon[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}

#assign variable to store those runs that don't converge
uncon2 <- c()
#loops through to identify runs that have a low probability and don't reach 
#convergence
for(i in 1:100){
  #looks at the mean rates of the two ascending rates from the model to see if 
  #they are less than 0
  if(mean(result.uncon[[i]]$asc2[451:500])-mean(result.uncon[[i]]$asc1[451:500]) < 0){
    #stores the current run that didn't meet convergence
    uncon2 <- c(uncon2, i)
  }
}

#load in the ttree depths to transform back into units of millions of years (MY)
tree.depths <- read.csv("../data/treedepths.csv")

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
write.csv(post.burn, file="../results/iucn.csv")
