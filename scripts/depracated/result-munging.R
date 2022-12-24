# processing and choosing among MCMC runs

load("../results/rangesize/trial1.Rdata")
load("../results/rangesize/trial2.Rdata")
load("../results/rangesize/trial3.Rdata")
load("../results/rangesize/trial4.Rdata")
load("../results/rangesize/trial5.Rdata")
load("../results/rangesize/trial6.Rdata")

results <- list()

trials <- list(x1,x2,x3,x4,x5,x6)
for(i in 1:length(x1)){
  liks <- c(trials[[1]][[i]]$p[500],
            trials[[2]][[i]]$p[500],
            trials[[3]][[i]]$p[500],
            trials[[4]][[i]]$p[500],
            trials[[5]][[i]]$p[500],
            trials[[6]][[i]]$p[500])
  results[[i]] <- trials[[which.max(liks)]][[i]]
}

tree.depths <- read.csv("../data/treedepths.csv")

##### Processing results #########
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
write.csv(post.burn, file="../results/rangesize.csv")

plot(results[[1]]$p, type="l", ylim=c(-450, -370))

drdesc <- drasc <- c()
drdesc <- results[[1]]$desc2 - results[[1]]$desc1
drp <- results[[1]]$p
for(i in 2:100){
  drdesc <- c(drdesc, results[[i]]$desc2 - results[[i]]$desc1)
  drp <- c(drp, results[[i]]$p)
}
plot(drdesc~drp, pch=16, cex=.5,col=rgb(0,0,0,.12))

plot(results[[1]]$desc2 - results[[1]]$desc1, type="l", cex=.4, col=rgb(1,0,0,.2),ylim=c(-20,30))
for(i in 2:100){
  lines(results[[i]]$desc2 - results[[i]]$desc1, pch=16,cex=.4, col=rainbow(100, alpha=.2)[i])
}

plot(x1[[1]]$asc1~x1[[1]]$asc2)

plot(x2[[1]]$asc2 - x2[[1]]$asc1, pch=16, cex=.4, col=rgb(1,0,0,.2),ylim=c(-20,20))
for(i in 2:100){
  lines(x2[[i]]$asc2 - x2[[i]]$asc1, pch=16,cex=.4, col=rainbow(100, alpha=.2)[i])
}

for(i in 2:100){
  lines(results[[i]]$p, col=rainbow(100,alpha=.5)[i])
}






