# this script plots a comparison of delta R stats for small/large range size

#load in required libraries
library(coda)

#load in data
carn.all <- read.csv("results_25small75large/carn_data_all.csv")

#load in tree depths
depths <- read.csv("results_25small75large/tree_depths.csv")
#name depths data frame with informative column names
colnames(depths) <- c("tree", "depth")
#remove row 52, tree could not run with model 
#(polytomy; branches not distinguishable)
depths <- depths[-c(52),]


# 500 observations for 99 trees (had to drop tree 52)
x  <- seq(1, 245001, by=2500)
x+2499
counter <- 1
#transform rates back into units of millions of years by the tree depths
#for each respective tree
for(i in x){
        carn.all[i:(i+2499), ] <- carn.all[i:(i+2499), ]/depths[counter,2]
        counter <- counter + 1
}

#calculate the delta R statistic for ascending rates
fission <- carn.all$asc1-carn.all$asc2
#calculate the delta R statistic for descending rates
fusion <- carn.all$desc1-carn.all$desc2

#calculate the jhighest posterior density interval for fissions 
hpdfis <- HPDinterval(as.mcmc(fission))
#calculate the jhighest posterior density interval for fusions
hpdfus <- HPDinterval(as.mcmc(fusion))


#state 1 = small range size, state 2 = large range size
HPDinterval(as.mcmc(carn.all$asc1))
HPDinterval(as.mcmc(carn.all$asc2))

HPDinterval(as.mcmc(carn.all$desc1))
HPDinterval(as.mcmc(carn.all$desc2))

colMeans(carn.all)

cols <- c(rgb(1, 0, 0, .5), rgb(0, 1, 0, .5))
plot(0,0,col="white",
     ylim=c(-1,10),
     xlim=c(-0.6,0.6),
     xlab=expression(paste(Delta, R[x])),
     ylab="density")
abline(v=0, lty=3,col="black")
polygon(density(fission), col = cols[1])
polygon(density(fusion),col = cols[2])


foo <- cbind(carn.all, rep(1:99, each=2500))
for(i in 1:99){
        hist(foo$asc1[foo[,9]==i]-foo$asc2[foo[,9]==i],
             main=paste("tree ", i))
}
plot(foo$p[foo[,9]==6])

points(pch = 22, bg = cols,
       x = rep(-.02, 2), y = c(50, 40))
text(x = rep(-.02, 2), y = c(50, 40), pos = 4, cex = .7,
     labels=c("fission", "fusion"))
     
lines(y=rep(-1.7, 2), x=hpdfis, col=cols[1], lwd=4)
lines(y=rep(-4, 2), x=hpdfus, col=cols[2], lwd=4)
