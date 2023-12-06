###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(doMC)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/pop_dens/datalists_popdens.RData")
#0 = small; 1 = large pop size

#load in tree data
trees <- read.nexus("../data/pop_dens/carnivora_pd_pruned.nex")

smap.trees <- list()
s1 <- foreach(i = 1:100) %dopar%{
  x <- datalist[[i]]$pop.dens
  names(x) <- datalist[[i]]$species
  smap.trees[[i]] <- make.simmap(tree = trees[[i]], 
              x = x, 
              model="ER", 
              nsim = 100)}

summary(s1[[1]])
as.factor(x)
cols<-setNames(c("blue","red"), unique(x))
plot(summary(s1[[100]]),colors= cols,fsize = 0.5, ftype="i", type = "fan")
legend("topleft",c("small range size","large range size"),
       pch=21,pt.bg=cols,pt.cex=2)
par(mar=c(5.1,4.1,4.1,2.1))

obj<-densityMap(s1[[100]],states=c("small RS","large RS"),plot=FALSE)
plot(obj,lwd=4,outline=TRUE,fsize=c(0.7,0.9),legend=50)

plot(s1[[1]])

