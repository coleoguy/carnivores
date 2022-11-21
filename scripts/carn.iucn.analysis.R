#testing carnivore chromosome data and range size in chromePlus

###LOAD IN PACKAGES###----------------------------------------------------------
#load in packages needed to run analysis
library(phytools)
library(chromePlus)
library(diversitree)
library(doMC)

###LOAD IN DATA NEEDED###-------------------------------------------------------

#load in chromosome number and binary trait data
load("../data/datalists_iucn.RData")

#load in tree data
trees <- read.nexus("../data/carnivorapruned.nex")

#load in tree depths
tree.depths <- read.csv("../data/treedepths.csv")

#assign prior from exponential distribution
prior <- make.prior.exponential(2)

# from primary analysis we can get our w
w <- c(14.44759, 12.13979, 16.06275, 14.23117, 2.477789, 3.194005)

# register cores to use in parallel
registerDoMC(detectCores(all.tests = T) - 25)

#create empty list to store results
result <- list()
#set iter to 500 for the number of steps to take in the model
iter <- 500

x1 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x1, file="trial1_iucn.Rdata")

x2 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x2, file="trial2_iucn.Rdata")

x3 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x3, file="trial3_iucn.Rdata")

x4 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x4, file="trial4_iucn.Rdata")

x5 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x5, file="trial5_iucn.Rdata")

x6 <- foreach(i = 1:100) %dopar%{
  #store chrom.range
  chrom.range <- range(datalist[[52]]$hap.chrom) + c(-1, 1)
  #run datatoMatrix function necessary for diversitree likelihood function
  datalist[[i]] <- datatoMatrix(x=datalist[[i]], range = chrom.range, hyper = T)
  
  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees[[i]], states = datalist[[i]],
                    k = ncol(datalist[[i]]), strict = F,
                    control = list(method = "ode"))
  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = datalist[[i]], lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))
  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = runif(6, 0, 10),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 50,
                      lower = 0)
}
save(x6, file="trial6_iucn.Rdata")

##### Processing results #########
#only process post burn-in results
iucn <- x[[1]][451:500,]
iucn[,2:7] <- iucn[,2:7]/tree.depths[1,2] 
for(i in 2:100){
  foo <- x[[i]][451:500,]
  foo[,2:7] <- foo[,2:7]/tree.depths[i,2] 
  iucn <- rbind(iucn, foo)
}
hist(iucn$asc1-iucn$asc2)
hist(iucn$desc1-iucn$desc2)

#save the results output
write.csv(iucn,file="../results/carn.iucn.csv")



#16.66098
#22.29159
#-19.81789
#24.81584
#23.93914
#18.72708
#17.01756
#24.27797
#21.50781
#20.99101
#19.27124
#18.21687
#24.32063
#17.81695
#26.69269
#19.1164
#20.69437
#18.25843
#18.57138
#24.60471
#18.60038
#24.48164


