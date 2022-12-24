library(ape)
library(chromePlus)

tree <- rcoal(100)
tips <- simChrom(tree, pars=c(1,.5,0,0,20),
                 model="2010", limits = c(3,50))
tips <- tips - (min(tips)-1)
Q <- matrix(0, max(tips), max(tips))
colnames(Q) <- c(1:max(tips))
for(i in 1:4){ # rows
  for(j in 1:4){ # cols
    if(i == (j+1)){
      Q[i,j] <- .5
    }
    if(i == (j-1)){
      Q[i,j] <- 1
    } 
  } 
}

diag(Q) <- -rowSums(Q)
tip.rates <- GetTipRates(tree = tree, 
                         Q = Q, 
                         tip.states = tips, 
                         hyper = FALSE)

