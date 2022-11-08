carn <- read.csv("Carn_data_incomplete.csv")
carn_prune <- carn[,-c(2:3,6:9,12:36)]


library(beeswarm)
beeswarm(log(carn_prune$final.pop)~carn_prune$rangesize,
         pch = 16,
         pwcol = factor(carn_prune$Family),
         xlab = "Range Size",
         ylab = "Estimated Population Size")
legend("topright", 
       legend = factor(unique(carn_prune$Family)),
       col = 1:14,
       pch = 16)
