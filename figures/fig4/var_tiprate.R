#Michelle Jonika
#This script will plot a visual that shows the varitation in tip rates between
#our two binary states



###LOAD IN AND PRUNE DATA###----------------------------------------------------
dat <- read.csv("../results/final.tiprates.csv")[,-1]

dat[dat$range.size == 0, ]$range.size <- "Small"
dat[dat$range.size == 1, ]$range.size <- "Large"

library(beeswarm)
beeswarm(dat$tip.rate~dat$range.size)

boxplot(log(dat$tip.rate + 1) ~ dat$range.size,
                xlab = "Range Size",
                ylab = "Tip Rate")
stripchart(log(dat$tip.rate) ~ dat$range.size,
           method = "jitter",
           pch = 19, 
           col = c("#39568CFF", "#FDE725FF"), 
           vertical = T, 
           add = T)


plot(dat$tip.rate~dat$hap.chrom,xlim=c(0,41))
