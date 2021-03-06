# Heath Blackmon
# 10 November
# This code makes the figure for the analysis of chromosome number in
# carnivores
library(coda)
library(ggplot2)
library(ggpubr)
library(viridis)
# this sets up a theme for ggplot
ggtheme <- theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              panel.border=element_blank(),
                              axis.line = element_line(colour="grey30"),
                              axis.title = element_text(colour="grey20"),
                              axis.text = (element_text(colour="grey30")),
                              legend.title = element_text(colour="grey20"),
                              legend.text = element_text(colour="grey30"))

#read in chromplus data
dat <- read.csv("../results/carn.med.hb.csv")
#format data in long format
longdat <- data.frame(c(dat$asc2-dat$asc1,
                        dat$desc2-dat$desc1),
                      rep(c("fission","fusion"), each=10000))
#assign column names based on the new format
colnames(longdat) <- c("rate","type")
#store data in a different data structure with a better name
dat <- longdat
#remove old data structure to clean environment
rm(longdat)
#calculate the highest posterior density internvals for each rate
hpd <- data.frame(X = c(HPDinterval(as.mcmc(dat$rate[dat$type=="fission"]))[1,],
                        HPDinterval(as.mcmc(dat$rate[dat$type=="fusion"]))[1,]),
                   Y = c(-1, -1, -2, -2),
                   types = rep(c("fission", "fusion"), each = 2))

#plot the rates along with the HPD for fission and fusion
ggplot(dat, aes(x=rate)) +
  geom_density(aes(fill=as.factor(type),y=..density..),
               stat="density", position="identity", alpha=0.35) +
  geom_line(data=hpd, aes(x=X, y=Y, color=types),
            alpha=0.45, size=1.4, lineend="round") +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  scale_fill_viridis_d(option="G", end=.6)+
  scale_color_viridis_d(option="G", end=.6)+
  guides(fill=guide_legend(title="parameter"),
         color="none") +
  xlab("rate difference (per MY)\n small - large range size") +
  ggtheme

#save as PDF 6"x6"
