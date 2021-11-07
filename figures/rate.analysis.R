# Heath Blackmon
# 10 November
# This code makes the figure for the analysis of chromosome number in
# carabidae
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

# This reads in the data without polyploidy

dat <- read.csv("../results/carn.med.hb.csv")
longdat <- data.frame(c(dat$asc2-dat$asc1,
                        dat$desc2-dat$desc1),
                      rep(c("fission","fusion"), each=10000))
colnames(longdat) <- c("rate","type")
dat <- longdat
rm(longdat)
hpd <- data.frame(X = c(HPDinterval(as.mcmc(dat$rate[dat$type=="fission"]))[1,],
                        HPDinterval(as.mcmc(dat$rate[dat$type=="fusion"]))[1,]),
                   Y = c(-1, -1, -2, -2),
                   types = rep(c("fission", "fusion"), each = 2))


ggplot(dat, aes(x=rate)) +
  geom_density(aes(fill=as.factor(type),y=..density..),
               stat="density", position="identity", alpha=0.35) +
  geom_line(data=hpd, aes(x=X, y=Y, color=types),
            alpha=0.45, size=1.4, lineend="round") +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  scale_fill_viridis_d(option="D", end=.6)+
  scale_color_viridis_d(option="D", end=.6)+
  guides(fill=guide_legend(title="parameter"),
         color="none") +
  xlab("Rate difference (per my)\n small - large range size") +
  ggtheme

