# Michelle Jonika
# 10 November
# This code makes the figure for the analysis of chromosome number in
# carnivores with range size as the binary trait

#load in libraries needed
library(viridis)
library(chromePlus)

#read in chromplus data for range size
data <- read.csv("../results/rangesize.csv")

plot.chromeplus(data = data,
                colors = c("#FDE725FF", "#39568CFF"),
                x_title = "rate difference (per MY)\n small - large range size",
                alpha_geom = 0.75,
                alpha_line = 0.45)

#save as PDF 6"x6"
