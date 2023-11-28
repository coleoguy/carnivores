# Michelle Jonika
# 10 November
# This code makes the figure for the analysis of chromosome number in
# carnivores with range size as the binary trait

#load in libraries needed
library(chromePlus)

#read in chromplus data for range size
data <- read.csv("../results/body_mass/initial_bm.csv")
fission <- data[,4] - data[,2]
fusion <- data[,5] - data[,3]
data_munge <- data.frame(fission, fusion)

plotChromeplus(data = data_munge,
                colors = c("#FDE725FF", "#39568CFF"),
                main_title = "",
               x_title = "rate difference (per MY)\n small - large body size",
                alpha_geom = 0.75,
                alpha_line = 0.45)

#save as PDF 6"x6"
